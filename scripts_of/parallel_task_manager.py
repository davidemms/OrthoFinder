# -*- coding: utf-8 -*-
#
# Copyright 2017 David Emms
#
# This program (OrthoFinder) is distributed under the terms of the GNU General Public License v3
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#  When publishing work that uses OrthoFinder please cite:
#      Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
#      improves orthogroup inference accuracy, Genome Biology 16:157
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.com
import os
import sys
import time
import types
import datetime
import traceback
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, wait, as_completed

try:
    import queue
except ImportError:
    import Queue as queue     

# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f)


PY2 = sys.version_info <= (3,)

if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# use user's executables by preference over OrthoFinder's
my_env = os.environ.copy()
my_env['PATH'] = my_env['PATH'] + ":" + os.path.join(__location__, 'bin')
# Fix LD_LIBRARY_PATH when using pyinstaller 
if getattr(sys, 'frozen', False):
    if 'LD_LIBRARY_PATH_ORIG' in my_env:
        my_env['LD_LIBRARY_PATH'] = my_env['LD_LIBRARY_PATH_ORIG']  
    else:
        my_env['LD_LIBRARY_PATH'] = ''  
    if 'DYLD_LIBRARY_PATH_ORIG' in my_env:
        my_env['DYLD_LIBRARY_PATH'] = my_env['DYLD_LIBRARY_PATH_ORIG']  
    else:
        my_env['DYLD_LIBRARY_PATH'] = ''    


def print_traceback(e):
    if PY2:
        traceback.print_exc()
    else:
        traceback.print_tb(e.__traceback__)


def stderr_exempt(stderr):
    ok_line_starts = {"diamond v", "Licensed under the GNU GPL", "Check http://github.com/"}
    try:
        stderr = stderr.decode()
    except (UnicodeDecodeError, AttributeError):
        stderr = stderr.encode()
    lines = stderr.split("\n")
    for line in lines:
        if line.rstrip() == "": continue
        if any(line.startswith(x) for x in ok_line_starts): continue
        return False
    return True

def PrintTime(message):
    print((str(datetime.datetime.now()).rsplit(".", 1)[0] + " : " + message))
    sys.stdout.flush()

def PrintNoNewLine(text):
    sys.stdout.write(text)

def ManageQueue(runningProcesses, cmd_queue):
    """Manage a set of runningProcesses working through cmd_queue.
    If there is an error the exit all processes as quickly as possible and 
    exit via Fail() methods. Otherwise return when all work is complete
    """            
    # set all completed processes to None
    qError = False
#    dones = [False for _ in runningProcesses]
    nProcesses = len(runningProcesses)
    nProcesses_list = list(range(nProcesses))
    while True:
        if runningProcesses.count(None) == len(runningProcesses): break
        time.sleep(.1)
#        for proc in runningProcesses:
        for i in nProcesses_list:
            proc = runningProcesses[i]
            if proc == None: continue
            if not proc.is_alive():
                if proc.exitcode != 0:
                    qError = True
                    while True:
                        try:
                            cmd_queue.get(True, .1)
                        except queue.Empty:
                            break
                runningProcesses[i] = None
    if qError:
        Fail()


def RunCommand_Simple(command):
    subprocess.call(command, env=my_env, shell=True)


def RunCommand(command, qPrintOnError=False, qPrintStderr=True):
    """ Run a single command """
    popen = subprocess.Popen(command, env=my_env, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if qPrintOnError:
        stdout, stderr = popen.communicate()
        if popen.returncode != 0:
            print(("\nERROR: external program called by OrthoFinder returned an error code: %d" % popen.returncode))
            print(("\nCommand: %s" % command))
            print(("\nstdout:\n%s" % stdout))
            print(("stderr:\n%s" % stderr))
        elif qPrintStderr and len(stderr) > 0 and not stderr_exempt(stderr):
            print("\nWARNING: program called by OrthoFinder produced output to stderr")
            print(("\nCommand: %s" % command))
            print(("\nstdout:\n%s" % stdout))
            print(("stderr:\n%s" % stderr))
        return popen.returncode
    else:
        popen.communicate()
        return popen.returncode

def CanRunCommand(command, qAllowStderr = False, qPrint = True, qRequireStdout=True, qCheckReturnCode=False):
    if qPrint: PrintNoNewLine("Test can run \"%s\"" % command)       # print without newline
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
    capture.wait()
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]
    if qCheckReturnCode:
        return_code_check = (capture.returncode == 0)
    else:
        return_code_check = True
    if (len(stdout) > 0 or not qRequireStdout) and (qAllowStderr or len(stderr) == 0) and return_code_check:
        if qPrint: print(" - ok")
        return True
    else:
        if qPrint: print(" - failed")
        if not return_code_check:
            print("Returned a non-zero code: %d" % capture.returncode)
        print("\nstdout:")        
        for l in stdout: print(l)
        print("\nstderr:")        
        for l in stderr: print(l)
        return False

q_print_first_traceback_0 = False
def Worker_RunCommands_And_Move(cmd_and_filename_queue, nProcesses, nToDo, qListOfLists, q_print_on_error, q_always_print_stderr):
    """
    Continuously takes commands that need to be run from the cmd_and_filename_queue until the queue is empty. If required, moves 
    the output filename produced by the cmd to a specified filename. The elements of the queue can be single cmd_filename tuples
    or an ordered list of tuples that must be run in the provided order.
  
    Args:
        cmd_and_filename_queue - queue containing (cmd, actual_target_fn) tuples (if qListOfLists is False) or a list of such 
            tuples (if qListOfLists is True). Alternatively, 'cmd' can be a python fn and actual_target_fn the fn to call it on.
        nProcesses - the number of processes that are working on the queue.
        nToDo - The total number of elements in the original queue
        qListOfLists - Boolean, whether each element of the queue corresponds to a single command or a list of ordered commands
        qShell - Boolean, should a shell be used to run the command.
        
    Implementation:
        nProcesses and nToDo are used to print out the progress.
    """
    while True:
        try:
            i, command_fns_list = cmd_and_filename_queue.get(True, 1)
            nDone = i - nProcesses + 1
            if nDone >= 0 and divmod(nDone, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
                PrintTime("Done %d of %d" % (nDone, nToDo))
            if not qListOfLists:
                command_fns_list = [command_fns_list]
            for command, fns in command_fns_list:
                if isinstance(command, types.FunctionType):
                    # This will block the process, but it is ok for trimming, it takes minimal time
                    fn = command
                    fn(fns)
                else:
                    if not isinstance(command, str):
                        print("ERROR: Cannot run command: " + str(command))
                        print("Please report this issue.")
                    else:
                        RunCommand(command, qPrintOnError=q_print_on_error, qPrintStderr=q_always_print_stderr)
                        if fns != None:
                            actual, target = fns
                            if os.path.exists(actual):
                                os.rename(actual, target)
        except queue.Empty:
            return     
        except Exception as e:
            print("WARNING: ")
            print(str(e))
            global q_print_first_traceback_0
            if not q_print_first_traceback_0:
                print_traceback(e)
                q_print_first_traceback_0 = True
        except:
            print("WARNING: Unknown caught unknown exception")

q_print_first_traceback_1 = False
def Worker_RunMethod(Function, args_queue):
    while True:
        try:
            args = args_queue.get(True, .1)
            Function(*args)
        except queue.Empty:
            return 
        except Exception as e:
            print("Error in function: " + str(Function))
            global q_print_first_traceback_1
            if not q_print_first_traceback_1:
                print_traceback(e)
                q_print_first_traceback_1 = True
            return

def RunMethodParallel(Function, args_queue, nProcesses):
    runningProcesses = [mp.Process(target=Worker_RunMethod, args=(Function, args_queue)) for i_ in range(nProcesses)]
    for proc in runningProcesses:
        proc.start()
    ManageQueue(runningProcesses, args_queue)


def _I_Spawn_Processes(message_to_spawner, message_to_PTM):
    """
    Args:
        message_queue - for passing messages that a new queue of tasks should be started (PTM -> I_Space_Processes) or that the tasks are complete
        cmds_queue - queue containing tasks that should be done
    Use:
        A process should be started as early as possible (while RAM usage is low) with this method as its target.
        This is now a separate process with low RAM usage.
        Each time some parallel work is required then the queue for that is placed in the message_queue by the PTM.
        _I_Spawn_Processes - will spawn parallel processes when instructed by the message_queue in the message_queue and get them
        working on the queue. When the queue is empty it will wait for the next one. It can receive a special signal to exit - the None
        object
    """
    while True:
        try:
            # peak in qoq - it is the only method that tried to remove things from the queue
            message = message_to_spawner.get(timeout=.1)
            if message is None:
                # Respond to request to terminate
                return
            # In which case, thread has been informed that there are tasks in the queue.
            func, args_list, n_parallel = message
            futures = []
            n_to_do = len(args_list)
            # Version 1: n worker threads for executing the method and a list of N arguments for calling the method
            with ProcessPoolExecutor(n_parallel) as pool:
                for args in args_list:
                    futures.append(pool.submit(func, *args))
                # for i, _ in as_completed(futures):
                #     n_done = i+1
                #     if n_done >= 0 and divmod(n_done, 10 if n_done <= 200 else 100 if n_done <= 2000 else 1000)[1] == 0:
                #         PrintTime("Done %d of %d" % (n_done, n_to_do))
            # Version 2: launch n worker threads each executing a worker method that takes tasks from a queue
            with ProcessPoolExecutor(n_parallel) as pool:
                for args in args_list:
                    futures.append(pool.submit(func, *args))
            wait(futures)
            message_to_PTM.put("Done")
            time.sleep(1)
        except queue.Empty:
            time.sleep(1)  # there wasn't anything this time, sleep then try again
    pass

class ParallelTaskManager_singleton:
    """
    Creating new process requires forking parent process and can lea to very high RAM usage. One way to mitigate this is
    to create the pool of processes as early in execution as possible so that the memory footprint is low. The
    ParallelTaskManager takes care of that, and can be used by calling `RunParallelOrderedCommandLists` above.
    Apr 2023 Update:
    When running external programs there is no need to use multiprocessing, multithreading is sufficient since new process
    will be created anyway, so the SIL is no longer an issue.
    """
    class __Singleton(object):
        def __init__(self):
            """Implementation:
            Allocate a thread that will perform all the tasks
            Communicate with it using a queue.
            When provided with a list of commands it should fire up some workers and get them to run the commands and then exit.
            An alternative would be they should always stay alive - but then they could die for some reason? And I'd have to check how many there are.
            """
            self.message_to_spawner = mp.Queue()
            self.message_to_PTM = mp.Queue()
            # Orders/Messages:
            # None (PTM -> spawn_thread) - thread should return (i.e. exit)
            # 'Done' (spawn_thread -> PTM) - the cmds from the cmd queue have completed
            # Anything else = (nParallel, nTasks) (PTM -> spawn_thread) - cmds (nTasks of them) have been placed in the cmd queue,
            #   they should be executed using nParallel threads
            self.manager_process = mp.Process(target=_I_Spawn_Processes, args=(self.message_to_spawner, self.message_to_PTM))
            self.manager_process.start()
    instance = None

    def __init__(self):
        if not ParallelTaskManager_singleton.instance:
            ParallelTaskManager_singleton.instance = ParallelTaskManager_singleton.__Singleton()

    def RunParallel(self, func, args_list, nParallel):
        """
        Args:
            cmd_list - list of commands or list of lists of commands (in which elements in inner list must be run in order)
            nParallel - number of parallel threads to use
            qShell - should the tasks be run in a shell
        """
        self.instance.message_to_spawner.put((func, args_list, nParallel))
        while True:
            try:
                signal = self.instance.message_to_PTM.get()
                if signal == "Done":
                    return
            except queue.Empty:
                pass
            time.sleep(1)

    def Stop(self):
        """Warning, cannot be restarted"""
        self.instance.message_to_spawner.put(None)
        self.instance.manager_process.join()


def RunParallelMethods(func, args_list, nProcesses):
    """nProcesss - the number of processes to run in parallel
    commands - list of lists of commands where the commands in the inner list are completed in order (the i_th won't run until
    the i-1_th has finished).
    """
    ptm = ParallelTaskManager_singleton()
    ptm.RunParallel(func, args_list, nProcesses)

def Success():
    ptm = ParallelTaskManager_singleton()
    ptm.Stop()
    sys.exit()

def Fail():
    sys.stderr.flush()
    ptm = ParallelTaskManager_singleton()
    ptm.Stop()
    print("ERROR: An error occurred, ***please review the error messages*** they may contain useful information about the problem.")
    sys.exit(1)


