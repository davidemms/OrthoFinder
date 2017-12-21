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

import time
import multiprocessing as mp
import Queue

import util

def _I_Spawn_Processes(message_queue, cmds_queue):
    """
    Args:
        message_queue - for passing messages that a new queue of tasks should be started (PTM -> I_Space_Processes) or that the tasks are complete
        cmds_queue - queue containing tasks that should be done
    Use:
        A process should be started as early as possible (while RAM usage is low) with this method as its target.
        This is now a separate process with low RAM usage.
        Each time some parallel work is required then the queue for that is placed in the message_queue by the PTM.
        _I_Spawn_Processes - will spawn parallel processes when instructed by the message_queue in the message_queue and get them 
        working on the queue. When the queue is empty it will wait for the next one. It can recieve a special signal to exit - the None
        object
    """
    while True:
        try:
            # peak in qoq - it is the only mehtod that tried to remove things from the queue
            message = message_queue.get()
            if message == 'Done': 
                # for saftey's sake put this check in. This message was intended for the PTM so put it back on the queue
                message_queue.put('Done')
                raise Queue.Empty
            if message == None: return
            # In which case, thread has been informed that there are tasks in the queue.
            nParallel, nTasks, qShell, qListOfLists, qHideStdout = message
            if qListOfLists:
                runningProcesses = [mp.Process(target=util.Worker_RunCommand, args = (cmds_queue, nParallel, nTasks, qShell, qHideStdout)) for i_ in xrange(nParallel)]           
            else:
                runningProcesses = [mp.Process(target=util.Worker_RunOrderedCommandList, args = (cmds_queue, nParallel, nTasks, qShell, qHideStdout)) for i_ in xrange(nParallel)] 
            for proc in runningProcesses:
                proc.start()
            for proc in runningProcesses:
                while proc.is_alive():
                    proc.join() 
            message_queue.put("Done")
            time.sleep(2)
        except Queue.Empty:
            time.sleep(4) # there wasn't anything this time, sleep then try again
    pass
    

class ParallelTaskManager_singleton:
    class __Singleton(object):
        def __init__(self):
            """Implementation:
            Allocate a thread that will perform all the tasks
            Communicate with it using a queue. 
            When provided with a list of commands it should fire up some workers and get them to run the commands and then exit.
            An alternative would be they should always stay alive - but then they could die for some reason? And I'd have to check how many there are.
            """
            self.message_queue = mp.Queue()   
            # Orders/Messages:
            # None (PTM -> spawn_thread) - thread should return (i.e. exit)
            # 'Done' (spawn_thread -> PTM) - the cmds from the cmd queue have completed
            # Anything else = (nParallel, nTasks) (PTM -> spawn_thread) - cmds (nTasks of them) have been placed in the cmd queue, 
            #   they should be executed using nParallel threads
            self.cmds_queue = mp.Queue()
            self.manager_process = mp.Process(target=_I_Spawn_Processes, args=(self.message_queue, self.cmds_queue))
            self.manager_process.start()
    instance = None
    
    def __init__(self):
        if not ParallelTaskManager_singleton.instance:
            ParallelTaskManager_singleton.instance = ParallelTaskManager_singleton.__Singleton()
        
    def RunParallel(self, cmd_list, qListOfLists, nParallel, qShell, qHideStdout):
        """
        Args:
            cmd_list - list of commands or list of lists of commands (in which elelments in inner list must be run in order)
            qListOfLists - is cmd_lists a list of lists
            nParallel - number of parallel threads to use
            qShell - should the tasks be run in a shell
        """    
        nTasks = len(cmd_list)
        for i, x in enumerate(cmd_list):
            self.instance.cmds_queue.put((i, x))
        self.instance.message_queue.put((nParallel, nTasks, qShell, qListOfLists, qHideStdout))
        while True:
            try:
                signal = self.instance.message_queue.get()
                if signal == "Done": 
                    return
            except Queue.Empty:
                time.sleep(1)
                print("PTM waiting for work to complete")
                
    def Stop(self):
        """Warning, cannot be restarted"""
        self.instance.message_queue.put(None)
        self.instance.manager_process.join()
        

# Create PTM right at start
ptm_initialised = ParallelTaskManager_singleton()
