# -*- coding: utf-8 -*-
#
# Copyright 2014 David Emms
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

nThreadsDefault = 16
nAlgDefault = 1

import os
import sys
import glob
import time
import numpy as np
import subprocess
import datetime
import Queue
import multiprocessing as mp
from collections import namedtuple

import tree, parallel_task_manager

"""
Utilities
-------------------------------------------------------------------------------
"""
SequencesInfo = namedtuple("SequencesInfo", "nSeqs nSpecies speciesToUse seqStartingIndices nSeqsPerSpecies")

picProtocol = 1
version = "2.3.0"

# Fix LD_LIBRARY_PATH when using pyinstaller 
my_env = os.environ.copy()
if getattr(sys, 'frozen', False):
    if 'LD_LIBRARY_PATH_ORIG' in my_env:
        my_env['LD_LIBRARY_PATH'] = my_env['LD_LIBRARY_PATH_ORIG']  
    else:
        my_env['LD_LIBRARY_PATH'] = ''  
    if 'DYLD_LIBRARY_PATH_ORIG' in my_env:
        my_env['DYLD_LIBRARY_PATH'] = my_env['DYLD_LIBRARY_PATH_ORIG']  
    else:
        my_env['DYLD_LIBRARY_PATH'] = ''    
    
def PrintNoNewLine(text):
    sys.stdout.write(text)

def PrintTime(message):
    print(str(datetime.datetime.now()).rsplit(".", 1)[0] + " : " + message)      

"""
Command & parallel command management
-------------------------------------------------------------------------------
"""

def RunCommand(command, qShell=True, qHideOutput = True):
    """ Run a single command """
    if qHideOutput:
        subprocess.call(command, env=my_env, shell=qShell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        subprocess.call(command, env=my_env, shell=qShell)
            
def RunOrderedCommandList(commandList, qShell=True, qHideStdout = True):
    """ Run a list of commands """
    FNULL = open(os.devnull, 'w')
    if qHideStdout:
        for cmd in commandList:
            subprocess.call(cmd, shell=qShell, stdout=subprocess.PIPE, stderr=FNULL, close_fds=True, env=my_env)
    else:
        for cmd in commandList:
            subprocess.call(cmd, shell=qShell, stderr=FNULL, close_fds=True, env=my_env)
    
def CanRunCommand(command, qAllowStderr = False, qPrint = True):
    if qPrint: PrintNoNewLine("Test can run \"%s\"" % command)       # print without newline
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]
#    print(stdout)
#    print(stderr)
    if len(stdout) > 0 and (qAllowStderr or len(stderr) == 0):
        if qPrint: print(" - ok")
        return True
    else:
        if qPrint: print(" - failed")
        return False
        
def Worker_RunCommand(cmd_queue, nProcesses, nToDo, qShell=True, qHideStdout=True):
    """ Run commands from queue until the queue is empty """
    while True:
        try:
            i, command = cmd_queue.get(True, 1)
            nDone = i - nProcesses + 1
            if nDone >= 0 and divmod(nDone, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
                PrintTime("Done %d of %d" % (nDone, nToDo))
            RunCommand(command, qShell, qHideStdout)
        except Queue.Empty:
            return   
            
def Worker_RunCommands_And_Move(cmd_and_filename_queue, nProcesses, nToDo, qListOfLists):
    """
    Continuously takes commands that need to be run from the cmd_and_filename_queue until the queue is empty. If required, moves 
    the output filename produced by the cmd to a specified filename. The elements of the queue can be single cmd_filename tuples
    or an ordered list of tuples that must be run in the provided order.
  
    Args:
        cmd_and_filename_queue - queue containing (cmd, actual_target_fn) tuples (if qListOfLists is False) of a list of such 
            tuples (if qListOfLists is True).
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
                subprocess.call(command, env=my_env, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if fns != None:
                    actual, target = fns
                    if os.path.exists(actual):
                        os.rename(actual, target)
        except Queue.Empty:
            return               
                            
def Worker_RunOrderedCommandList(cmd_queue, nProcesses, nToDo, qShell=True, qHideStdout=True):
    """ repeatedly takes items to process from the queue until it is empty at which point it returns. Does not take a new task
        if it can't acquire queueLock as this indicates the queue is being rearranged.
        
        Writes each commands output and stderr to a file
    """
    while True:
        try:
            i, commandSet = cmd_queue.get(True, 1)
            nDone = i - nProcesses + 1
            if nDone >= 0 and divmod(nDone, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
                PrintTime("Done %d of %d" % (nDone, nToDo))
            RunOrderedCommandList(commandSet, qShell, qHideStdout)
        except Queue.Empty:
            return   
        
def RunParallelOrderedCommandLists(nProcesses, commands, qHideStdout = False):
    """nProcesss - the number of processes to run in parallel
    commands - list of lists of commands where the commands in the inner list are completed in order (the i_th won't run until
    the i-1_th has finished).
    """
    ptm = parallel_task_manager.ParallelTaskManager_singleton()
    ptm.RunParallel(commands, True, nProcesses, qShell=True, qHideStdout = qHideStdout)              
    
def ManageQueue(runningProcesses, cmd_queue):
    """Manage a set of runningProcesses working through cmd_queue.
    If there is an error the exit all processes as quickly as possible and 
    exit via Fail() methods. Otherwise return when all work is complete
    """            
    # set all completed processes to None
    qError = False
#    dones = [False for _ in runningProcesses]
    nProcesses = len(runningProcesses)
    while True:
        if runningProcesses.count(None) == len(runningProcesses): break
        time.sleep(2)
#        for proc in runningProcesses:
        for i in xrange(nProcesses):
            proc = runningProcesses[i]
            if proc == None: continue
            if not proc.is_alive():
                if proc.exitcode != 0:
                    qError = True
                    while True:
                        try:
                            cmd_queue.get(True, 1)
                        except Queue.Empty:
                            break
                runningProcesses[i] = None
    if qError:
        Fail()              

""" 
Run a method in parallel
"""      
              
def Worker_RunMethod(Function, args_queue):
    while True:
        try:
            args = args_queue.get(True, 1)
            Function(*args)
        except Queue.Empty:
            return 

def RunMethodParallel(Function, args_queue, nProcesses):
    runningProcesses = [mp.Process(target=Worker_RunMethod, args=(Function, args_queue)) for i_ in xrange(nProcesses)]
    for proc in runningProcesses:
        proc.start()
    ManageQueue(runningProcesses, args_queue)
    
def ExampleRunMethodParallel():
    F = lambda x, y: x**2
    args_queue = mp.Queue()
    for i in xrange(100): args_queue.put((3,i))
    RunMethodParallel(F, args_queue, 16)
       
"""
Directory and file management
-------------------------------------------------------------------------------
"""               
               
def GetDirectoryName(baseDirName, i):
    if i == 0:
        return baseDirName + os.sep
    else:
        return baseDirName + ("_%d" % i) + os.sep

"""Call GetNameForNewWorkingDirectory before a call to CreateNewWorkingDirectory to find out what directory will be created"""
def CreateNewWorkingDirectory(baseDirectoryName):
    dateStr = datetime.date.today().strftime("%b%d") 
    iAppend = 0
    newDirectoryName = GetDirectoryName(baseDirectoryName + dateStr, iAppend)
    while os.path.exists(newDirectoryName):
        iAppend += 1
        newDirectoryName = GetDirectoryName(baseDirectoryName + dateStr, iAppend)
    os.mkdir(newDirectoryName)
    return newDirectoryName

def GetUnusedFilename(baseFilename, ext):
    iAppend = 0
    newFilename = baseFilename + ext
    while os.path.exists(newFilename):
        iAppend += 1
        newFilename = baseFilename + ("_%d" % iAppend) + ext
    return newFilename, iAppend
       
def SortArrayPairByFirst(useForSortAr, keepAlignedAr, qLargestFirst=False):
    sortedTuples = sorted(zip(useForSortAr, keepAlignedAr), reverse=qLargestFirst)
    useForSortAr = [i for i, j in sortedTuples]
    keepAlignedAr = [j for i, j in sortedTuples]
    return useForSortAr, keepAlignedAr      

# Get Info from seqs IDs file?
def GetSeqsInfo(inputDirectory, speciesToUse, nSpAll):
    seqStartingIndices = [0]
    nSeqs = 0
    nSeqsPerSpecies = dict()
    for iFasta in xrange(nSpAll):
        fastaFilename = inputDirectory + "Species%d.fa" % iFasta
        n = 0
        with open(fastaFilename) as infile:
            for line in infile:
                if len(line) > 1 and line[0] == ">":
                    n+=1
        nSeqsPerSpecies[iFasta] = n
        if iFasta in speciesToUse:
            nSeqs += n
            seqStartingIndices.append(nSeqs)
    seqStartingIndices = seqStartingIndices[:-1]
    nSpecies = len(speciesToUse)
    return SequencesInfo(nSeqs=nSeqs, nSpecies=nSpecies, speciesToUse=speciesToUse, seqStartingIndices=seqStartingIndices, nSeqsPerSpecies=nSeqsPerSpecies)
 
def GetSpeciesToUse(speciesIDsFN):
    """Returns species indices to use and total number of species available """
    speciesToUse = []
    speciesToUse_names = []
    nSkipped = 0
    with open(speciesIDsFN, 'rb') as speciesF:
        for line in speciesF:
            line = line.rstrip()
            if not line: continue
            if line.startswith("#"): nSkipped += 1
            else: 
                iSp, spName = line.split(": ")
                speciesToUse.append(int(iSp))
                speciesToUse_names.append(spName)
    return speciesToUse, len(speciesToUse) + nSkipped, speciesToUse_names
 
def Success():
    ptm = parallel_task_manager.ParallelTaskManager_singleton()
    ptm.Stop()  
    sys.exit()
   
def Fail():
    sys.stderr.flush()
    ptm = parallel_task_manager.ParallelTaskManager_singleton()
    ptm.Stop()
    print("ERROR: An error occurred, please review error messages for more information.")
    sys.exit()
    
"""
IDExtractor
-------------------------------------------------------------------------------
"""

def GetIDPairFromString(line):
    return map(int, line.split("_"))

class IDExtractor(object):
    """IDExtractor deals with the fact that for different datasets a user will
    want to extract a unique sequence ID from the fasta file accessions uin different 
    ways."""
    def GetIDToNameDict(self):
        raise NotImplementedError("Should not be implemented")
    def GetNameToIDDict(self):
        raise NotImplementedError("Should not be implemented")

class FullAccession(IDExtractor):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        self.nameToIDDict = dict()
        with open(idsFilename, 'rb') as idsFile:
            for line in idsFile:
                line = line.rstrip()
                if not line: continue
#                if line.startswith("#"): continue
                id, accession = line.split(": ", 1)
                id = id.replace("#", "")
                # Replace problematic characters
                accession = accession.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_") #.replace(".", "_")
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession                
                self.nameToIDDict[accession] = id 
                
    def GetIDToNameDict(self):
        return self.idToNameDict
        
    def GetNameToIDDict(self):
        return self.nameToIDDict
                
class FirstWordExtractor(IDExtractor):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        self.nameToIDDict = dict()
        with open(idsFilename, 'rb') as idsFile:
            for line in idsFile:
                id, rest = line.split(": ", 1)
                accession = rest.split(None, 1)[0]
                # Replace problematic characters
                accession = accession.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_") #.replace(".", "_")
                if accession in self.nameToIDDict:
                    raise RuntimeError("A duplicate accession was found using just first part: % s" % accession)
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession                
                self.nameToIDDict[accession] = id   
                
    def GetIDToNameDict(self):
        return self.idToNameDict
        
    def GetNameToIDDict(self):
        return self.nameToIDDict    

def HaveSupportValues(speciesTreeFN_ids):
    qHaveSupport = False
    try:
        tree.Tree(speciesTreeFN_ids, format=2)
        qHaveSupport = True
    except:
        pass
    return qHaveSupport

def RenameTreeTaxa(treeFN_or_tree, newTreeFilename, idsMap, qSupport, qFixNegatives=False, inFormat=None, label=None):
    if label != None: qSupport = False
    try:
        if type(treeFN_or_tree) == tree.TreeNode:
            t = treeFN_or_tree
        else:
            qHaveSupport = False
            if inFormat == None:
                try:
                    t = tree.Tree(treeFN_or_tree, format=2)
                    qHaveSupport = True
                except:
                    t = tree.Tree(treeFN_or_tree)
            else:
                t = tree.Tree(treeFN_or_tree, format=inFormat)
        for node in t.get_leaves():
            node.name = idsMap[node.name]
        if qFixNegatives:
            tree_length = sum([n.dist for n in t.traverse() if n != t])
            sliver = tree_length * 1e-6
        iNode = 1
        for n in t.traverse():
            if qFixNegatives and n.dist < 0.0: n.dist = sliver
            if label != None:
                if (not n.is_leaf()) and (not n.is_root()):
                    n.name = label + ("%d" % iNode)
                    iNode += 1
        if label != None: 
            with open(newTreeFilename, 'wb') as outfile:
                outfile.write(t.write(format=3)[:-1] + label + "0;")  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        else:
            if qSupport or qHaveSupport:
                t.write(outfile = newTreeFilename, format=2)  
            else:
                t.write(outfile = newTreeFilename, format=5)  
    except:
        pass

def IsWorkingDirectory(orthofinderWorkingDir):
    ok = True
    ok = ok and len(glob.glob(orthofinderWorkingDir + "clusters_OrthoFinder_*.txt_id_pairs.txt")) > 0
    ok = ok and len(glob.glob(orthofinderWorkingDir + "Species*.fa")) > 0
    return ok
    
"""
Find results of previous run    
-------------------------------------------------------------------------------
"""

def GetSpeciesDirectory():
    # Confirms all required Sequence files and BLAST etc are present
    pass

def GetOGsFile(userArg):
    """returns the WorkingDirectory, ResultsDirectory and clusters_id_pairs filename"""
    qSpecifiedResultsFile = False
    if userArg == None:
        print("ERROR: orthofinder_results_directory has not been specified")
        Fail()
    if os.path.isfile(userArg):
        fn = os.path.split(userArg)[1]
        if ("clusters_OrthoFinder_" not in fn) or ("txt_id_pairs.txt" not in fn):
            print("ERROR:\n    %s\nis neither a directory or a clusters_OrthoFinder_*.txt_id_pairs.txt file." % userArg)
            Fail()
        qSpecifiedResultsFile = True
        # user has specified specific results file
    elif userArg[-1] != os.path.sep: 
        userArg += os.path.sep
    
    # find required files
    if qSpecifiedResultsFile:
        orthofinderWorkingDir = os.path.split(userArg)[0] + os.sep
        if not IsWorkingDirectory(orthofinderWorkingDir):
            print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s" % orthofinderWorkingDir)
            Fail()
    else:
        orthofinderWorkingDir = os.path.split(userArg)[0] if qSpecifiedResultsFile else userArg
        if not IsWorkingDirectory(orthofinderWorkingDir):
            orthofinderWorkingDir = userArg + "WorkingDirectory" + os.sep   
            if not IsWorkingDirectory(orthofinderWorkingDir):
                print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s\nor\n   %s\n" % (userArg, orthofinderWorkingDir))
                Fail()
            
    if qSpecifiedResultsFile:
        print("\nUsing orthogroups in file:\n    %s" % userArg)
        return orthofinderWorkingDir, orthofinderWorkingDir, userArg
    else:     
        # identify orthogroups file
        clustersFiles = glob.glob(orthofinderWorkingDir + "clusters_OrthoFinder_*.txt_id_pairs.txt")
        orthogroupFiles = glob.glob(orthofinderWorkingDir + "OrthologousGroups*.txt") + glob.glob(orthofinderWorkingDir + "Orthogroups*.txt")
        if orthofinderWorkingDir != userArg:
            orthogroupFiles += glob.glob(userArg + "OrthologousGroups*.txt")
            orthogroupFiles += glob.glob(userArg + "Orthogroups*.txt")
        # User may have specified a WorkingDirectory and results could be in directory above
        if len(orthogroupFiles) < len(clustersFiles):
            orthogroupFiles += glob.glob(userArg + ".." + os.sep + "OrthologousGroups*.txt")
            orthogroupFiles += glob.glob(userArg + ".." + os.sep + "Orthogroups*.txt")
        clustersFiles = sorted(clustersFiles)
        orthogroupFiles = sorted(orthogroupFiles)
        if len(clustersFiles) > 1 or len(orthogroupFiles) > 1:
            print("ERROR: Results from multiple OrthoFinder runs found\n")
            print("Tab-delimiter Orthogroups*.txt/OrthologousGroups*.txt files:")
            for fn in orthogroupFiles:
                print("    " + fn)
            print("With corresponding cluster files:")
            for fn in clustersFiles:
                print("    " + fn)
            print("\nPlease run with only one set of results in directories or specifiy the specific clusters_OrthoFinder_*.txt_id_pairs.txt file on the command line")
            Fail()        
            
        if len(clustersFiles) != 1 or len(orthogroupFiles) != 1:
            print("ERROR: Results not found in <orthofinder_results_directory> or <orthofinder_results_directory>/WorkingDirectory")
            print("\nCould not find:\n    Orthogroups*.txt/OrthologousGroups*.txt\nor\n    clusters_OrthoFinder_*.txt_id_pairs.txt")
            Fail()
            
        print("\nUsing orthogroups in file:\n    %s" % orthogroupFiles[0])
        print("and corresponding clusters file:\n    %s" % clustersFiles[0])
        return orthofinderWorkingDir, userArg, clustersFiles[0]

def PrintCitation():  
    print ("\nCITATION:")  
    print (" When publishing work that uses OrthoFinder please cite:")
    print (" Emms D.M. & Kelly S. (2015), Genome Biology 16:157\n")   

    print (" If you use the species tree in your work then please also cite:")
    print (" Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278")
    print (" Emms D.M. & Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914")

def PrintUnderline(text, qHeavy=False):
    print("\n" + text)
    n = len(text)
    if text.startswith("\n"): n -= 1
    print(("=" if qHeavy else "-") * n)

def FlowText(text, n=60):
    """Split text onto lines of no more that n characters long
    """
    lines = ""
    while len(text) > 0:
        if len(lines) > 0: lines += "\n"
        if len(text) > n:
            # split at no more than 60
            iEnd = n
            while iEnd > 0 and text[iEnd] != " ": iEnd-=1
            if iEnd == 0:
                # there was nowhere to split it at a blank, just have to split at 60
                lines += text[:n]
                text = text[n:]
            else:
                # split at blank
                lines += text[:iEnd]
                text = text[iEnd+1:]  # skip blank
        else:
            lines += text
            text = ""
    return lines
    
"""
-------------------------------------------------------------------------------
"""

class nOrtho_sp(object):
    """ matrix of number of genes in species i that have orthologues/an orthologue in species j"""
    def __init__(self, nSp):
        self.n = np.zeros((nSp, nSp))
        self.n_121 = np.zeros((nSp, nSp))  # genes in i that have one orthologue in j
        self.n_12m = np.zeros((nSp, nSp))  # genes in i that have many orthologues in j
        self.n_m21 = np.zeros((nSp, nSp))  # genes in i that are in a many-to-one orthology relationship with genes in j
        self.n_m2m = np.zeros((nSp, nSp))  # genes in i that are in a many-to-many orthology relationship with genes in j
        
    def __iadd__(self, other):
        self.n += other.n
        self.n_121 += other.n_121
        self.n_12m += other.n_12m
        self.n_m21 += other.n_m21
        self.n_m2m += other.n_m2m
        return self
        
class Finalise(object):
    def __enter__(self):
        pass
    def __exit__(self, type, value, traceback):
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()