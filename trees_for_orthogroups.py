#!/usr/bin/env python
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

"""
Created on Thu Sep 25 13:15:22 2014

@author: david
"""
import os
import sys
import multiprocessing as mp
import time
import subprocess
import Queue
import glob

import orthofinder    

version = "0.4.0"
nProcessesDefault = 16
    

def RunCommandSet(commandSet):
    orthofinder.util.PrintTime("Runing command: %s" % commandSet[-1])
    for cmd in commandSet:
        subprocess.call(cmd, shell=True)
    orthofinder.util.PrintTime("Finshed command: %s" % commandSet[-1])
    
class FastaWriter(object):
    def __init__(self, fastaFileDir):
        self.SeqLists = dict()
        qFirst = True
        accession = ""
        sequence = ""
        for fn in glob.glob(fastaFileDir + "Species*.fa"):
            with open(fn, 'rb') as fastaFile:
                for line in fastaFile:
                    if line[0] == ">":
                        # deal with old sequence
                        if not qFirst:
                            self.SeqLists[accession] = sequence
                            sequence = ""
                        qFirst = False
                        # get id for new sequence
                        accession = line[1:].rstrip()
                    else:
                        sequence += line
                self.SeqLists[accession] = sequence
    
    def WriteSeqsToFasta(self, seqs, outFilename):
        with open(outFilename, 'wb') as outFile:
            for seq in seqs:
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % seq)
                    outFile.write(self.SeqLists[seq])
                else:
                    print("ERROR: %s not found" % seq)
                                
    def WriteSeqsToFasta_withNewAccessions(self, seqs, outFilename, idDict):
        with open(outFilename, 'wb') as outFile:
            for seq in seqs:
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % idDict[seq])
                    outFile.write(self.SeqLists[seq])
                    
                    
def Worker_RunCommand(cmd_queue):
    """ repeatedly takes items to process from the queue until it is empty at which point it returns. Does not take a new task
        if it can't acquire queueLock as this indicates the queue is being rearranged.
        
        Writes each commands output and stderr to a file
    """
    while True:
        try:
            commandSet = cmd_queue.get(True, 10)
            RunCommandSet(commandSet)
        except Queue.Empty:
            return   
    
def RunParallelCommandSets(nProcesses, commands):
    
    # Setup the workers and run
    cmd_queue = mp.Queue()
    for cmd in commands:
        cmd_queue.put(cmd)
    runningProcesses = [mp.Process(target=Worker_RunCommand, args=(cmd_queue,)) for i_ in xrange(nProcesses)]
    for proc in runningProcesses:
        proc.start()
    
    for proc in runningProcesses:
        while proc.is_alive():
            proc.join(60.)
            time.sleep(2)    

def WriteTestFile(workingDir):
    testFN = workingDir + "SimpleTest.fa"
    with open(testFN, 'wb') as outfile:
        outfile.write(">a\nA\n>b\nA")
    return testFN

def IsWorkingDirectory(orthofinderWorkingDir):
    ok = True
    ok = ok and len(glob.glob(orthofinderWorkingDir + "clusters_OrthoFinder_*.txt_id_pairs.txt")) > 0
    ok = ok and len(glob.glob(orthofinderWorkingDir + "Species*.fa")) > 0
    return ok
      
class TreesForOrthogroups(object):
    def __init__(self, baseOutputDir, orthofinderWorkingDir):
        self.baseOgFormat = "OG%07d"
        self.baseOutputDir = baseOutputDir
        self.orthofinderWorkingDir = orthofinderWorkingDir
    
    def Align_linsi(self, fasta, alignedFasta, alignmentReport, nThreads=1):
        return "mafft --localpair --maxiterate 1000 --anysymbol --thread %d %s > %s 2> %s" % (nThreads, fasta, alignedFasta, alignmentReport)   
        
    def Align_mafft(self, fasta, alignedFasta, alignmentReport, nThreads=1):
        """ For larger numbers of sequences (>500 perhaps)"""
        return "mafft --anysymbol --thread %d %s > %s 2> %s" % (nThreads, fasta, alignedFasta, alignmentReport)   
    
    def GetFastaFilename(self, iOG):
        return self.baseOutputDir + "Sequences/" + (self.baseOgFormat % iOG) + ".fa"
    def GetAlignmentFilename(self, iOG):
        return self.baseOutputDir + "Alignments/" + (self.baseOgFormat % iOG) + ".fa"
    def GetTreeFilename(self, iOG):
        return self.baseOutputDir + "Trees/" + (self.baseOgFormat % iOG) + "_tree.txt"
        
    def WriteFastaFiles(self, fastaWriter, ogs, idDict):
        for iOg, og in enumerate(ogs):
            filename = self.GetFastaFilename(iOg)
            fastaWriter.WriteSeqsToFasta_withNewAccessions(og, filename, idDict)
      
    def OGsStillToDo(self, ogs):
        retOGs = []
        nDone = 0
        for i, og in enumerate(ogs):
            treeFilename = self.GetTreeFilename(i)
            if os.path.isfile(treeFilename) and os.path.getsize(treeFilename) != 0:
                nDone +=1
                pass
            elif len(og) > 1:
                retOGs.append((i, og))
        return retOGs, nDone
              
    def GetAlignmentCommands(self, IandOGs_toDo, nSwitchToMafft, nThreads=1):
        commands = []
        for i, og in IandOGs_toDo:
            ogFastaFilename = self.GetFastaFilename(i)
            alignedFilename = self.GetAlignmentFilename(i)
            reportFilename = "/dev/null"
            if len(og) < nSwitchToMafft:
                commands.append(self.Align_linsi(ogFastaFilename, alignedFilename, reportFilename, nThreads=nThreads))
            else:
                commands.append(self.Align_mafft(ogFastaFilename, alignedFilename, reportFilename, nThreads=nThreads))
        return commands
        
    def GetTreeCommands(self, alignmenstsForTree, IandOGs_toDo):
        commands = []
        for (i, og), alignFN in zip(IandOGs_toDo, alignmenstsForTree):
            treeFilename = self.GetTreeFilename(i)
            commands.append("FastTree %s > %s 2> /dev/null" % (alignFN, treeFilename))
        return commands
               
    def DoTrees(self, ogs, idDict, nProcesses, nSwitchToMafft=500):
        
        testFN = WriteTestFile(self.orthofinderWorkingDir)
        if not orthofinder.CanRunCommand("mafft %s" % testFN, qAllowStderr=True):
            print("ERROR: Cannot run mafft")
            print("Please check MAFFT is installed and that the executables are in the system path\n")
            return False
        if not orthofinder.CanRunCommand("mafft %s" % testFN, qAllowStderr=True):
            print("ERROR: Cannot run mafft")
            print("Please check mafft is installed and that the executables are in the system path\n")
            return False
        if not orthofinder.CanRunCommand("FastTree %s" % testFN, qAllowStderr=True):
            print("ERROR: Cannot run FastTree")
            print("Please check FastTree is installed and that the executables are in the system path\n")
            return False
        os.remove(testFN)
        
        # 0
        dirs = ['Sequences', 'Alignments', 'Trees']
        for d in dirs:
            if not os.path.exists(self.baseOutputDir + d):
                os.mkdir(self.baseOutputDir + d)
        
        # 1.
        fastaWriter = FastaWriter(self.orthofinderWorkingDir)
        self.WriteFastaFiles(fastaWriter, ogs, idDict)
        print("\nFasta files for orthogroups have been written to:\n   %s" % self.baseOutputDir + "Sequences/")
        
        # 2
        IandOGs_toDo, nDone = self.OGsStillToDo(ogs)
        if nDone != 0: print("\nAlignments and trees have already been generated for %d orthogroups" % nDone)
        print("\nAlignments and trees will be generated for %d orthogroups" % len(IandOGs_toDo)) 
        
        # 3
        alignCommands = self.GetAlignmentCommands(IandOGs_toDo, nSwitchToMafft)
        alignmentFilesToUse = [self.GetAlignmentFilename(i) for i, og in IandOGs_toDo]
        treeCommands = self.GetTreeCommands(alignmentFilesToUse, IandOGs_toDo)
        commandsSet = [(alignCmd, treeCms) for alignCmd, treeCms in zip(alignCommands, treeCommands)]
            
        # 4
        if len(commandsSet) > 0:
            print("\nExample commands that will be run:")
            for cmdSet in commandsSet[:10]:
                for cmd in cmdSet:
                    print(cmd)
            print("")
            
        RunParallelCommandSets(nProcesses, commandsSet)
        
        orthofinder.PrintCitation()
        print("\nFasta files for orthogroups have been written to:\n   %s\n" % (self.baseOutputDir + "Sequences/"))
        print("Multiple sequences alignments have been written to:\n   %s\n" % (self.baseOutputDir + "Alignments/"))
        print("Gene trees have been written to:\n   %s\n" % (self.baseOutputDir + "Trees/"))
 
def PrintHelp():
    print("Usage")    
    print("-----")
    print("python trees_for_orthogroups.py orthofinder_results_directory [-t max_number_of_threads]")
    print("python trees_for_orthogroups.py -h")
    print("\n")
    
    print("Arguments")
    print("---------")
    print("""orthofinder_results_directory
    Generate multiple sequence alignments and trees for the orthogroups in orthofinder_results_directory.\n""")
    
    print("""-t max_number_of_threads, --threads max_number_of_threads
    The maximum number of processes to be run simultaneously. The deafult is %d but this 
    should be increased by the user to the maximum number of cores available.\n""" % nProcessesDefault)
        
    print("""-h, --help
   Print this help text""")
    orthofinder.PrintCitation()   

def GetIDsDict(orthofinderWorkingDir):
    # sequence IDs
    idsFilename = orthofinderWorkingDir + "SequenceIDs.txt"
    try:
        idExtract = orthofinder.FirstWordExtractor(idsFilename)
        idDict = idExtract.GetIDToNameDict()
    except RuntimeError as error:
        print(error.message)
        if error.message.startswith("ERROR"):
            print("ERROR: %s contains a duplicate ID. If %s was prepared manually then please check the IDs are correct. " % (idsFilename, idsFilename))
            orthofinder.Fail()
        else:
            print("Tried to use only the first part of the accession in order to list the sequences in each orthologous group more concisely but these were not unique. Will use the full accession line instead.")     
            try:
                idExtract = orthofinder.FullAccession(idsFilename)
                idDict = idExtract.GetIDToNameDict()
            except:
                print("ERROR: %s contains a duplicate ID. If %s was prepared manually then please check the IDs are correct. " % (idsFilename, idsFilename))
                orthofinder.Fail()
    
    # species names
    speciesDict = dict()
    with open(orthofinderWorkingDir + "SpeciesIDs.txt", 'rb') as idsFile:
        for line in idsFile:
            iSpecies, filename = line.rstrip().split(": ", 1)
            speciesName = os.path.splitext(os.path.split(filename)[1])[0]
            speciesDict[iSpecies] = speciesName   
    idDict = {seqID:speciesDict[seqID.split("_")[0]] + "_" + name for seqID, name in idDict.items()}
    return idDict    

def GetOGsFile(userArg):
    """returns the WorkingDirectory, ResultsDirectory and clusters_id_pairs filename"""
    qSpecifiedResultsFile = False
    if userArg == None:
        print("ERROR: orthofinder_results_directory has not been specified")
        orthofinder.Fail()
    if os.path.isfile(userArg):
        fn = os.path.split(userArg)[1]
        if ("clusters_OrthoFinder_" not in fn) or ("txt_id_pairs.txt" not in fn):
            print("ERROR:\n    %s\nis neither a directory or a clusters_OrthoFinder_*.txt_id_pairs.txt file." % userArg)
            orthofinder.Fail()
        qSpecifiedResultsFile = True
        # user has specified specific results file
    elif userArg[-1] != os.path.sep: 
        userArg += os.path.sep
    
    # find required files
    if qSpecifiedResultsFile:
        orthofinderWorkingDir = os.path.split(userArg)[0] + os.sep
        if not IsWorkingDirectory(orthofinderWorkingDir):
            print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s" % orthofinderWorkingDir)
            orthofinder.Fail()
    else:
        orthofinderWorkingDir = os.path.split(userArg)[0] if qSpecifiedResultsFile else userArg
        if not IsWorkingDirectory(orthofinderWorkingDir):
            orthofinderWorkingDir = userArg + "WorkingDirectory" + os.sep   
            if not IsWorkingDirectory(orthofinderWorkingDir):
                print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s\nor\n   %s\n" % (userArg, orthofinderWorkingDir))
                orthofinder.Fail()
            
    if qSpecifiedResultsFile:
        print("Generating trees for orthogroups in file:\n    %s" % userArg)
        return orthofinderWorkingDir, orthofinderWorkingDir, userArg
    else:     
        # identify orthogroups file
        clustersFiles = glob.glob(orthofinderWorkingDir + "clusters_OrthoFinder_*.txt_id_pairs.txt")
        orthogroupFiles = glob.glob(orthofinderWorkingDir + "OrthologousGroups*.txt") 
        if orthofinderWorkingDir != userArg:
            orthogroupFiles += glob.glob(userArg + "OrthologousGroups*.txt")
        # User may have specified a WorkingDirectory and results could be in directory above
        if len(orthogroupFiles) < len(clustersFiles):
            orthogroupFiles += glob.glob(userArg + ".." + os.sep + "OrthologousGroups*.txt")
        clustersFiles = sorted(clustersFiles)
        orthogroupFiles = sorted(orthogroupFiles)
        if len(clustersFiles) > 1 or len(orthogroupFiles) > 1:
            print("ERROR: Results from multiple OrthoFinder runs found\n")
            print("Tab-delimiter OrthologousGroups*.txt files:")
            for fn in orthogroupFiles:
                print("    " + fn)
            print("With corresponding cluster files:")
            for fn in clustersFiles:
                print("    " + fn)
            print("\nPlease run with only one set of results in directories or specifiy the specific clusters_OrthoFinder_*.txt_id_pairs.txt file on the command line")
            orthofinder.Fail()        
            
        if len(clustersFiles) != 1 or len(orthogroupFiles) != 1:
            print("ERROR: Results not found in <orthofinder_results_directory> or <orthofinder_results_directory>/WorkingDirectory")
            print("\nCould not find:\n    OrthologousGroups*.txt\nor\n    clusters_OrthoFinder_*.txt_id_pairs.txt")
            orthofinder.Fail()
            
        print("Generating trees for orthogroups in file:\n    %s" % orthogroupFiles[0])
        print("and corresponding clusters file:\n    %s" % clustersFiles[0])
        return orthofinderWorkingDir, userArg, clustersFiles[0]

if __name__ == "__main__":
    print("\nOrthoFinder Alignments and Trees version %s Copyright (C) 2015 David Emms\n" % version)
    print("""    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software, and you are welcome to redistribute it under certain conditions.
    For details please see the License.md that came with this software.\n""")
    if len(sys.argv) == 1 or sys.argv[1] == "--help" or sys.argv[1] == "help" or sys.argv[1] == "-h":
        PrintHelp()
        sys.exit()
        
    v = map(int, orthofinder.version.split("."))
    v = 100 * v[0] + 10*v[1] + v[2] 
    if v < 28: 
        print("ERROR: OrthoFinder program has not been updated, please update 'orthofinder.py' to the version %s\n" % version)
        orthofinder.Fail()

    # Get arguments    
    userDir = None
    nProcesses = None
    
    args = sys.argv[1:]    
    while len(args) != 0:
        arg = args.pop(0)
        if arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument -t")
                orthofinder.Fail()
            arg = args.pop(0)
            try:
                nProcesses = int(arg)
            except:
                print("Incorrect argument for number of threads: %s" % arg)
                orthofinder.Fail()   
        else:
            userDir = arg
    
    # Check arguments
    orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs = GetOGsFile(userDir)

    if nProcesses == None:
        print("""Number of parallel processes has not been specified, will use the default value.  
   Number of parallel processes can be specified using the -t option\n""")
        nProcesses = nProcessesDefault
    print("Using %d threads for alignments and trees\n" % nProcesses)
    
    ogs = orthofinder.MCL.GetPredictedOGs(clustersFilename_pairs)     
    idDict = GetIDsDict(orthofinderWorkingDir)
    
    treeGen = TreesForOrthogroups(orthofinderResultsDir, orthofinderWorkingDir)
    treeGen.DoTrees(ogs, idDict, nProcesses, nSwitchToMafft=500)

