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
import glob

import util 

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
            for seq in self.SortSeqs([s.ToString() for s in seqs]):
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % seq)
                    outFile.write(self.SeqLists[seq])
                else:
                    print("ERROR: %s not found" % seq)
                                
    def WriteSeqsToFasta_withNewAccessions(self, seqs, outFilename, idDict):
        with open(outFilename, 'wb') as outFile:
            for seq in self.SortSeqs([s.ToString() for s in seqs]):
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % idDict[seq])
                    outFile.write(self.SeqLists[seq])
                    
    def SortSeqs(self, seqs):
        return sorted(seqs, key=lambda x: map(int, x.split("_")))

def WriteTestFile(workingDir):
    d = workingDir + "/_dependencies_check/"
    if not os.path.exists(d):
        os.mkdir(d)
    testFN = d + "SimpleTest.fa"
    with open(testFN, 'wb') as outfile:
        outfile.write(">a\nA\n>b\nA")
    return testFN, d
      
class TreesForOrthogroups(object):
    def __init__(self, program_caller, msa_program, tree_program, resultsDir, ogsWorkingDir):
        self.program_caller = program_caller
        self.msa_program = msa_program
        self.tree_program = tree_program
        self.baseOgFormat = "OG%07d"
        self.resultsDir = resultsDir
        self.workingDir = resultsDir + "WorkingDirectory/"
        if not os.path.exists(self.workingDir): os.mkdir(self.workingDir)
        self.ogsWorkingDir = ogsWorkingDir
    
    def Align_linsi(self, fasta, alignedFasta, alignmentReport):
        return "mafft --localpair --maxiterate 1000 --anysymbol %s > %s 2> %s" % (fasta, alignedFasta, alignmentReport)   
        
    def Align_mafft(self, fasta, alignedFasta, alignmentReport):
        """ For larger numbers of sequences (>500 perhaps)"""
        return "mafft --anysymbol %s > %s 2> %s" % (fasta, alignedFasta, alignmentReport)   
    
    def GetFastaFilename(self, iOG, qResults=False):
        if qResults:
            return self.resultsDir + "Sequences/" + (self.baseOgFormat % iOG) + ".fa"
        else:
            return self.workingDir + "Sequences_ids/" + (self.baseOgFormat % iOG) + ".fa"
            
    def GetAlignmentFilename(self, iOG, qResults=False):
        if qResults:
            return self.resultsDir + "Alignments/" + (self.baseOgFormat % iOG) + ".fa"
        else:
            return self.workingDir + "Alignments_ids/" + (self.baseOgFormat % iOG) + ".fa"
            
    def GetTreeFilename(self, iOG, qResults=False):
        if qResults:
            return self.resultsDir + "Gene_Trees/" + (self.baseOgFormat % iOG) + "_tree.txt"
        else:
            return self.workingDir + "Trees_ids/" + (self.baseOgFormat % iOG) + "_tree_id.txt"
        
    def WriteFastaFiles(self, fastaWriter, ogs, idDict):
        for iOg, og in enumerate(ogs):
            fastaWriter.WriteSeqsToFasta_withNewAccessions(og, self.GetFastaFilename(iOg, True), idDict)
            fastaWriter.WriteSeqsToFasta(og, self.GetFastaFilename(iOg))
              
    def GetAlignmentCommands(self, ogs, nSwitchToMafft):
        if self.msa_program != "mafft":
            infn_list = [self.GetFastaFilename(i) for i, og in enumerate(ogs) if len(og) >= 2]
            outfn_list = [self.GetAlignmentFilename(i) for i, og in enumerate(ogs) if len(og) >= 2]
            id_list = ["OG%07d" % i for i, og in enumerate(ogs) if len(og) >= 2]
            print "GetMSACommmand params:", self, "\n", self.msa_program, "\n", infn_list[:5], "\n", outfn_list[:5], "\n", id_list[:5]
            return self.program_caller.GetMSACommands(self.msa_program, infn_list, outfn_list, id_list) 
        else:
            commands = []
            for i, og in enumerate(ogs):
                if len(og) < 2: break
                ogFastaFilename = self.GetFastaFilename(i)
                alignedFilename = self.GetAlignmentFilename(i)
                reportFilename = "/dev/null"
                if len(og) < nSwitchToMafft:
                    commands.append(self.Align_linsi(ogFastaFilename, alignedFilename, reportFilename))
                else:
                    commands.append(self.Align_mafft(ogFastaFilename, alignedFilename, reportFilename))
            return commands
        
    def GetTreeCommands(self, alignmentsForTree, ogs):
        if self.tree_program != "fasttree":
            outfn_list = [self.GetTreeFilename(i) for i, og in enumerate(ogs) if len(og) >= 2]
            id_list = ["OG%07d" % i for i, og in enumerate(ogs) if len(og) >= 2]
            return self.program_caller.GetTreeCommands(self.tree_program, alignmentsForTree, outfn_list, id_list) 
        else:
            commands = []
            for i, (alignFN, og) in enumerate(zip(alignmentsForTree, ogs)):
                if len(og) < 4: break
                treeFilename = self.GetTreeFilename(i)
                commands.append("FastTree %s > %s 2> /dev/null" % (alignFN, treeFilename))
            return commands
               
    def DoTrees(self, ogs, idDict, nProcesses, qStopAfterSeqs, qStopAfterAlignments, nSwitchToMafft=500):
        # 0       
        resultsDirsFullPath = []
        for fn in [self.GetFastaFilename, self.GetAlignmentFilename, self.GetTreeFilename]:
            for qIDs in [True, False]:
                d = os.path.split(fn(0, not qIDs))[0]
                if not os.path.exists(d): os.mkdir(d)
                if not qIDs: resultsDirsFullPath.append(d)
            if qStopAfterSeqs: break
            if qStopAfterAlignments and fn == self.GetAlignmentFilename: break
        
        # 1.
        fastaWriter = FastaWriter(self.ogsWorkingDir)
        self.WriteFastaFiles(fastaWriter, ogs, idDict)
        if qStopAfterSeqs: return resultsDirsFullPath
        
        # 2
        if qStopAfterAlignments:
            util.PrintUnderline("Inferring multiple sequence alignments") 
        else:
            util.PrintUnderline("Inferring multiple sequence alignments and gene trees") 
        
        # 3
        alignCommands = self.GetAlignmentCommands(ogs, nSwitchToMafft)
        if qStopAfterAlignments:
            util.RunParallelCommands(nProcesses, alignCommands, qShell=True)
            return resultsDirsFullPath[:2]
        alignmentFilesToUse = [self.GetAlignmentFilename(i) for i, _ in enumerate(alignCommands)]
        treeCommands = self.GetTreeCommands(alignmentFilesToUse, ogs)
        commandsSet = []
        for i in xrange(len(treeCommands)):
            commandsSet.append([alignCommands[i], treeCommands[i]])
        for i in xrange(len(treeCommands), len(alignCommands)):
            commandsSet.append([alignCommands[i]])
        util.RunParallelOrderedCommandLists(nProcesses, commandsSet)
        
        # Convert ids to accessions
        for i, alignFN in enumerate(alignmentFilesToUse):
            with open(alignFN, 'rb') as infile, open(self.GetAlignmentFilename(i, True), 'wb') as outfile:
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(">" + idDict[line[1:].rstrip()] + "\n")
                    else:
                        outfile.write(line)
            if os.path.exists(self.GetTreeFilename(i)):
                util.RenameTreeTaxa(self.GetTreeFilename(i), self.GetTreeFilename(i, True), idDict, qFixNegatives=True)
        
        return resultsDirsFullPath[:2]
