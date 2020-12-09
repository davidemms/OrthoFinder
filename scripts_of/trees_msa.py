#!/usr/bin/env python3
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
import numpy as np
from collections import Counter, defaultdict

from . import util, program_caller as pc
from . import files
from . import trim

class FastaWriter(object):
    def __init__(self, fastaFileDir_list, speciesToUse):
        """
        The ist of previous directories could incude species that are no longer used, 
        but not ones that are conflicting
        """
        self.SeqLists = dict()
        qFirst = True
        accession = ""
        sequence = ""
        required_files = set(["Species%d.fa" %i for i in speciesToUse])
        for d in fastaFileDir_list:
            for fn in glob.glob(d + "Species*.fa"):
                if os.path.basename(fn) not in required_files: continue
                with open(fn, 'r') as fastaFile:
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
        with open(outFilename, 'w') as outFile:
            for seq in self.SortSeqs([s.ToString() for s in seqs]):
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % seq)
                    outFile.write(self.SeqLists[seq])
                else:
                    print(("ERROR: %s not found" % seq))
                                
    def WriteSeqsToFasta_withNewAccessions(self, seqs, outFilename, idDict):
        with open(outFilename, 'w') as outFile:
            for seq in self.SortSeqs([s.ToString() for s in seqs]):
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % idDict[seq])
                    outFile.write(self.SeqLists[seq])
                    
    def SortSeqs(self, seqs):
        return sorted(seqs, key=lambda x: list(map(int, x.split("_"))))

def WriteTestFile(workingDir):
    d = workingDir + "/_dependencies_check/"
    if not os.path.exists(d):
        os.mkdir(d)
    testFN = d + "SimpleTest.fa"
    with open(testFN, 'w') as outfile:
        outfile.write(">a\nA\n>b\nA")
    return testFN, d
 
""" 
-----------------------------------------------------------------------------
                    Identify orthogroups for species tree            
-----------------------------------------------------------------------------
"""

def GetMulticopyCutoff(nSpecies, factor = 0.25, pMax = 0.25):
    """ Get the maximum number of species with multicopy genes for each number of non-single copy genes. A species is non-single
    copy if it is absent or if it is multicopy
    Args:
        nSpecies - number of species in the analysis
        factor - probability that a single copy gene is actually non-orthologous relative to the proportion of multi-copy species
        pMax - maximum allowed probability that one of the genes is non-orthologous
    Returns:
        maxMulticopy - array, i-th element is maximum number of allowed multicopy species if ispecies are non-single copy
    """        
    allowed_multicopy = []
    for iNonSingle in range(nSpecies):
        nSingle = nSpecies - iNonSingle
        qFoundMax = False
        for nMultiple in range(iNonSingle+1):
            pFalse = factor*nMultiple/float(nSpecies)
            pAnyFalse = 1.-(1-pFalse)**nSingle
            if pAnyFalse > pMax:
                allowed_multicopy.append(nMultiple - 1)
                qFoundMax = True
                break
        if not qFoundMax:
            allowed_multicopy.append(iNonSingle)
    return allowed_multicopy  
  
def SingleCopy_WithProbabilityTest(fraction, ogMatrix):
    """ X fraction of species have a single copy of the gene, all other species can have whatever number of copies including 0"""
    nSpecies = ogMatrix.shape[1]
    allowed_multicopy = GetMulticopyCutoff(nSpecies)
    multi = (ogMatrix > 1 * np.ones((1, nSpecies))).sum(1) # multicopy
    excluded = nSpecies - (ogMatrix == 1 * np.ones((1, nSpecies))).sum(1) # single-copy
    temp2 = (ogMatrix == np.ones((1, nSpecies))).sum(1) / float(nSpecies)   # proportion of species with gene present
    c2 = temp2 > fraction 
    passFracOGs = c2.nonzero()[0] 
    ok_ogs = [iog for iog in passFracOGs if multi[iog] <= allowed_multicopy[excluded[iog]]]
    return ok_ogs 
    
def GetOrthogroupOccupancyInfo(m):
    """
    Args:
        m - orthogroup matrix
    """
    N = m.shape[1]
    f = 1./N
    fractions = []
    nOrtho = []
    for n in range(N):
        F = 1.-n*f
        fractions.append(F)
        nOrtho.append(len(SingleCopy_WithProbabilityTest(F-1e-5, m)))
    nOrtho = list(map(float, nOrtho))
    return fractions, nOrtho
    
def DetermineOrthogroupsForSpeciesTree(m, nOGsMin=100, nSufficient=1000, increase_required=2.):
    """Orthogroups can be used if at least a fraction f of the species in the orthogroup are single copy, f is determined as described 
    below. Species that aren't single copy are allowed in the orthogroup as long as there aren't too many (criterion is quite strict). 
    The presence of species with multiple copies suggests that the single copy species might actually have arisen from duplication 
    and loss. Therefore, we use a probability model to determine how many of the excluded species can be multicopy before it is likely
    that the single copy genes are actually hidden paralogues.
    Args:
        m - the orthogroup matrix shape=(nOGs, nSpecies), each entry is number of genes from that species
        nOGsMin - the minimum number of orthogroups. Need to have at least this many before considering the next criteria
        p - proportionalIncreaseMin: Each decrease in the proportion of single copy species should bring about a relative increase 
        in the number of orthogroups that will be used of 'p', otherwise it is not worth lowering the bar
    """
    fractions, nOrtho = GetOrthogroupOccupancyInfo(m)
    if nOrtho[0] > nSufficient:
        ogsToUse = SingleCopy_WithProbabilityTest(1.0-1e-5, m)
        return ogsToUse, 1.0
    nSpecies = m.shape[1]
    for i in range(1, nSpecies):
        if nOrtho[i-1] < nOGsMin: continue
        if nOrtho[i-1] > nSufficient: break
        p = ((nOrtho[i] - nOrtho[i-1])/nOrtho[i-1]) /  (-(fractions[i] - fractions[i-1])/fractions[i-1])
        if fractions[i] > 0.5:
            if p < increase_required: break
        else:
            # if fewer than half the species are in any orthogroup then set highr bar for reducing fraction further
            if p < 2*increase_required: break
    f = fractions[i-1]
    ogsToUse = SingleCopy_WithProbabilityTest(f-1e-5, m)
    return ogsToUse, f   

class MSA(object):
    def __init__(self, msa_dict):
        self.seqs = msa_dict
        self.length = len(list(msa_dict.values())[0])

def ReadAlignment(fn):
    msa = dict()
    accession = None
    length = None
    seq = ""
    with open(fn, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith(">"):
                if accession != None:
                    if length != None and len(seq) != length:
                        text  = "ERROR: Sequence length mismatch in MSA: %s & %d" % (length, len(seq))
                        files.FileHandler.LogFailAndExit(text)
                    msa[accession] = seq
                accession = line[1:]
                seq = ""
            else:
                seq += line
        if accession != None:
            if length != None and len(seq) != length:
                text = "Error: Sequence length mismatch in MSA: %s & %d" % (length, len(seq))
                files.FileHandler.LogFailAndExit(text)
            msa[accession] = seq
    return MSA(msa)

def CreateConcatenatedAlignment(ogsToUse_ids, ogs, alignment_filename_function, output_filename, fSingleCopy, fMaxGap=0.5):
    allSpecies = {str(gene.iSp) for og in ogs for gene in og}
    concatentaedAlignments = defaultdict(str)
    for iOg in ogsToUse_ids:
        try:
            speciesCounts = Counter([gene.iSp for gene in ogs[iOg]])
            selectedSeqs = {gene.ToString() for gene in ogs[iOg] if speciesCounts[gene.iSp] == 1}
            alignment = ReadAlignment(alignment_filename_function(iOg))
            speciesInThisOg = set()
            for name, al in alignment.seqs.items():
                if name.split()[0] in selectedSeqs:
                    iSp = name.split("_")[0]
                    speciesInThisOg.add(iSp)  # this allows for the MSA method to have failed to put the sequence in the MSA
                    al = al.replace('*', '-')
                    concatentaedAlignments[iSp] += al
            # now put blanks for the missing species
            for iSp in allSpecies.difference(speciesInThisOg):
                concatentaedAlignments[iSp] += "-"*(alignment.length)
        except IndexError:
            # allow empty MSA (could fail for unknown reason)
            print("WARNING: An MSA failed for an unknown reason: %s" % alignment_filename_function(iOg))
            print("No tree or orthologs will be inferred for this orthogroup. To correct the issue, identify & correct the problematic gene sequence and rerun.")
            pass
    # Trim the completed alignment: to 50% of fraction of species present
    maxGap = (1.-fMaxGap*fSingleCopy)*len(allSpecies)
    # vectorise this as the for-loop method took too long (45mins on 3M length alignment vs 1.5s)
    # The write method doesn't require further optimisation. Took 0.15s on the same alignment
    if len(concatentaedAlignments) == 0 or len(list(concatentaedAlignments.values())[0]) == 0:
        print("All MSAs for the concatenated multiple sequence alignment were empty.")
        print("Please correct the error and re-run.")
        util.Fail()
    names = list(concatentaedAlignments.keys())
    M = np.array([list(concatentaedAlignments[name]) for name in names])
    gap_counts = sum(M == "-")
    i_keep = np.where(gap_counts <= maxGap)
    M = M[:, i_keep]
    nChar = 80
    with open(output_filename, 'w') as outfile:
        for iSeq, name in enumerate(names):
            outfile.write(">%s\n" % name)
            seq = M[iSeq,:].tolist()[0]
            for i in range(0, len(seq), nChar):
                outfile.write("".join(seq[i:i+nChar]) + "\n")
            
def trim_fn(fn):
    trim.main(fn, fn, 0.1, 500, 0.75)

""" 
-----------------------------------------------------------------------------
                             TreesForOrthogroups            
-----------------------------------------------------------------------------
"""    
     
class TreesForOrthogroups(object):
    def __init__(self, program_caller, msa_program, tree_program):
        self.program_caller = program_caller
        self.msa_program = msa_program
        self.tree_program = tree_program
    
    def GetFastaFilename(self, iOG, qResults=False):
        return files.FileHandler.GetOGsSeqFN(iOG, qResults)
            
    def GetAlignmentFilename(self, iOG, qResults=False):
        return files.FileHandler.GetOGsAlignFN(iOG, qResults)
            
    def GetTreeFilename(self, iOG, qResults=False):
        return files.FileHandler.GetOGsTreeFN(iOG, qResults)
        
    def WriteFastaFiles(self, fastaWriter, ogs, idDict, qBoth):
        # The results ones are now written by default after orthogroups, check they're not already there
        if not os.path.exists(self.GetFastaFilename(0, True)):
            for iOg, og in enumerate(ogs):
                fastaWriter.WriteSeqsToFasta_withNewAccessions(og, self.GetFastaFilename(iOg, True), idDict)
        if qBoth: 
            for iOg, og in enumerate(ogs):
                fastaWriter.WriteSeqsToFasta(og, self.GetFastaFilename(iOg))
              
    def GetAlignmentCommandsAndNewFilenames(self, ogs):
#        if self.msa_program != "mafft":
        infn_list = [self.GetFastaFilename(i) for i, og in enumerate(ogs) if len(og) >= 2]
        outfn_list = [self.GetAlignmentFilename(i) for i, og in enumerate(ogs) if len(og) >= 2]
        id_list = ["OG%07d" % i for i, og in enumerate(ogs) if len(og) >= 2]
        nSeqs = [len(og) for og in ogs if len(og) >= 2]
        return self.program_caller.GetMSACommands(self.msa_program, infn_list, outfn_list, id_list, nSeqs) 
        
    def GetTreeCommands(self, alignmentsForTree, ogs):
        outfn_list = [self.GetTreeFilename(i) for i, og in enumerate(ogs) if len(og) >= 3]
        id_list = ["OG%07d" % i for i, og in enumerate(ogs) if len(og) >= 3]
        nSeqs = [len(og) for og in ogs if len(og) >= 3]
        return self.program_caller.GetTreeCommands(self.tree_program, alignmentsForTree, outfn_list, id_list, nSeqs) 
     
    def RenameAlignmentTaxa(self, idsAlignFNS, accAlignFNs, idsDict):
        for i, (alignFN, outAlignFN) in enumerate(zip(idsAlignFNS, accAlignFNs)):
            if not os.path.exists(alignFN): continue
            with open(alignFN, 'r') as infile, open(outAlignFN, 'w') as outfile:
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(">" + idsDict[line[1:].rstrip()] + "\n")
                    else:
                        outfile.write(line)
          
    def DoTrees(self, ogs, ogMatrix, idDict, speciesIdDict, speciesToUse, nProcesses, qStopAfterSeqs, qStopAfterAlignments, qDoSpeciesTree, qTrim):
        idDict.update(speciesIdDict) # same code will then also convert concatenated alignment for species tree
        # 0       
        resultsDirsFullPath = [files.FileHandler.GetResultsSeqsDir(), files.FileHandler.GetResultsAlignDir(), files.FileHandler.GetResultsTreesDir()]
        
        # 1.
        fastaWriter = FastaWriter(files.FileHandler.GetSpeciesSeqsDir(), speciesToUse)
        self.WriteFastaFiles(fastaWriter, ogs, idDict, True)
        if qStopAfterSeqs: return resultsDirsFullPath

        # 3
        # Get OGs to use for species tree
        if qDoSpeciesTree:
            iOgsForSpeciesTree, fSingleCopy = DetermineOrthogroupsForSpeciesTree(ogMatrix)            
            concatenated_algn_fn = files.FileHandler.GetSpeciesTreeConcatAlignFN()
        else:
            iOgsForSpeciesTree = []
        alignCommands_and_filenames = self.GetAlignmentCommandsAndNewFilenames(ogs)
        if qStopAfterAlignments:
            util.PrintUnderline("Inferring multiple sequence alignments")
            pc.RunParallelCommandsAndMoveResultsFile(nProcesses, alignCommands_and_filenames, False)
            if qDoSpeciesTree:
                CreateConcatenatedAlignment(iOgsForSpeciesTree, ogs, self.GetAlignmentFilename, concatenated_algn_fn, fSingleCopy)
                # write OGs used to file
                dSpeciesTree = os.path.split(files.FileHandler.GetSpeciesTreeResultsFN(0, True))[0] + "/"
                with open(dSpeciesTree + "Orthogroups_for_concatenated_alignment.txt", 'w') as outfile:
                    for iog in iOgsForSpeciesTree: outfile.write("OG%07d\n" % iog)
            # ids -> accessions
            alignmentFilesToUse = [self.GetAlignmentFilename(i) for i, _ in enumerate(alignCommands_and_filenames)]        
            accessionAlignmentFNs = [self.GetAlignmentFilename(i, True) for i in range(len(alignmentFilesToUse))]
            if qDoSpeciesTree: 
                alignmentFilesToUse.append(concatenated_algn_fn)
                accessionAlignmentFNs.append(files.FileHandler.GetSpeciesTreeConcatAlignFN(True))
            self.RenameAlignmentTaxa(alignmentFilesToUse, accessionAlignmentFNs, idDict)
            return resultsDirsFullPath[:2]
        
        # Otherwise, alignments and trees
        # Strategy is
        # 1. Do alignments (and trees) require for species tree
        # 2. Create concatenated alignment
        # 3. Create second list of commands [speciestree] + [remaining alignments and trees]
        alignmentFilesToUse = [self.GetAlignmentFilename(i) for i, _ in enumerate(alignCommands_and_filenames)]
        treeCommands_and_filenames = self.GetTreeCommands(alignmentFilesToUse, ogs)
        commands_and_filenames = []
        if qDoSpeciesTree:
            print(("Species tree: Using %d orthogroups with minimum of %0.1f%% of species having single-copy genes in any orthogroup" % (len(iOgsForSpeciesTree), 100.*fSingleCopy)))
            util.PrintUnderline("Inferring multiple sequence alignments for species tree") 
            # Do required alignments and trees
            speciesTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            if qTrim:
                for i in iOgsForSpeciesTree:
                    commands_and_filenames.append([alignCommands_and_filenames[i], 
                                                (trim_fn, alignmentFilesToUse[i]), 
                                                treeCommands_and_filenames[i]])
            else:
                for i in iOgsForSpeciesTree:
                    commands_and_filenames.append([alignCommands_and_filenames[i],
                                                   treeCommands_and_filenames[i]])
            pc.RunParallelCommandsAndMoveResultsFile(nProcesses, commands_and_filenames, True)
            CreateConcatenatedAlignment(iOgsForSpeciesTree, ogs, self.GetAlignmentFilename, concatenated_algn_fn, fSingleCopy)
            # write OGs used to file
            dSpeciesTree = os.path.split(files.FileHandler.GetSpeciesTreeResultsFN(0, True))[0] + "/"
            with open(dSpeciesTree + "Orthogroups_for_concatenated_alignment.txt", 'w') as outfile:
                for iog in iOgsForSpeciesTree: outfile.write("OG%07d\n" % iog)
            # Add species tree to list of commands to run
            commands_and_filenames = [self.program_caller.GetTreeCommands(self.tree_program, [concatenated_algn_fn], [speciesTreeFN_ids], ["SpeciesTree"])]
            util.PrintUnderline("Inferring remaining multiple sequence alignments and gene trees") 
        else:
            util.PrintUnderline("Inferring multiple sequence alignments and gene trees") 

        # Now continue as before
        iOgsForSpeciesTree = set(iOgsForSpeciesTree)                         
        for i in range(len(treeCommands_and_filenames)):
            if i in iOgsForSpeciesTree: continue
            if qTrim:
                commands_and_filenames.append([alignCommands_and_filenames[i], 
                                            (trim_fn, alignmentFilesToUse[i]),
                                            treeCommands_and_filenames[i]])
            else:
                commands_and_filenames.append([alignCommands_and_filenames[i], 
                                            treeCommands_and_filenames[i]])
        for i in range(len(treeCommands_and_filenames), len(alignCommands_and_filenames)):
            if i in iOgsForSpeciesTree: continue
            commands_and_filenames.append([alignCommands_and_filenames[i]])
        pc.RunParallelCommandsAndMoveResultsFile(nProcesses, commands_and_filenames, True)
        
        # Convert ids to accessions for MSA
        accessionAlignmentFNs = [self.GetAlignmentFilename(i, True) for i in range(len(alignmentFilesToUse))]
        # Add concatenated Alignment
        if qDoSpeciesTree:
            alignmentFilesToUse.append(concatenated_algn_fn)
            accessionAlignmentFNs.append(files.FileHandler.GetSpeciesTreeConcatAlignFN(True))
            qHaveSupport = util.HaveSupportValues(speciesTreeFN_ids)
            if os.path.exists(speciesTreeFN_ids):
                util.RenameTreeTaxa(speciesTreeFN_ids, files.FileHandler.GetSpeciesTreeUnrootedFN(True), idDict, qSupport=qHaveSupport, qFixNegatives=True)
            else:
                text = "ERROR: Species tree inference failed"
                files.FileHandler.LogFailAndExit(text)
        self.RenameAlignmentTaxa(alignmentFilesToUse, accessionAlignmentFNs, idDict)
        return resultsDirsFullPath[:2]
