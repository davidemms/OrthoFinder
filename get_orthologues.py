#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 16:00:33 2016

@author: david
"""

import sys
import os
import subprocess
import numpy as np
from collections import Counter, defaultdict
import cPickle as pic
import itertools
import multiprocessing as mp
import Queue

import tree
import trees_for_orthogroups as tfo
import orthofinder
import root_from_duplications as rfd
import process_trees as pt

#print(orthofinder.__file__)

nThreads = orthofinder.nThreadsDefault

class Seq(object):
    def __init__(self, seqInput):
        """ Constructor takes sequence in any format and returns generators the 
        Seq object accordingly. If performance is really important then can write 
        individual an @classmethod to do that without the checks"""
        if type(seqInput) is str:
            self.iSp, self.iSeq = map(int, seqInput.split("_"))
        elif len(seqInput) == 2:
            if seqInput[0] is str:
                self.iSp, self.iSeq = map(int, seqInput)
            else:
                self.iSp= seqInput[0]
                self.iSeq = seqInput[1]
        else:
            raise NotImplemented
    
    def __eq__(self, other):
        return (isinstance(other, self.__class__)
            and self.__dict__ == other.__dict__)

    def __ne__(self, other):
        return not self.__eq__(other)         
        
    def __repr__(self):
        return self.ToString()
#  
#    def __hash__(self):
#        return 100000 * self.iSp + self.iSeq
    
    def ToString(self):
        return "%d_%d" % (self.iSp, self.iSeq)

# ==============================================================================================================================
        
class OrthoGroupsSet(object):
    def __init__(self, orthofinderWorkingDir, clustersFilename_pairs, idExtractor = orthofinder.FirstWordExtractor):
        self.workingDirOF = orthofinderWorkingDir
        self.seqIDsFN = orthofinderWorkingDir + "SequenceIDs.txt"
        self.speciesIDsFN = orthofinderWorkingDir + "SpeciesIDs.txt"
        self.speciesIDsEx = orthofinder.FullAccession(self.speciesIDsFN)
        self._Spec_SeqIDs = None
        self._extractor = idExtractor
        self.clustersFN = clustersFilename_pairs
        self.seqIDsEx = None
        self.ogs = self.OGs()
        self.speciesToUse = orthofinder.GetSpeciesToUse(self.speciesIDsFN)
        self.seqsInfo = orthofinder.GetSeqsInfo(orthofinderWorkingDir, self.speciesToUse)
        self.fileInfo = orthofinder.FileInfo(inputDir=orthofinderWorkingDir, outputDir = orthofinderWorkingDir, graphFilename="")
#        nSeqs, nSpecies, speciesStartingIndices = orthofinder.BlastFileProcessor.GetNumberOfSequencesInFileFromDir(workingDir, self.speciesToUse)
#        self.bfp = orthofinder.BlastFileProcessor(workingDir, self.speciesToUse, nSeqs, nSpecies, speciesStartingIndices)

    def SequenceDict(self):
        if self.seqIDsEx == None:
            self.seqIDsEx = self._extractor(self.seqIDsFN)
        return self.seqIDsEx.GetIDToNameDict()
        
    def SpeciesDict(self):
        d = self.speciesIDsEx.GetIDToNameDict()
        return {k:v.rsplit(".",1)[0] for k,v in d.items()}
        
    def Spec_SeqDict(self, qSplitBar=True):
        if self._Spec_SeqIDs != None:
            return self._Spec_SeqIDs
#        try:
        seqs = self.SequenceDict()
        specs = self.SpeciesDict()
        self._Spec_SeqIDs = dict()
        for seq in seqs:
            iSpecies = seq.split("_")[0]
            if iSpecies not in specs: continue
            if qSplitBar:
                self._Spec_SeqIDs[seq] = specs[iSpecies].replace(". ", "_").replace(" ", "_") + "_" + seqs[seq].split("|")[0]
            else:
                self._Spec_SeqIDs[seq] = specs[iSpecies].replace(". ", "_").replace(" ", "_") + "_" + seqs[seq]
        return self._Spec_SeqIDs
    
    def OGs(self):
        ogs = orthofinder.MCL.GetPredictedOGs(self.clustersFN)     
        ogs = [[Seq(g) for g in og] for og in ogs if len(og) >= 4]   
        return ogs
        
#    def BitScores(self, iSp, jSp):
##        return self.bfp.GetBLAST6Scores(self.speciesToUse.index(iSp), self.speciesToUse.index(jSp), False)
#        return orthofinder.BlastFileProcessor.GetBLAST6Scores(self.seqsInfo, self.fileInfo, iSp, jSp, False)
##        with open(self.workingDir + "B%d_%d.pic" % (iSp, jSp), 'rb') as infile:
##            return pic.load(infile)

# ==============================================================================================================================

def lil_min(M):
    n = M.shape[0]
    mins = np.ones((n, 1), dtype = np.float64) * 9e99
    for kRow in xrange(n):
        values=M.getrowview(kRow)
        if values.nnz == 0:
            continue
        mins[kRow] = min(values.data[0])
    return mins 

def lil_max(M):
    n = M.shape[0]
    maxes = np.zeros((n, 1), dtype = np.float64)
    for kRow in xrange(n):
        values=M.getrowview(kRow)
        if values.nnz == 0:
            continue
        maxes[kRow] = max(values.data[0])
    return maxes
     
# ==============================================================================================================================      
# ASTRAL

def GetOGsToUse(ogSet):
    return range(100, min(10000, len(ogSet.ogs)))

def CreateTaxaMapFile(ogSet, i_ogs_to_use, outputFN):
    """Get max number of sequences per species"""
    sp_max = defaultdict(int)
    for iog in i_ogs_to_use:
        thisCount = defaultdict(int)
        for seq in ogSet.ogs[iog]:
            thisCount[seq.iSp] += 1
        for iSp in thisCount:
            sp_max[iSp] = max(sp_max[iSp], thisCount[iSp])
    with open(outputFN, 'wb') as outfile:
        for iSp, m in sp_max.items():
            outfile.write(("%d:" % iSp) + ",".join(["%d_%d" % (iSp, j) for j in xrange(m)]) + "\n")
        
def ConvertTree(treeString):
    """for trees with sequence names iSp_jSeq replaces the jSeq with 0, 1,..."""
    t = tree.Tree(treeString)
    sp_counts = defaultdict(int)
    for seq in t:
        iSp, jSeq = seq.name.split("_")
        kSeq = sp_counts[iSp]
        sp_counts[iSp] += 1
        seq.name = "%s_%d" % (iSp, kSeq)
    return (t.write() + "\n")

def ConcatenateTrees(i_ogs_to_use, treesPat, outputFN):
    with open(outputFN, 'wb') as outfile:
        for iog in i_ogs_to_use:
            with open(treesPat % iog, 'rb') as infile:
                for line in infile: 
                  if ";" in line: outfile.write(ConvertTree(line))
    
def RunAstral(ogSet, treesPat, workingDir):
    dir_astral = workingDir + "ASTRAL/"
    os.mkdir(dir_astral)
    i_ogs_to_use = GetOGsToUse(ogSet)
    tmFN = dir_astral + "Taxa_map.txt"
    CreateTaxaMapFile(ogSet, i_ogs_to_use, tmFN)
    treesFN = dir_astral + "TreesFile.txt"
    ConcatenateTrees(i_ogs_to_use, treesPat, treesFN)
    speciesTreeFN = workingDir + "SpeciesTree_astral.txt"
    subprocess.call(" ".join(["java", "-Xmx6000M", "-jar", "~/software/ASTRAL-multiind/Astral/astral.4.8.0.jar", "-a", tmFN, "-i", treesFN, "-o", speciesTreeFN]), shell=True)
    return speciesTreeFN

# ==============================================================================================================================      
# DendroBlast   

def Worker_BlastScores(cmd_queue, seqsInfo, fileInfo, nProcesses, nToDo):
    while True:
        try:
            i, args = cmd_queue.get(True, 1)
            nDone = i - nProcesses + 1
            if nDone >= 0 and divmod(nDone, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
                orthofinder.util.PrintTime("Done %d of %d" % (nDone, nToDo))
            B = orthofinder.BlastFileProcessor.GetBLAST6Scores(seqsInfo, fileInfo, *args, qExcludeSelfHits = False)
            with open(fileInfo.outputDir + "Bit%d_%d.pic" % args, 'wb') as outfile:
                pic.dump(B, outfile, protocol = orthofinder.picProtocol)
        except Queue.Empty:
            return 
                
class DendroBLASTTrees(object):
    def __init__(self, ogSet, outD, nProcesses):
        self.outD = outD
        self.ogSet = ogSet
        self.nProcesses = nProcesses
        self.species = sorted(map(int, self.ogSet.SpeciesDict().keys()))
        treesDir = outD + "Trees/"
        self.workingDir = outD + "WorkingDirectory/"
        treesIDsDir = self.workingDir + "Trees_ids/"
        distancesDir = self.workingDir + "Distances/"
        dirs = [self.workingDir, treesDir, distancesDir, treesIDsDir]
        for d in dirs:
            if not os.path.exists(d):
                os.mkdir(d)
        self.treesPatIDs = treesIDsDir + "OG%07d_tree_id.txt"
        self.treesPat = treesDir + "OG%07d_tree.txt"
        self.distPat = distancesDir + "OG%07d.phy"
        # Check files exist
        
    def ReadAndPickle(self): 
        cmd_queue = mp.Queue()
        i = 0
        for iSp in xrange(len(self.ogSet.seqsInfo.speciesToUse)):
            for jSp in xrange(len(self.ogSet.seqsInfo.speciesToUse)):
                cmd_queue.put((i, (iSp, jSp)))           
                i+=1
        runningProcesses = [mp.Process(target=Worker_BlastScores, args=(cmd_queue, self.ogSet.seqsInfo, self.ogSet.fileInfo, nThreads, i)) for i_ in xrange(nThreads)]
        for proc in runningProcesses:
            proc.start()
        for proc in runningProcesses:
            while proc.is_alive():
                proc.join() 
                
    def NumberOfSequences(self, species):
        ids = self.ogSet.SequenceDict()
        counts = Counter([g.split("_")[0] for g in ids.keys()])
        counts = {int(k):v for k,v in counts.items() if int(k) in species}
        return counts
           
    def GetOGMatrices(self):
        """
        ogMatrices contains matrix M for each OG where:
            Mij = 0.5*max(Bij, Bmin_i)/Bmax_i
        """
        ogs = self.ogSet.OGs()
        ogsPerSpecies = [[[(g, i) for i, g in enumerate(og) if g.iSp == iSp] for iSp in self.species] for og in ogs]
        nGenes = [len(og) for og in ogs]
        nSeqs = self.NumberOfSequences(self.species)
        ogMatrices = [np.zeros((n, n)) for n in nGenes]
        for iiSp, sp1 in enumerate(self.species):
            orthofinder.util.PrintTime("Processing species %d" % sp1)
            Bs = [orthofinder.LoadMatrix("Bit", self.ogSet.fileInfo, iiSp, jjSp) for jjSp in xrange(len(self.species))]
            mins = np.ones((nSeqs[sp1], 1), dtype=np.float64)*9e99 
            maxes = np.zeros((nSeqs[sp1], 1), dtype=np.float64)
            for B, sp2 in zip(Bs, self.species):
                mins = np.minimum(mins, lil_min(B))
                maxes = np.maximum(maxes, lil_max(B))
            for jjSp, B  in enumerate(Bs):
                for og, m in zip(ogsPerSpecies, ogMatrices):
                    for gi, i in og[iiSp]:
                        for gj, j in og[jjSp]:
                                m[i, j] = 0.5*max(B[gi.iSeq, gj.iSeq], mins[gi.iSeq]) /  maxes[gi.iSeq]
        return ogs, ogMatrices
    
    def WriteOGMatrices(self, ogs, ogMatrices):
        newMatrices = []
        for iog, (og, m) in enumerate(zip(ogs, ogMatrices)):
            # dendroblast scores
            n = m.shape[0]
            m2 = np.zeros(m.shape)
            for i in xrange(n):
                for j in xrange(i):
                    m2[i, j] = -np.log(m[i,j] + m[j, i])
                    m2[j, i] = m2[i, j]
            self.WritePhylipMatrix(m2, [g.ToString() for g in og], self.distPat % iog)
            newMatrices.append(m2)
        return newMatrices
    
    def WritePhylipMatrix(self, m, names, outFN):
        with open(outFN, 'wb') as outfile:
            n = m.shape[0]
            outfile.write("%d\n" % n)
            for i in xrange(n):
                outfile.write(names[i] + " ")
                values = " ".join(["%.6g" % (0. + m[i,j]) for j in range(n)])   # hack to avoid printing out "-0"
                outfile.write(values + "\n")
    
    def SpeciesTreeDistances(self, ogs, ogMatrices, method = 0):
        spPairs = list(itertools.combinations(self.species, 2))
        D = [[] for _ in spPairs]
        if method == 0:
            """ closest distance for each species pair in each orthogroup"""
            for og, m in zip(ogs, ogMatrices):
                spDict = defaultdict(list)
                for i, g in enumerate(og):
                    spDict[g.iSp].append(i)
                for (sp1, sp2), d_list in zip(spPairs, D):
                    distances = [m[i,j] for i in spDict[sp1] for j in spDict[sp2]]
                    if len(distances) > 0: d_list.append(min(distances))
#                    d_list.append(min(distances) if len(distances) > 0 else None)
        return D, spPairs
    
    def PrepareSpeciesTreeCommand(self, D, spPairs):
        n = len(self.species)
        M = np.zeros((n, n))
        for (sp1, sp2), d in zip(spPairs, D):
            sp1 = self.species.index(sp1)
            sp2 = self.species.index(sp2)
            x = np.median(d)
            M[sp1, sp2] = x
            M[sp2, sp1] = x
        speciesMatrixFN = os.path.split(self.distPat)[0] + "/SpeciesMatrix.phy"
        with open(speciesMatrixFN, 'wb') as outfile:
            outfile.write("%d\n" % n)
#            speciesDict = self.ogSet.SpeciesDict()
            for i in xrange(n):
#                outfile.write(speciesDict[str(self.species[i])] + " ")
                outfile.write(str(self.species[i]) + " ")
                values = " ".join(["%.6g" % (0. + M[i,j]) for j in range(n)])   # hack to avoid printing out "-0"
                outfile.write(values + "\n")       
        treeFN = os.path.split(self.treesPatIDs)[0] + "/SpeciesTree_ids.txt"
        cmd = " ".join(["fastme", "-i", speciesMatrixFN, "-o", treeFN, "-w", "O"] + (["-s"] if n < 1000 else []))
        return cmd, treeFN
                
    def PrepareGeneTreeCommand(self):
        cmds = []
        for iog in xrange(len(self.ogSet.ogs)):
            nTaxa = len(self.ogSet.ogs[iog])
            cmds.append([" ".join(["fastme", "-i", self.distPat % iog, "-o", self.treesPatIDs % iog, "-w", "O"] + (["-s"] if nTaxa < 1000 else []))])
        return cmds

    def RenameTreeTaxa(self, treeFN, newTreeFilename, idsMap, qFixNegatives=False):     
#        with open(treeFN, "rb") as inputTree: treeString = inputTree.next()
        try:
            t = tree.Tree(treeFN)
            for node in t.get_leaves():
                node.name = idsMap[node.name]
            if qFixNegatives:
                for n in t.traverse():
                    if n.dist < 0.0: n.dist = 0.0
            t.write(outfile = newTreeFilename, format=4)  
        except:
            pass
    
    def RunAnalysis(self):
        ogs, ogMatrices_partial = self.GetOGMatrices()
        ogMatrices = self.WriteOGMatrices(ogs, ogMatrices_partial)
        D, spPairs = self.SpeciesTreeDistances(ogs, ogMatrices)
        cmd_spTree, spTreeFN_ids = self.PrepareSpeciesTreeCommand(D, spPairs)
        cmds_geneTrees = self.PrepareGeneTreeCommand()
        print("\n3. Inferring gene and species trees")
        print(  "-----------------------------------")
        tfo.RunParallelCommandSets(self.nProcesses, [[cmd_spTree]] + cmds_geneTrees, qHideStdout = True)
        seqDict = self.ogSet.Spec_SeqDict()
        for iog in xrange(len(self.ogSet.ogs)):
            self.RenameTreeTaxa(self.treesPatIDs % iog, self.treesPat % iog, seqDict, qFixNegatives=True)
        self.RenameTreeTaxa(spTreeFN_ids, self.workingDir + "SpeciesTree_unrooted.txt", self.ogSet.SpeciesDict(), qFixNegatives=True)       
        #spTreeFN_ids = RunAstral(self.ogSet, self.treesPatIDs, self.workingDir)   
        #self.RenameTreeTaxa(spTreeFN_ids, self.workingDir + "SpeciesTree_astral_unrooted.txt", self.ogSet.SpeciesDict(), qFixNegatives=True)     
        return len(ogs), D, spPairs, spTreeFN_ids

# ==============================================================================================================================      
# DLCPar

def GetTotalLength(t):
    return sum([node.dist for node in t])
  
def AllEqualBranchLengths(t):
    lengths = [node.dist for node in t]
    return (len(lengths) > 1 and len(set(lengths)) == 1)

def RootGeneTreesArbitrarily(treesPat, nOGs, outputDir):
    filenames = [treesPat % i for i in xrange(nOGs)]
    outFilenames = [outputDir + os.path.split(treesPat % i)[1] for i in xrange(nOGs)]
    treeFilenames = [fn for fn in filenames if fn.endswith(".txt")]
    nErrors = 0
    with open(outputDir + 'root_errors.txt', 'wb') as errorfile:
        for treeFN, outFN in zip(treeFilenames, outFilenames):
            try:    
                t = tree.Tree(treeFN)
                if len(t.get_children()) != 2:
                    R = t.get_midpoint_outgroup()
                    # if it's a tree with 3 genes all with zero length branches then root arbitrarily (it's possible this could happen with more than 3 nodes)
                    if GetTotalLength(t) == 0.0:
                      for leaf in t:
                        R = leaf
                        break
                    elif AllEqualBranchLengths(t):
                      # more generally, for any branch length all branches could have that same length
                      for leaf in t:
                        R = leaf
                        break
                    t.set_outgroup(R)
                t.resolve_polytomy()
                t.write(outfile = outFN)
            except Exception as err:
                try:
                    t = tree.Tree(treeFN)
                    for leaf in t:
                       R = leaf
                       break
                    t.set_outgroup(R)
                    t.resolve_polytomy()
                    t.write(outfile = outFN)
                except:
                    errorfile.write(treeFN + ": " + str(err) + '\n')
                    nErrors += 1    
    if nErrors != 0:
      print("WARNING: Some trees could not be rooted")
      print("Usually this is because the tree contains genes from a single species.")    

def WriteGeneSpeciesMap(d, ogSet):
    fn = d + "GeneMap.smap"
    iSpecies = ogSet.SpeciesDict().keys()
    with open(fn, 'wb') as outfile:
        for iSp in iSpecies:
            outfile.write("%s_*\t%s\n" % (iSp, iSp))
    return fn

def RunDlcpar(treesPat, ogSet, nOGs, speciesTreeFN, workingDir):
    """
    
    Implementation:
    - (skip: label species tree)
    - sort out trees (midpoint root, resolve plytomies etc)
    - run
    
    """
    rootedTreeDir = workingDir + "Trees_ids_arbitraryRoot/"
    if not os.path.exists(rootedTreeDir): os.mkdir(rootedTreeDir)
    RootGeneTreesArbitrarily(treesPat, nOGs, rootedTreeDir)
    geneMapFN = WriteGeneSpeciesMap(rootedTreeDir, ogSet)

    dlcparResultsDir = workingDir + 'dlcpar/'
    if not os.path.exists(dlcparResultsDir): os.mkdir(dlcparResultsDir)
    filenames = [rootedTreeDir + os.path.split(treesPat % i)[1] for i in xrange(nOGs)]
    
    dlcCommands = ['dlcpar_search -s %s -S %s -D 1 -C 0.125 %s -O %s' % (speciesTreeFN, geneMapFN, fn, dlcparResultsDir + os.path.splitext(os.path.split(fn)[1])[0]) for fn in filenames]
#    print(dlcCommands[0])
    # use this to run in parallel
    tfo.RunParallelCommandSets(nThreads, [[c] for c in dlcCommands], qHideStdout = True)
    return dlcparResultsDir

# ==============================================================================================================================      
# Main

def WriteTestDistancesFile(workingDir):
    testFN = workingDir + "SimpleTest.phy"
    with open(testFN, 'wb') as outfile:
        outfile.write("4\n1_1 0 0 0.2 0.25\n0_2 0 0 0.21 0.28\n3_1 0.21 0.21 0 0\n4_1 0.25 0.28 0 0")
    return testFN

def CanRunDependencies(workingDir):
    # FastME
    testFN = WriteTestDistancesFile(workingDir)
    outFN = workingDir + "SimpleTest.tre"
    if os.path.exists(outFN): os.remove(outFN)        
    if not orthofinder.CanRunCommand("fastme -i %s -o %s" % (testFN, outFN), qAllowStderr=False):
        print("ERROR: Cannot run fastme")
        print("Please check FastME is installed and that the executables are in the system path\n")
        return False
    os.remove(testFN)
    os.remove(outFN)
    # DLCPar
    if not orthofinder.CanRunCommand("dlcpar_search --version", qAllowStderr=False):
        print("ERROR: Cannot run dlcpar_search")
        print("Please check DLCpar is installed and that the executables are in the system path\n")
        return False
    return True    
        
def PrintHelp():
    print("Usage")    
    print("-----")
    print("get_orthologues.py orthofinder_results_directory [-t max_number_of_threads]")
    print("get_orthologues.py -h")
    print("\n")
    
    print("Arguments")
    print("---------")
    print("""orthofinder_results_directory
    Generate gene trees for the orthogroups, generated rooted species tree and infer ortholgues.\n""")
    
    print("""-t max_number_of_threads, --threads max_number_of_threads
    The maximum number of processes to be run simultaneously. The deafult is %d but this 
    should be increased by the user to the maximum number of cores available.\n""" % orthofinder.nThreadsDefault)
        
    print("""-h, --help
   Print this help text""")
    orthofinder.PrintCitation()   

def GetResultsFilesString(rootedSpeciesTreeFN):
    st = ""
    baseResultsDir = os.path.abspath(os.path.split(rootedSpeciesTreeFN[0])[0] + "./../Trees/")
    st += "Gene trees:\n   %s\n" % baseResultsDir
    if len(rootedSpeciesTreeFN) == 1:
        resultsDir = os.path.split(rootedSpeciesTreeFN[0])[0]
        st += "Rooted species tree:\n   %s\n" % rootedSpeciesTreeFN[0]
        st += "Species-by-species orthologues:\n   %s\n" % resultsDir
    else:
        st += "\nMultiple potential outgroups were identified for the species tree. Each case has been analysed separately.\n" 
        st+=  "Please review the rooted species trees and use the results corresponding to the correct one.\n\n"        
        for tFN in rootedSpeciesTreeFN:
            resultsDir = os.path.split(tFN)[0] + "/"
            st += "Rooted species tree:\n   %s\n" % tFN
            st += "Species-by-species orthologues directory:\n   %s\n\n" % resultsDir
    return st
            

def GetOrthologues(orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs, nProcesses):
    ogSet = OrthoGroupsSet(orthofinderWorkingDir, clustersFilename_pairs, idExtractor = orthofinder.FirstWordExtractor)
    if len(ogSet.speciesToUse) < 4: 
        print("ERROR: Not enough species to infer species tree")
        orthofinder.Fail()

    print("\n1. Checking required programs are installed")
    print(  "-------------------------------------------")
    if not CanRunDependencies(orthofinderWorkingDir): orthofinder.Fail()   
    
    print("\n2. Reading sequence similarity scores")
    print(  "-------------------------------------")
    resultsDir = orthofinder.util.CreateNewWorkingDirectory(orthofinderResultsDir + "Orthologues_")
    
    db = DendroBLASTTrees(ogSet, resultsDir, nProcesses)
    db.ReadAndPickle()
    nOGs, D, spPairs, spTreeFN_ids = db.RunAnalysis()
     
    print("\n4. Best outgroup(s) for species tree")
    print(  "------------------------------------")   
    spDict = ogSet.SpeciesDict()
    roots, clusters, rootedSpeciesTreeFN, nSupport = rfd.GetRoot(spTreeFN_ids, os.path.split(db.treesPatIDs)[0] + "/", rfd.GeneToSpecies_dash, nProcesses, treeFmt = 1)
    if len(roots) > 1:
        print("Observed %d duplications. %d support the best roots and %d contradict them." % (len(clusters), nSupport, len(clusters) - nSupport))
        print("Best outgroups for species tree:")  
    else:
        print("Observed %d duplications. %d support the best root and %d contradict it." % (len(clusters), nSupport, len(clusters) - nSupport))
        print("Best outgroup for species tree:")  
    for r in roots: print("  " + (", ".join([spDict[s] for s in r]))  )
    
    qMultiple = len(roots) > 1
    if qMultiple: print("\nAnalysing each of the potential species tree roots.")
    resultsSpeciesTrees = []
    for i, (r, speciesTree_fn) in enumerate(zip(roots, rootedSpeciesTreeFN)):
        if qMultiple: 
            resultsDir_new = resultsDir + "Orthologues_for_potential_outgroup_%d/" % i
        else:
            resultsDir_new = resultsDir + "Orthologues/"
        os.mkdir(resultsDir_new)
        resultsSpeciesTrees.append(resultsDir_new + "SpeciesTree_rooted.txt")
        db.RenameTreeTaxa(speciesTree_fn, resultsSpeciesTrees[-1], db.ogSet.SpeciesDict(), qFixNegatives=True)

        print("\n5%s. Reconciling gene and species trees" % ("-%d"%i if qMultiple else "")) 
        print(  "-------------------------------------" + ("--" if qMultiple else ""))   
        print("Root: " + (", ".join([spDict[s] for s in r])))
        dlcparResultsDir = RunDlcpar(db.treesPatIDs, ogSet, nOGs, speciesTree_fn, db.workingDir)

        # Orthologue lists
        print("\n6%s. Inferring orthologues from gene trees" % ("-%d"%i if qMultiple else ""))
        print(  "----------------------------------------" + ("--" if qMultiple else ""))      
        pt.get_orthologue_lists(ogSet, resultsDir_new, dlcparResultsDir, db.workingDir)     
     
    print("\n7. Writing results files")
    print(  "------------------------")   
    
    return GetResultsFilesString(resultsSpeciesTrees)
     
if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "--help" or sys.argv[1] == "-h":
        PrintHelp()
        sys.exit()
        
    # get arguments
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
    print("0. Getting Orthologues")
    print("----------------------")
    if nProcesses == None:
        print("""\nNumber of parallel processes has not been specified, will use the default value.  
Number of parallel processes can be specified using the -t option.""")
        nProcesses = orthofinder.nThreadsDefault
    print("Using %d threads for alignments and trees" % nProcesses)
    
    orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs = tfo.GetOGsFile(userDir)
    resultsString = GetOrthologues(orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs, nProcesses)
    print(resultsString)
    orthofinder.PrintCitation()
    
    
    
