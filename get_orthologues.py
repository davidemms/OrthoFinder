# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 16:00:33 2016

@author: david
"""

import sys
import os
import ete2
import glob
import numpy as np
from collections import Counter, defaultdict
import cPickle as pic
#import sys
import itertools
import subprocess
import multiprocessing as mp
import Queue

sys.path = [p for p in sys.path if "trunk" not in p]
import trees_for_orthogroups as tfo
import idextractor
import orthofinder

# temp
#import matplotlib.pyplot as plt
#plt.interactive(True)

nThreads = 8

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
        
#    def __repr__(self):
#        return self.ToString()
#  
#    def __hash__(self):
#        return 100000 * self.iSp + self.iSeq
    
    def ToString(self):
        return "%d_%d" % (self.iSp, self.iSeq)

# ==============================================================================================================================
        
class OrthoGroupsSet(object):
    def __init__(self, orthofinderWorkingDir, clustersFilename_pairs, idExtractor = idextractor.FirstWordExtractor):
        self.workingDirOF = orthofinderWorkingDir
        self.seqIDsFN = orthofinderWorkingDir + "SequenceIDs.txt"
        self.speciesIDsFN = orthofinderWorkingDir + "SpeciesIDs.txt"
        self.speciesIDsEx = idextractor.FullAccession(self.speciesIDsFN)
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
        return self.speciesIDsEx.GetIDToNameDict()
        
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

#def CreateTaxaMapFile(directory):

def GetOGsToUse(ogSet):
    return range(50, min(1000, len(ogSet.ogs)))

def CreateTaxaMapFile(ogSet, i_ogs_to_use, outputFN):
    d = defaultdict(list)
    for iog in i_ogs_to_use:
        for seq in ogSet.ogs[iog]:
            d[seq.iSp].append(seq.ToString())
    with open(outputFN, 'wb') as outfile:
        for k, v in d.items():
            outfile.write(("%d:" % k) + ",".join(v) + "\n")

def ConcatenateTrees(i_ogs_to_use, treesPat, outputFN):
    with open(outputFN, 'wb') as outfile:
        for iog in i_ogs_to_use:
            with open(treesPat % iog, 'rb') as infile:
                for line in infile: outfile.write(line)
    
def RunAstral(ogSet, treesPat, workingDir):
    dir_astral = workingDir + "ASTRAL/"
    os.mkdir(dir_astral)
    i_ogs_to_use = GetOGsToUse(ogSet)
    tmFN = dir_astral + "Taxa_map.txt"
    CreateTaxaMapFile(ogSet, i_ogs_to_use, tmFN)
    treesFN = dir_astral + "TreesFile.txt"
    ConcatenateTrees(i_ogs_to_use, treesPat, treesFN)
#    subprocess.call(["java", "-jar", "/home/david/Astral/astral.4.10.6.jar", "-a", tmFN, "-i", treesFN])
    print(" ".join(["java", "-Xmx6000M", "-jar", "/home/david/software/ASTRAL-multiind/Astral/astral.4.8.0.jar", "-a", tmFN, "-i", treesFN, "-o", workingDir+"SpeciesTree.txt"]))

# ==============================================================================================================================      
        
def Worker_BlastScores(cmd_queue, seqsInfo, fileInfo):
    while True:
        try:
            args = cmd_queue.get(True, 1)
            B = orthofinder.BlastFileProcessor.GetBLAST6Scores(seqsInfo, fileInfo, *args)
            with open(fileInfo.outputDir + "Bit%d_%d.pic" % args, 'wb') as outfile:
                pic.dump(B, outfile, protocol = orthofinder.picProtocol)
        except Queue.Empty:
            return 
                
class DendroBLASTTrees(object):
    def __init__(self, ogSet, outD):
        self.outD = outD
        self.ogSet = ogSet
        self.species = sorted(map(int, self.ogSet.SpeciesDict().keys()))
        # Check files exist
        
    def ReadAndPickle(self): 
        cmd_queue = mp.Queue()
        for iSp in xrange(len(self.ogSet.seqsInfo.speciesToUse)):
            for jSp in xrange(len(self.ogSet.seqsInfo.speciesToUse)):
                cmd_queue.put((iSp, jSp))           
        runningProcesses = [mp.Process(target=Worker_BlastScores, args=(cmd_queue, self.ogSet.seqsInfo, self.ogSet.fileInfo)) for i_ in xrange(nThreads)]
        for proc in runningProcesses:
            proc.start()
        for proc in runningProcesses:
            while proc.is_alive():
                proc.join() 
                
    def NumberOfSequences(self, species):
        ids = ogSet.SequenceDict()
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
            orthofinder.util.PrintTime(str(sp1))
#            Bs = [self.ogSet.BitScores(iiSp, jjSp) for jjSp in xrange(len(self.species))]
            Bs = [orthofinder.LoadMatrix("Bit", self.ogSet.fileInfo, iiSp, jjSp) for jjSp in xrange(len(self.species))]
            mins = np.ones((nSeqs[sp1], 1), dtype=np.float64)*9e99 
            maxes = np.zeros((nSeqs[sp1], 1), dtype=np.float64)
            for B, sp2 in zip(Bs, self.species):
                mins = np.minimum(mins, lil_min(B))
                maxes = np.maximum(maxes, lil_max(B))
            orthofinder.util.PrintTime("got mins and maxes")
#            print(mins)
#            print(maxes)
#            sys.exit()
            for jjSp, B  in enumerate(Bs):
                for og, m in zip(ogsPerSpecies, ogMatrices):
                    for gi, i in og[iiSp]:
                        for gj, j in og[jjSp]:
#                            try:
                                m[i, j] = 0.5*max(B[gi.iSeq, gj.iSeq], mins[gi.iSeq]) /  maxes[gi.iSeq]
#                                if gi.ToString() == "3_201" and gj.ToString() == "0_817": print((B[gi.iSeq, gj.iSeq], mins[gi.iSeq], maxes[gi.iSeq]))
#                                if gi.ToString() == "0_817" and gj.ToString() == "3_201": 
#                                    print((B[gi.iSeq, gj.iSeq], mins[gi.iSeq], maxes[gi.iSeq]))
#                                    print(m[i, j])
#                            except RuntimeWarning:
#                                print(m.shape)
#                                print(B.shape)
#                                print((B[gi.iSeq, gj.iSeq], mins[gi.iSeq]))
#                                print((i, j, gi.iSeq, gj.iSeq))
#                                raise
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
            self.WritePhylipMatrix(m2, [g.ToString() for g in og], self.outD + "Distances_OG%07d.phy" % iog)
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
    
    def SpeciesTree(self, D, spPairs, speciesMatrixFN):
        n = len(self.species)
        M = np.zeros((n, n))
        for (sp1, sp2), d in zip(spPairs, D):
            sp1 = self.species.index(sp1)
            sp2 = self.species.index(sp2)
            x = np.median(d)
            M[sp1, sp2] = x
            M[sp2, sp1] = x
        with open(speciesMatrixFN, 'wb') as outfile:
            outfile.write("%d\n" % n)
#            speciesDict = self.ogSet.SpeciesDict()
            for i in xrange(n):
#                outfile.write(speciesDict[str(self.species[i])] + " ")
                outfile.write(str(self.species[i]) + " ")
                values = " ".join(["%.6g" % (0. + M[i,j]) for j in range(n)])   # hack to avoid printing out "-0"
                outfile.write(values + "\n")           
                
    def InferTrees(self, treesPat):
        seqDict = self.ogSet.Spec_SeqDict()
        for iog in xrange(len(self.ogSet.ogs)):
            fn = self.outD + "Distances_OG%07d.phy" % iog
            treeFN = treesPat % iog
            subprocess.call(["fastme", "-i", fn, "-o", treeFN])
            self.RenameTreeTaxa(treeFN, os.path.splitext(treeFN)[0] + "_names.txt", seqDict)
        treeFN = self.outD + "SpeciesTree.txt"
        subprocess.call(["fastme", "-i", self.outD + "SpeciesMatrix.phy", "-s", "-o", treeFN])
        self.RenameTreeTaxa(treeFN, os.path.splitext(treeFN)[0] + "_names.txt", self.ogSet.SpeciesDict())

    def RenameTreeTaxa(self, treeFN, newTreeFilename, idsMap):     
#        with open(treeFN, "rb") as inputTree: treeString = inputTree.next()
        try:
            tree = ete2.Tree(treeFN, format=1)
            for node in tree.get_leaves():
                node.name = idsMap[node.name]
#            outputTreeString = tree.write(format=1)
#            with open(newTreeFilename, 'wb') as outFile: outFile.write(outputTreeString)     
            tree.write(format=1, outfile = newTreeFilename)  
        except:
            pass
    
    def RunAnalysis(self, treesPat):
        ogs, ogMatrices_partial = self.GetOGMatrices()
        ogMatrices = self.WriteOGMatrices(ogs, ogMatrices_partial)
        D, spPairs = self.SpeciesTreeDistances(ogs, ogMatrices)
        self.SpeciesTree(D, spPairs, self.outD + "SpeciesMatrix.phy")
        self.InferTrees(treesPat)
        return D, spPairs
        
def PrintHelp():
    pass
    
if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "--help" or sys.argv[1] == "-h":
        PrintHelp()
        sys.exit()
        
    resultsDir = sys.argv[1]
    resultsDir = "/home/david/projects/Orthology/OrthoFinder/Trees/Tests/Vertebrata/"
    orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs = tfo.GetOGsFile(resultsDir)
#    resultsDir = "/home/david/vtemp/ExampleDataset/Results_Jun23/"
#    workingDir_of = resultsDir + "WorkingDirectory/"
#    if not os.path.exists(workingDir_of): workingDir_of = resultsDir
    
    workingDir = orthofinder.util.CreateNewWorkingDirectory(orthofinderWorkingDir + "WorkingDir_DB")
    ogSet = OrthoGroupsSet(orthofinderWorkingDir, clustersFilename_pairs, idExtractor = idextractor.FirstWordExtractor)
    
    db = DendroBLASTTrees(ogSet, workingDir)
    db.ReadAndPickle()
    treesPat = workingDir + "Tree_OG%07d.txt"
    D, spPairs = db.RunAnalysis(treesPat)
    RunAstral(ogSet, treesPat, workingDir)
    
    # Next, infer orthologues - dlcpar