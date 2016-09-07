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
# david_emms@hotmail.comhor: david

import os
import sys
import glob
import shutil
import subprocess
import numpy as np
from collections import Counter, defaultdict
import cPickle as pic
import itertools
import multiprocessing as mp
import Queue
import warnings

import util
import tree
import matrices
import mcl as MCL
import root_from_duplications as rfd
import orthologues_from_recon_trees as pt
import blast_file_processor as BlastFileProcessor

nThreads = util.nThreadsDefault

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
    def __init__(self, orthofinderWorkingDir, clustersFilename_pairs, idExtractor = util.FirstWordExtractor):
        self.workingDirOF = orthofinderWorkingDir
        self.seqIDsFN = orthofinderWorkingDir + "SequenceIDs.txt"
        self.speciesIDsFN = orthofinderWorkingDir + "SpeciesIDs.txt"
        self.speciesIDsEx = util.FullAccession(self.speciesIDsFN)
        self._Spec_SeqIDs = None
        self._extractor = idExtractor
        self.clustersFN = clustersFilename_pairs
        self.seqIDsEx = None
        self.ogs = None
        self.speciesToUse = util.GetSpeciesToUse(self.speciesIDsFN)
        self.seqsInfo = util.GetSeqsInfo(orthofinderWorkingDir, self.speciesToUse)
        self.fileInfo = util.FileInfo(inputDir=orthofinderWorkingDir, outputDir = orthofinderWorkingDir, graphFilename="")
        self.id_to_og = None

    def SequenceDict(self):
        if self.seqIDsEx == None:
            try:
                self.seqIDsEx = self._extractor(self.seqIDsFN)
            except RuntimeError as error:
                print(error.message)
                if error.message.startswith("ERROR"): 
                    util.Fail()
                else:
                    print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup\nmore concisely but these were not unique. The full accession line will be used instead.\n")     
                    self.seqIDsEx = util.FullAccession(self.seqIDsFN)
        return self.seqIDsEx.GetIDToNameDict()
        
    def SpeciesDict(self):
        d = self.speciesIDsEx.GetIDToNameDict()
        return {k:v.rsplit(".",1)[0] for k,v in d.items()}
        
    def Spec_SeqDict(self):
        if self._Spec_SeqIDs != None:
            return self._Spec_SeqIDs
        seqs = self.SequenceDict()
        specs = self.SpeciesDict()
        specs_ed = {k:v.replace(".", "_").replace(" ", "_") for k,v in specs.items()}
        self._Spec_SeqIDs = {seqID:specs_ed[seqID.split("_")[0]] + "_" + name for seqID, name in seqs.items()}
        return self._Spec_SeqIDs
    
    def OGs(self):
        if self.ogs != None:
            return self.ogs
        self.ogs = MCL.GetPredictedOGs(self.clustersFN)     
        self.ogs = [[Seq(g) for g in og] for og in self.ogs if len(og) >= 4]   
        return self.ogs
        
    def ID_to_OG_Dict(self):
        if self.id_to_og != None:
            return self.id_to_og
        self.id_to_og = {g.ToString():iog for iog, og in enumerate(self.OGs()) for g in og}
        return self.id_to_og
        

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
    return range(100, min(10000, len(ogSet.OGs())))

def CreateTaxaMapFile(ogSet, i_ogs_to_use, outputFN):
    """Get max number of sequences per species"""
    sp_max = defaultdict(int)
    ogs = ogSet.OGs()
    for iog in i_ogs_to_use:
        thisCount = defaultdict(int)
        for seq in ogs[iog]:
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
                util.PrintTime("Done %d of %d" % (nDone, nToDo))
            B = BlastFileProcessor.GetBLAST6Scores(seqsInfo, fileInfo, *args, qExcludeSelfHits = False)
            with open(fileInfo.outputDir + "Bit%d_%d.pic" % args, 'wb') as outfile:
                pic.dump(B, outfile, protocol = util.picProtocol)
        except Queue.Empty:
            return 
                
class DendroBLASTTrees(object):
    def __init__(self, ogSet, outD, nProcesses):
        self.outD = outD
        self.ogSet = ogSet
        self.nProcesses = nProcesses
        self.species = sorted(map(int, self.ogSet.SpeciesDict().keys()))
        treesDir = outD + "Gene_Trees/"
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
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
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
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            ogs = self.ogSet.OGs()
            ogsPerSpecies = [[[(g, i) for i, g in enumerate(og) if g.iSp == iSp] for iSp in self.species] for og in ogs]
            nGenes = [len(og) for og in ogs]
            nSeqs = self.NumberOfSequences(self.species)
            ogMatrices = [np.zeros((n, n)) for n in nGenes]
            for iiSp, sp1 in enumerate(self.species):
                util.PrintTime("Processing species %d" % sp1)
                Bs = [matrices.LoadMatrix("Bit", self.ogSet.fileInfo, iiSp, jjSp) for jjSp in xrange(len(self.species))]
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
    
    def DeleteBlastMatrices(self):
        for f in glob.glob(self.ogSet.fileInfo.outputDir + "Bit*_*.pic"):
            if os.path.exists(f): os.remove(f)
        
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
        ogs = self.ogSet.OGs()
        for iog in xrange(len(ogs)):
            nTaxa = len(ogs[iog])
            cmds.append([" ".join(["fastme", "-i", self.distPat % iog, "-o", self.treesPatIDs % iog, "-w", "O"] + (["-s"] if nTaxa < 1000 else []))])
        return cmds
        
    def RunAnalysis(self):
        ogs, ogMatrices_partial = self.GetOGMatrices()
        ogMatrices = self.WriteOGMatrices(ogs, ogMatrices_partial)
        
        D, spPairs = self.SpeciesTreeDistances(ogs, ogMatrices)
        cmd_spTree, spTreeFN_ids = self.PrepareSpeciesTreeCommand(D, spPairs)
        cmds_geneTrees = self.PrepareGeneTreeCommand()
        print("\n3. Inferring gene and species trees")
        print(  "-----------------------------------")
        util.RunParallelOrderedCommandLists(self.nProcesses, [[cmd_spTree]] + cmds_geneTrees, qHideStdout = True)
        seqDict = self.ogSet.Spec_SeqDict()
        for iog in xrange(len(self.ogSet.OGs())):
            util.RenameTreeTaxa(self.treesPatIDs % iog, self.treesPat % iog, seqDict, qFixNegatives=True)
#        util.RenameTreeTaxa(spTreeFN_ids, self.workingDir + "SpeciesTree_unrooted.txt", self.ogSet.SpeciesDict(), qFixNegatives=True)        
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
    util.RunParallelOrderedCommandLists(nThreads, [[c] for c in dlcCommands], qHideStdout = True)
    return dlcparResultsDir

# ==============================================================================================================================      
# Main

def WriteTestDistancesFile(testFN):
    with open(testFN, 'wb') as outfile:
        outfile.write("4\n1_1 0 0 0.2 0.25\n0_2 0 0 0.21 0.28\n3_1 0.21 0.21 0 0\n4_1 0.25 0.28 0 0")
    return testFN

def CanRunDependencies(workingDir):
    # FastME
    testFN = workingDir + "SimpleTest.phy"
    WriteTestDistancesFile(testFN)
    outFN = workingDir + "SimpleTest.tre"
    if os.path.exists(outFN): os.remove(outFN)        
    if not util.CanRunCommand("fastme -i %s -o %s" % (testFN, outFN), qAllowStderr=False):
        print("ERROR: Cannot run fastme")
        print("Please check FastME is installed and that the executables are in the system path.\n")
        return False
    os.remove(testFN)
    os.remove(outFN)
    fastme_stat_fn = workingDir + "SimpleTest.phy_fastme_stat.txt"
    if os.path.exists(fastme_stat_fn): os.remove(fastme_stat_fn)
    # DLCPar
    if not util.CanRunCommand("dlcpar_search --version", qAllowStderr=False):
        print("ERROR: Cannot run dlcpar_search")
        print("Please check DLCpar is installed and that the executables are in the system path.\n")
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
    should be increased by the user to the maximum number of cores available.\n""" % util.nThreadsDefault)
        
    print("""-h, --help
   Print this help text""")
    util.PrintCitation()   

def GetResultsFilesString(rootedSpeciesTreeFN):
    st = ""
    baseResultsDir = os.path.abspath(os.path.split(rootedSpeciesTreeFN[0])[0] + "./../Trees/")
    st += "\nGene trees:\n   %s\n" % baseResultsDir
    if len(rootedSpeciesTreeFN) == 1:
        resultsDir = os.path.split(rootedSpeciesTreeFN[0])[0]
        st += "\nRooted species tree:\n   %s\n" % rootedSpeciesTreeFN[0]
        st += "\nSpecies-by-species orthologues:\n   %s\n" % resultsDir
    else:
        st += "\nMultiple potential outgroups were identified for the species tree. Each case has been analysed separately.\n" 
        st+=  "Please review the rooted species trees and use the results corresponding to the correct one.\n\n"        
        for tFN in rootedSpeciesTreeFN:
            resultsDir = os.path.split(tFN)[0] + "/"
            st += "Rooted species tree:\n   %s\n" % tFN
            st += "Species-by-species orthologues directory:\n   %s\n\n" % resultsDir
    return st
            

def CleanWorkingDir(dendroBlast):
    dendroBlast.DeleteBlastMatrices()
    dirs = ['Distances/', "matrices_orthologues/", "Trees_ids_arbitraryRoot/"]
    for d in dirs:
        dFull = dendroBlast.workingDir + d
        if os.path.exists(dFull): 
            shutil.rmtree(dFull)
            
def GetOrthologues(orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs, nProcesses):
    ogSet = OrthoGroupsSet(orthofinderWorkingDir, clustersFilename_pairs, idExtractor = util.FirstWordExtractor)
    if len(ogSet.speciesToUse) < 4: 
        print("ERROR: Not enough species to infer species tree")
        util.Fail()

    print("\n1. Checking required programs are installed")
    print(  "-------------------------------------------")
    if not CanRunDependencies(orthofinderWorkingDir): 
        print("Orthogroups have been inferred but the dependencies for inferring gene trees and\northologues have not been met. Please review previous messages for more information.")
        sys.exit()
        
    
    print("\n2. Calculating gene distances")
    print(  "-----------------------------")
    resultsDir = util.CreateNewWorkingDirectory(orthofinderResultsDir + "Orthologues_")
    
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
            resultsDir_new = resultsDir + "Orthologues_using_outgroup_%d/" % i
            resultsSpeciesTrees.append(resultsDir_new + "SpeciesTree_rooted_at_outgroup_%d.txt" % i)
        else:
            resultsDir_new = resultsDir + "Orthologues/"
            resultsSpeciesTrees.append(resultsDir + "SpeciesTree_rooted.txt")
        os.mkdir(resultsDir_new)
        util.RenameTreeTaxa(speciesTree_fn, resultsSpeciesTrees[-1], db.ogSet.SpeciesDict(), qFixNegatives=True)

        print("\n5%s. Reconciling gene and species trees" % ("-%d"%i if qMultiple else "")) 
        print(  "-------------------------------------" + ("--" if qMultiple else ""))   
        print("Outgroup: " + (", ".join([spDict[s] for s in r])))
        dlcparResultsDir = RunDlcpar(db.treesPatIDs, ogSet, nOGs, speciesTree_fn, db.workingDir)
        reconTreesRenamedDir = db.workingDir + "Recon_Gene_Trees/"
        os.mkdir(reconTreesRenamedDir)
        for iog in xrange(len(db.ogSet.OGs())):
            util.RenameTreeTaxa(dlcparResultsDir + "OG%07d_tree_id.locus.tree" % iog, reconTreesRenamedDir + "OG%07d_tree.txt" % iog, db.ogSet.Spec_SeqDict(), qFixNegatives=False, inFormat=8)

        # Orthologue lists
        print("\n6%s. Inferring orthologues from gene trees" % ("-%d"%i if qMultiple else ""))
        print(  "----------------------------------------" + ("--" if qMultiple else ""))      
        pt.get_orthologue_lists(ogSet, resultsDir_new, dlcparResultsDir, db.workingDir)     
     
    CleanWorkingDir(db)
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
                util.Fail()
            arg = args.pop(0)
            try:
                nProcesses = int(arg)
            except:
                print("Incorrect argument for number of threads: %s" % arg)
                util.Fail()   
        else:
            userDir = arg
        
    # Check arguments
    print("0. Getting Orthologues")
    print("----------------------")
    if nProcesses == None:
        print("""\nNumber of parallel processes has not been specified, will use the default value.  
Number of parallel processes can be specified using the -t option.""")
        nProcesses = util.nThreadsDefault
    print("Using %d threads for alignments and trees" % nProcesses)
    
    orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs = util.GetOGsFile(userDir)
    resultsString = GetOrthologues(orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs, nProcesses)
    print(resultsString)
    util.PrintCitation()
    
    
    
