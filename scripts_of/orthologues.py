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
# david_emms@hotmail.comhor: david
import os
import sys
import csv
import time
import shutil
import numpy as np
import subprocess
from collections import Counter, defaultdict
import itertools
import multiprocessing as mp
import warnings
try: 
    import queue
except ImportError:
    import Queue as queue  

PY2 = sys.version_info <= (3,)       
csv_write_mode = 'wb' if PY2 else 'wt'

from . import util
from . import tree
from . import mcl as MCL
from . import stride
from . import trees2ologs_dlcpar
from . import trees2ologs_of
from . import blast_file_processor as BlastFileProcessor
from . import trees_msa
from . import wrapper_phyldog
from . import stag
from . import files
from . import parallel_task_manager

nThreads = util.nThreadsDefault

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
    
class Seq(object):
    def __init__(self, seqInput):
        """ Constructor takes sequence in any format and returns generators the 
        Seq object accordingly. If performance is really important then can write 
        individual an @classmethod to do that without the checks"""
        if type(seqInput) is str:
            a,b = seqInput.split("_")
            self.iSp = int(a)
            self.iSeq = int(b)
        elif len(seqInput) == 2:
            if seqInput[0] is str:
                self.iSp, self.iSeq = list(map(int, seqInput))
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
    
    def ToString(self):
        return "%d_%d" % (self.iSp, self.iSeq)

# ==============================================================================================================================
        
class OrthoGroupsSet(object):
    def __init__(self, orthofinderWorkingDir_list, speciesToUse, nSpAll, qAddSpeciesToIDs, idExtractor = util.FirstWordExtractor):
        self.speciesIDsEx = util.FullAccession(files.FileHandler.GetSpeciesIDsFN())
        self._Spec_SeqIDs = None
        self._extractor = idExtractor
        self.seqIDsEx = None
        self.ogs_all = None
        self.iOgs4 = 0
        self.speciesToUse = speciesToUse     # list of ints
        self.seqsInfo = util.GetSeqsInfo(orthofinderWorkingDir_list, self.speciesToUse, nSpAll)
        self.id_to_og = None
        self.qAddSpeciesToIDs = qAddSpeciesToIDs

    def SequenceDict(self):
        if self.seqIDsEx == None:
            try:
                self.seqIDsEx = self._extractor(files.FileHandler.GetSequenceIDsFN())
            except RuntimeError as error:
                print(str(error))
                if str(error).startswith("ERROR"): 
                    files.FileHandler.LogFailAndExit()
                else:
                    print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup\nmore concisely but these were not unique. The full accession line will be used instead.\n")     
                    self.seqIDsEx = util.FullAccession(files.FileHandler.GetSequenceIDsFN())
        return self.seqIDsEx.GetIDToNameDict()
        
    def SpeciesDict(self):
        d = self.speciesIDsEx.GetIDToNameDict()
        return {k:v.rsplit(".",1)[0] for k,v in d.items()}
        
    def Spec_SeqDict(self):
        if self._Spec_SeqIDs != None:
            return self._Spec_SeqIDs
        seqs = self.SequenceDict()
        seqs = {k:v for k,v in seqs.items() if int(k.split("_")[0]) in self.speciesToUse}
        if not self.qAddSpeciesToIDs:
            self._Spec_SeqIDs = seqs
            return seqs
        specs = self.SpeciesDict()
        specs_ed = {k:v.replace(".", "_").replace(" ", "_") for k,v in specs.items()}
        self._Spec_SeqIDs = {seqID:specs_ed[seqID.split("_")[0]] + "_" + name for seqID, name in seqs.items()}
        return self._Spec_SeqIDs
    
    def OGs(self, qInclAll=False):
        if self.ogs_all != None:
            if qInclAll:
                return self.ogs_all
            else:
                return self.ogs_all[:self.iOgs4]
        ogs = MCL.GetPredictedOGs(files.FileHandler.GetClustersFN())     
        self.ogs_all = [[Seq(g) for g in og] for og in ogs]   
        self.iOgs4 = len(self.ogs_all) if len(self.ogs_all[-1]) >= 4 else next(i for i, og in enumerate(self.ogs_all) if len(og) < 4) 
        if qInclAll:
            return self.ogs_all
        else:
            return self.ogs_all[:self.iOgs4]
        
    def OrthogroupMatrix(self):
        """ qReduce give a matrix with only as many columns as species for cases when 
        clustering has been performed on a subset of species"""
        ogs = self.OGs()
        iSpecies = sorted(set([gene.iSp for og in ogs for gene in og]))
        speciesIndexDict = {iSp:iCol for iCol, iSp in enumerate(iSpecies)}
        nSpecies = len(iSpecies)
        nGroups = len(ogs)
        # (i, j)-th entry of ogMatrix gives the number of genes from i in orthologous group j
        ogMatrix = np.zeros((nGroups, nSpecies)) 
        for i_og, og in enumerate(ogs):
            for gene in og:
                ogMatrix[i_og, speciesIndexDict[gene.iSp]] += 1
        return ogMatrix
        
    def ID_to_OG_Dict(self):
        if self.id_to_og != None:
            return self.id_to_og
        self.id_to_og = {g.ToString():iog for iog, og in enumerate(self.OGs()) for g in og}
        return self.id_to_og
        
# ==============================================================================================================================

def lil_min(M):
    n = M.shape[0]
    mins = np.ones((n, 1), dtype = np.float64) * 9e99
    for kRow in range(n):
        values=M.getrowview(kRow)
        if values.nnz == 0:
            continue
        mins[kRow] = min(values.data[0])
    return mins 

def lil_max(M):
    n = M.shape[0]
    maxes = np.zeros((n, 1), dtype = np.float64)
    for kRow in range(n):
        values=M.getrowview(kRow)
        if values.nnz == 0:
            continue
        maxes[kRow] = max(values.data[0])
    return maxes

def lil_minmax(M):
    n = M.shape[0]
    mins = np.ones((n, 1), dtype = np.float64) * 9e99
    maxes = np.zeros((n, 1), dtype = np.float64)
    for kRow in range(n):
        values=M.getrowview(kRow)
        if values.nnz == 0:
            continue
        mins[kRow] = min(values.data[0])
        maxes[kRow] = max(values.data[0])
    return mins, maxes
 
   
# ==============================================================================================================================    
# Species trees for two- & three-species analyses

def WriteSpeciesTreeIDs_TwoThree(taxa, outFN):
    """
    Get the unrooted species tree for two or three species
    Args:
        taxa - list of species names
    Returns:
    
    """
    t = tree.Tree()
    for s in taxa:
        t.add_child(tree.TreeNode(name=s))
    t.write(outfile=outFN)
    
def GetSpeciesTreeRoot_TwoTaxa(taxa):
    speciesTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
    t = tree.Tree("(%s,%s);" % (taxa[0], taxa[1]))  
    t.write(outfile=speciesTreeFN_ids)
    return speciesTreeFN_ids
    
# ==============================================================================================================================      
# DendroBlast   

def Worker_OGMatrices_ReadBLASTAndUpdateDistances(cmd_queue, worker_status_queue, iWorker, ogMatrices, nGenes, seqsInfo, blastDir_list, ogsPerSpecies, qDoubleBlast):
    speciesToUse = seqsInfo.speciesToUse
    with np.errstate(divide='ignore'):
        while True:
            try:
                iiSp, sp1, nSeqs_sp1 = cmd_queue.get(True, 1)
                worker_status_queue.put(("start", iWorker, iiSp))
                Bs = [BlastFileProcessor.GetBLAST6Scores(seqsInfo, blastDir_list, sp1, sp2, qExcludeSelfHits = False, qDoubleBlast=qDoubleBlast) for sp2 in speciesToUse]    
                mins = np.ones((nSeqs_sp1, 1), dtype=np.float64)*9e99 
                maxes = np.zeros((nSeqs_sp1, 1), dtype=np.float64)
                for B in Bs:
                    m0, m1 = lil_minmax(B)
                    mins = np.minimum(mins, m0)
                    maxes = np.maximum(maxes, m1)
                maxes_inv = 1./maxes
                for jjSp, B  in enumerate(Bs):
                    for og, m in zip(ogsPerSpecies, ogMatrices):
                        for gi, i in og[iiSp]:
                            for gj, j in og[jjSp]:
                                    m[i][j] = 0.5*max(B[gi.iSeq, gj.iSeq], mins[gi.iSeq]) * maxes_inv[gi.iSeq]
                del Bs, B, mins, maxes, m0, m1, maxes_inv, m    # significantly reduces RAM usage
                worker_status_queue.put(("finish", iWorker, iiSp))
            except queue.Empty:
                worker_status_queue.put(("empty", iWorker, None))
                return 

def GetRAMErrorText():
    text = "ERROR: The computer ran out of RAM and killed OrthoFinder processes\n"
    text += "Try using a computer with more RAM. If you used the '-a' option\n"
    text += "it may be possible to complete the run by removing this option."
    return text
                
class DendroBLASTTrees(object):
    def __init__(self, ogSet, nProcesses_alg, nProcess_std, qDoubleBlast):
        self.ogSet = ogSet
        self.nProcesses = nProcesses_alg
        self.nProcess_std = nProcess_std
        self.qDoubleBlast = qDoubleBlast
        # Check files exist
    
    def TreeFilename_IDs(self, iog):
        return files.FileHandler.GetOGsTreeFN(iog)
        
    def GetOGMatrices_FullParallel(self):
        """
        read the blast files as well, remove need for intermediate pickle and unpickle
        ogMatrices contains matrix M for each OG where:
            Mij = 0.5*max(Bij, Bmin_i)/Bmax_i
        """
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            ogs = self.ogSet.OGs()
            ogsPerSpecies = [[[(g, i) for i, g in enumerate(og) if g.iSp == iSp] for iSp in self.ogSet.seqsInfo.speciesToUse] for og in ogs]
            nGenes = [len(og) for og in ogs]
            nSeqs = self.ogSet.seqsInfo.nSeqsPerSpecies
            ogMatrices = [[mp.Array('d', n, lock=False) for _ in range(n)] for n in nGenes]
            blastDir_list = files.FileHandler.GetBlastResultsDir()
            cmd_queue = mp.Queue()
            for iiSp, sp1 in enumerate(self.ogSet.seqsInfo.speciesToUse):
                cmd_queue.put((iiSp, sp1, nSeqs[sp1]))
            worker_status_queue = mp.Queue()
            runningProcesses = [mp.Process(target=Worker_OGMatrices_ReadBLASTAndUpdateDistances, args=(cmd_queue, worker_status_queue, iWorker, ogMatrices, nGenes, self.ogSet.seqsInfo, blastDir_list, ogsPerSpecies, self.qDoubleBlast)) for iWorker in range(self.nProcesses)]
            for proc in runningProcesses:
                proc.start()
            rota = [None for iWorker in range(self.nProcesses)]
            unfinished = []
            while True:
                # get process alive/dead
                time.sleep(1)
                alive = [proc.is_alive() for proc in runningProcesses]
                # read latest updates from queue, update rota
                try:
                    while True:
                        status, iWorker, iTask = worker_status_queue.get(True, 0.1)
                        if status == "start":
                            rota[iWorker] = iTask
                        elif status == "finish":
                            rota[iWorker] = None
                        elif status == "empty":
                            rota[iWorker] = "empty"
                except queue.Empty:
                    pass
                # if worker is dead but didn't finish task, issue warning
                for al, r in zip(alive, rota):
                    if (not al) and (r != "empty"):
                        text = GetRAMErrorText()
                        files.FileHandler.LogFailAndExit(text)
                        unfinished.append(r)
                if not any(alive):
                    break
                
            if len(unfinished) != 0:
                files.FileHandler.LogFailAndExit()
#                print("WARNING: Computer ran out of RAM and killed OrthoFinder processes")
#                print("OrthoFinder will attempt to run these processes once more. If it is")
#                print("unsuccessful again then it will have to exit. Consider using")
#                print("the option '-a 1' or running on a machine with more RAM")
            #ogMatrices = [np.matrix(m) for m in ogMatrices]
            return ogs, ogMatrices      
                   
    def CompleteOGMatrices(self, ogs, ogMatrices):
        newMatrices = []
        for iog, (og, m) in enumerate(zip(ogs, ogMatrices)):
            # dendroblast scores
            n = m.shape[0]
            m2 = np.zeros(m.shape)
            max_og = -9e99
            for i in range(n):
                for j in range(i):
                    m2[i, j] = -np.log(m[i,j] + m[j,i])  
                    m2[j, i] = m2[i, j]  
                    max_og = max(max_og, m2[i,j])
            newMatrices.append(m2)
        return newMatrices
        
    def CompleteAndWriteOGMatrices(self, ogs, ogMatrices):
        """
        ogMatrices - each matrix is a list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        for iog, (og, m) in enumerate(zip(ogs, ogMatrices)):
            # dendroblast scores
            n = len(m)
            max_og = -9e99
            # Careful not to over-write a value and then attempt to try to use the old value
            for i in range(n):
                for j in range(i):
                    m[i][j] = -np.log(m[i][j] + m[j][i])  
                    m[j][i] = m[i][j]  
                    max_og = max(max_og, m[i][j])
            self.WritePhylipMatrix(m, [g.ToString() for g in og], files.FileHandler.GetOGsDistMatFN(iog), max_og)
        return ogMatrices
    
    @staticmethod
    def WritePhylipMatrix(m, names, outFN, max_og):
        """
        m - list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        max_og = 1.1*max_og
        sliver = 1e-6
        with open(outFN, 'w') as outfile:
            n = len(m)
            outfile.write("%d\n" % n)
            for i in range(n):
                outfile.write(names[i] + " ")
                # values could be -inf, these are the most distantly related so replace with max_og
                V = [0. + (0. if i==j else m[i][j] if m[i][j] > -9e99 else max_og) for j in range(n)] # "0. +": hack to avoid printing out "-0"
                V = [sliver if 0 < v < sliver  else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
                values = " ".join(["%.6f" % v for v in V])   
                outfile.write(values + "\n")
    
    def SpeciesTreeDistances(self, ogs, ogMatrices, method = 0):
        """
        ogMatrices - each matrix is a list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        spPairs = list(itertools.combinations(self.ogSet.seqsInfo.speciesToUse, 2))
        D = [[] for _ in spPairs]
        if method == 0:
            """ closest distance for each species pair in each orthogroup"""
            for og, m in zip(ogs, ogMatrices):
                spDict = defaultdict(list)
                for i, g in enumerate(og):
                    spDict[g.iSp].append(i)
                for (sp1, sp2), d_list in zip(spPairs, D):
                    distances = [m[i][j] for i in spDict[sp1] for j in spDict[sp2]]
                    if len(distances) > 0: d_list.append(min(distances))
#                    d_list.append(min(distances) if len(distances) > 0 else None)
        return D, spPairs
    
    def PrepareSpeciesTreeCommand(self, D, spPairs, qPutInWorkingDir=False):
        n = len(self.ogSet.seqsInfo.speciesToUse)
        M = np.zeros((n, n))
        for (sp1, sp2), d in zip(spPairs, D):
            sp1 = self.ogSet.seqsInfo.speciesToUse.index(sp1)
            sp2 = self.ogSet.seqsInfo.speciesToUse.index(sp2)
            x = np.median(d)
            M[sp1, sp2] = x
            M[sp2, sp1] = x
        speciesMatrixFN = files.FileHandler.GetSpeciesTreeMatrixFN(qPutInWorkingDir)  
        sliver = 1e-6
        with open(speciesMatrixFN, 'w') as outfile:
            outfile.write("%d\n" % n)
            for i in range(n):
                outfile.write(str(self.ogSet.seqsInfo.speciesToUse[i]) + " ")
                V = [(0. + M[i,j]) for j in range(n)]  # hack to avoid printing out "-0"
                V = [sliver if 0 < v < sliver else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
                values = " ".join(["%.6f" % v for v in V])   
                outfile.write(values + "\n")       
        treeFN = files.FileHandler.GetSpeciesTreeUnrootedFN()
        cmd = " ".join(["fastme", "-i", speciesMatrixFN, "-o", treeFN, "-N", "-w", "O"] + (["-s"] if n < 1000 else []))
        return cmd, treeFN
                
    def PrepareGeneTreeCommand(self):
        cmds = []
        ogs = self.ogSet.OGs()
        for iog in range(len(ogs)):
            nTaxa = len(ogs[iog])
            cmds.append([" ".join(["fastme", "-i", files.FileHandler.GetOGsDistMatFN(iog), "-o", files.FileHandler.GetOGsTreeFN(iog), "-N", "-w", "O"] + (["-s"] if nTaxa < 1000 else []))])
        return cmds

    @staticmethod    
    def EnoughOGsForSTAG(ogs, speciesToUse):
        nSp = len(speciesToUse)
        nSp_perOG = [len(set([g.iSp for g in og])) for og in ogs]
        return (nSp_perOG.count(nSp) >= 100)
    
    def RunAnalysis(self, qSpeciesTree=True):
        """
        Args:
            qSpeciesTree - Bool: infer the species tree
        """
        util.PrintUnderline("Calculating gene distances")
        ogs, ogMatrices_partial = self.GetOGMatrices_FullParallel()
        ogMatrices = self.CompleteAndWriteOGMatrices(ogs, ogMatrices_partial)
        del ogMatrices_partial
        util.PrintTime("Done")
        cmds_trees = self.PrepareGeneTreeCommand()
        qLessThanFourSpecies = len(self.ogSet.seqsInfo.speciesToUse) < 4
        if not qSpeciesTree:
            qSTAG = False
        elif qLessThanFourSpecies:
            qSTAG = False
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            WriteSpeciesTreeIDs_TwoThree(self.ogSet.seqsInfo.speciesToUse, spTreeFN_ids)
        else:
            qSTAG = self.EnoughOGsForSTAG(ogs, self.ogSet.seqsInfo.speciesToUse)
            if not qSTAG:
                print("Using fallback species tree inference method")
                D, spPairs = self.SpeciesTreeDistances(ogs, ogMatrices)
                cmd_spTree, spTreeFN_ids = self.PrepareSpeciesTreeCommand(D, spPairs)
                cmds_trees = [[cmd_spTree]] + cmds_trees
        del ogMatrices
        util.PrintUnderline("Inferring gene and species trees" if qSpeciesTree else "Inferring gene trees")
        parallel_task_manager.RunParallelOrderedCommandLists(self.nProcess_std, cmds_trees)
        if qSTAG:
            # Trees must have been completed
            print("")
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            stag.Run_ForOrthoFinder(files.FileHandler.GetOGsTreeDir(), files.FileHandler.GetWorkingDirectory_Write(), self.ogSet.seqsInfo.speciesToUse, spTreeFN_ids)
        if qSpeciesTree:
            util.RenameTreeTaxa(spTreeFN_ids, files.FileHandler.GetSpeciesTreeUnrootedFN(True), self.ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True)        
            return spTreeFN_ids, qSTAG
        else:      
            return None, qSTAG
            
    def SpeciesTreeOnly(self):
        qLessThanFourSpecies = len(self.ogSet.seqsInfo.speciesToUse) < 4
        if qLessThanFourSpecies:
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            WriteSpeciesTreeIDs_TwoThree(self.ogSet.seqsInfo.speciesToUse, spTreeFN_ids)
        else:
            ogs, ogMatrices_partial = self.GetOGMatrices_FullParallel()
            ogMatrices = self.CompleteOGMatrices(ogs, ogMatrices_partial)
            del ogMatrices_partial
            D, spPairs = self.SpeciesTreeDistances(ogs, ogMatrices)
            del ogMatrices
            cmd_spTree, spTreeFN_ids = self.PrepareSpeciesTreeCommand(D, spPairs, True)
            parallel_task_manager.RunOrderedCommandList([cmd_spTree])
        spTreeUnrootedFN = files.FileHandler.GetSpeciesTreeUnrootedFN(True) 
        util.RenameTreeTaxa(spTreeFN_ids, spTreeUnrootedFN, self.ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True)  
        return spTreeFN_ids

# ==============================================================================================================================      
# Main
            
def CheckUserSpeciesTree(speciesTreeFN, expSpecies):
    # File exists
    if not os.path.exists(speciesTreeFN):
        print(("Species tree file does not exist: %s" % speciesTreeFN))
        util.Fail()
    # Species in tree are unique
    try:
        t = tree.Tree(speciesTreeFN, format=1)
    except Exception as e:
        print("\nERROR: Incorrectly formated user-supplied species tree")
        print(str(e))
        util.Fail()
    actSpecies = (t.get_leaf_names())
    c = Counter(actSpecies)
    if 1 != c.most_common()[0][1]:
        print("ERROR: Species names in species tree are not unique")
        for sp, n in c.most_common():
            if 1 != n:
                print(("Species '%s' appears %d times" % (sp, n)))
        util.Fail()
    # All required species are present
    actSpecies = set(actSpecies)
    ok = True
    for sp in expSpecies:
        if sp not in actSpecies:
            print(("ERROR: '%s' is missing from species tree" % sp))
            ok = False
    # expected species are unique
    c = Counter(expSpecies)
    if 1 != c.most_common()[0][1]:
        print("ERROR: Species names are not unique")
        for sp, n in c.most_common():
            if 1 != n:
                print(("Species '%s' appears %d times" % (sp, n)))
        util.Fail()
    expSpecies = set(expSpecies)
    for sp in actSpecies:
        if sp not in expSpecies:
            print(("ERROR: Additional species '%s' in species tree" % sp))
            ok = False
    if not ok: util.Fail()
    # Tree is rooted
    if len(t.get_children()) != 2:
        print("ERROR: Species tree is not rooted")
        util.Fail()

def ConvertUserSpeciesTree(speciesTreeFN_in, speciesDict, speciesTreeFN_out):
    t = tree.Tree(speciesTreeFN_in, format=1)  
    t.prune(t.get_leaf_names())
    revDict = {v:k for k,v in speciesDict.items()}
    for sp in t:
        sp.name = revDict[sp.name]       
    t.write(outfile=speciesTreeFN_out)
    
def WriteTestDistancesFile(testFN):
    with open(testFN, 'w') as outfile:
        outfile.write("4\n1_1 0 0 0.2 0.25\n0_2 0 0 0.21 0.28\n3_1 0.2 0.21 0 0\n4_1 0.25 0.28 0 0")
    return testFN

def CanRunOrthologueDependencies(workingDir, qMSAGeneTrees, qPhyldog, qStopAfterTrees, msa_method, tree_method, recon_method, program_caller, qStopAfterAlignments):  
    # FastME
    if (not qMSAGeneTrees):
        testFN = workingDir + "SimpleTest.phy"
        WriteTestDistancesFile(testFN)
        outFN = workingDir + "SimpleTest.tre"
        if os.path.exists(outFN): os.remove(outFN)        
        if not parallel_task_manager.CanRunCommand("fastme -i %s -o %s" % (testFN, outFN), qAllowStderr=False):
            print("ERROR: Cannot run fastme")
            print("Please check FastME is installed and that the executables are in the system path.\n")
            return False
        os.remove(testFN)
        os.remove(outFN)
        fastme_stat_fn = workingDir + "SimpleTest.phy_fastme_stat.txt"
        if os.path.exists(fastme_stat_fn): os.remove(fastme_stat_fn)
    # DLCPar
    if ("dlcpar" in recon_method) and not (qStopAfterTrees or qStopAfterAlignments):
        if not parallel_task_manager.CanRunCommand("dlcpar_search --version", qAllowStderr=False):
            print("ERROR: Cannot run dlcpar_search")
            print("Please check DLCpar is installed and that the executables are in the system path.\n")
            return False
        if recon_method == "dlcpar_convergedsearch":
            capture = subprocess.Popen("dlcpar_search --version", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
            stdout = [x for x in capture.stdout]
            try:
                stdout = "".join([x.decode() for x in stdout])
            except (UnicodeDecodeError, AttributeError):
                stdout = "".join([x.encode() for x in stdout])
            version = stdout.split()[-1]
            tokens = list(map(int, version.split(".")))
            major, minor = tokens[:2]
            release = tokens[2] if len(tokens) > 2 else 0
            # require 1.0.1 or above            
            actual = (major, minor, release)
            required = [1,0,1]
            versionOK = True
            for r, a in zip(required, actual):
                if a > r:
                    versionOK = True
                    break
                elif a < r:
                    versionOK = False
                    break
                else:
                    pass
                    # need to check next level down
            if not versionOK:
                print("ERROR: dlcpar_convergedsearch requires dlcpar_search version 1.0.1 or above")
                return False                   
    
    # FastTree & MAFFT
    if qMSAGeneTrees or qPhyldog:
        testFN, temp_dir = trees_msa.WriteTestFile(workingDir)
        if msa_method == "mafft":
            if not parallel_task_manager.CanRunCommand("mafft %s" % testFN, qAllowStderr=True):
                print("ERROR: Cannot run mafft")
                print("Please check MAFFT is installed and that the executables are in the system path\n")
                return False
        elif msa_method != None:
            if not program_caller.TestMSAMethod(temp_dir, msa_method):
                print(("ERROR: Cannot run user-configured MSA method '%s'" % msa_method))
                print("Please check program is installed and that it is correctly configured in the orthofinder/config.json file\n")
                return False
        if tree_method == "fasttree":
            if qMSAGeneTrees and (not qStopAfterAlignments) and not parallel_task_manager.CanRunCommand("FastTree %s" % testFN, qAllowStderr=True):
                print("ERROR: Cannot run FastTree")
                print("Please check FastTree is installed and that the executables are in the system path\n")
                return False      
        elif tree_method != None:
            if not program_caller.TestTreeMethod(temp_dir, tree_method):
                print(("ERROR: Cannot run user-configured tree method '%s'" % tree_method))
                print("Please check program is installed and that it is correctly configured in the orthofinder/config.json file\n")
                return False
        try:
            shutil.rmtree(temp_dir)
        except OSError:
            time.sleep(1)
            shutil.rmtree(temp_dir, True)  # shutil / NFS bug - ignore errors, it's less crucial that the files are deleted
            
    if qPhyldog:
        if not parallel_task_manager.CanRunCommand("mpirun -np 1 phyldog", qAllowStderr=False):
            print("ERROR: Cannot run mpirun -np 1 phyldog")
            print("Please check phyldog is installed and that the executable is in the system path\n")
            return False
        
    return True    
        
def PrintHelp():
    print("Usage")    
    print("-----")
    print("orthologues.py orthofinder_results_directory [-t max_number_of_threads]")
    print("orthologues.py -h")
    print("\n")
    
    print("Arguments")
    print("---------")
    print("""orthofinder_results_directory
    Generate gene trees for the orthogroups, generated rooted species tree and infer ortholgues.\n""")
    
    print(("""-t max_number_of_threads, --threads max_number_of_threads
    The maximum number of processes to be run simultaneously. The deafult is %d but this 
    should be increased by the user to the maximum number of cores available.\n""" % util.nThreadsDefault))
        
    print("""-h, --help
   Print this help text""")
    util.PrintCitation()       
            
def WriteOrthologuesMatrix(fn, matrix, speciesToUse, speciesDict):
    with open(fn, csv_write_mode) as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow([""] + [speciesDict[str(index)] for index in speciesToUse])
        for ii, iSp in enumerate(speciesToUse):
            overlap = [matrix[ii, jj] for jj, jSp in enumerate(speciesToUse)]
            writer.writerow([speciesDict[str(iSp)]] + overlap)   
    

def WriteOrthologuesStats(ogSet, nOrtho_sp):
    """
    nOrtho_sp is a util.nOrtho_sp object
    """
    speciesToUse = ogSet.speciesToUse
    speciesDict = ogSet.SpeciesDict()
    d = files.FileHandler.GetOGsStatsResultsDirectory()
    WriteOrthologuesMatrix(d + "OrthologuesStats_Totals.tsv", nOrtho_sp.n, speciesToUse, speciesDict)
    WriteOrthologuesMatrix(d + "OrthologuesStats_one-to-one.tsv", nOrtho_sp.n_121, speciesToUse, speciesDict)
    WriteOrthologuesMatrix(d + "OrthologuesStats_one-to-many.tsv", nOrtho_sp.n_12m, speciesToUse, speciesDict)
    WriteOrthologuesMatrix(d + "OrthologuesStats_many-to-one.tsv", nOrtho_sp.n_m21, speciesToUse, speciesDict)
    WriteOrthologuesMatrix(d + "OrthologuesStats_many-to-many.tsv", nOrtho_sp.n_m2m, speciesToUse, speciesDict)
    # Duplications
    nodeCount = defaultdict(int)
    nodeCount_50 = defaultdict(int)
    ogCount = defaultdict(int)
    ogCount_50 = defaultdict(int)
    if not os.path.exists(files.FileHandler.GetDuplicationsFN()): return
    with open(files.FileHandler.GetDuplicationsFN(), 'rb' if PY2 else 'rt') as infile:
        reader = csv.reader(infile, delimiter="\t")
        next(reader)
        # for line in reader:
        #     try:
        #         og, node, _, support, _, _, _ = line
        #     except:
        #         print(line)
        #         raise
        for og, node, _, support, _, _, _ in reader:
            support = float(support)
            nodeCount[node] += 1
            ogCount[og] += 1
            if support >= 0.5:
                nodeCount_50[node] += 1
                ogCount_50[og] += 1
    with open(d + "Duplications_per_Species_Tree_Node.tsv", csv_write_mode) as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Species Tree Node", "Duplications (all)", "Duplications (50% support)"])
#        max_node = max([int(s[1:]) for s in nodeCount.keys()])    # Get largest node number
        for node in nodeCount:
            writer.writerow([node, nodeCount[node], nodeCount_50[node]])
    # Write on species tree
    in_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
    out_tree_fn = os.path.split(files.FileHandler.GetDuplicationsFN())[0] + "/SpeciesTree_Gene_Duplications_0.5_Support.txt"
    t = tree.Tree(in_tree_fn, format=1)
    for n in t.traverse():
        n.name = n.name + "_" + str(nodeCount_50[n.name])
    with open(out_tree_fn, 'w') as outfile:
        outfile.write(t.write(format=1)[:-1] + t.name + ";")
    with open(d + "Duplications_per_Orthogroup.tsv", csv_write_mode) as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Orthogroup", "Duplications (all)", "Duplications (50% support)"])
        if len(ogCount) > 0:
            max_og = max([int(s[2:]) for s in ogCount.keys()]) 
            pat = files.FileHandler.baseOgFormat 
            for i in range(max_og + 1):
                og = pat % i
                writer.writerow([og, ogCount[og], ogCount_50[og]])

def TwoAndThreeGeneHOGs(ogSet, st_rooted_labelled, hog_writer):
    ogs = ogSet.OGs(qInclAll=True)
    for iog, og in enumerate(ogs):
        n = len(og) 
        if n < 2 or n > 3: continue
        og_name = "OG%07d" % iog
        sp_present = set([str(g.iSp) for g in og])
        stNode = trees2ologs_of.MRCA_node(st_rooted_labelled, sp_present)
        hogs_to_write = hog_writer.get_skipped_nodes(stNode, None)  
        if len(sp_present) > 1:
            # We don't create files for 'species specific HOGs'
            st_node = trees2ologs_of.MRCA_node(st_rooted_labelled, sp_present)
            hogs_to_write = hogs_to_write + [st_node.name]
        genes = [g.ToString() for g in og] # Inefficient as will convert back again, but trivial cost I think
        hog_writer.write_hog_genes(genes, hogs_to_write, og_name)

def TwoAndThreeGeneOrthogroups(ogSet, resultsDir):
    speciesDict = ogSet.SpeciesDict()
    sequenceDict = ogSet.SequenceDict()
    ogs = ogSet.OGs(qInclAll=True)
    nOrthologues_SpPair = util.nOrtho_sp(len(ogSet.speciesToUse))
    all_orthologues = []
    d_empty = defaultdict(list)
    for iog, og in enumerate(ogs):
        n = len(og) 
        if n == 1: break
        elif n == 2:
            if og[0].iSp == og[1].iSp: continue
            # orthologues is a list of tuples of dictionaries
            # each dictionary is sp->list of genes in species
            d0 = defaultdict(list)
            d0[str(og[0].iSp)].append(str(og[0].iSeq))
            d1 = defaultdict(list)
            d1[str(og[1].iSp)].append(str(og[1].iSeq))
            orthologues = [(d0, d1, d_empty, d_empty)]  
        elif n == 3:
            sp = [g.iSp for g in og]
            c = Counter(sp) 
            nSp = len(c)
            if nSp == 3:
                g = [(str(g.iSp), str(g.iSeq)) for g in og]
                d0 = defaultdict(list)
                d0[g[0][0]].append(g[0][1])
                d1 = defaultdict(list)
                d1[g[1][0]].append(g[1][1])
                d1[g[2][0]].append(g[2][1])
                orthologues = [(d0, d1, d_empty, d_empty)]  
                d0 = defaultdict(list)
                d0[g[1][0]].append(g[1][1])
                d1 = defaultdict(list)
                d1[g[2][0]].append(g[2][1])
                orthologues.append((d0,d1, d_empty, d_empty))
            elif nSp == 2:             
                sp0, sp1 = list(c.keys())
                d0 = defaultdict(list)
                d0[str(sp0)] = [str(g.iSeq) for g in og if g.iSp == sp0]
                d1 = defaultdict(list)
                d1[str(sp1)] = [str(g.iSeq) for g in og if g.iSp == sp1]
                orthologues = [(d0, d1, d_empty, d_empty)]
            else: 
                continue # no orthologues
        elif n >= 4:
            continue
        all_orthologues.append((iog, orthologues))
    nspecies = len(ogSet.speciesToUse)
    sp_to_index = {str(sp):i for i, sp in enumerate(ogSet.speciesToUse)}
    with trees2ologs_of.OrthologsFiles(resultsDir, speciesDict, ogSet.speciesToUse, nspecies, sp_to_index) as (olog_files_handles, suspect_genes_file_handles):    
        olog_lines_tot = [["" for j in range(nspecies)] for i in range(nspecies)]
        olog_sus_lines_tot = ["" for i in range(nspecies)]
        nOrthologues_SpPair += trees2ologs_of.GetLinesForOlogFiles(all_orthologues, speciesDict, ogSet.speciesToUse, sequenceDict, 
                                                                      False, olog_lines_tot, olog_sus_lines_tot)
        # olog_sus_lines_tot will be empty
        lock_dummy = mp.Lock()
        for i in range(nspecies):
            for j in range(nspecies):
                trees2ologs_of.WriteOlogLinesToFile(olog_files_handles[i][j], olog_lines_tot[i][j], lock_dummy)
    return nOrthologues_SpPair
    
def ReconciliationAndOrthologues(recon_method, ogSet, nHighParallel, nLowParallel, iSpeciesTree=None, stride_dups=None, q_split_para_clades=False):
    """
    ogSet - info about the orthogroups, species etc
    resultsDir - where the Orthologues top level results directory will go (should exist already)
    reconTreesRenamedDir - where to put the reconcilled trees that use the gene accessions
    iSpeciesTree - which of the potential roots of the species tree is this
    method - can be dlcpar, dlcpar_deep, of_recon
    """
    speciesTree_ids_fn = files.FileHandler.GetSpeciesTreeIDsRootedFN()
    labeled_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
    util.RenameTreeTaxa(speciesTree_ids_fn, labeled_tree_fn, ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N')
    workingDir = files.FileHandler.GetWorkingDirectory_Write()    # workingDir - Orthologues working dir
    resultsDir_ologs = files.FileHandler.GetOrthologuesDirectory()
    reconTreesRenamedDir = files.FileHandler.GetOGsReconTreeDir(True)
    if "dlcpar" in recon_method:
        qDeepSearch = (recon_method == "dlcpar_convergedsearch")
        util.PrintTime("Starting DLCpar")
        dlcparResultsDir, dlcparLocusTreePat = trees2ologs_dlcpar.RunDlcpar(ogSet, speciesTree_ids_fn, workingDir, nHighParallel, qDeepSearch)
        util.PrintTime("Done DLCpar")
        spec_seq_dict = ogSet.Spec_SeqDict()
        for iog in range(len(ogSet.OGs())):
            util.RenameTreeTaxa(dlcparResultsDir + dlcparLocusTreePat % iog, files.FileHandler.GetOGsReconTreeFN(iog), spec_seq_dict, qSupport=False, qFixNegatives=False, inFormat=8, label='n')
    
        # Orthologue lists
        util.PrintUnderline("Inferring orthologues from gene trees" + (" (root %d)"%iSpeciesTree if iSpeciesTree != None else ""))
        pickleDir = files.FileHandler.GetPickleDir()
        nOrthologues_SpPair = trees2ologs_dlcpar.create_orthologue_lists(ogSet, resultsDir_ologs, dlcparResultsDir, pickleDir)  

    elif "phyldog" == recon_method:
        util.PrintTime("Starting Orthologues from Phyldog")
        nOrthologues_SpPair = trees2ologs_of.DoOrthologuesForOrthoFinder_Phyldog(ogSet, workingDir, trees2ologs_of.GeneToSpecies_dash, resultsDir_ologs, reconTreesRenamedDir)
        util.PrintTime("Done Orthologues from Phyldog")
    else:
        start = time.time()
        util.PrintTime("Starting OF Orthologues")
        qNoRecon = ("only_overlap" == recon_method)
        # The next function should not create the HOG writer and label the species tree. This should be done here and passed as arguments
        species_tree_rooted_labelled = tree.Tree(speciesTree_ids_fn)
        # Label nodes of species tree
        species_tree_rooted_labelled.name = "N0"    
        iNode = 1
        node_names = [species_tree_rooted_labelled.name]
        for n in species_tree_rooted_labelled.traverse():
            if (not n.is_leaf()) and (not n.is_root()):
                n.name = "N%d" % iNode
                node_names.append(n.name)
                iNode += 1
        # HOG Writer
        speciesDict = ogSet.SpeciesDict()
        SequenceDict = ogSet.SequenceDict()
        hog_writer = trees2ologs_of.HogWriter(species_tree_rooted_labelled, node_names, SequenceDict, speciesDict, ogSet.speciesToUse)
        nOrthologues_SpPair = trees2ologs_of.DoOrthologuesForOrthoFinder(ogSet, species_tree_rooted_labelled, trees2ologs_of.GeneToSpecies_dash, 
                                                                         stride_dups, qNoRecon, hog_writer, q_split_para_clades, nLowParallel)
        util.PrintTime("Done OF Orthologues")
        TwoAndThreeGeneHOGs(ogSet, species_tree_rooted_labelled, hog_writer)
        hog_writer.close_files()
    nOrthologues_SpPair += TwoAndThreeGeneOrthogroups(ogSet, resultsDir_ologs)
    if nLowParallel > 1 and "phyldog" != recon_method and "dlcpar" not in recon_method:
        trees2ologs_of.SortParallelFiles(nLowParallel, ogSet.speciesToUse, speciesDict)
    stop = time.time()
    # print("%fs for orthologs etc" % (stop-start))
    WriteOrthologuesStats(ogSet, nOrthologues_SpPair)
#    print("Identified %d orthologues" % nOrthologues)
        
                
def OrthologuesFromTrees(recon_method, nHighParallel, nLowParallel, userSpeciesTree_fn, qAddSpeciesToIDs, q_split_para_clades):
    """
    userSpeciesTree_fn - None if not supplied otherwise rooted tree using user species names (not orthofinder IDs)
    qUserSpTree - is the speciesTree_fn user-supplied
    
    Just infer orthologues from trees, don't do any of the preceeding steps.
    """
    speciesToUse, nSpAll, _ = util.GetSpeciesToUse(files.FileHandler.GetSpeciesIDsFN())    
    ogSet = OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesToUse, nSpAll, qAddSpeciesToIDs, idExtractor = util.FirstWordExtractor)
    if userSpeciesTree_fn != None:
        speciesDict = files.FileHandler.GetSpeciesDict()
        speciesToUseNames = [speciesDict[str(iSp)] for iSp in ogSet.speciesToUse]
        CheckUserSpeciesTree(userSpeciesTree_fn, speciesToUseNames)
        speciesTreeFN_ids = files.FileHandler.GetSpeciesTreeIDsRootedFN()
        ConvertUserSpeciesTree(userSpeciesTree_fn, speciesDict, speciesTreeFN_ids)
    util.PrintUnderline("Running Orthologue Prediction", True)
    util.PrintUnderline("Reconciling gene and species trees") 
    ReconciliationAndOrthologues(recon_method, ogSet, nHighParallel, nLowParallel, q_split_para_clades=q_split_para_clades)
    util.PrintUnderline("Writing results files")
    util.PrintTime("Writing results files")
    files.FileHandler.CleanWorkingDir2()
    
def OrthologuesWorkflow(speciesToUse, nSpAll, 
                       tree_options,
                       msa_method,
                       tree_method,
                       recon_method,
                       nHighParallel,
                       nLowParallel,
                       qDoubleBlast,
                       qAddSpeciesToIDs,
                       qTrim,
                       userSpeciesTree = None, 
                       qStopAfterSeqs = False,
                       qStopAfterAlign = False,
                       qStopAfterTrees = False, 
                       qMSA = False,
                       qPhyldog = False,
                       results_name = "",
                       q_split_para_clades=False):
    """
    1. Setup:
        - ogSet, directories
        - DendroBLASTTress - object
    2. DendrobBLAST:
        - read scores
        - RunAnalysis: Get distance matrices, do trees
    3. Root species tree
    4. Reconciliation/Orthologues
    5. Clean up
    
    Variables:
    - ogSet - all the relevant information about the orthogroups, species etc.
    """
    ogSet = OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesToUse, nSpAll, qAddSpeciesToIDs, idExtractor = util.FirstWordExtractor)
    
    tree_generation_method = "msa" if qMSA or qPhyldog else "dendroblast"
    stop_after = "seqs" if qStopAfterSeqs else "align" if qStopAfterAlign else ""
    files.FileHandler.MakeResultsDirectory2(tree_generation_method, stop_after, results_name)    
    """ === 1 === ust = UserSpeciesTree
    MSA:               Sequences    Alignments                        GeneTrees    db    SpeciesTree
    Phyldog:           Sequences    Alignments                        GeneTrees    db    SpeciesTree  
    Dendroblast:                                  DistanceMatrices    GeneTrees    db    SpeciesTree
    MSA (ust):         Sequences    Alignments                        GeneTrees    db
    Phyldog (ust):     Sequences    Alignments                        GeneTrees    db      
    Dendroblast (ust):                            DistanceMatrices    GeneTrees    db        
    """
    qDB_SpeciesTree = False
    if userSpeciesTree:
        util.PrintUnderline("Using user-supplied species tree") 
        spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
        ConvertUserSpeciesTree(userSpeciesTree, ogSet.SpeciesDict(), spTreeFN_ids)
    
    if qMSA or qPhyldog:
        qLessThanFourSpecies = len(ogSet.seqsInfo.speciesToUse) < 4
        treeGen = trees_msa.TreesForOrthogroups(tree_options, msa_method, tree_method)       
        if (not userSpeciesTree) and qLessThanFourSpecies:
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            WriteSpeciesTreeIDs_TwoThree(ogSet.seqsInfo.speciesToUse, spTreeFN_ids)
            util.RenameTreeTaxa(spTreeFN_ids, files.FileHandler.GetSpeciesTreeUnrootedFN(True), ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True)
        qDoMSASpeciesTree = (not qLessThanFourSpecies) and (not userSpeciesTree)
        util.PrintTime("Starting MSA/Trees")
        seqs_alignments_dirs = treeGen.DoTrees(ogSet.OGs(qInclAll=True), 
                                               ogSet.OrthogroupMatrix(), 
                                               ogSet.Spec_SeqDict(), 
                                               ogSet.SpeciesDict(), 
                                               ogSet.speciesToUse, 
                                               nHighParallel, 
                                               qStopAfterSeqs, 
                                               qStopAfterAlign or qPhyldog, 
                                               qDoSpeciesTree=qDoMSASpeciesTree,
                                               qTrim = qTrim) 
        util.PrintTime("Done MSA/Trees")
        if qDoMSASpeciesTree:
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
        if qStopAfterSeqs:
            print("")
            return
        elif qStopAfterAlign:
            print("")
            return 
        db = DendroBLASTTrees(ogSet, nLowParallel, nHighParallel, qDoubleBlast)
        if qDB_SpeciesTree and not userSpeciesTree and not qLessThanFourSpecies:
            util.PrintUnderline("Inferring species tree (calculating gene distances)")
            print("Loading BLAST scores")
            spTreeFN_ids = db.SpeciesTreeOnly()
        if qPhyldog:
#            util.PrintTime("Do species tree for phyldog")
#            spTreeFN_ids, spTreeUnrootedFN = db.SpeciesTreeOnly()
            if userSpeciesTree: 
                userSpeciesTree = ConvertUserSpeciesTree(userSpeciesTree, ogSet.SpeciesDict(), files.FileHandler.GetSpeciesTreeUnrootedFN())
            util.PrintTime("Starting phyldog")
            species_tree_ids_labelled_phyldog = wrapper_phyldog.RunPhyldogAnalysis(files.FileHandler.GetPhyldogWorkingDirectory(), ogSet.OGs(), speciesToUse, nHighParallel)
    else:
        db = DendroBLASTTrees(ogSet, nLowParallel, nHighParallel, qDoubleBlast)
        spTreeFN_ids, qSTAG = db.RunAnalysis(userSpeciesTree == None)
        if userSpeciesTree != None:
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
    files.FileHandler.LogWorkingDirectoryTrees()
    qSpeciesTreeSupports = False if (userSpeciesTree or qMSA or qPhyldog) else qSTAG
    """
    SpeciesTree
    spTreeFN_ids, or equivalently FileHandler.GetSpeciesTreeUnrootedFN() in all cases (user, inferred etc)
    Thus, we always have the species tree ids format
    
    With phyldog, we also have species_tree_ids_labelled_phyldog - with the node labels given by phyldog
    """    
         
    """ === 3 ===
    MSA:               RootSpeciesTree
    Phyldog:           RootSpeciesTree    
    Dendroblast:       RootSpeciesTree  
    MSA (ust):         ConvertSpeciesTreeIDs
    Phyldog (ust):     ConvertSpeciesTreeIDs
    Dendroblast (ust): ConvertSpeciesTreeIDs
    """    
    if qPhyldog:
        rootedSpeciesTreeFN = [species_tree_ids_labelled_phyldog]
        roots = [None]
        qMultiple = False
        stride_dups = None
    elif userSpeciesTree:
        rootedSpeciesTreeFN = [spTreeFN_ids]
        roots = [None]
        qMultiple = False
        stride_dups = None
    elif len(ogSet.seqsInfo.speciesToUse) == 2:
        hardcodeSpeciesTree = GetSpeciesTreeRoot_TwoTaxa(ogSet.seqsInfo.speciesToUse)
        rootedSpeciesTreeFN = [hardcodeSpeciesTree]
        roots = [None]
        qMultiple = False
        stride_dups = None
    else:
        util.PrintUnderline("Best outgroup(s) for species tree") 
        util.PrintTime("Starting STRIDE")
        roots, clusters_counter, rootedSpeciesTreeFN, nSupport, _, _, stride_dups = stride.GetRoot(spTreeFN_ids, files.FileHandler.GetOGsTreeDir(), stride.GeneToSpecies_dash, nHighParallel, qWriteRootedTree=True)
        util.PrintTime("Done STRIDE")
        nAll = sum(clusters_counter.values())
        nFP_mp = nAll - nSupport
        n_non_trivial = sum([v for k, v in clusters_counter.items() if len(k) > 1])
        if len(roots) > 1:
            print(("Observed %d well-supported, non-terminal duplications. %d support the best roots and %d contradict them." % (n_non_trivial, n_non_trivial-nFP_mp, nFP_mp)))
            print("Best outgroups for species tree:")  
        else:
            print(("Observed %d well-supported, non-terminal duplications. %d support the best root and %d contradict it." % (n_non_trivial, n_non_trivial-nFP_mp, nFP_mp)))
            print("Best outgroup for species tree:")  
        spDict = ogSet.SpeciesDict()
        for r in roots: print(("  " + (", ".join([spDict[s] for s in r]))  ))
        qMultiple = len(roots) > 1
    shutil.copy(rootedSpeciesTreeFN[0], files.FileHandler.GetSpeciesTreeIDsRootedFN())

    """
    SpeciesTree:
    We now have a list of rooted species trees: rootedSpeciesTreeFN (this should be recorded by the file handler)
    """
    if qStopAfterTrees:
        # root the gene trees using the species tree and write out their accessions - really I could remove the whole '-ot, -os, -oa' options, they are probably rarely used if ever.
        if userSpeciesTree:
            return
        # otherwise, root species tree
        resultsSpeciesTrees = []
        for i, (r, speciesTree_fn) in enumerate(zip(roots, rootedSpeciesTreeFN)):
            resultsSpeciesTrees.append(files.FileHandler.GetSpeciesTreeResultsFN(i, not qMultiple))
            util.RenameTreeTaxa(speciesTree_fn, resultsSpeciesTrees[-1], db.ogSet.SpeciesDict(), qSupport=qSpeciesTreeSupports, qFixNegatives=True)
            labeled_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
            util.RenameTreeTaxa(speciesTree_fn, labeled_tree_fn, db.ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N')
        idDict = ogSet.Spec_SeqDict()
        qHaveSupport = None 
        for iog in range(len(ogSet.OGs())):
            infn = files.FileHandler.GetOGsTreeFN(iog)
            if os.path.exists(infn):
                if qHaveSupport is None: qHaveSupport = util.HaveSupportValues(infn)
                util.RenameTreeTaxa(infn, files.FileHandler.GetOGsTreeFN(iog, True), idDict, qSupport=qHaveSupport, qFixNegatives=True)       
        files.FileHandler.CleanWorkingDir2()
    
    if qMultiple: print("\nWARNING: Multiple potential species tree roots were identified, only one will be analyed.")
    resultsSpeciesTrees = []
    i = 0
    r = roots[0]
    speciesTree_fn = rootedSpeciesTreeFN[0]
    util.PrintUnderline("Reconciling gene trees and species tree")         
    resultsSpeciesTrees.append(files.FileHandler.GetSpeciesTreeResultsFN(0, True))
    if (not userSpeciesTree) and (not qPhyldog) and len(ogSet.seqsInfo.speciesToUse) != 2:
        print(("Outgroup: " + (", ".join([spDict[s] for s in r]))))
    util.RenameTreeTaxa(speciesTree_fn, resultsSpeciesTrees[-1], db.ogSet.SpeciesDict(), qSupport=qSpeciesTreeSupports, qFixNegatives=True)
    util.PrintTime("Starting Recon and orthologues")
    ReconciliationAndOrthologues(recon_method, db.ogSet, nHighParallel, nLowParallel, i if qMultiple else None, stride_dups=stride_dups, q_split_para_clades=q_split_para_clades) 
    # util.PrintTime("Done Recon")
    
    if qMultiple:
        for i, (r, speciesTree_fn) in enumerate(zip(roots, rootedSpeciesTreeFN)):
            unanalysedSpeciesTree = files.FileHandler.GetSpeciesTreeResultsFN(i, False)
            util.RenameTreeTaxa(speciesTree_fn, unanalysedSpeciesTree, db.ogSet.SpeciesDict(), qSupport=qSpeciesTreeSupports, qFixNegatives=True, label='N')
    
    """
    SpeciesTree: If it's been inferred, there is now at least one rooted results species trees: GetSpeciesTreeResultsFN()
    """
    
    files.FileHandler.CleanWorkingDir2()
    util.PrintUnderline("Writing results files", True)
