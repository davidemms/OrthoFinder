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

import re
import os
import sys
import csv
import glob
import itertools
import numpy as np
try:
    import cPickle as pic
except ImportError:
    import pickle as pic
from scipy import sparse
from collections import defaultdict

from . import tree
from . import util
from . import files
from . import parallel_task_manager

PY2 = sys.version_info <= (3,)
csv_write_mode = 'wb' if PY2 else 'wt'

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]


# ==============================================================================================================================      
# Run DLCPar

def GetTotalLength(t):
    return sum([node.dist for node in t])
  
def AllEqualBranchLengths(t):
    lengths = [node.dist for node in t]
    return (len(lengths) > 1 and len(set(lengths)) == 1)

def RootGeneTreesArbitrarily(nOGs, outputDir):
    filenames = [files.FileHandler.GetOGsTreeFN(i) for i in range(nOGs)]
    outFilenames = [outputDir + os.path.split(files.FileHandler.GetOGsTreeFN(i))[1] for i in range(nOGs)]
    treeFilenames = [fn for fn in filenames if fn.endswith(".txt")]
    nErrors = 0
    with open(outputDir + 'root_errors.txt', 'w') as errorfile:
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

def WriteGeneSpeciesMap(d, speciesDict):
    fn = d + "GeneMap.smap"
    iSpecies = list(speciesDict.keys())
    with open(fn, 'w') as outfile:
        for iSp in iSpecies:
            outfile.write("%s_*\t%s\n" % (iSp, iSp))
    return fn

def RunDlcpar(ogSet, speciesTreeFN, workingDir, nParallel, qDeepSearch):
    """
    
    Implementation:
    - (skip: label species tree)
    - sort out trees (midpoint root, resolve plytomies etc)
    - run
    
    """
    ogs = ogSet.OGs()
    nOGs = len(ogs)
    dlcparResultsDir = workingDir + 'dlcpar/'
    if not os.path.exists(dlcparResultsDir): os.mkdir(dlcparResultsDir)
    RootGeneTreesArbitrarily(nOGs, dlcparResultsDir)
    spec_seq_dict = ogSet.Spec_SeqDict()
    for iog in range(len(ogs)):
        util.RenameTreeTaxa(files.FileHandler.GetOGsTreeFN(iog), files.FileHandler.GetOGsTreeFN(iog, True), spec_seq_dict, qSupport=False, qFixNegatives=True, qViaCopy=False)
    geneMapFN = WriteGeneSpeciesMap(dlcparResultsDir, ogSet.SpeciesDict())
    filenames = [dlcparResultsDir + os.path.split(files.FileHandler.GetOGsTreeFN(i))[1] for i in range(nOGs)]
    if qDeepSearch:
        nTaxa = [len(og) for og in ogs[:nOGs]]
        nIter =     [1000 if n < 25 else 25000 if n < 200 else 50000 for n in nTaxa]
        nNoImprov = [ 100 if n < 25 else  1000 if n < 200 else  2000 for n in nTaxa]
        dlcCommands = ['dlcpar_search -s %s -S %s -D 1 -C 0.125 %s -I .txt -i %d --nprescreen 100 --nconverge %d' % (speciesTreeFN, geneMapFN, fn, i, n) for (fn, i, n) in zip(filenames, nIter, nNoImprov)]
    else:
        dlcCommands = ['dlcpar_search -s %s -S %s -D 1 -C 0.125 %s -I .txt -x 1' % (speciesTreeFN, geneMapFN, fn) for fn in filenames]
    parallel_task_manager.RunParallelOrderedCommandLists(nParallel, [[c] for c in dlcCommands])
    return dlcparResultsDir, "OG%07d_tree_id.dlcpar.locus.tree"


# ==============================================================================================================================      
# Analyse DLCPar results

#make_dicts checks every leaf against every other leaf to find the ancestor node and checks this node against speclist and duplist to see which dictionary the gene-pair should be placed in.
def make_dicts(dlcparResultsDir):
    treeFNs = glob.glob(dlcparResultsDir + '/*.locus.tree')
    recons = glob.glob(dlcparResultsDir + '/*.locus.recon')
    treeFNs.sort(key = natural_sort_key)
    recons.sort(key = natural_sort_key)
    Orthologs = defaultdict(list)
    for treeFN, recon in zip(treeFNs, recons):
        t = tree.Tree(treeFN, format = 8)
        reconlist = []
        with open(recon, 'r') as rec:
            for line in rec:
                reconlist.append(line.split())
        speclist = [x[0] for x in reconlist if 'spec' == x[2]]
        speclist = [t&n if n in t else t for n in speclist]
        for node in speclist:
            c1, c2 = node.get_children()
            for leaf1 in c1.get_leaf_names():
                for leaf2 in c2.get_leaf_names(): 
                  Orthologs[leaf1].append(leaf2)
                  Orthologs[leaf2].append(leaf1)     
    return Orthologs
    
def GetSpeciesGenesInfo():
    speciesLabels, nSpAll, _ = util.GetSpeciesToUse(files.FileHandler.GetSpeciesIDsFN()) 
    seqsInfo = util.GetSeqsInfo(files.FileHandler.GetSpeciesSeqsDir(), speciesLabels, nSpAll)
    genenumbers = list(np.diff(seqsInfo.seqStartingIndices))
    genenumbers.append(seqsInfo.nSeqs - seqsInfo.seqStartingIndices[-1])
    return speciesLabels, genenumbers
        
def one_to_one_efficient(orthodict, genenumbers, speciesLabels, iSpecies, pickleDir):
    """ speciesLabels is an ordered list of the speciesIDs
        try to mostly deal with iSpecies which is the ordinal number not the label it is given
    """
    #Creates all matrices and appends them to matrixlist.
    util.PrintTime("Processing orthologues for species %d" % iSpecies)
    matrixlist = []
    numspecies = len(speciesLabels)
    speciesLabelsReverse = {label:i for i, label in enumerate(speciesLabels)}
    for j in range(numspecies):
        if iSpecies > j:
            matrixlist.append(sparse.lil_matrix((genenumbers[iSpecies], genenumbers[j]), dtype=np.dtype(np.int8)))
        else:
            matrixlist.append(None)
    #Fill matrices with orthodata
    iSpecieslist = [x for x in orthodict if x.startswith('%d_' % speciesLabels[iSpecies])]
    for count, queryGene in enumerate(iSpecieslist):
        _,iGene = list(map(int, queryGene.split('_')))
        for Gene in orthodict[queryGene]:
            jSpLabel,jGene = list(map(int,Gene.split('_')))
            jSp = speciesLabelsReverse[jSpLabel]
            if iSpecies > jSp:
                matrixlist[jSp][iGene, jGene] = 1
    for j, m in enumerate(matrixlist):    
        with open(pickleDir + 'ortholog_%d_%d_matrix.pic' % (iSpecies, j), 'wb') as file:
            pic.dump(m, file)
    return matrixlist   
    
def multiply(specmin, specmax, pickleDir):
    with open(pickleDir + 'ortholog_%d_%d_matrix.pic' % (specmax, specmin), 'rb') as F:
        M = pic.load(F)    
    M = M.tocsr()
    product = M.dot(M.transpose())
    return product, M

def WriteOrthologues(resultsDir, spec1, spec2, orthologues, ogSet, nOrtho_sp, i, j): 
    """
    spec1 (int) - the ID for the first species
    spec2 (int) - the ID for the second species
    i (int) - the ordinal for the first species (after excluded species have been removed)
    j (int) - the ordinal for the second species (after excluded species have been removed)
    """
    speciesDict = ogSet.SpeciesDict()
    id_to_og = ogSet.ID_to_OG_Dict()
    sequenceDict = ogSet.SequenceDict()
    d1 = resultsDir + "Orthologues_" + speciesDict[str(spec1)] + "/"
    d2 = resultsDir + "Orthologues_" + speciesDict[str(spec2)] + "/"
    with open(d1 + '%s__v__%s.tsv' % (speciesDict[str(spec1)], speciesDict[str(spec2)]), csv_write_mode) as outfile1, open(d2 + '%s__v__%s.tsv' % (speciesDict[str(spec2)], speciesDict[str(spec1)]), csv_write_mode) as outfile2:
        writer1 = csv.writer(outfile1, delimiter="\t")
        writer2 = csv.writer(outfile2, delimiter="\t")
        writer1.writerow(("Orthogroup", speciesDict[str(spec1)], speciesDict[str(spec2)]))
        writer2.writerow(("Orthogroup", speciesDict[str(spec2)], speciesDict[str(spec1)]))
        for genes1, genes2 in orthologues:
            n1 = len(genes1)
            n2 = len(genes2)
            nOrtho_sp.n[i, j] += n1 
            nOrtho_sp.n[j, i] += n2 
            if n1 == 1 and n2 == 1:
                nOrtho_sp.n_121[i, j] += 1
                nOrtho_sp.n_121[j, i] += 1
            elif n1 == 1:
                nOrtho_sp.n_12m[i, j] += 1
                nOrtho_sp.n_m21[j, i] += n2
            elif n2 == 1:
                nOrtho_sp.n_m21[i, j] += n1
                nOrtho_sp.n_12m[j, i] += 1
            else:
                nOrtho_sp.n_m2m[i, j] += n1
                nOrtho_sp.n_m2m[j, i] += n2
            og = "OG%07d" % id_to_og["%d_%d" % (spec1, genes1[0])]
            writer1.writerow((og, ", ".join([sequenceDict["%d_%d" % (spec1, o)] for o in genes1]), ", ".join([sequenceDict["%d_%d" % (spec2, o)] for o in genes2])))
            writer2.writerow((og, ", ".join([sequenceDict["%d_%d" % (spec2, o)] for o in genes2]), ", ".join([sequenceDict["%d_%d" % (spec1, o)] for o in genes1])))

def GetOrthologues(orig_matrix, orig_matrix_csc, index):
    orthologues = orig_matrix.getrowview(index).nonzero()[1]
    index = orthologues[0]
    originalSpeciesGenes = orig_matrix_csc.getcol(index).nonzero()[0]  # SparseEfficiencyWarning: changing the sparsity structure of a csc_matrix is expensive. lil_matrix is more efficient
    return (originalSpeciesGenes, orthologues)
                
#Takes in the output of multiply, finds all of the orthology relationships which it writes to textfiles and returns the number of each type of orthology relationship.
def find_all(matrix, orig_matrix):
    orig_matrix_csc = orig_matrix.tocsc()   # most efficient for column access
    orig_matrix = orig_matrix.tolil()       # most efficient for rowaccess
    orthologues = []
    done = set()
    for iIndex in range(matrix.shape[0]):
        if matrix[iIndex, iIndex] == 0 or iIndex in done: continue
        orthologuesSp1, orthologuesSp2 = GetOrthologues(orig_matrix, orig_matrix_csc, iIndex)
        done.update(orthologuesSp1)
        orthologues.append((orthologuesSp1, orthologuesSp2))
    return orthologues
        
def species_write_all(ogSet, pickleDir, resultsDir):
    speciesDict = ogSet.SpeciesDict()
    # Calls multiply and find_all on each species pair, and appends the numbers from find_all's output to the relevant csv lists.
    speciesIDs = ogSet.speciesToUse
    nspecies = len(speciesIDs)           
    nOrthologues_SpPair = util.nOrtho_sp(nspecies)
    for index1 in range(nspecies):
        d = resultsDir + "Orthologues_" + speciesDict[str(speciesIDs[index1])]
        if not os.path.exists(d): os.mkdir(d)
    for index1, index2 in itertools.product(range(nspecies), range(nspecies)):      
        if index1 >= index2: continue
        product, M = multiply(index1, index2, pickleDir)
        orthologues = find_all(product, M)
        WriteOrthologues(resultsDir, speciesIDs[index2], speciesIDs[index1], orthologues, ogSet, nOrthologues_SpPair, index2 ,index1)   
    return nOrthologues_SpPair
    
def create_orthologue_lists(ogSet, resultsDir, dlcparResultsDir, pickleDir):
    # -> Matrices
#    matrixDir = workingDir + "matrices_orthologues/" 
    orthodict = make_dicts(dlcparResultsDir)
    
    # -> dictionary
    speciesLabels, genenumbers = GetSpeciesGenesInfo()
    for iSpecies in range(len(speciesLabels)):
        one_to_one_efficient(orthodict, genenumbers, speciesLabels, iSpecies, pickleDir)
        
    # -> csv files
    nOrthologues_SpPair = species_write_all(ogSet, pickleDir, resultsDir)
    for fn in glob.glob(pickleDir + "ortholog_*.pic"):
        if os.path.exists(fn): os.remove(fn)
    return nOrthologues_SpPair
    
    
    
        
