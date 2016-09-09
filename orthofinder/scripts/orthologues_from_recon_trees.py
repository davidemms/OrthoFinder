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
import csv
import glob
import itertools
import numpy as np
import cPickle as pic
import tree
from scipy import sparse
from collections import defaultdict

import util

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]

#make_dicts checks every leaf against every other leaf to find the ancestor node and checks this node against speclist and duplist to see which dictionary the gene-pair should be placed in.
def make_dicts(dlcparResultsDir, outputDir):
    if not os.path.exists(outputDir): os.mkdir(outputDir) 
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
        speclist = [t&n for n in speclist]
        for node in speclist:
            c1, c2 = node.get_children()
            for leaf1 in c1.get_leaf_names():
                for leaf2 in c2.get_leaf_names(): 
                  Orthologs[leaf1].append(leaf2)
                  Orthologs[leaf2].append(leaf1)     
    picFN = outputDir + 'totalorthodict_fast.pic'
    with open(picFN, 'wb') as F:
        pic.dump(Orthologs, F)
    return Orthologs, picFN
    
def GetSpeciesGenesInfo(ogSet):
    speciesLabels, nSpAll = util.GetSpeciesToUse(ogSet.speciesIDsFN) 
    seqsInfo = util.GetSeqsInfo(ogSet.workingDirOF, speciesLabels, nSpAll)
    genenumbers = list(np.diff(seqsInfo.seqStartingIndices))
    genenumbers.append(seqsInfo.nSeqs - seqsInfo.seqStartingIndices[-1])
    return speciesLabels, genenumbers
        
def one_to_one_efficient(orthodict, genenumbers, speciesLabels, iSpecies, outputDir):
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
        _,iGene = map(int, queryGene.split('_'))
        for Gene in orthodict[queryGene]:
            jSpLabel,jGene = map(int,Gene.split('_'))
            jSp = speciesLabelsReverse[jSpLabel]
            if iSpecies > jSp:
                matrixlist[jSp][iGene, jGene] = 1
    for j, m in enumerate(matrixlist):    
        with open(outputDir + 'ortholog_%d_%d_matrix.pic' % (iSpecies, j), 'wb') as file:
            pic.dump(m, file)
    return matrixlist   
    
def multiply(specmin, specmax, matrixDir):
    with open(matrixDir + 'ortholog_%d_%d_matrix.pic' % (specmax, specmin), 'rb') as F:
        M = pic.load(F)    
    M = M.tocsr()
    product = M.dot(M.transpose())
    return product, M

def WriteOrthologues(resultsDir, spec1, spec2, orthologues, ogSet):
    speciesDict = ogSet.SpeciesDict()
    id_to_og = ogSet.ID_to_OG_Dict()
    sequenceDict = ogSet.SequenceDict()
    d1 = resultsDir + "Orthologues_" + speciesDict[str(spec1)] + "/"
    d2 = resultsDir + "Orthologues_" + speciesDict[str(spec2)] + "/"
    with open(d1 + '%s__v__%s.csv' % (speciesDict[str(spec1)], speciesDict[str(spec2)]), 'wb') as outfile1, open(d2 + '%s__v__%s.csv' % (speciesDict[str(spec2)], speciesDict[str(spec1)]), 'wb') as outfile2:
        writer1 = csv.writer(outfile1)
        writer2 = csv.writer(outfile2)
        writer1.writerow(("Orthogroup", speciesDict[str(spec1)], speciesDict[str(spec2)]))
        writer2.writerow(("Orthogroup", speciesDict[str(spec2)], speciesDict[str(spec1)]))
        for genes1, genes2 in orthologues:
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
    for iIndex in xrange(matrix.shape[0]):
        if matrix[iIndex, iIndex] == 0 or iIndex in done: continue
        orthologuesSp1, orthologuesSp2 = GetOrthologues(orig_matrix, orig_matrix_csc, iIndex)
        done.update(orthologuesSp1)
        orthologues.append((orthologuesSp1, orthologuesSp2))
    return orthologues
        
def species_find_all(ogSet, matrixDir, resultsDir):
    speciesDict = ogSet.SpeciesDict()
    # Calls multiply and find_all on each species pair, and appends the numbers from find_all's output to the relevant csv lists.
    speciesIDs = ogSet.speciesToUse
    nspecies = len(speciesIDs)           
    for index1 in xrange(nspecies):
        d = resultsDir + "Orthologues_" + speciesDict[str(speciesIDs[index1])]
        if not os.path.exists(d): os.mkdir(d)
    for index1, index2 in itertools.product(xrange(nspecies), xrange(nspecies)):      
        if index1 >= index2: continue
        product, M = multiply(index1, index2, matrixDir)
        orthologues = find_all(product, M)
        WriteOrthologues(resultsDir, speciesIDs[index2], speciesIDs[index1], orthologues, ogSet)    
    
def get_orthologue_lists(ogSet, resultsDir, dlcparResultsDir, workingDir):
    # -> Matrices
    matrixDir = workingDir + "matrices_orthologues/"
    orthodict, orthodictFN = make_dicts(dlcparResultsDir, matrixDir)
    
    # -> dictionary
    speciesLabels, genenumbers = GetSpeciesGenesInfo(ogSet)
    for iSpecies in xrange(len(speciesLabels)):
        one_to_one_efficient(orthodict, genenumbers, speciesLabels, iSpecies, matrixDir)
        
    # -> csv files
    species_find_all(ogSet, matrixDir, resultsDir)
    
        
