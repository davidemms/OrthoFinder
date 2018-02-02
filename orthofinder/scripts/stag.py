#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018 David Emms
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
# For any enquiries send an email to David Emms
# david_emms@hotmail.com

import os
import sys
import glob
import datetime
import argparse
import subprocess
import numpy as np
from itertools import combinations
try:
    import tree
except ImportError:
    try:
        import ete3 as tree
    except ImportError:
        try:
            import ete2 as tree
        except ImportError:
            print("ERROR: Could not import tree library ete3 or ete2. Please install one of these.")
            sys.exit()
        
import consensus_tree as cons
        
def CanRunCommand(command, qAllowStderr = False, qPrint = True):
    if qPrint: sys.stdout.write("Test can run \"%s\"" % command)      
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]

    if len(stdout) > 0 and (qAllowStderr or len(stderr) == 0):
        if qPrint: print(" - ok")
        return True
    else:
        if qPrint: print(" - failed")
        return False
        
def CheckFastME(workingDir):
    testFN = workingDir + "SimpleTest.phy"
    with open(testFN, 'wb') as outfile:
        outfile.write("4\n1_1 0 0 0.2 0.25\n0_2 0 0 0.21 0.28\n3_1 0.2 0.21 0 0\n4_1 0.25 0.28 0 0")
    outFN = workingDir + "SimpleTest.tre"
    if os.path.exists(outFN): os.remove(outFN)  
    print("")      
    if not CanRunCommand("fastme -i %s -o %s" % (testFN, outFN), qAllowStderr=False):
        print("ERROR: Cannot run fastme")
        print('Please check "fastme" is installed and that the executable is in the system path.\n')
        return False
    os.remove(testFN)
    os.remove(outFN)
    fastme_stat_fn = workingDir + "SimpleTest.phy_fastme_stat.txt"
    if os.path.exists(fastme_stat_fn): os.remove(fastme_stat_fn)

def WritePhylipMatrix(m, names, outFN, max_og=1e6):
    """
    m - list of mp.Array  (so that each represents an nSeq x nSeq matrix
    """
    sliver = 1e-6
    with open(outFN, 'wb') as outfile:
        n = len(m)
        outfile.write("%d\n" % n)
        for i in xrange(n):
            outfile.write(names[i] + " ")
            # values could be -inf, these are the most distantly related so replace with max_og
            V = [0. + (m[i][j] if m[i][j] > -9e99 else max_og) for j in range(n)] # "0. +": hack to avoid printing out "-0"
            V = [sliver if 0 < v < sliver  else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
            values = " ".join(["%.6f" % v for v in V])   
            outfile.write(values + "\n")

#def GeneToSpecies_dash(g):
#  return g.split("_", 1)[0]
#  
#def GeneToSpecies_secondDash(g):
#  return "_".join(g.split("_", 2)[:2])
#  
#def GeneToSpecies_3rdDash(g):
#  return "_".join(g.split("_", 3)[:3])
#  
#def GeneToSpecies_dot(g):
#  return g.split(".", 1)[0]
#  
#def GeneToSpecies_hyphen(g):
#  return g.split("-", 1)[0]

class UnrecognisedGene(Exception):
    pass

class GeneToSpecies(object):
    def __init__(self, gene_map_fn):
        if not os.path.isfile(gene_map_fn):
            print("ERROR: Could not open %s" % gene_map_fn)
            sys.exit()
        self.exact = dict()
        self.startswith = dict()
        if gene_map_fn.endswith("SpeciesIDs.txt"):
            with open(gene_map_fn, 'rb') as infile:
                for line in infile:
                    _, sp = line.rstrip().rsplit(".", 1)[0].split()
                    sp = sp.replace(".", "_")
                    self.startswith[sp] = sp
        else:
            with open(gene_map_fn, 'rb') as infile:
                for i, line in enumerate(infile):
                    t = line.rstrip().split()
                    if len(t) != 2: 
                        print("ERROR: Invalid format in gene to species mapping file line %d" % (i+1)) 
                        print(line)
                        sys.exit()
                    g, sp = t
                    if g.endswith("*"):
                        self.startswith[g[:-1]] = sp
                    else:
                        self.exact[g] = sp
        self.species = list(set(self.exact.values() + self.startswith.values()))
        print("%d species in mapping file:" % len(self.species))
        for s in self.species:
            print(s)
        print("\nSTAG will infer a species tree from each gene tree with all species present") 
        self.sp_to_i = {s:i for i,s in enumerate(self.species)}
                
    def ToSpecies(self, gene):
        if gene in self.exact: return self.exact[gene]
        else:
            for k, v in self.startswith.items():
                if gene.startswith(k): return v
            # if not found then raise an Exceptions
            raise UnrecognisedGene(gene)
                
    def SpeciesToIndexDict(self):
        return self.sp_to_i
    
    def NumberOfSpecies(self):
        return len(self.species)

class GeneToSpecies_OrthoFinder(GeneToSpecies):
    def __init__(self, speciesToUse):
        # no exact conversions
        self.exact = dict()
        self.startswith = {("%d_" % sp):str(sp) for sp in speciesToUse}
        self.species = map(str,speciesToUse)
        self.sp_to_i = {s:i for i,s in enumerate(self.species)}
            

def GetDirectoryName(baseDirName, i):
    if i == 0:
        return baseDirName + os.sep
    else:
        return baseDirName + ("_%d" % i) + os.sep
        
def CreateNewWorkingDirectory(baseDirectoryName):
    dateStr = datetime.date.today().strftime("%b%d") 
    iAppend = 0
    newDirectoryName = GetDirectoryName(baseDirectoryName + dateStr, iAppend)
    while os.path.exists(newDirectoryName):
        iAppend += 1
        newDirectoryName = GetDirectoryName(baseDirectoryName + dateStr, iAppend)
    os.mkdir(newDirectoryName)
    return newDirectoryName

def GetDistances_fast(t, nSp, g_to_i):
    D = np.ones((nSp, nSp)) * 9e99
    for n in t.traverse('postorder'):
        if n.is_leaf():
            n.add_feature('d', {g_to_i[n.name]:n.dist})
#            print(n.name)
#            print(d)
        else:
            children = n.get_children()
            for ch0, ch1 in combinations(children,2):
                for sp0,dist0 in ch0.d.items():
                    for sp1,dist1 in ch1.d.items():
                        if sp0==sp1: continue
                        i = sp0 if sp0<sp1 else sp1
                        j = sp1 if sp0<sp1 else sp0
                        D[i,j] = min(D[i,j], dist0+dist1)
                spp = {k for ch in children for k in ch.d.keys()}
                d = {k:(min([ch.d[k] for ch in children if k in ch.d])+n.dist) for k in spp}
                n.add_feature('d', d)
#                print(d)
    for i in xrange(nSp):
        for j in xrange(i):
            D[i,j] = D[j, i]
        D[i,i]=0.
    return D
    
def GetDistances(t, nSp, g_to_i):
    D = np.ones((nSp, nSp)) * 9e99
    for g0 in t:
        s0 = g_to_i[g0.name]
        for g1 in t:
            s1 = g_to_i[g1.name]
            if s0 == s1: continue
            d = g0.get_distance(g1)
            if D[s0,s1] > d: 
                D[s0,s1] = d
                D[s1,s0] = d
    np.fill_diagonal(D, 0)
    return D

def ProcessTrees(dir_in, dir_matrices, dir_trees_out, GeneToSpecies, qVerbose=True):
    nSp = GeneToSpecies.NumberOfSpecies()
    s_to_i = GeneToSpecies.SpeciesToIndexDict()
    nSuccess = 0
    nNotAllPresent = 0
    nFail = 0
    if qVerbose: print("\nProcessing gene trees:")
    for fn in glob.glob(dir_in + "/*"):
        try:
            t = tree.Tree(fn)
        except tree.parser.newick.NewickError:
            print(os.path.split(fn)[1] + " - WARNING: ETE could not interpret tree file, it will be ignored")
            nFail += 1
            continue
        try:
            genes = t.get_leaf_names()
            species = map(GeneToSpecies.ToSpecies, genes)
        except UnrecognisedGene, e:
            print(os.path.split(fn)[1] + " - WARNING: unrecognised gene, %s" % e.message)
            nFail += 1
            continue
        nThis = len(set(species))
        if nThis != nSp:
#            print(os.path.split(fn)[1] + " - Only %d species, skipping" % nThis)
            nNotAllPresent += 1
            continue
        g_to_i = {g:s_to_i[s] for g,s in zip(genes, species)}
        D = GetDistances_fast(t, nSp, g_to_i)
        species_names_fastme = map(str,xrange(nSp))
        matrixFN = dir_matrices + os.path.split(fn)[1] + ".dist.phylip"
        treeOutFN = dir_trees_out + os.path.split(fn)[1] + ".tre"
        WritePhylipMatrix(D, species_names_fastme, matrixFN, max_og=1e6)
        subprocess.call("fastme -i %s -o %s -w O -s -N" % (matrixFN, treeOutFN), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        nSuccess += 1
        if qVerbose: print(os.path.split(fn)[1] + " - Processed")
    if qVerbose: print("\nExamined %d trees" % (nSuccess + nNotAllPresent + nFail))
    print("%d trees had all species present and will be used by STAG to infer the species tree\n" % nSuccess)
      
def InferSpeciesTree(tree_dir, species, outputFN):
    t = cons.ConsensusTree(tree_dir)
    for n in t:
        n.name = species[int(n.name)]
    t.write(outfile=outputFN)

def Run_ForOrthoFinder(dir_in, d_out, speciesToUse):
    dir_matrices = d_out + "Distances_SpeciesTree/"
    os.mkdir(dir_matrices)
    dir_trees_out = d_out + "SpeciesTrees_ids/"
    os.mkdir(dir_trees_out)
    gene_to_species = GeneToSpecies_OrthoFinder(speciesToUse)
    ProcessTrees(dir_in, dir_matrices, dir_trees_out, gene_to_species, qVerbose=False)
    outputFN = d_out + "STAG_SpeciesTree_ids.tre"
    InferSpeciesTree(dir_trees_out, gene_to_species.species, outputFN)
    return outputFN
    
def main(args):
    dir_in = args.gene_trees
    gene_to_species = GeneToSpecies(args.species_map)
    dir_out = CreateNewWorkingDirectory(dir_in + "/../STAG_Results")
    CheckFastME(dir_out)
    dir_matrices = dir_out + "DistanceMatrices/"
    os.mkdir(dir_matrices)
    dir_trees_out = dir_out + "Trees/"
    os.mkdir(dir_trees_out)
    ProcessTrees(dir_in, dir_matrices, dir_trees_out, gene_to_species)
    outputFN = dir_out + "SpeciesTree.tre"
    InferSpeciesTree(dir_trees_out, gene_to_species.species, outputFN)
    print("STAG species tree: " + os.path.abspath(outputFN) + "\n")

def TestTimings():
    import time
    gene_to_species = GeneToSpecies("/home/david/projects/SpeciesTreeFromAllGenes/speed_test/mapping.txt")
    nSp = gene_to_species.NumberOfSpecies()
    s_to_i = gene_to_species.SpeciesToIndexDict()
    t = tree.Tree("/home/david/projects/SpeciesTreeFromAllGenes/speed_test/OG0000117_tree_id.txt")
    genes = t.get_leaf_names()
    species = map(gene_to_species.ToSpecies, genes)
    nThis = len(set(species))
    if nThis != nSp:
        print("Not all species present: %d/%d" % (nThis, nSp))
        return
    g_to_i = {g:s_to_i[s] for g,s in zip(genes, species)}
    start = time.time()
    D2 = GetDistances_fast(t, nSp, g_to_i)
    stop = time.time()
    tnew = stop-start 
    print(tnew)
    D = GetDistances(t, nSp, g_to_i)
    torig = time.time() - stop
    print(D[0,:])
    print(D2[0,:])
    n = D.shape[0]
    for i in xrange(n):
        for j in xrange(n):
            if D[i,j] == 0:
                assert(D2[i,j] < 1e-6)
            else:
                assert(abs(D2[i,j] - D[i,j])/abs(D[i,j]) < 1e-6)
    print(torig)
    print(torig/tnew)    

if __name__ == "__main__":
    text = """
*********************************************************
*                                                       *
*      STAG: Species Tree inference from All Genes      *
*                                                       *
*********************************************************"""
    print(text[1:] + "\n")
    parser = argparse.ArgumentParser()
    parser.add_argument("species_map", help = "Map file from gene names to species names, or SpeciesIDs.txt file from OrthoFinder")
    parser.add_argument("gene_trees", help = "Directory conaining gene trees")
    args = parser.parse_args()
    main(args)
#    TestTimings()