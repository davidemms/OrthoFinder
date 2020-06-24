#!/usr/bin/env python3
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

from . import tree, newick, util, parallel_task_manager
from . import consensus_tree as cons

# import tree, newick, util, parallel_task_manager
# import consensus_tree as cons
        
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
    with open(testFN, 'w') as outfile:
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
    with open(outFN, 'w') as outfile:
        n = len(m)
        outfile.write("%d\n" % n)
        for i in range(n):
            outfile.write(names[i] + " ")
            # values could be -inf, these are the most distantly related so replace with max_og
            V = [0. + (m[i][j] if m[i][j] > -9e99 else max_og) for j in range(n)] # "0. +": hack to avoid printing out "-0"
            V = [sliver if v < sliver  else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
            values = " ".join(["%.6f" % v for v in V])   
            outfile.write(values + "\n")

class UnrecognisedGene(Exception):
    pass

class GeneToSpecies(object):
    def __init__(self, gene_map_fn):
        if not os.path.isfile(gene_map_fn):
            print(("ERROR: Could not open %s" % gene_map_fn))
            sys.exit()
        self.exact = dict()
        self.startswith = dict()
        if gene_map_fn.endswith("SpeciesIDs.txt"):
            with open(gene_map_fn, 'r') as infile:
                for line in infile:
                    _, sp = line.rstrip().rsplit(".", 1)[0].split()
                    sp = sp.replace(".", "_")
                    self.startswith[sp] = sp
        else:
            with open(gene_map_fn, 'r') as infile:
                for i, line in enumerate(infile):
                    t = line.rstrip().split()
                    if len(t) != 2: 
                        print(("ERROR: Invalid format in gene to species mapping file line %d" % (i+1))) 
                        print(line)
                        sys.exit()
                    g, sp = t
                    if g.endswith("*"):
                        self.startswith[g[:-1]] = sp
                    else:
                        self.exact[g] = sp
        self.species = sorted(list(set(list(self.exact.values()) + list(self.startswith.values()))))
        print(("%d species in mapping file:" % len(self.species)))
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
        self.species = list(map(str,speciesToUse))
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
            n.add_feature('d', {g_to_i[n.name]:max(0.0, n.dist)})
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
                d = {k:(min([ch.d[k] for ch in children if k in ch.d])+max(0.0, n.dist)) for k in spp}
                n.add_feature('d', d)
    for i in range(nSp):
        for j in range(i):
            D[i,j] = D[j, i]
        D[i,i]=0.
    return D

def ProcessTrees(dir_in, dir_matrices, dir_trees_out, GeneToSpecies, qVerbose=True, qSkipSingleCopy=False, qForOF=False):
    nSp = GeneToSpecies.NumberOfSpecies()
    s_to_i = GeneToSpecies.SpeciesToIndexDict()
    nSuccess = 0
    nNotAllPresent = 0
    nFail = 0
    if qVerbose: print("\nProcessing gene trees:")
    for fn in glob.glob(dir_in + "/*"):
        try:
            t = tree.Tree(fn)
        except newick.NewickError:
            print((os.path.split(fn)[1] + " - WARNING: ETE could not interpret tree file, it will be ignored"))
            nFail += 1
            continue
        try:
            genes = t.get_leaf_names()
            species = list(map(GeneToSpecies.ToSpecies, genes))
        except UnrecognisedGene as e:
            print((os.path.split(fn)[1] + " - WARNING: unrecognised gene, %s" % str(e)))
            nFail += 1
            continue
        nThis = len(set(species))
        if nThis != nSp:
#            print(os.path.split(fn)[1] + " - Only %d species, skipping" % nThis)
            nNotAllPresent += 1
            continue
        if qSkipSingleCopy and nThis == len(genes):
            # Single copy - don't recalculate the tree
            treeOutFN = dir_trees_out + os.path.split(fn)[1] + ".tre"
            for n in t:
                n.name = s_to_i[GeneToSpecies.ToSpecies(n.name)]
            t.write(outfile = treeOutFN, format=5)
            if qVerbose: print((os.path.split(fn)[1] + " - Processed"))
            continue
        g_to_i = {g:s_to_i[s] for g,s in zip(genes, species)}
        D = GetDistances_fast(t, nSp, g_to_i)
        species_names_fastme = list(map(str,range(nSp)))
        matrixFN = dir_matrices + os.path.split(fn)[1] + ".dist.phylip"
        treeOutFN = dir_trees_out + os.path.split(fn)[1] + ".tre"
        WritePhylipMatrix(D, species_names_fastme, matrixFN, max_og=1e6)
        command = "fastme -i %s -o %s -w O -s -n" % (matrixFN, treeOutFN)
        if qForOF:
            parallel_task_manager.RunCommand(command, qPrintOnError=True, qPrintStderr=True)
        else:
            popen = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            popen.communicate()
        nSuccess += 1
        if qVerbose: print((os.path.split(fn)[1] + " - Processed"))
    if qVerbose: print(("\nExamined %d trees" % (nSuccess + nNotAllPresent + nFail)))
    print(("%d trees had all species present and will be used by STAG to infer the species tree" % nSuccess))
      
def Astral(tree_dir, astral_jar_file, qForOF=False):
    treesFN = tree_dir + "../TreesFile.txt"
    with open(treesFN, 'w') as outfile:
        for fn in glob.glob(tree_dir + "/*"):
            t = tree.Tree(fn)    
            outfile.write(t.write(format=9) + "\n")
    speciesTreeFN = tree_dir + "../SpeciesTree_ids.txt"
    command = " ".join(["java", "-Xmx6000M", "-jar", astral_jar_file, "-i", treesFN, "-o", speciesTreeFN])
    if qForOF:
        parallel_task_manager.RunCommand(command, qPrintOnError=True, qPrintStderr=True)
    else:
        subprocess.call(command, shell=True)
    return tree.Tree(speciesTreeFN)     
      
def InferSpeciesTree(tree_dir, species, outputFN, astral_jar=None):
    if astral_jar == None:
        t = cons.ConsensusTree(tree_dir)
    else:
        t = Astral(tree_dir, astral_jar)
    for n in t:
        n.name = species[int(n.name)]
    t.write(outfile=outputFN)

def Run_ForOrthoFinder(dir_in, d_working, speciesToUse, speciesTreeIds_FN_out):
    dir_matrices = d_working + "Distances_SpeciesTree/"
    os.mkdir(dir_matrices)
    dir_trees_out = d_working + "SpeciesTrees_ids/"
    os.mkdir(dir_trees_out)
    gene_to_species = GeneToSpecies_OrthoFinder(speciesToUse)
    ProcessTrees(dir_in, dir_matrices, dir_trees_out, gene_to_species, qVerbose=False, qForOF=True)
    InferSpeciesTree(dir_trees_out, gene_to_species.species, speciesTreeIds_FN_out)
    
def main(args):
    dir_in = args.gene_trees
    astral_jar = None
    gene_to_species = GeneToSpecies(args.species_map)
    dir_out = CreateNewWorkingDirectory(dir_in + "/../STAG_Results")
    CheckFastME(dir_out)
    dir_matrices = dir_out + "DistanceMatrices/"
    os.mkdir(dir_matrices)
    dir_trees_out = dir_out + "Trees/"
    os.mkdir(dir_trees_out)
    ProcessTrees(dir_in, dir_matrices, dir_trees_out, gene_to_species, qVerbose=(not args.quiet))
    outputFN = dir_out + "SpeciesTree.tre"
#    if args.astral_jar == None:
    if astral_jar == None:
        InferSpeciesTree(dir_trees_out, gene_to_species.species, outputFN)
    else:
        InferSpeciesTree(dir_trees_out, gene_to_species.species, outputFN, astral_jar=astral_jar)
    print(("STAG species tree: " + os.path.abspath(outputFN) + "\n"))

if __name__ == "__main__":
    text = """
*********************************************************
*                                                       *
*      STAG: Species Tree inference from All Genes      *
*                                                       *
*********************************************************"""
    print((text[1:] + "\n"))
    parser = argparse.ArgumentParser()
    parser.add_argument("species_map", help = "Map file from gene names to species names, or SpeciesIDs.txt file from OrthoFinder")
    parser.add_argument("gene_trees", help = "Directory conaining gene trees")
    parser.add_argument("-q", "--quiet", help = "Only print sparse output", action="store_true")
#    parser.add_argument("-a", "--astral_jar", help = "ASTRAL jar file. Use ASTRAL to combine STAG species tree estimates instead of greedy consensus tree.")
    args = parser.parse_args()
    main(args)