# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:11:11 2017

@author: david

Perform directed 'reconciliation' first and then apply EggNOG method

1 - root gene trees on outgroup: unique one this time
2 - infer orthologues
"""
import os
import csv
import glob
import argparse
import operator
import itertools
import numpy as np
import multiprocessing as mp
from collections import Counter, defaultdict

import tree as tree_lib
import resolve, util

def GeneToSpecies_dash(g):
  return g.split("_", 1)[0]
  
OrthoFinderIDs = GeneToSpecies_dash
  
def GeneToSpecies_secondDash(g):
  return "_".join(g.split("_", 2)[:2])
  
def GeneToSpecies_3rdDash(g):
  return "_".join(g.split("_", 3)[:3])
  
def GeneToSpecies_dot(g):
  return g.split(".", 1)[0]
  
def GeneToSpecies_hyphen(g):
  return g.split("-", 1)[0]  
    
class RootMap(object):
    def __init__(self, setA, setB, GeneToSpecies):
        self.setA = setA
        self.setB = setB
        self.GeneToSpecies = GeneToSpecies
        
    def GeneMap(self, gene_name):
        sp = self.GeneToSpecies(gene_name)
        if sp in self.setA: return True 
        elif sp in self.setB: return False
        else: 
            print(gene_name)
            print(sp)
            raise Exception
        
def StoreSpeciesSets(t, GeneMap, tag="sp_"):
    tag_up = tag + "up"
    tag_down = tag + "down"  
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node.add_feature(tag_down, {GeneMap(node.name)})
        elif node.is_root():
            continue
        else:
            node.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in node.get_children()]))
    for node in t.traverse('preorder'):
        if node.is_root():
            node.add_feature(tag_up, set())
        else:
            parent = node.up
            if parent.is_root():
                others = [ch for ch in parent.get_children() if ch != node]
                node.add_feature(tag_up, set.union(*[other.__getattribute__(tag_down) for other in others]))
            else:
                others = [ch for ch in parent.get_children() if ch != node]
                sp_downs = set.union(*[other.__getattribute__(tag_down) for other in others])
                node.add_feature(tag_up, parent.__getattribute__(tag_up).union(sp_downs))
    t.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in t.get_children()]))

def OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, n1, n2):
    f_dup = len(sp_up.intersection(sett1)) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * len(sp_down.intersection(sett2)) * N_recip
    f_a = len(sp_up.intersection(sett1)) * (n2-len(sp_up.intersection(sett2))) * (n1-len(sp_down.intersection(sett1))) * len(sp_down.intersection(sett2)) * N_recip
    f_b = (n1-len(sp_up.intersection(sett1))) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * (n2-len(sp_down.intersection(sett2))) * N_recip
    choice = (f_dup, f_a, f_b)
#    print(choice)
    return max(choice)
    
def GetRoots(tree, species_tree_rooted, GeneToSpecies):
    """
    Allow non-binary gene or species trees.
    (A,B,C) => consider splits A|BC, B|AC, C|AB - this applies to gene and species tree
    If a clean ingroup/outgroup split cannot be found then score root by geometric mean of fraction of expected species actually 
    observed for the two splits
    """
    speciesObserved = set([GeneToSpecies(g) for g in tree.get_leaf_names()])
    if len(speciesObserved) == 1:
        return []
    
    # use species tree to find correct outgroup according to what species are present in the gene tree
    n = species_tree_rooted
    children = n.get_children()
    leaves = [set(ch.get_leaf_names()) for ch in children]
    have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]
    while sum(have) < 2:
        n = children[have.index(True)]
        children = n.get_children()
        leaves = [set(ch.get_leaf_names()) for ch in children]
        have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]

    # Get splits to look for
    roots_list = []
    scores_list = []   # the fraction completeness of the two clades
#    roots_set = set()
    for i in xrange(len(leaves)):
        t1 = leaves[i]
        t2 = set.union(*[l for j,l in enumerate(leaves) if j!=i])
        # G - set of species in gene tree
        # First relevant split in species tree is (A,B), such that A \cap G \neq \emptyset and A \cap G \neq \emptyset
        # label all nodes in gene tree according the whether subsets of A, B or both lie below node
        StoreSpeciesSets(tree, GeneToSpecies)   # sets of species
        root_mapper = RootMap(t1, t2, GeneToSpecies)    
        sett1 = set(t1)
        sett2 = set(t2)
        nt1 = float(len(t1))
        nt2 = float(len(t2))
        N_recip = 1./(nt1*nt1*nt2*nt2)
        GeneMap = root_mapper.GeneMap
        StoreSpeciesSets(tree, GeneMap, "inout_") # ingroup/outgroup identification
        # find all possible locations in the gene tree at which the root should be

        T = {True,}
        F = {False,}
        TF = set([True, False])
        for m in tree.traverse('postorder'):
            if m.is_leaf(): 
                if len(m.inout_up) == 1 and m.inout_up != m.inout_down:
                    # this is the unique root
                    return [m]
            else:
                if len(m.inout_up) == 1 and len(m.inout_down) == 1 and m.inout_up != m.inout_down:
                    # this is the unique root
                    return [m]
                nodes = m.get_children() if m.is_root() else [m] + m.get_children()
                clades = [ch.inout_down for ch in nodes] if n.is_root() else [m.inout_up] + [ch.inout_down for ch in m.get_children()]               
                # do we have the situation A | B or (A,B),S?
                if len(nodes) == 3:
                    if all([len(c) == 1 for c in clades]) and T in clades and F in clades:
                        # unique root
                        if clades.count(T) == 1:
                            return [nodes[clades.index(T)]]
                        else:
                            return [nodes[clades.index(F)]]
                    elif T in clades and F in clades:
                        #AB-(A,B) or B-(AB,A)
                        ab = [c == TF for c in clades]
                        i = ab.index(True)
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(np.sqrt(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2)))
                    elif clades.count(TF) >= 2:  
                        # (A,A,A)-excluded, (A,A,AB)-ignore as want A to be bigest without including B, (A,AB,AB), (AB,AB,AB) 
                        i = 0
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(np.sqrt(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2)))
                elif T in clades and F in clades:
                    roots_list.append(m)
                    scores_list.append(0)  # last choice
    # If we haven't found a unique root then use the scores for completeness of ingroup/outgroup to root
    if len(roots_list) == 0: 
        return [] # This shouldn't occur
    return [sorted(zip(scores_list, roots_list), reverse=True)[0][1]]
                
def WriteQfO2(orthologues_list_pairs_list, outfilename, qAppend = True):
    """ takes a list where each entry is a pair, (genes1, genes2), which are orthologues of one another
    """
    with open(outfilename, 'ab' if qAppend else 'wb') as outfile:
        for gs1, gs2 in orthologues_list_pairs_list:
            for g1 in gs1:
                g1 = g1.split()[0].split("_")[-1]
                for g2 in gs2:
                    g2 = g2.split()[0].split("_")[-1]
                    outfile.write("%s\t%s\n" % (g1, g2))
    
def GetGeneToSpeciesMap(args):
    GeneToSpecies = GeneToSpecies_dash
    if args.separator and args.separator == "dot":
        GeneToSpecies = GeneToSpecies_dot  
    elif args.separator and args.separator == "second_dash":
        GeneToSpecies = GeneToSpecies_secondDash  
    elif args.separator and args.separator == "3rd_dash":
        GeneToSpecies = GeneToSpecies_3rdDash  
    elif args.separator and args.separator == "hyphen":
        GeneToSpecies = GeneToSpecies_hyphen  
    return GeneToSpecies
  
def OverlapSize(node, GeneToSpecies):  
    descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in node.get_children()]
    return len(descendents[0].intersection(descendents[1])), descendents[0], descendents[1]
          
def GetRoot(tree, species_tree_rooted, GeneToSpecies):
        roots = GetRoots(tree, species_tree_rooted, GeneToSpecies)
#        print("%d roots" % len(roots))
        # 1. Store distance of each root from each leaf
#        for root in roots:
#            print(root)
#        return roots[0]
        if len(roots) > 0:
            root_dists = [r.get_closest_leaf()[1] for r in roots]
            i, _ = max(enumerate(root_dists), key=operator.itemgetter(1))
            return roots[i]
        else:
            return None # single species tree
    
def GetOrthologues_for_tree(iog, treeFN, species_tree_rooted, GeneToSpecies, qWrite=False, dupsWriter=None, seqIDs=None, spIDs=None, all_stride_dup_genes=None):
    """ if dupsWriter != None then seqIDs and spIDs must also be provided"""
    qPrune=True
    orthologues = []
    if (not os.path.exists(treeFN)) or os.stat(treeFN).st_size == 0: return set(orthologues), treeFN
    try:
        tree = tree_lib.Tree(treeFN)
    except:
        tree = tree_lib.Tree(treeFN, format=3)
#    if qPrune: tree.prune(tree.get_leaf_names())
    if len(tree) == 1: return set(orthologues), tree
    root = GetRoot(tree, species_tree_rooted, GeneToSpecies)
    if root == None: return set(orthologues), tree
    # Pick the first root for now
    if root != tree:
        tree.set_outgroup(root)

    tree = Resolve(tree, GeneToSpecies)
    if qPrune: tree.prune(tree.get_leaf_names())
    if len(tree) == 1: return set(orthologues), tree
    """ At this point need to label the tree nodes """
    iNode = 1
    tree.name = "n0"
    for n in tree.traverse():
        if n.is_leaf(): continue
        if not n.is_root():
            n.name = "n%d" % iNode
            iNode += 1
        ch = n.get_children()
        if len(ch) == 2: 
            o, sp0, sp1 = OverlapSize(n, GeneToSpecies)
            if o != 0:
                if dupsWriter != None:
                    sp_present = sp0.union(sp1)
                    if len(sp_present) == 1:
                        stNode = species_tree_rooted & next(sp for sp in sp_present)
                        isSTRIDE = "Terminal"
                    else:
                        stNode = species_tree_rooted.get_common_ancestor(sp_present)
                        isSTRIDE = "Shared" if all_stride_dup_genes == None else "STRIDE" if frozenset(n.get_leaf_names()) in all_stride_dup_genes else ""
                    dupsWriter.writerow(["OG%07d" % iog, spIDs[stNode.name] if len(stNode) == 1 else stNode.name, n.name, float(o)/(len(stNode)), isSTRIDE, ", ".join([seqIDs[g] for g in ch[0].get_leaf_names()]), ", ".join([seqIDs[g] for g in ch[1].get_leaf_names()])]) 
            else:
                d0 = defaultdict(list)
                for g in ch[0].get_leaf_names():
                    sp, seq = g.split("_")
                    d0[sp].append(seq)
                d1 = defaultdict(list)
                for g in ch[1].get_leaf_names():
                    sp, seq = g.split("_")
                    d1[sp].append(seq)
                orthologues.append((d0, d1))
        elif len(ch) > 2:
            species = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in ch]
            for (n0, s0), (n1, s1) in itertools.combinations(zip(ch, species), 2):
                if len(s0.intersection(s1)) == 0:
                    d0 = defaultdict(list)
                    for g in n0.get_leaf_names():
                        sp, seq = g.split("_")
                        d0[sp].append(seq)
                    d1 = defaultdict(list)
                    for g in n1.get_leaf_names():
                        sp, seq = g.split("_")
                        d1[sp].append(seq)
                    orthologues.append((d0, d1))
#    raise Exception("WriteQfO2")
    if qWrite:
        directory = os.path.split(treeFN)[0]
        WriteQfO2(orthologues, directory + "/../Orthologues_M3/" + os.path.split(treeFN)[1], qAppend=False)
    return orthologues, tree

def AppendOrthologuesToFiles(orthologues_alltrees, speciesDict, iSpeciesToUse, sequenceDict, resultsDir):
    # Sort the orthologues according to speices pairs
    sp_to_index = {str(sp):i for i, sp in enumerate(iSpeciesToUse)}
    nOrtho = util.nOrtho_sp(len(iSpeciesToUse))    
    species = speciesDict.keys()
#    left = [[] for sp in species]  
#    right = [[] for sp in species]
    # reorder orthologues on a per-species basis
    nSpecies = len(species)
    for i in xrange(nSpecies):
        sp0 = species[i]
        strsp0 = sp0 + "_"
        isp0 = sp_to_index[sp0]
        d0 = resultsDir + "Orthologues_" + speciesDict[sp0] + "/"
        for j in xrange(i, nSpecies):
            sp1 = species[j]
            if sp1 == sp0: continue
            strsp1 = sp1 + "_"
            isp1 = sp_to_index[sp1]
            d1 = resultsDir + "Orthologues_" + speciesDict[sp1] + "/"
            with open(d0 + '%s__v__%s.csv' % (speciesDict[sp0], speciesDict[sp1]), 'ab') as outfile1, open(d1 + '%s__v__%s.csv' % (speciesDict[sp1], speciesDict[sp0]), 'ab') as outfile2:
                writer1 = csv.writer(outfile1, delimiter="\t")
                writer2 = csv.writer(outfile2, delimiter="\t")
                for iog, ortholouges_onetree in orthologues_alltrees:                   
                    og = "OG%07d" % iog
                    for leavesL, leavesR in ortholouges_onetree:
                        nL0 = len(leavesL[sp0])
                        nR0 = len(leavesR[sp0])
                        nL1 = len(leavesL[sp1])
                        nR1 = len(leavesR[sp1])
                        if nL0*nR1 + nL1*nR0 == 0: continue # no orthologues
                        # each species can be in only one of L and R at most: they might both be in the same half
                        if nL0 > 0:
                            # then nR0 == 0 so nR1 > 0 since checked (nL0*nR1 + nL1*nR0 != 0)
                            n0 = nL0
                            n1 = nR1
                            text0 = ", ".join([sequenceDict[strsp0 + g] for g in leavesL[sp0]])
                            text1 = ", ".join([sequenceDict[strsp1 + g] for g in leavesR[sp1]])
                        else:
                            n0 = nR0
                            n1 = nL1
                            text0 = ", ".join([sequenceDict[strsp0 + g] for g in leavesR[sp0]])
                            text1 = ", ".join([sequenceDict[strsp1 + g] for g in leavesL[sp1]])
                        writer1.writerow((og, text0, text1))
                        writer2.writerow((og, text1, text0))
                        nOrtho.n[isp0, isp1] += n0
                        nOrtho.n[isp1, isp0] += n1
                        if n0 == 1 and n1 == 1:
                            nOrtho.n_121[isp0, isp1] += 1
                            nOrtho.n_121[isp1, isp0] += 1
                        elif n0 == 1:
                            nOrtho.n_12m[isp0, isp1] += 1
                            nOrtho.n_m21[isp1, isp0] += n1
                        elif n1 == 1:
                            nOrtho.n_m21[isp0, isp1] += n0
                            nOrtho.n_12m[isp1, isp0] += 1
                        else:
                            nOrtho.n_m2m[isp0, isp1] += n0
                            nOrtho.n_m2m[isp1, isp0] += n1
    return nOrtho   
                                      
def Resolve(tree, GeneToSpecies):
    StoreSpeciesSets(tree, GeneToSpecies)
    for n in tree.traverse("postorder"):
        tree = resolve.resolve(n, GeneToSpecies)
    return tree

def GetOrthologuesStandalone_Parallel(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree):
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    args_queue = mp.Queue()
    for treeFn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")): args_queue.put((0, treeFn, species_tree_rooted, GeneToSpecies, True))
    util.RunMethodParallel(GetOrthologues_for_tree, args_queue, 8)
    
def GetOrthologuesStandalone_Serial(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree):
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
#    args_queue = mp.Queue()
    for treeFn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")): 
        GetOrthologues_for_tree(0, treeFn, species_tree_rooted, GeneToSpecies, True)        
        
def DoOrthologuesForOrthoFinder(ogSet, treesIDsPatFn, species_tree_rooted_fn, GeneToSpecies, workingDir, output_dir, reconTreesRenamedDir, all_stride_dup_genes):    # Create directory structure
    speciesDict = ogSet.SpeciesDict()
    SequenceDict = ogSet.SequenceDict()
    # Write directory and file structure
    speciesIDs = ogSet.speciesToUse
    nspecies = len(speciesIDs)           
    for index1 in xrange(nspecies):
        d = output_dir + "Orthologues_" + speciesDict[str(speciesIDs[index1])] + "/"
        if not os.path.exists(d): os.mkdir(d)     
        for index2 in xrange(nspecies):
            if index2 == index1: continue
            with open(d + '%s__v__%s.csv' % (speciesDict[str(speciesIDs[index1])], speciesDict[str(speciesIDs[index2])]), 'wb') as outfile:
                writer1 = csv.writer(outfile, delimiter="\t")
                writer1.writerow(("Orthogroup", speciesDict[str(speciesIDs[index1])], speciesDict[str(speciesIDs[index2])]))
    # Infer orthologues and write them to file           
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    # Label nodes of species tree
    species_tree_rooted.name = "N0"    
    iNode = 1
    for n in species_tree_rooted.traverse():
        if (not n.is_leaf()) and (not n.is_root()):
            n.name = "N%d" % iNode
            iNode += 1
    nOgs = len(ogSet.OGs())
    nOrthologues_SpPair = util.nOrtho_sp(nspecies)
    with open(reconTreesRenamedDir + "../Duplications.csv", 'wb') as outfile:
        dupWriter = csv.writer(outfile, delimiter="\t")
        dupWriter.writerow(["Orthogroup", "Species Tree Node", "Gene Tree Node", "Support", "Type",	"Genes 1", "Genes 2"])
        for iog in xrange(nOgs):
            allOrthologues = []
            orthologues, recon_tree = GetOrthologues_for_tree(iog, treesIDsPatFn(iog), species_tree_rooted, GeneToSpecies, dupsWriter=dupWriter, seqIDs=ogSet.Spec_SeqDict(), spIDs=ogSet.SpeciesDict(), all_stride_dup_genes=all_stride_dup_genes)
            allOrthologues.append((iog, orthologues))
            util.RenameTreeTaxa(recon_tree, reconTreesRenamedDir + "OG%07d_tree.txt" % iog, ogSet.Spec_SeqDict(), qSupport=False, qFixNegatives=True, label='n') 
            if iog >= 0 and divmod(iog, 10 if nOgs <= 200 else 100 if nOgs <= 2000 else 1000)[1] == 0:
                util.PrintTime("Done %d of %d" % (iog, nOgs))
            nOrthologues_SpPair += AppendOrthologuesToFiles(allOrthologues, speciesDict, ogSet.speciesToUse, SequenceDict, output_dir)
    return nOrthologues_SpPair

def RootAllTrees():
    import tree
    speciesIDs = util.FirstWordExtractor("SpeciesIDs.txt").GetIDToNameDict()
    species_tree_rooted = tree.Tree("SpeciesTree_ids_0_rooted_unresolved.txt")
    GeneToSpecies = GeneToSpecies_dash
    for fn in glob.glob("Trees_ids/OG*txt"):
        print("*** " + fn + " ***")
        t = tree.Tree(fn)
        root = GetRoot(t, species_tree_rooted, GeneToSpecies)
        if root == None: 
            print("Fail: " + fn)
        else:
            if root != t:
                t.set_outgroup(root)
            for n in t:
                n.name = speciesIDs[n.name.split("_")[0]]
            t.write(outfile="Trees_ids_rooted/" + os.path.split(fn)[1])
    util.Success()

if __name__ == "__main__":
#    RootAllTrees()
    parser = argparse.ArgumentParser()
    parser.add_argument("trees_dir")
    parser.add_argument("rooted_species_tree")
#    parser.add_argument("-p", "--prune", action='store_true')
    parser.add_argument("-s", "--separator", choices=("dot", "dash", "second_dash", "3rd_dash", "hyphen"), help="Separator been species name and gene name in gene tree taxa")
    args = parser.parse_args()
    output_dir = os.path.split(args.trees_dir)[0]
    qSingleTree = False
    try:
        tree_lib.Tree(args.trees_dir)
        qSingleTree = True
        print("Analysing single tree")
    except:
        try:
            tree = tree_lib.Tree(args.trees_dir)
            qSingleTree = True
        except:
            pass
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    GeneToSpecies = GetGeneToSpeciesMap(args)
    output_dir = output_dir + "/../Orthologues_M3/"
    print(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
#    GetOrthologuesStandalone_Parallel(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree)
    GetOrthologuesStandalone_Serial(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree)

