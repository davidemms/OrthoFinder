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
                        scores_list.append(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2))
                    elif clades.count(TF) >= 2:  
                        # (A,A,A)-excluded, (A,A,AB)-ignore as want A to be bigest without including B, (A,AB,AB), (AB,AB,AB) 
                        i = 0
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2))
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
        for gs1, gs2, _, _ in orthologues_list_pairs_list:
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
    intersection = descendents[0].intersection(descendents[1])
    return len(intersection), intersection, descendents[0], descendents[1]

nSuccess = 0
nFail = 0
nOrthoNew = 0

def ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours, relOverlapCutoff=4):
    """
    Is an overlap suspicious and if so can ift be resolved by identifying genes that are out of place?
    Args:
        overlap - the species with genes in both clades
        sp0 - the species below ch[0]
        sp1 - the species below ch[1]
        ch - the two child nodes
        tree - the gene tree
        neighbours - dictionary species->neighbours, where neighbours is a list of the sets of species observed at successive topological distances from the species
    Returns:
        qSuccess - has the overlap been resolved
        genes_removed - the out-of-place genes that have been removed so as to resolve the overlap
    
    Implementation:
        - The number of species in the overlap must be a 5th or less of the number of species in each clade
        - for each species with genes in both clades: the genes in one clade must all be more out of place (according to the 
          species tree) than all the gene from that species in the other tree
    """
    global nSuccess
    global nFail
    global nOrthoNew
    oSize = len(overlap)
    if relOverlapCutoff*oSize >= len(sp0) and relOverlapCutoff*oSize >= len(sp1): return False, []
    # The overlpa looks suspect, misplaced genes?
    # for each species, we'd need to be able to determine that all genes from A or all genes from B are misplaced
    genes_removed = []
    nA_removed = 0
    nB_removed = 0
    qResolved = True
    for sp in overlap:
        A = [g for g in ch[0].get_leaf_names() if g.split("_")[0] in overlap]
        B = [g for g in ch[1].get_leaf_names() if g.split("_")[0] in overlap]
#        print("Overlap: " + str((A,B)))
        sp_set = {sp}
        A_levels = []
        B_levels = []
        for X, level in zip((A,B),(A_levels, B_levels)):
            for g in X:
                r = (tree & g).up
                nextSpecies = set([gg.split("_")[0] for gg in r.get_leaf_names()])
                while len(nextSpecies) == 1:
                    r = r.up
                    nextSpecies = set([gg.split("_")[0] for gg in r.get_leaf_names()])
                nextSpecies = nextSpecies.difference(sp_set)
                # get the level
                observed = [bool(neigh & nextSpecies) for neigh in neighbours[sp]]
                N = len(observed)
#                                level.append((observed.index(True), next(i for i, oo in zip(xrange(N-1,-1,-1), reversed(observed)) if oo)))
                # the sum of the closest and furthest expected distance toplological distance for the closest genes in the gene tree (based on species tree topology)
                level.append(observed.index(True) + next(i for i, oo in zip(xrange(N-1,-1,-1), reversed(observed)) if oo))  
#        print(sp)
#        print(A)
#        print(A_levels)
#        print(B)
#        print(B_levels)
#        print("")
        qRemoveA = max(B_levels) < min(A_levels)                           
        qRemoveB = max(A_levels) < min(B_levels)                            
        if qRemoveA and relOverlapCutoff*oSize < len(sp0):
            nA_removed += len(A_levels)
            genes_removed.extend(A)
        elif qRemoveB and relOverlapCutoff*oSize < len(sp1):
            nB_removed += len(B_levels)
            genes_removed.extend(B)
        else:
            qResolved = False
            break
    if qResolved:
        nSuccess +=1
        nOrthoNew += (len(ch[0]) - nA_removed)  * (len(ch[1])-nB_removed)
#        print((nSuccess, nFail, nOrthoNew))
        return True, set(genes_removed)
    else:
        nFail += 1
        return False, set()
          
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

def GetOrthologues_for_tree(iog, treeFN, species_tree_rooted, GeneToSpecies, neighbours, qWrite=False, dupsWriter=None, seqIDs=None, spIDs=None, all_stride_dup_genes=None):
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
    if root == None: return set(orthologues), tree, set()
    # Pick the first root for now
    if root != tree:
        tree.set_outgroup(root)

    tree = Resolve(tree, GeneToSpecies)
    if qPrune: tree.prune(tree.get_leaf_names())
#    tree.write(outfile="/home/david/projects/OrthoFinder/Development/Orthologues/ReplacingDLCpar/Overlaps/Orthologues_M3/" + os.path.split(treeFN)[1] + ".rec.txt")
    if len(tree) == 1: return set(orthologues), tree, set()
    """ At this point need to label the tree nodes """
    iNode = 1
    tree.name = "n0"
    suspect_genes = set()
    empty_set = set()
    # preorder traverse so that suspect genes can be identified first, before their closer ortholgoues are proposed
    for n in tree.traverse('preorder'):
        if n.is_leaf(): continue
        if not n.is_root():
            n.name = "n%d" % iNode
            iNode += 1
        ch = n.get_children()
        if len(ch) == 2: 
            oSize, overlap, sp0, sp1 = OverlapSize(n, GeneToSpecies)
            if oSize != 0:
                qResolved, misplaced_genes = ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours)
            else:
                misplaced_genes = empty_set
            if oSize != 0 and not qResolved:
                if dupsWriter != None:
                    sp_present = sp0.union(sp1)
                    if len(sp_present) == 1:
                        stNode = species_tree_rooted & next(sp for sp in sp_present)
                        isSTRIDE = "Terminal"
                    else:
                        stNode = species_tree_rooted.get_common_ancestor(sp_present)
                        isSTRIDE = "Shared" if all_stride_dup_genes == None else "STRIDE" if frozenset(n.get_leaf_names()) in all_stride_dup_genes else ""
                    dupsWriter.writerow(["OG%07d" % iog, spIDs[stNode.name] if len(stNode) == 1 else stNode.name, n.name, float(oSize)/(len(stNode)), isSTRIDE, ", ".join([seqIDs[g] for g in ch[0].get_leaf_names()]), ", ".join([seqIDs[g] for g in ch[1].get_leaf_names()])]) 
            else:
                # sort out bad genes - no orthology for all the misplaced genes at this level (misplaced_genes). 
                # For previous levels, (suspect_genes) have their orthologues written to suspect orthologues file
                d0 = defaultdict(list)
                d0_sus = defaultdict(list)
                for g in [g for g in ch[0].get_leaf_names() if g not in misplaced_genes]:
                    sp, seq = g.split("_")
                    if g in suspect_genes:
                        d0_sus[sp].append(seq)
                    else:
                        d0[sp].append(seq)
#                if len(d0_sus) > 0: print(d0_sus)
                d1 = defaultdict(list)
                d1_sus = defaultdict(list)
                for g in [g for g in ch[1].get_leaf_names() if g not in misplaced_genes]:
                    sp, seq = g.split("_")
                    if g in suspect_genes:
                        d1_sus[sp].append(seq)
                    else:
                        d1[sp].append(seq)
                orthologues.append((d0, d1, d0_sus, d1_sus))
                suspect_genes.update(misplaced_genes)
        elif len(ch) > 2:
            species = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in ch]
            for (n0, s0), (n1, s1) in itertools.combinations(zip(ch, species), 2):
                if len(s0.intersection(s1)) == 0:
                    d0 = defaultdict(list)
                    d0_sus = defaultdict(list)
                    for g in n0.get_leaf_names():
                        sp, seq = g.split("_")
                        if g in suspect_genes:
                            d0_sus[sp].append(seq)
                        else:
                            d0[sp].append(seq)
                    d1 = defaultdict(list)
                    d1_sus = defaultdict(list)
                    for g in n1.get_leaf_names():
                        sp, seq = g.split("_")
                        if g in suspect_genes:
                            d1_sus[sp].append(seq)
                        else:
                            d1[sp].append(seq)
                    orthologues.append((d0, d1, d0_sus, d1_sus))
#    raise Exception("WriteQfO2")
    if qWrite:
        directory = os.path.split(treeFN)[0]
        WriteQfO2(orthologues, directory + "/../Orthologues_M3/" + os.path.split(treeFN)[1], qAppend=False)
    return orthologues, tree, suspect_genes

def AppendOrthologuesToFiles(orthologues_alltrees, speciesDict, iSpeciesToUse, sequenceDict, resultsDir, qContainsSuspectOlogs):
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
        if qContainsSuspectOlogs: 
            outfile1_sus = open(resultsDir + "Suspect_Orthologues/%s.csv" % speciesDict[sp0], 'ab')
            writer1_sus = csv.writer(outfile1_sus, delimiter="\t")
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
                if qContainsSuspectOlogs:
                    outfile2_sus = open(resultsDir + "Suspect_Orthologues/%s.csv" % speciesDict[sp1], 'ab')
                    writer2_sus = csv.writer(outfile2_sus, delimiter="\t")
                writer1 = csv.writer(outfile1, delimiter="\t")
                writer2 = csv.writer(outfile2, delimiter="\t")
                for iog, ortholouges_onetree in orthologues_alltrees:                   
                    og = "OG%07d" % iog
                    for leavesL, leavesR, leavesL_sus, leavesR_sus  in ortholouges_onetree:
                        # suspect_genes are the genes which, for this level, the orthologues should be considered suspect as the gene appears misplaced (at this level)
                        nL0 = len(leavesL[sp0])
                        nR0 = len(leavesR[sp0])
                        nL1 = len(leavesL[sp1])
                        nR1 = len(leavesR[sp1])
                        if nL0*nR1 + nL1*nR0 != 0: 
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
                        # Write suspect orthologues
                        if not qContainsSuspectOlogs: continue
                        nL0s = len(leavesL_sus[sp0])
                        nR0s = len(leavesR_sus[sp0])
                        nL1s = len(leavesL_sus[sp1])
                        nR1s = len(leavesR_sus[sp1])
                        if nL0s*(nR1+nR1s) + (nL1+nL1s)*nR0s != 0: 
                            # each species can be in only one of L and R at most: they might both be in the same half
                            if nL0s > 0:
                                # then nR0 == 0 so nR1 > 0 since checked (nL0*nR1 + nL1*nR0 != 0)
                                text0 = ", ".join([sequenceDict[strsp0 + g] for g in leavesL_sus[sp0]])
                                text1 = ", ".join([sequenceDict[strsp1 + g] for g in leavesR[sp1]+leavesR_sus[sp1]])
                            else:
                                text0 = ", ".join([sequenceDict[strsp0 + g] for g in leavesR_sus[sp0]])
                                text1 = ", ".join([sequenceDict[strsp1 + g] for g in leavesL[sp1]+leavesL_sus[sp1]])
                            writer1_sus.writerow((og, text0, text1))
                            writer2_sus.writerow((og, text1, text0))
                        outfile2_sus.close()
        if qContainsSuspectOlogs:
            outfile1_sus.close()
    return nOrtho   
                                      
def Resolve(tree, GeneToSpecies):
    StoreSpeciesSets(tree, GeneToSpecies)
    for n in tree.traverse("postorder"):
        tree = resolve.resolve(n, GeneToSpecies)
    return tree

def GetSpeciesNeighbours(t):
    """
    Args: t = rooted species tree
    """
    species = t.get_leaf_names()
    levels = {s:[] for s in species}
    for n in t.traverse('postorder'):
        if n.is_leaf(): continue
        children = n.get_children()
        leaf_sets = [set(ch.get_leaf_names()) for ch in children]
        not_i = [set.union(*[l for j, l in enumerate(leaf_sets) if j != i]) for i in xrange(len(children))]
        for l,n in zip(leaf_sets, not_i):
            for ll in l:
                levels[ll].append(n)
    return levels

def GetOrthologuesStandalone_Parallel(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree):
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    neighbours = GetSpeciesNeighbours(species_tree_rooted)
    args_queue = mp.Queue()
    for treeFn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")): args_queue.put((0, treeFn, species_tree_rooted, GeneToSpecies, neighbours, True))
    util.RunMethodParallel(GetOrthologues_for_tree, args_queue, 8)

def RootTreeStandalone_Serial(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree):
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
#    args_queue = mp.Queue()
    for treeFN in glob.glob(trees_dir + ("*" if qSingleTree else "/*")): 
        if (not os.path.exists(treeFN)) or os.stat(treeFN).st_size == 0: return 
        try:
            tree = tree_lib.Tree(treeFN)
        except:
            tree = tree_lib.Tree(treeFN, format=3)
        if len(tree) == 1: return 
        root = GetRoot(tree, species_tree_rooted, GeneToSpecies)
        if root == None: return 
        # Pick the first root for now
        if root != tree:
            tree.set_outgroup(root)
        tree.write(outfile=treeFN + ".rooted.txt")
    
def GetOrthologuesStandalone_Serial(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree):
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    neighbours = GetSpeciesNeighbours(species_tree_rooted)
#    args_queue = mp.Queue()
    for treeFn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")): 
        print(treeFn)
        GetOrthologues_for_tree(0, treeFn, species_tree_rooted, GeneToSpecies, neighbours, True)        
        
def DoOrthologuesForOrthoFinder(ogSet, treesIDsPatFn, species_tree_rooted_fn, GeneToSpecies, workingDir, output_dir, reconTreesRenamedDir, all_stride_dup_genes):    # Create directory structure
    speciesDict = ogSet.SpeciesDict()
    SequenceDict = ogSet.SequenceDict()
    # Write directory and file structure
    speciesIDs = ogSet.speciesToUse
    nspecies = len(speciesIDs)      
    dSuspect = output_dir + "Suspect_Orthologues/"
    if not os.path.exists(dSuspect): os.mkdir(dSuspect)     
    for index1 in xrange(nspecies):
        with open(dSuspect + '%s.csv' % speciesDict[str(speciesIDs[index1])], 'wb') as outfile:
            writer1 = csv.writer(outfile, delimiter="\t")
            writer1.writerow(("Orthogroup", speciesDict[str(speciesIDs[index1])], "Other"))
        d = output_dir + "Orthologues_" + speciesDict[str(speciesIDs[index1])] + "/"
        if not os.path.exists(d): os.mkdir(d)     
        for index2 in xrange(nspecies):
            if index2 == index1: continue
            with open(d + '%s__v__%s.csv' % (speciesDict[str(speciesIDs[index1])], speciesDict[str(speciesIDs[index2])]), 'wb') as outfile:
                writer1 = csv.writer(outfile, delimiter="\t")
                writer1.writerow(("Orthogroup", speciesDict[str(speciesIDs[index1])], speciesDict[str(speciesIDs[index2])]))
    # Infer orthologues and write them to file           
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    neighbours = GetSpeciesNeighbours(species_tree_rooted)
    # Label nodes of species tree
    species_tree_rooted.name = "N0"    
    iNode = 1
    for n in species_tree_rooted.traverse():
        if (not n.is_leaf()) and (not n.is_root()):
            n.name = "N%d" % iNode
            iNode += 1
    nOgs = len(ogSet.OGs())
    nOrthologues_SpPair = util.nOrtho_sp(nspecies) 
    species = speciesDict.keys()
    with open(reconTreesRenamedDir + "../Duplications.csv", 'wb') as outfile:
        dupWriter = csv.writer(outfile, delimiter="\t")
        dupWriter.writerow(["Orthogroup", "Species Tree Node", "Gene Tree Node", "Support", "Type",	"Genes 1", "Genes 2"])
        for iog in xrange(nOgs):
            orthologues, recon_tree, suspect_genes = GetOrthologues_for_tree(iog, treesIDsPatFn(iog), species_tree_rooted, GeneToSpecies, neighbours, dupsWriter=dupWriter, seqIDs=ogSet.Spec_SeqDict(), spIDs=ogSet.SpeciesDict(), all_stride_dup_genes=all_stride_dup_genes)
            for index0 in xrange(nspecies):
                strsp0 = species[index0]
                strsp0_ = strsp0+"_"
                these_genes = [g for g in suspect_genes if g.startswith(strsp0_)]
                if len(these_genes) > 0:
                    with open(output_dir + "Orthologues_" + speciesDict[strsp0] + "/MisplacedGenes.txt", 'ab') as outfile:
                        outfile.write("\n".join([SequenceDict[g]]))
            allOrthologues = [(iog, orthologues)]
            util.RenameTreeTaxa(recon_tree, reconTreesRenamedDir + "OG%07d_tree.txt" % iog, ogSet.Spec_SeqDict(), qSupport=False, qFixNegatives=True, label='n') 
            if iog >= 0 and divmod(iog, 10 if nOgs <= 200 else 100 if nOgs <= 2000 else 1000)[1] == 0:
                util.PrintTime("Done %d of %d" % (iog, nOgs))
            nOrthologues_SpPair += AppendOrthologuesToFiles(allOrthologues, speciesDict, ogSet.speciesToUse, SequenceDict, output_dir, True)
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
#    RootTreeStandalone_Serial(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree)
    util.Success()

