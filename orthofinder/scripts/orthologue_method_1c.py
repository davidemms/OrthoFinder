# -*- coding: utf-8 -*-
"""
Created on Wed May 24 14:43:42 2017

@author: david

EggNOG-like approach

1 - root gene trees on outgroup
2 - infer orthologues

Version 1c: 
- single root
- make allowances for misplaced genes causing an overlap

"""

import os
import sys
import time
import ete3 as ete
import glob
import argparse
import itertools
import operator
#import random
#import cProfile
from collections import Counter
from Bio import Phylo

import orthologue_method_1 as om1
  
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
        
def IsDup(node, node_identities, GeneToSpecies, qPrint = False):
    """
    Is node a duplication node
    Args:
        node - is this node a duplication node
        node_identities - dictionary of the identities of nodes already calculated
        GeneToSpecies - Function taking genes names and returning species names
    """
    parent = None if node.is_root() else node.up
    if node in node_identities:
        parent_c, type_c = node_identities[node]
        if parent_c == parent:
            return type_c
    descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in node.get_children()]
    if qPrint: print(descendents)
    # *** What if non-binary? Should really compare the relevant two branches
    dup = len(descendents[0].intersection(descendents[1])) != 0
    node_identities[node] = (parent, dup)
    return dup
    
def IsDup_Simple(node, GeneToSpecies):
    descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in node.get_children()]
    try:
        return len(descendents[0].intersection(descendents[1])) != 0
    except:
        print(node)
        raise
            
    
def IsDup_Overlaps(node, GeneToSpecies):
    """
    Determine whether a node is a duplication based on the overlapping species. Also return genes that should be ignored if it is a speciation
    Args:
        node - the node to detmine duplication/speciation identity of
        GeneToSpecies - function taking a gene name and returning the species name
    Returns:
        isDup - boolean
        ignore - list of genes to be ignored (genes are from the species that overlap and the genes appear to be incorrectly placed
            and so orthology shouldn't be inferred for them across the node)
    """
    """
    Is node a duplication node
    Args:
        node - is this node a duplication node
        node_identities - dictionary of the identities of nodes already calculated
        GeneToSpecies - Function taking genes names and returning species names
        
    6861133 v 7084653 or 6856703
    """
    NULL = set()
    ch = node.get_children()
    if len(ch) < 2:
        return True, NULL, NULL
    elif len(ch) >2:
        l = map(len, ch)
        x = sorted(zip(l,ch), reverse=True)
        ch = [x[0][1], x[1][1]]
    g0 = ch[0].get_leaf_names()
    g1 = ch[1].get_leaf_names()
    s0 = map(GeneToSpecies, g0)
    s1 = map(GeneToSpecies, g1)
    d0 = set(s0)
    d1 = set(s1)
#    descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in node.get_children()] # this could be cached or iteratively determined
    # *** What if non-binary? Should really compare the relevant two branches
    overlap = d0.intersection(d1)
    ns0 = len(d0)
    ns1 = len(d1)
    o = len(overlap)
    genes_overlap0 = len([True for s in s0 if s in overlap])
    genes_overlap1 = len([True for s in s1 if s in overlap])
    if o > 0 and (genes_overlap0 != len(g0) and genes_overlap1 != len(g1)):
        sys.stderr.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (o, ns0, ns1, genes_overlap0, genes_overlap1, len(g0), len(g1)))
#        if (0 < genes_overlap0 < 3 or 0 < genes_overlap1 < 3) and len(g0)*len(g1) > 5 and (genes_overlap0 != len(g0) and genes_overlap1 != len(g1)):
#            sys.stderr.write(node.get_ascii(compact=False, show_internal=False) + "\n")
    if o == 0:
        return False, NULL, NULL
#    else:
#        return True, None, None
    # Can't have more than 1/3 of either species sets in the overlap w/o it being a duplication
    if 4*o > ns0 or 4*o > ns1:
        return True, None, None
    # How many genes would have to be eliminated compared to the size of the clade? Can we eliminate fewer than a tenth?
    e0 = set([g for g,s in zip(g0,s0) if s in overlap])
    e1 = set([g for g,s in zip(g0,s0) if s in overlap])
#    print(("Elim", len(e0), len(e1)))
    qElim0 = 4*len(e0) <= len(g0)
    qElim1 = 4*len(e1) <= len(g1)
    if qElim0 and qElim1:
        print(1)
        return False, e0 if len(e0) <= len(e1) else NULL, NULL if len(e0) <= len(e1) else e1
    elif qElim0:
        print(1)
        return False, e0, NULL
    elif qElim1:
        print(1)
        return False, NULL, e1
    else: 
        return True, None, None

def GetRoot(tree, species_tree_rooted, GeneToSpecies):
        roots, j, GeneMap = om1.GetRoots(tree, species_tree_rooted, GeneToSpecies)
        print("%d roots" % len(roots))
        # 1. Store distance of each root from each leaf
#        for root in roots:
#            print(root)
#        return roots[0]
        if len(roots) > 0:
            root_dists = [r.get_closest_leaf()[1] for r in roots]
            i, _ = max(enumerate(root_dists), key=operator.itemgetter(1))
            return roots[i]
    
def GetOrthologues_for_tree(tree, species_tree_rooted, GeneToSpecies, outfn, qPrune, qRoot=True):
        orthologues2 = []
        if qPrune: tree.prune(tree.get_leaf_names())
        if len(tree) == 1: return
        if qRoot:
            root = GetRoot(tree, species_tree_rooted, GeneToSpecies)
    #        if root == None: 
    #            print("No root")
    #            return
            # Walk through tree
            # for each node, look through all pairs for which this is the divergence
            if root != tree:
                print(root)
                tree.set_outgroup(root)
#        ch = tree.get_children()
#        print(ch[0] if len(ch[0]) < len(ch[1]) else ch[1])
#        print(tree.write())
        for n in tree.traverse('postorder'):
            ch = n.get_children()
            if len(ch) != 2: continue
            dup, elim0, elim1 = IsDup_Overlaps(n, GeneToSpecies)
            if not dup:
                orthologues2 += [(l0,l1) for l0 in ch[0].get_leaf_names() if l0 not in elim0 for l1 in ch[1].get_leaf_names() if l1 not in elim1]
#            dup = IsDup_Simple(n, GeneToSpecies)
#            if not dup:
#                orthologues2 += [(l0,l1) for l0 in ch[0].get_leaf_names() for l1 in ch[1].get_leaf_names()]
        orthologues2 = set(orthologues2)
        print("%d orthologues" % len(orthologues2))
        om1.WriteQfO(orthologues2, outfn, False)

    
def GetOrthologues(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree, qPrune=False):
    species_tree_rooted = ete.Tree(species_tree_rooted_fn)
    for fn in reversed(glob.glob(trees_dir + ("*" if qSingleTree else "/*"))):
        outfn = output_dir + os.path.split(fn)[1]
        if (not os.path.exists(fn)) or os.stat(fn).st_size == 0: continue
        try:
            tree = ete.Tree(fn)
        except:
            tree = ete.Tree(fn, format=3)
        GetOrthologues_for_tree(tree, species_tree_rooted, GeneToSpecies, outfn, qPrune)

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
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("trees_dir")
    parser.add_argument("rooted_species_tree")
    parser.add_argument("-p", "--prune", action='store_true')
    parser.add_argument("-s", "--separator", choices=("dot", "dash", "second_dash", "3rd_dash", "hyphen"), help="Separator been species name and gene name in gene tree taxa")
    args = parser.parse_args()
    output_dir = os.path.split(args.trees_dir)[0]
    qSingleTree = False
    try:
        ete.Tree(args.trees_dir)
        qSingleTree = True
    except:
        try:
            tree = ete.Tree(args.trees_dir, format=3)
            qSingleTree = True
        except:
            pass
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    GeneToSpecies = GetGeneToSpeciesMap(args)
    output_dir = output_dir + "/../Orthologues_M1c/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    GetOrthologues(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree, qPrune=args.prune)
