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
import ete3 as ete
import glob
import argparse
import operator
from collections import Counter
 
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
  
#SpeciesName = OrthoFinderIDs
#SpeciesName = GeneToSpecies_3rdDash
#SpeciesName = GeneToSpecies_dash
  
class RootMap(object):
    def __init__(self, setA, setB, GeneToSpecies):
        self.setA = setA
        self.setB = setB
        self.GeneToSpecies = GeneToSpecies
        
    def GeneMap(self, gene_name):
        sp = self.GeneToSpecies(gene_name)
        if sp in self.setA: return True 
        elif sp in self.setB: return False
        else: raise Exception
        
def StoreSpeciesSets(t, GeneMap):
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node.add_feature('sp_down', {GeneMap(node.name)})
        elif node.is_root():
            continue
        else:
            node.add_feature('sp_down', set.union(*[ch.sp_down for ch in node.get_children()]))
    for node in t.traverse('preorder'):
        if node.is_root():
            node.add_feature('sp_up', set())
        else:
            parent = node.up
            if parent.is_root():
                others = [ch for ch in parent.get_children() if ch != node]
                node.add_feature('sp_up', set.union(*[other.sp_down for other in others]))
            else:
                others = [ch for ch in parent.get_children() if ch != node]
                sp_downs = set.union(*[other.sp_down for other in others])
                node.add_feature('sp_up', parent.sp_up.union(sp_downs))
    t.add_feature('sp_down', set.union(*[ch.sp_down for ch in t.get_children()]))

def GetRoots(tree, species_tree_rooted, GeneToSpecies):
    species = set([GeneToSpecies(g) for g in tree.get_leaf_names()])
    if len(species) == 1:
        return [], 0, None
    n1, n2 = species_tree_rooted.get_children()
    t1 = set(n1.get_leaf_names())
    t2 = set(n2.get_leaf_names())
    have1 = len(species.intersection(t1)) != 0
    have2 = len(species.intersection(t2)) != 0
#    print(tree)
    while not (have1 and have2):
        # Doesn't contain outgroup, step down in species tree until it does
        if have1:
            n = n1
#            print(1)
        else:
            n = n2
#            print(2)
#        print(n1)
#        print(n2)
        n1, n2 = n.get_children()
        t1 = n1.get_leaf_names()
        t2 = n2.get_leaf_names()
        have1 = len(species.intersection(t1)) != 0
        have2 = len(species.intersection(t2)) != 0
        
    n1, n2 = species_tree_rooted.get_children()
    root_mapper = RootMap(t1, t2, GeneToSpecies)    
    GeneMap = root_mapper.GeneMap
    StoreSpeciesSets(tree, GeneMap)
    found = set()
    TF = set([True, False])
    TFfr = frozenset([True, False])
    Tfr = frozenset([True])
    Ffr = frozenset([False])
    fail = 0
    for m in tree:
        n = m.up
#        print("")
#        print(m)
#        print(n)
#        print(n.sp_down)
        while not n.is_root() and n.sp_down != TF:
            m = n
            n = m.up
#            print(n)
#            print(n.sp_down)
        if n.sp_down == TF:
            children = n.get_children()
            if n.is_root():
                colour = m.sp_down
                if any([x.sp_down != colour and len(x.sp_down) == 1 for x in children]):
                    comb = Counter([frozenset(x.sp_down) for x in children])
                    # case 0
                    if comb[TFfr] == 0:
                        # case 0A - one of the branches is the root
                        for c in children:
                            if sum([c.sp_down == x.sp_down for x in children]) == 1:
                                found.add(c) # only holds for one of them
                                break
                    elif comb[TFfr] == 1 and (comb[Tfr] == 2 or comb[Ffr] == 2):
                        # case 0B - one mixed branch, two identical True/False branches
                        # we don't know this is the division, stepping down in the mixed branch might still be all same as the single state ones
                        # we'll find this division while walking up the tree
                        pass
                    elif comb[TFfr] == 1 and comb[Tfr] == 1:
                        # case 0C - one mixed branch, one True & one False
                        found.add([c for c in children if c.sp_down == TF][0])
                    else:
                        # case 0D - two mixed branches
                        # while find the transition while walking up the tree
                        pass 
#                    found.add(n)
#                    print("*** Root1 ***")
            elif len(children) == 2:
#                found.add(n)
                c1, c2 = children
                single_state = c1.sp_down if len(c1.sp_down) == 1 else c2.sp_down
                if len(c1.sp_down) == 1 and len(c2.sp_down) == 1:
                    # Case 1 - Both single state
                    if len(n.sp_up) == 1:
                        # Case 1A - 3rd clade also single state
                        # Root is the bipartition separating True from False
                        found.add(c1 if n.sp_up == c2.sp_down else c2)
                    else:
                        # Case 1B - 3rd clade is mixed
                        found.add(n)
#                    print("*** Root2 ***")
#                    print(c1.sp_down)
#                    print(c2.sp_down)
#                    print(n.sp_up)
#                    print("------")
                else:
                    # Case 2 - only one is single state and it's not the same as the 3rd clade
                    if single_state != n.sp_up:
                        # Case 2A - only one is single state and it's not the same as the 3rd clade
#                        print("*** Root3 ***")
                        found.add(c1 if len(c1.sp_down) == 1 else c2)
#                    else:
#                        # Case 2A - only one is single state and it's the same as the 3rd clade
#                        # root is in the mixed clade and will be found while walking up that
#                        pass
            else:
                fail += 1
#    for f in found:
#        print(f)
#        print(len(f.get_leaf_names()))
    return list(found), fail, GeneMap  

def WriteQfO(ortho_pairs, outfilename, qForQFO = True):
    with open(outfilename, 'ab' if qForQFO else 'wb') as outfile:
        for p in ortho_pairs:
            g1, g2 = list(p)
            if qForQFO:
                g1 = g1.split("_")[-1]
                g2 = g2.split("_")[-1]
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
  
def GeneToSpecies_dash(g):
  return g.split("_", 1)[0]
  
OrthoFinderIDs = GeneToSpecies_dash
    
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
        roots, j, GeneMap = GetRoots(tree, species_tree_rooted, GeneToSpecies)
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
        WriteQfO(orthologues2, outfn, False)

    
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
