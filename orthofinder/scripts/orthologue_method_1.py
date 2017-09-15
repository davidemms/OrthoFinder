# -*- coding: utf-8 -*-
"""
Created on Wed May 24 14:43:42 2017

@author: david

EggNOG-like approach

1 - root gene trees on outgroup
2 - infer orthologues

"""

import os
import ete3 as ete
import glob
import argparse
import itertools
#import random
#import cProfile
from collections import Counter
from Bio import Phylo

#def RootAtClade(t, accs_in_clade):
#    if len(accs_in_clade) == 1:
#        t.set_outgroup(list(accs_in_clade)[0])
#        return t
#    accs = set(t.get_leaf_names())
#    dummy = list(accs.difference(accs_in_clade))[0]
#    t.set_outgroup(dummy)
#    node = t.get_common_ancestor(accs_in_clade)
#    t.set_outgroup(node)
#    return t

c0 = 0  
c1 = 0
c2 = 0
  
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

def IsDup(node, node_identities, GeneToSpecies, qPrint = False):
    """
    Is node a duplication node
    Args:
        node - is this node a duplication node
        node_identities - distionary of the identities of nodes already calculated
        GeneToSpecies - Function taking genes names and returning species names
    """
#    global c0, c1,c2
    parent = None if node.is_root() else node.up
    if node in node_identities:
#        c0+=1
        parent_c, type_c = node_identities[node]
        if parent_c == parent:
#            c1+=1
            return type_c
    descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in node.get_children()]
    if qPrint: print(descendents)
    # *** What if non-binary? Should really compare the relevant two branches
    dup = len(descendents[0].intersection(descendents[1])) != 0
    node_identities[node] = (parent, dup)
#    c2+=1
    return dup

def _dev_PrintNumberOfRoots(trees_dir, species_tree_rooted_fn, output_dir, qSingleTree):
    species_tree_rooted = ete.Tree(species_tree_rooted_fn)
    for fn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")):
#        print("")
#        print(fn)
        try:
            tree = ete.Tree(fn)
        except:
            tree = ete.Tree(fn, format=3)
        if len(tree) == 1: continue
        roots, j, _ = GetRoots(tree, species_tree_rooted)
        width = len(fn)
        formatString = "%d %" + ("%ds" % (width+3))
        outputLine = formatString % (len(roots), fn)
        formatString2 = "%" + ("%ds" % (width+3+6))
        print(formatString2 % outputLine)
        if j != 0: print(j)
        print(roots)

def WriteQfO(ortho_pairs, outfilename, qForQFO = True):
    with open(outfilename, 'ab' if qForQFO else 'wb') as outfile:
        for p in ortho_pairs:
            g1, g2 = list(p)
            if qForQFO:
                g1 = g1.split("_")[-1]
                g2 = g2.split("_")[-1]
            outfile.write("%s\t%s\n" % (g1, g2))

def GetOrthologues_for_tree(tree, species_tree_rooted, GeneToSpecies, outfn, qPrune):
#        print(tree)
#        tree_phylo = Phylo.read(fn, 'newick')
        if qPrune: tree.prune(tree.get_leaf_names())
        if len(tree) == 1: return
        roots, j, GeneMap = GetRoots(tree, species_tree_rooted, GeneToSpecies)
        print("%d roots" % len(roots))
#        print("Root:")
#        print(roots[0])
        # 1. Store distance of each root from each leaf
        orthologues2 = []
        if len(roots) > 0:
            leaves = tree.get_leaves()
            for l in leaves: 
                l.add_feature('root_dists', [l.get_distance(r) for r in roots])
        
            # 2. For each pair of leaves, find the root that maximisest eh closest of the two distances
            pairs = [frozenset(p) for p in itertools.combinations(leaves, 2)]
            chosen = []
            for l1, l2 in pairs:
                dists = [min(d1,d2) for d1,d2 in zip(l1.root_dists, l2.root_dists)]
                # What if multiple roots give the same value? Actually very unlikely, it's a continuous value
                iRoot = dists.index(max(dists))
    #            print((roots[iRoot].name, l1.name, l2.name))
                chosen.append(iRoot)
            
            # node_identities,  node -> (parent, specdup). If the cached parent node is the same as this parent node then don't need to recalculate identity
            node_identities = dict() 
    #        orthologues1 = []
            
    #        # Approach 1 - for each root:
    #        # find pairs which use this root
    #        # find their divergence node, decide if spciation or duplication
    #        results = []
    #        for i, root in enumerate(roots):
    #            tree.set_outgroup(root)
    #            for iThis, (l1, l2) in zip(chosen, pairs):
    #                if iThis != i: continue
    #                # 3. Identify divergence node for the pair of leaves
    ##                    GetDivergenceNode(l1.name, l2.name, roots[iRoot].name, tree_phylo)
    #                n = l1.get_common_ancestor([l2])
    #                
    #                # 4. If node is duplication => paralogue, else orthologue
    #                dup = IsDup(n, node_identities, GeneMap)
    #                results.append(({l1.name, l2.name}, n.name if n != tree else root.name, not dup))
    #                if not dup: orthologues1.append(frozenset((l1.name, l2.name)))
    #        orthologues1 = set(orthologues1)
            
            # Approach 2 - for each root
            # Walk through tree
            # for each node, look through all pairs for which this is the divergence
            root_dict = {p:iRoot for p,iRoot in zip(pairs, chosen)}
            results = {p:None for p in pairs}
            for i, root in enumerate(roots):
    #            print(i)
                if root != tree:
                    tree.set_outgroup(root)
                for iThis, p in zip(chosen, pairs):
                    if iThis != i: continue
                    if results[p] != None: continue
                    l1, l2 = p
                    # 3. Identify divergence node for the pair of leaves
    #                    GetDivergenceNode(l1.name, l2.name, roots[iRoot].name, tree_phylo)
                    n = l1.get_common_ancestor([l2])
                    
                    # 4. If node is duplication => paralogue, else orthologue
                    dup = IsDup(n, node_identities, GeneToSpecies, qPrint=False)
                    if not dup: orthologues2.append(frozenset((l1.name, l2.name)))
                    results[p] = not dup
                        
                    # Now, do all pairs below this node which use the corresponding root
                    children = n.get_children()
                    if len(children) != 2: continue
                    d1 = children[0].get_leaves()
                    d2 = children[1].get_leaves()
                    for l1 in d1:
                        for l2 in d2:
                            this_pair = frozenset([l1,l2])
                            if root_dict[this_pair] == i: 
                                results[this_pair] = not dup
#                                print((n.get_leaf_names(), this_pair, not dup))
#                                print((this_pair, not dup))
                                o = frozenset((l1.name, l2.name))
                                if not dup: orthologues2.append(o)
    #                            if o not in orthologues1:
    #                                pass
    #                                print((o, n.name, root.name))
            orthologues2 = set(orthologues2)
#        WriteQfO(orthologues2, "/home/david/projects/OrthoFinder/benchmarks/QuestForOrthologues/Dev/1_HuertaCepas/X_Euk_Fungi/" + os.path.split(fn)[1] + ".txt")
#        WriteQfO(orthologues2, "/home/david/projects/OrthoFinder/benchmarks/Orthologues-TreeInterpretaion/Fungi/Sample_0.30_Jun01/M_1_" + os.path.split(fn)[1] + ".txt")
        print("%d orthologues" % len(orthologues2))
        WriteQfO(orthologues2, outfn, False)
#        for _ in xrange(100):
#            i = random.randint(0, len(results))
#            print(results[i])
#        print(Counter(chosen).most_common())
#        print([r for r in results if {'16_13022', '17_20334'} == r[0]])
#        print("Identified %d orthologues" % len(orthologues2))
#        print("Identified %d orthologues (original method)" % len(orthologues1))
#        print("%d pairs in overlap" % len(orthologues1.intersection(orthologues2)))
#        print(set(orthologues2).symmetric_difference(orthologues1))
#        nSpecies = len(set([GeneToSpecies(l) for l in tree.get_leaf_names()]))
#        print("%f orthologues per gene per species" % (2*float(len(orthologues2)) / len(tree) / nSpecies))
#        return orthologues1, orthologues2
        
#        print("Identified %d orthologues" % len(orthologues2))
#        nSpecies = len(set([GeneToSpecies(l) for l in tree.get_leaf_names()]))
#        print("%f orthologues per gene per species" % (2*float(len(orthologues2)) / len(tree) / nSpecies))
    
def GetOrthologues(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree, qPrune=False):
#    outfn = output_dir + "Orthologues_M1.txt"
#    print(outfn)
#    return
#    if os.path.exists(outfn): return
    species_tree_rooted = ete.Tree(species_tree_rooted_fn)
    for fn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")):
        outfn = output_dir + os.path.split(fn)[1]
#        print("")
#        print(fn)
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
#    parser.add_argument("-t", "--tree", action='store_true')
    args = parser.parse_args()
#    print(args.tree)
    output_dir = os.path.split(args.trees_dir)[0]
    qSingleTree = False
    try:
        ete.Tree(args.trees_dir)
        qSingleTree = True
#        print("Analysing single tree")
    except:
        try:
            tree = ete.Tree(args.trees_dir, format=3)
            qSingleTree = True
        except:
            pass
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
#    _dev_PrintNumberOfRoots(args.trees_dir, args.rooted_species_tree, output_dir, qSingleTree)
    
    GeneToSpecies = GetGeneToSpeciesMap(args)
    output_dir = output_dir + "/../Orthologues_M1/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    GetOrthologues(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree, qPrune=args.prune)
#    print((c0, c1,c2))
