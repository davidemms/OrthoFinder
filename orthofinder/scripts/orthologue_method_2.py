# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:11:11 2017

@author: david

Perform directed 'reconciliation' first and then apply EggNOG method

1 - root gene trees on outgroup: unique one this time
2 - infer orthologues
"""
import os
import glob
import ete2 as ete
import argparse

import rearrange
import orthologue_method_1c as om1

def Reconcile(tree, GeneToSpecies):
    om1.om1.StoreSpeciesSets(tree, GeneToSpecies)
    for n in tree.traverse("postorder"):
#        print("Next")
#        print("")        
#        print(n)
#        if 'recon' in n.features: 
#            print("Now marked")
#            for nn in n.recon:
#                print(nn)
        rearrange.rearrange(n, GeneToSpecies)
#        if n.up != None and 'recon' in n.up.features: 
#            print("Marked")
#            for nn in n.up.recon:
#                print(nn)
#        for nn in tree.traverse():
#            if 'sp_down' not in nn.features:
#                print("No sp_down")
#                print(tree)
#                print(n)
#                raise Exception
#        print(tree)
    # Do a second iteration. Does it improve the results
#    om1.StoreSpeciesSets(tree, GeneToSpecies)
#    for n in tree.traverse("postorder"):
#        rearrange.rearrange(n, GeneToSpecies)

def GetOrthologues(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree, qPrune=False):
    species_tree_rooted = ete.Tree(species_tree_rooted_fn)
    for fn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")):
        print("")
        print(fn)
        outfn = output_dir + os.path.split(fn)[1]
        if (not os.path.exists(fn)) or os.stat(fn).st_size == 0: continue
        try:
            tree = ete.Tree(fn)
        except:
            tree = ete.Tree(fn, format=3)
#        print(tree)
#        tree_phylo = Phylo.read(fn, 'newick')
        if qPrune: tree.prune(tree.get_leaf_names())
        if len(tree) == 1: continue
#        print(tree)
        root = om1.GetRoot(tree, species_tree_rooted, GeneToSpecies)
#        print(root)
        if root == None: continue
        # Pick the first root for now
        if root != tree:
            tree.set_outgroup(root)

#        print(tree.write(format=9))
        Reconcile(tree, GeneToSpecies)
#        print(tree)
#        print(tree.write())
        om1.GetOrthologues_for_tree(tree, species_tree_rooted, GeneToSpecies, outfn, qPrune=False, qRoot=False)


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
        print("Analysing single tree")
    except:
        try:
            tree = ete.Tree(args.trees_dir)
            qSingleTree = True
        except:
            pass
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
#    _dev_PrintNumberOfRoots(args.trees_dir, args.rooted_species_tree, output_dir, qSingleTree)
    
    GeneToSpecies = om1.GetGeneToSpeciesMap(args)
    output_dir = output_dir + "/../Orthologues_M2c/"
    print(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    GetOrthologues(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree, qPrune=args.prune)