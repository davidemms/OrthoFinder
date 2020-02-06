#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts_of import tree, util

def AveDist(node):
    return np.average([node.get_distance(l) for l in node.get_leaf_names()]) 

def Fail():
    print("ERROR: An error occurred, please review error messages for more information.")
    sys.exit(1)

def CheckTree(t):
    if not 2 == len(t.get_children()):
        print("Input tree must be rooted")
        Fail()

def main():
    parser = argparse.ArgumentParser(description="Modify branch lengths on a rooted tree so that it is ultrametric")
    parser.add_argument("tree_fn", help="File containing a rooted tree in newick format")
    parser.add_argument("-r", "--root_age", type=float, help="Rescale branch lengths by a multiplicative factor to achieve requested root age")
    args = parser.parse_args()
    tree_fn = args.tree_fn
    r = args.root_age
    if not os.path.exists(tree_fn):
        print("Input tree file does not exist: %s" % tree_fn)
    t = tree.Tree(tree_fn, format=1)
    CheckTree(t)
    d = AveDist(t)
    print("Average distance from root to leaves: %f" % d)
    for n in t.traverse('preorder'):
        if n.is_root():
            n.dist = 0
            continue
        # work downwards, setting the branch distances from the top down
        x = t.get_distance(n) - n.dist
        y = n.dist
        print("\nTaxa:")
        print(", ".join(n.get_leaf_names()))
        if n.is_leaf():
            z = 0.
        else:
            z = AveDist(n)
        print("Distance of parent node from root: %f" % x)
        print("Current branch length: %f" % y)
        print("Average distance to leaves: %f" % z)
        if (y+z) == 0.:
            n.dist = 0
        else:
            f = (d-x)/(y + z)
            n.dist = f * n.dist
        print("Branch length for ultrametric tree: %f" % n.dist)
    if r != None:
        x = r/d
        print("\nRescaling branch lengths by factor of %0.2f so that root age is %f" % (x, r))
        for n in t.traverse():
            if n.is_root(): continue
            n.dist = x * n.dist
    outfn = tree_fn + ".ultrametric.tre"
    t.write(outfile=outfn, format=5)
    print("\nUltrametric tree written to: %s\n" % outfn)

if __name__ == "__main__":
    with util.Finalise():
        main()
