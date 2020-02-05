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

import sys
import time
import glob
import argparse
from collections import Counter, defaultdict

from . import tree, newick

if sys.version_info > (3,):
    long = int

class BitVector(object):
    """
    An efficient representation of splits of a tree.
    WARNING: By assumption, the taxa_index is the same between all BitVectors - this will not be checked!
    See: https://stackoverflow.com/questions/2147848/how-do-i-represent-and-work-with-n-bit-vectors-in-python
    """
    def __init__(self, taxa_index, taxon=None):
        """
        Creates a BitVector corresponding to the single-taxon split for taxon, or no taxa if taxon==None
        
        Args:
            taxa_index - the dictionary from taxa to their indices 
            taxon - the name of the taxon (string)
        """
        self.X = long()
        self.n = len(taxa_index)
        self.taxa_index = taxa_index
        if taxon != None:
            self.X |= 1<<taxa_index[taxon]
        
    def Add(self, other):
        self.X = self.X | other.X
        
    def AddList(self, others):
        for o in others:
            self.X = self.X | o.X
    
    def Invert(self):
        """
        Inverts the bits
        Args:
            None
        
        Returns:
            None        
        """
        mask = ((1<<self.n)-1)
        self.X ^= mask
        
    def Canonical(self):
        """
        Returns the split for this BitVector in the canonical form (arbitrarily, but consistently, defined) as a long int
        Args:
            None
        Returns:
            long int
        """
        if (self.X >> 0) & 1: 
            c = BitVector(self.taxa_index)
            c.X = self.X
            c.Invert()
            return c.X
        else:
            return self.X
            
    def Is(self, taxon):
        return (self.X >> self.taxa_index[taxon]) & 1

def UpdateSplits(all_splits, tree, taxa_index):
    """
    Adds the splits in tree to the list all_splits
    Args:
        all_splits - list of splits
        tree - Tree object
        taxa_index - dictionary from taxa to their index
        taxa_set - a set object of the taxa in taxa_index
        
    Returns:
        None
        
    Requirements:
        - The taxa in tree should be identical to the taxa in taxa_index
    """
    nRoot = len(tree.get_children())  
    if nRoot == 2:
        # don't double count at root. If there are two children then must only do once
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                if node.up.is_root():
                    nRoot -= 1
                    if nRoot == 0: continue
                s = BitVector(taxa_index, node.name)
                node.add_feature('split_down', s)
                all_splits.append((s.Canonical(), node.dist))
            elif node.is_root():
                continue
            else:
                if node.up.is_root():
                    nRoot -= 1
                    if nRoot == 0: continue
                s = BitVector(taxa_index)
                s.AddList([ch.split_down for ch in node.get_children()])    
                node.add_feature('split_down', s)
                all_splits.append((s.Canonical(), node.dist))
    else:
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                s = BitVector(taxa_index, node.name)
                node.add_feature('split_down', s)
                all_splits.append((s.Canonical(), node.dist))
            elif node.is_root():
                continue
            else:
                s = BitVector(taxa_index)
                s.AddList([ch.split_down for ch in node.get_children()])    
                node.add_feature('split_down', s)
                all_splits.append((s.Canonical(), node.dist))

def GetAllSplits(trees_dir):
    """
    Gets the set of splits in a set of trees
    
    Args:
        trees_dir - directory containing the input trees    
    
    Returns:
        list of (splits,length) tuples, taxa_index, taxa_ordered, number of trees
    """
    treesFNs = glob.glob(trees_dir + "/*")
    if len(treesFNs) == 0:
        print("ERROR: No trees found in directory")
        raise Exception()
    iFirst = 0
    first_tree = None
    while iFirst < len(treesFNs):
        try:
            first_tree = tree.Tree(treesFNs[iFirst])
            break
        except newick.NewickError:
            print("Could not read tree, skipping: %s" % treesFNs[iFirst])
            iFirst += 1
    if first_tree is None:
        print("ERROR: Could not read any of the trees for STAG species tree inference in %s" % trees_dir)
        util.Fail()
    taxa_ordered = first_tree.get_leaf_names()
    taxa_index = {taxon:i for i, taxon in enumerate(taxa_ordered)}
    splits = []
    UpdateSplits(splits, first_tree, taxa_index)
    for treeFN in treesFNs[iFirst+1:]:
        try:
            t = tree.Tree(treeFN)
            UpdateSplits(splits, t, taxa_index)
        except newick.NewickError:
            print("Could not read tree, skipping: %s" % treeFN)
    return splits, taxa_index, taxa_ordered, len(treesFNs)

def GetCompatibleSplits(splits_lengths):
    """
    Args:
        splits_lengths - list of (split, length) tuples
    Returns:
        compatible_splits, dict:split->list_of_lengths
    """
    splits = [s for s,l in splits_lengths]
    lengths = defaultdict(list)
    for s,l in splits_lengths:
        lengths[s].append(l)
    counter = Counter(splits)
    compatible_splits = []
    for x, n in counter.most_common():
        # Check if it is consistent with the tree 
        ok = True
        for y,_ in compatible_splits:
            z = x & y
            if not (z==x or z==y or z==0):
                ok = False
                break
        # Add the split
        if ok: compatible_splits.append((x,n))
    return compatible_splits, lengths
    
def ConstructTree(compatible_splits, split_lengths, taxa_index, taxa_ordered, nTrees):
    # Progressively build the tree - Inspired by Day's algorithm (but I've not checked how similar)
    # Start from the leaves 
#    t = ConstructTree(compatible_splits)
    nTrees = float(nTrees)
    t = tree.Tree()
    nodes_list = []
    for i, taxon in enumerate(taxa_ordered):
        n = tree.TreeNode()
        n.name = taxon
        x = BitVector(taxa_index, taxon).Canonical()
        n.dist = sum(split_lengths[x]) / float(len(split_lengths[x]))
        nodes_list.append(n)
        t.add_child(n)
    compatible_splits = sorted(compatible_splits) # ascending
    # now nodes are in preorder    
#    one = '1'
    for x, nSup in compatible_splits:    
        # have already put leaves on:
        if bin(x).count("1") == 1: continue 
        # 1's correspond to the cluster for the tree rooted on the first
        # Nodes are in the array: [0,1,2,3,4]
        # if cluster (2,3) is observed then -> [0, 1, (2,3), None, 4]
        # if cluster (2,3,4) is then observed then go to the first non-zero index for the node?
#        rep = bin(x)
        # the tree may be non-binary. We may add multiple child nodes at once
        iInsert = None
        for i, node in enumerate(nodes_list):
            if (x >> i) & 1 and node != None:
                # then this is one of the child nodes
                if iInsert == None:
                    iInsert = i
                    node_new = tree.TreeNode()
                    node_new.support = nSup/nTrees
                    node_new.dist = sum(split_lengths[x]) / float(len(split_lengths[x]))
                    node_new.add_child(node.detach())
                    nodes_list[i] = node_new
                    t.add_child(node_new)
                else:
                    node_new.add_child(node.detach())
                    nodes_list[i] = None
#        print(t)
    t.unroot()
    return t
    
    
def ConsensusTree(trees_dir):
    """
    Calcualtes a greedy consensus tree for the set of trees
    
    Args:
        trees_dir - directory containing the input trees    
    
    Returns:
        A greedy consensus tree
    
    Requires (will raise exception if not):
        - All input trees should have the same leaf set
        - All trees should be unrooted
    """
    # Get the set of splits in the trees
    splits_lengths, taxa_index, taxa_ordered, nTrees = GetAllSplits(trees_dir)
    compatible_splits, split_lengths = GetCompatibleSplits(splits_lengths)
    t = ConstructTree(compatible_splits, split_lengths, taxa_index, taxa_ordered, nTrees)
    return t

def main(args):
#    try:
    start = time.time()
    t = ConsensusTree(args.trees_dir)
    stop = time.time()    
    print(t)
    print("\n%f seconds" % (stop-start))
    # despite unrooting, t.write(outfile="greedy_consensus_tree.tre") still writes it as a rooted tree
    with open("greedy_consensus_tree.tre", 'wb') as outfile:
        outfile.write(t.write())
#    except:
#        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("trees_dir")
    args = parser.parse_args()
    main(args)



