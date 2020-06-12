#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2017 David Emms
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
#  When publishing work that uses OrthoFinder please cite:
#      Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
#      improves orthogroup inference accuracy, Genome Biology 16:157
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.com

import os
import sys
import csv
import glob
#import shelve
import datetime
import argparse
import itertools
import multiprocessing as mp
from collections import Counter, defaultdict

from . import probroot
from . import tree 

PY2 = sys.version_info <= (3,)
csv_write_mode = 'wb' if PY2 else 'wt'

def compare(exp, act):
    """exp - expected set of species
    act - actual set of species
    """
    return act.issubset(exp)

criteria = "crit_all"

nProcs = mp.cpu_count()
spTreeFormat = 1  

tooLarge = frozenset(["Too many genes"])

class Node(object):
    """ The class allows the user to get the 'child' nodes in any of the three directions
    rather than just for the two 'child' nodes in the tree model
    """
    def __init__(self, node):
        """
        node - tree node
        """
        self.node = node
        
    def get_child_species_clades(self, c):
        nodes = [c] if c.is_leaf() else c.get_children()  # what if non-binary?
        return [ch.sp_down for ch in nodes]
        
    def get_grandchild_species_clades(self, c):
        if c.is_root(): raise Exception
        if c.is_leaf():
            return [[c.sp_down]]
        else:
            return [self.get_child_species_clades(ch) for ch in c.get_children()]

    def get_up_grand_species_clades(self):
        """ 
        Exceptions:
        - if tree is non-binary
        """
        if self.node.is_root(): 
            children = self.node.get_children()
            if len(children) != 3: raise Exception("Expected binary tree")
            raise Exception("Don't know which of the 3 to return")
            # should it be self.get_grandchild_species_clades(children[-1])
            return [self.get_grandchild_species_clades(c) for c in children]            
        ancestors = self.node.get_ancestors() # ordered list starting with parent and progressing upwards
        c = ancestors[0]
        # must be at least one ancestor
        ch_temp = c.get_children()
        if c.is_root() and len(ch_temp) == 3: 
            # both easy
            c1, c2 = ch_temp[:2] if ch_temp[2] == self.node else ch_temp[1:3] if ch_temp[0] == self.node else (ch_temp[0], ch_temp[2])
            c1clades = self.get_child_species_clades(c1) 
            c2clades = self.get_child_species_clades(c2) 
            return [c1clades, c2clades]
        elif len(ch_temp) == 2:
            # easy one
            c1 = ch_temp[0] if ch_temp[0] != self.node else ch_temp[1]
            c1clades = self.get_child_species_clades(c1)
            # difficult one
            if len(ancestors) < 2: raise Exception("Expected binary tree")
            c2 = ancestors[1]
            ch_temp = c2.get_children()  
            if c2.is_root() and len(ch_temp) == 3:
                # easy - both are children
                c21, c22 = ch_temp[:2] if ch_temp[2] == c else ch_temp[1:3] if ch_temp[0] == c else (ch_temp[0], ch_temp[2])
                c21clade = c21.sp_down
                c22clade = c22.sp_down
            elif (not c2.is_root()) and len(ch_temp) == 2:
                c21 = ch_temp[0] if ch_temp[0] != c else ch_temp[1]
                c21clade = c21.sp_down
                c22clade = c2.sp_up     # c21clade = everything not in the others, work out at the end
            else:
                raise Exception("Expected binary tree")
            return [c1clades, [c21clade, c22clade]]
        else:
            raise Exception("Expected binary tree")
      
    def get_up_genes(self, nMax):
        if self.node.is_root(): 
            raise Exception("Error in duplicate gene identification, no 'up' node.")
        nUp = len(self.node.get_tree_root()) - len(self.node)
        if nUp > nMax: return tooLarge
        else:
            return frozenset(self.node.get_tree_root().get_leaf_names()).difference(frozenset(self.node.get_leaf_names()))        
      
    def get_gene_sets(self, i, j, nMax=2000):
        """
        Mirrors get_grandrelative_clades_stored. 
        i - 0,1,2: set of clades with respect to order in get_grandrelative_clades_stored
        j - 0,1,2: set of clades with respect to order in get_grandrelative_clades_stored
        
        Returns the actuall genes in the clades
        """
        if self.node.is_root():
            children = self.node.get_children()
            if len(children) != 3: return None
            return (tooLarge if len(children[i]) > nMax else frozenset(children[i].get_leaf_names()), tooLarge if len(children[j]) > nMax else frozenset(children[j].get_leaf_names()))
        else:
            children = self.node.get_children()
            if len(children) != 2: return frozenset([]), frozenset([])
            if i <= 1:
                iGenes = frozenset(children[i].get_leaf_names()) if len(children[i]) < nMax else tooLarge
            else:
                iGenes = frozenset(self.get_up_genes(nMax))
            if j <=1:
                jGenes = frozenset(children[j].get_leaf_names()) if len(children[j]) < nMax else tooLarge
            else:
                jGenes = frozenset(self.get_up_genes(nMax))
            return iGenes, jGenes
      
    def get_grandrelative_clades_stored(self):
        """
        returns the hierarchically stored sets of species 
        """
        if self.node.is_root():
            children = self.node.get_children()
            if len(children) != 3: return None
            return [self.get_grandchild_species_clades(c) for c in children]
        else:
            children = self.node.get_children()
            if len(children) != 2: return None
            down_clades = [self.get_grandchild_species_clades(c) for c in children]
            return down_clades + [self.get_up_grand_species_clades()]
            
    def GetSpeciesSets(self, allTaxa, GeneMap):
        if self.node.is_root():
            na, nb, nc = self.node.get_children()
            a = na.get_leaf_names()
            b = nb.get_leaf_names()
            c = nc.get_leaf_names()
        else:
            ch = self.node.get_children()
            if len(ch) != 2: return None
            na, nb = ch
            a = set(na.get_leaf_names())
            b = set(nb.get_leaf_names())
            # now get node with set of species, c, below it
            c = allTaxa.difference(a).difference(b)
        return [set(map(GeneMap, x)) for x in (a,b,c)]

    @staticmethod
    def ToSpecies(clades_of_clades, GeneMap):
        return [[set(map(GeneMap, [gene for gene in grandchild])) for grandchild in child] for child in clades_of_clades]
        
    @staticmethod            
    def FlatenSpeciesSet(clades_of_clades):
        return set([sp for child in clades_of_clades for grandchild in child for sp in grandchild])

               
def get_partitions(tree):
     all_leaves = frozenset(tree.get_leaf_names())
     all_partitions = set([all_leaves])
     for n in tree.iter_descendants():
        p1 = frozenset(n.get_leaf_names())
        p2 = frozenset(all_leaves - p1)
        all_partitions.add(p1)
        all_partitions.add(p2)
     return all_partitions      
   
def StoreSpeciesSets(t, GeneMap, allTaxa):
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
      
def SaveTree(tree, root_clade, cladeName, treeName, iExample):
    fn = outputDir + "/%s_%s_%d_%d.tre" % (cladeName, os.path.split(treeName)[1].split(".")[0], iExample, len(root_clade))
    t = RootAtClade(tree, root_clade)
    t.write(outfile = fn)
    print(fn)
    
def StoreGeneSets(t):
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node.add_feature('g_down', [node.name])
        elif node.is_root():
            continue
        else:
            node.add_feature('g_down', [g for ch in node.get_children() for g in ch.g_down])
    for node in t.traverse('preorder'):
        if node.is_root():
            node.add_feature('g_up', [])
        else:
            parent = node.up
            if parent.is_root():
                others = [ch for ch in parent.get_children() if ch != node]
                x = [g for other in others for g in other.g_down]
                node.add_feature('g_up', x)
            else:
                others = [ch for ch in parent.get_children() if ch != node]
                g_downs = [g for other in others for g in other.g_down]
                x = parent.g_up + g_downs
                node.add_feature('g_up', x)

def GetStoredSpeciesSets(node):
    children = node.get_children()
    if node.is_root():
        if len(children) != 3: return None
        return [ch.sp_down for ch in children]
    else:
        if len(children) != 2: return None
        return [ch.sp_down for ch in children] + [node.sp_up]
        
def GetStoredGeneSets(node):
    children = node.get_children()
    if node.is_root():
        if len(children) != 3: return None
        return [ch.g_down for ch in children]
    else:
        if len(children) != 2: return None
        return [ch.g_down for ch in children] + [node.g_up]

def GeneToSpecies_dash(g):
  return g.split("_", 1)[0]
  
def GeneToSpecies_secondDash(g):
  return "_".join(g.split("_", 2)[:2])
  
def GeneToSpecies_3rdDash(g):
  return "_".join(g.split("_", 3)[:3])
  
def GeneToSpecies_dot(g):
  return g.split(".", 1)[0]
  
def GeneToSpecies_hyphen(g):
  return g.split("-", 1)[0]
                   
def LocalCheck_clades(clade1, clade2, expClades, GeneToSpecies):
    """Expected clades are now in tree structure going down two levels: [[A,B], [C,D]]
    Observed clades should match up with expected clades in terms of containing no unexpected species
    """   
    X = set.union(*expClades[0])
    Y = set.union(*expClades[1])
    x0, x1 = expClades[0] if len(expClades[0]) == 2 else (None, None) if len(expClades[0]) == 1 else (Exception, Exception)
    y0, y1 = expClades[1] if len(expClades[1]) == 2 else (None, None) if len(expClades[1]) == 1 else (Exception, Exception)
    if x0 == Exception or y0 == Exception: 
        return False   # Can't get a well-supported duplication meeting topology criteria if tree is not fully resolved
    for actClades in [clade1, clade2]:
            iUsed = None
            for clade in actClades:
                if len(clade) == 1:
                    c = clade[0]
                    if compare(X, c) and iUsed != 0:
                        iUsed = 0
                        continue
                    elif compare(Y, c) and iUsed != 1:
                        iUsed = 1
                        continue
                    else:
                        return False
                else:
                    c0 = clade[0] 
                    c1 = clade[1]
                    if (iUsed != 0 and x0 != None) and ((compare(x0, c0) and compare(x1, c1)) or (compare(x1, c0) and compare(x0, c1))):
                        iUsed = 0
                        continue
                    elif iUsed != 0 and (compare(X, c0) and compare(X, c1)):
                        iUsed = 0
                        continue
                    elif (iUsed != 1 and y0 != None) and ((compare(y0, c0) and compare(y1, c1)) or (compare(y1, c0) and compare(y0, c1))):
                        iUsed = 1
                        continue
                    elif iUsed != 1 and (compare(Y, c0) and compare(Y, c1)):
                        iUsed = 1
                        continue
                    else: 
                        return False
    return True 

def Join(clades):
     return set([s for child in clades for grandchild in child for s in grandchild])
     
def SupportedHierachies(t, G, S, GeneToSpecies, species, dict_clades, clade_names, treeName, qWriteDupTrees=False):
    """
    Only get the species sets in the first instance as work out the clades as and when
    """
    qAncient = False
    supported = defaultdict(int)
    genesPostDup = set()
    # Pre-calcualte species sets on tree: traverse tree from leaves inwards
    StoreSpeciesSets(t, GeneToSpecies, G)
    if qWriteDupTrees:
        t_write = None
        StoreGeneSets(t)
        for n in t.traverse():
            if n.is_root(): continue
        iExample = 0
    range3 = list(range(3))
    for counter, n in enumerate(t.traverse()):
        if n.is_leaf(): continue
        # just get the list of putative descendant species at this point without worrying about grandchild clades
        spSets = GetStoredSpeciesSets(n)
        if spSets == None: continue # non-binary
        clades = None
        # check each of three directions to the root
        for i, j in itertools.combinations(range3, 2):
            s1 = spSets[i]
            s2 = spSets[j]
            # check for terminal duplications
            if len(s1) == 1 and s1 == s2:
                supported[frozenset(s1)] += 1
            if min(len(s1), len(s2)) < 2: continue
            k1 = frozenset(species)
            k2 = frozenset(species)
            for kk in dict_clades: 
                if  s1.issubset(kk) and len(kk) < len(k1): k1 = kk
                if  s2.issubset(kk) and len(kk) < len(k2): k2 = kk
            edges1 = set([kk for kk in dict_clades if s1.issubset(kk) and len(kk) == len(k1)])
            edges2 = set([kk for kk in dict_clades if s2.issubset(kk) and len(kk) == len(k2)])
            for k1 in edges1.intersection(edges2):
                if len(k1) == 1:
                    # terminal duplciations dealt with above (since no extra info about species tree topology is needed)
                    break
                elif k1 == species:
                    # putative ancient duplication
                    if not qAncient: 
                        #print(treeName)
                        qAncient = True
                    break
                elif all([not clade.isdisjoint(s1) for clade0 in dict_clades[k1] for clade in clade0]) and all([not clade.isdisjoint(s2) for clade0 in dict_clades[k1] for clade in clade0]):
                    # Passed the check that the required species are present but at this point don't know where in the tree
                    # Get grandchild clades as required (could still avoid the more costly up clade if it's not required)
                    if clades == None:
                        N = Node(n)
                        clades = N.get_grandrelative_clades_stored()
                        if Join(clades[0]) != spSets[0] or Join(clades[1]) != spSets[1] or Join(clades[2]) != spSets[2]:
                            print((clades[0]))
                            print((spSets[0]))
                            print("")
                            print((clades[1]))
                            print((spSets[1]))
                            print("")
                            print((clades[2]))
                            print((spSets[2]))
                            print("")
                            raise Exception("Mismatch")
                        if clades == None: break   # locally non-binary in vicinity of node, skip to next node
                    if not LocalCheck_clades(clades[i], clades[j], dict_clades[k1], GeneToSpecies): continue
                    supported[frozenset(k1)] +=1
                    genes = N.get_gene_sets(i, j)
                    genesPostDup.add(genes[0].union(genes[1]))
                    if qWriteDupTrees:
                        if t_write == None: 
                            try:
                                t_write = t.copy()
                            except:
                                t_write = tree.Tree(treeName, format=1)
                        ii = 0 if (0!= i and 0!=j) else 1 if (1!=i and 1!=j) else 2
                        gSets = GetStoredGeneSets(n)
                        SaveTree(t_write, gSets[ii], clade_names[k1], treeName, iExample)
                        iExample += 1       
    return supported, genesPostDup  
    
"""
Parallelisation wrappers
================================================================================================================================
"""

def SupportedHierachies_wrapper(treeName, GeneToSpecies, species, dict_clades, clade_names, qWriteDupTrees=False):
    if not os.path.exists(treeName): return [], []
    try:
        t = tree.Tree(treeName, format=1)
    except:
        return [], []
    G = set(t.get_leaf_names())
    S = set(map(GeneToSpecies, G))
    if not S.issubset(species):
        print(("ERROR in %s" % treeName))
        print("Some genes cannot be mapped to species in the species tree")
        print((S.difference(species)))
        return None
    if len(S) < 3:
        return defaultdict(int), []
    result = SupportedHierachies(t, G, S, GeneToSpecies, species, dict_clades, clade_names, treeName, qWriteDupTrees)
    return result
    
def SupportedHierachies_wrapper2(args):
    return SupportedHierachies_wrapper(*args)
    
"""
End of Parallelisation wrappers
================================================================================================================================
"""
   
def AnalyseSpeciesTree(speciesTree):
    species = frozenset(speciesTree.get_leaf_names())
    parts = list(get_partitions(speciesTree))
    nSpecies = len(species)
    dict_clades = dict() # dictionary of clades we require evidence of duplicates from for each partition
    clade_names = dict()
    for p in parts:
        if len(p) == nSpecies: continue
        speciesTree.set_outgroup(list(species.difference(p))[0])
        if len(p) == 1:
            # no use for identifying clades
            continue
        n = speciesTree.get_common_ancestor(p)
        clade_names[p] = n.name + "_" + str(hash("".join(p)))[-8:]
        # go down two steps
        ch = n.get_children()
        ch0 = [ch[0]] if ch[0].is_leaf() else ch[0].get_children() 
        ch1 = [ch[1]] if ch[1].is_leaf() else ch[1].get_children() 
        dict_clades[p] = [[set(c.get_leaf_names()) for c in ch0], [set(c.get_leaf_names()) for c in ch1]]
    return species, dict_clades, clade_names

def RootAtClade(t, accs_in_clade):
    if len(accs_in_clade) == 1:
        t.set_outgroup(list(accs_in_clade)[0])
        return t
    accs = set(t.get_leaf_names())
    dummy = list(accs.difference(accs_in_clade))[0]
    t.set_outgroup(dummy)
    node = t.get_common_ancestor(accs_in_clade)
    t.set_outgroup(node)
    return t
    
def ParsimonyRoot(allSpecies, clades, supported_clusters_counter):
    contradictions = dict()
    for clade in clades:
        clade_p = allSpecies.difference(clade)
        against = 0
        for observed, n in supported_clusters_counter.items():
            if (not observed.issubset(clade)) and (not observed.issubset(clade_p)):
                against += n
        contradictions[clade] = against
    m = min(contradictions.values())
    n = sum(supported_clusters_counter.values())
    nSupport = n-m
    roots = []
    for clade, score in contradictions.items():
        if score == m:
            if len(clade) > len(allSpecies)/2:
                roots.append(allSpecies.difference(clade))
            elif len(clade) == len(allSpecies)/ 2:
                if allSpecies.difference(clade) not in roots: roots.append(clade) 
            else:
                roots.append(clade)
    return roots, nSupport

def GetRoot(speciesTreeFN, treesDir, GeneToSpeciesMap, nProcessors, qWriteDupTrees=False, qWriteRootedTree=False):
    """ 
                    ******* The Main method ******* 
    """
    qHaveBranchSupport = False
    try:
        speciesTree = tree.Tree(speciesTreeFN, format=2)
        qHaveBranchSupport = True
    except:
        speciesTree = tree.Tree(speciesTreeFN, format=1)
    species, dict_clades, clade_names = AnalyseSpeciesTree(speciesTree)
    pool = mp.Pool(nProcessors, maxtasksperchild=1)       
    list_of_dicts = pool.map(SupportedHierachies_wrapper2, [(fn, GeneToSpeciesMap, species, dict_clades, clade_names, qWriteDupTrees) for fn in glob.glob(treesDir + "/*")])
    pool.close()
    clusters = Counter()
    all_stride_dup_genes = set()
    for l, stride_dup_genes in list_of_dicts:
        if l == None:
            sys.exit()
        clusters.update(l)
        all_stride_dup_genes.update(stride_dup_genes)
    roots, nSupport = ParsimonyRoot(species, list(dict_clades.keys()), clusters)
    roots = list(set(roots))
    speciesTrees_rootedFNs =[]
    # Get distance of each from a supported clade
    topDist = []
    branchDist = []
    if len(clusters) > 0 and len(roots) > 1:
        # Evaluate which is the best one
        for r in roots:
            # 1. Root at potential outgroup
            speciesTree = RootAtClade(speciesTree, r)
            # 2. Get minimum topological distance to any clade with gene duplication evidence
            topDist.append(min([speciesTree.get_distance(n, topology_only=True) for n in speciesTree.traverse('preorder') if frozenset(n.get_leaf_names()) in clusters]))
            # 3. Get minimum branch length to outgroups that are also the minimum topological distance
            branchDist.append(min([speciesTree.get_distance(n) for n in speciesTree.traverse('preorder') if (frozenset(n.get_leaf_names()) in clusters and speciesTree.get_distance(n, topology_only=True) == topDist[-1])]))        
        # 4. Get the distance for the root which is furthest from its closest outgroup
        maxTopDist = max(topDist)
        bestDist = -1
        for i, (t, d) in enumerate(zip(topDist, branchDist)):
            if t == maxTopDist and d > bestDist:
                bestDist = d
                iRootToUse = i
        rootToUse = roots.pop(iRootToUse)
        roots = [rootToUse] + roots
        
    if qWriteRootedTree:
        for i, r in enumerate(roots):
            speciesTree = RootAtClade(speciesTree, r) 
            speciesTree_rootedFN = os.path.splitext(speciesTreeFN)[0] + "_%d_rooted.txt" % i 
    #    speciesTree = LabelNodes()
            speciesTree.write(outfile=speciesTree_rootedFN, format = 2 if qHaveBranchSupport else 5)   # 5 With all branch lengths. No support or node names. 1 has node names
            speciesTrees_rootedFNs.append(speciesTree_rootedFN)
    return roots, clusters, speciesTrees_rootedFNs, nSupport, list(dict_clades.keys()), species, all_stride_dup_genes

def PrintRootingSummary(roots, clusters_counter, nSupport):
    nAll = sum(clusters_counter.values())
    nFP_mp = nAll - nSupport
    n_non_trivial = sum([v for k, v in clusters_counter.items() if len(k) > 1])
    if len(roots) > 1: print(("Identified %d non-terminal duplications.\n%d support the best roots and %d contradict them." % (n_non_trivial, n_non_trivial-nFP_mp, nFP_mp)))
    else: print(("Identified %d non-terminal duplications.\n%d support the best root and %d contradict it." % (n_non_trivial, n_non_trivial-nFP_mp, nFP_mp)))
    print("Most parsimonious outgroup(s) for species tree:")
    for r in roots[:5]: 
        print(("{" + ", ".join(r) + "}"))
    if len(roots) > 5:
        print("Etc...")
        print(("%d possible roots" % len(roots)))
    return nFP_mp, n_non_trivial 
    
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
        
def GetCluseterName(species_tree, S, cluster):
    if len(cluster) == 1: return list(cluster)[0]
    else:
        n = species_tree.get_common_ancestor(cluster) 
    if n == species_tree or (n.up == species_tree and len(S.difference(cluster)) < len(cluster)):
        complement = S.difference(cluster) # complement
#            print(complement)
        if len(complement) == 1:
#                continue
            n = species_tree & list(complement)[0]
        else:
            n = species_tree.get_common_ancestor(complement)
    return n.name + "_" + list(cluster)[0]
      
def WriteResults(species_tree_fn_or_text, roots, S, clades, clusters_counter, output_dir):
#    for c in clusters_counter:
#        print((clusters_counter[c], c))
    print(("\nResults written to:\n" + os.path.realpath(output_dir)))
    # Label species tree nodes
    species_tree = tree.Tree(species_tree_fn_or_text)
    thisRoot = roots[0]
    species_tree = RootAtClade(species_tree, thisRoot) 
    iNode = 0
    for n in species_tree.traverse():
        if not n.is_leaf():
            n.name = "N%d" % iNode
            iNode+=1
    species_tree.write(outfile=output_dir + "Species_tree_labelled.tre", format=1)
#    print(species_tree)
#    species_tree = tree.Tree(output_dir + "Species_tree_labelled.tre", format=1)
    # Calculate probabilities
    qBinary = True
    for n in species_tree.traverse():
        if len(n.get_children()) > 2:
            qBinary = False
    if qBinary:
        p_final = probroot.GetProbabilities(species_tree, S, clades, clusters_counter)
    else:
        print("Probability distribution for root location is not supported for non-binary trees")
        print("To get a probability distribution for the root, please supply a fully resolved input species tree")
    # Write numbers of duplications
    table = dict()
    new_tree = tree.Tree(output_dir + "Species_tree_labelled.tre", format=1)
    for clade in clades + [frozenset([s]) for s in S]:
        qAnti = False
        anticlade = S.difference(clade)
        if len(clade) == 1:
            node = new_tree & list(clade)[0]
        else:
            node = new_tree.get_common_ancestor(clade)
        if node == new_tree:
            node = new_tree.get_common_ancestor(anticlade)
            qAnti = True
        x = anticlade if qAnti else clade
        y = clade if qAnti else anticlade
        X = ("(%d)" % clusters_counter[x]) if len(clade) == 1 else clusters_counter[x] 
        if qBinary:
            p = p_final[clade] if clade in p_final else p_final[anticlade]
        else:
            p = 0.
        table[node.name] = [node.name, "X" if (clade in roots or anticlade in roots) else "", "%0.1f%%" % (100.*p) , X, clusters_counter[y]]
    with open(output_dir + "Duplication_counts.tsv", csv_write_mode) as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Branch", "MP Root", "Probability", "Duplications supporting clade", "Duplications supporting opposite clade"])
        qSingle = len(thisRoot) == 1
        root_branches = [n.name for n in new_tree.get_children()]
        writer.writerow([root_branches[0] + " (& " + root_branches[1] + ")"] + table[root_branches[0]][1:])
        for i in range(2 if qSingle else 3, iNode):  
            name = "N%d" % i
            if name in table:
                writer.writerow(table[name])
            else:
                print(("Skipping %s" % name))
        for sp in S:
            if sp in table:
                if qSingle and sp in thisRoot: continue
                writer.writerow(table[sp])
   
def Main_Full(args):
    text = """
****************************************************************************
*     STRIDE: Species Tree Root Inference from Gene Duplication Events     *
*                                                                          *
****************************************************************************"""
#    text = "STRIDE: Species Tree Root Inference from Gene Duplication Events"
    print((text[1:]))
#    print(text + "\n" + "="*len(text))
    GeneToSpecies = GeneToSpecies_dash
    if args.separator and args.separator == "dot":
        GeneToSpecies = GeneToSpecies_dot  
    elif args.separator and args.separator == "second_dash":
        GeneToSpecies = GeneToSpecies_secondDash  
    elif args.separator and args.separator == "3rd_dash":
        GeneToSpecies = GeneToSpecies_3rdDash  
    elif args.separator and args.separator == "hyphen":
        GeneToSpecies = GeneToSpecies_hyphen 
        
    if not args.directory:
        speciesTree = tree.Tree(args.Species_tree, format=spTreeFormat)
        species, dict_clades, clade_names = AnalyseSpeciesTree(speciesTree)
        c, stride_dup_genes = SupportedHierachies_wrapper(args.gene_trees, GeneToSpecies, species, dict_clades, clade_names)      
        for k, v in c.items(): print((k, v))
#    elif args.debug:
#        speciesTree = tree.Tree(args.Species_tree, format=spTreeFormat)
#        species, dict_clades, clade_names = AnalyseSpeciesTree(speciesTree)
#        clusters_counter = Counter()
#        for fn in glob.glob(args.gene_trees + "/*"):
#            c, stride_dup_genes = SupportedHierachies_wrapper(fn, GeneToSpecies, species, dict_clades, clade_names)
#            clusters_counter.update(c)
#        roots, nSupport = ParsimonyRoot(species, dict_clades.keys(), clusters_counter)
#        PrintRootingSummary(roots, clusters_counter, nSupport)
    else:
        nTrees = len(glob.glob(args.gene_trees + "/*"))
        if nTrees == 0:
            print(("No trees found in %s\nExiting" % args.gene_trees))
            sys.exit()
        print(("Analysing %d gene trees" % nTrees))
#        roots, clusters_counter, _, nSupport, clades, species = GetRoot(args.Species_tree, args.gene_trees, GeneToSpecies, nProcs, treeFmt = 1, qWriteDupTrees=args.output)
        roots, clusters_counter, _, nSupport, clades, species, all_stride_dup_genes = GetRoot(args.Species_tree, args.gene_trees, GeneToSpecies, nProcs)
        PrintRootingSummary(roots, clusters_counter, nSupport)
        outputDir = CreateNewWorkingDirectory(args.gene_trees + "/../STRIDE_Results")
#        shelveFN = outputDir + "STRIDE_data.shv"
#        d = shelve.open(shelveFN)
#        d['roots'] = roots
#        d['clusters_counter'] = clusters_counter
#        d['species'] = species
#        d['nSupport'] = nSupport
#        d['SpeciesTreeFN'] = os.path.abspath(args.Species_tree)
#        with open(args.Species_tree, 'r') as infile:
#            tree_text = "".join([l.rstrip() for l in infile.readlines()])
#        d['SpeciesTreeText'] = tree_text
#        d['TreesDir'] = os.path.abspath(args.gene_trees)
#        d['clades'] = clades
#        d.close()
#        outputFigFN = outputFN_base + ".pdf"
        #DrawDuplicationsTree(args.Species_tree, clusters_counter, outputFigFN)
        WriteResults(args.Species_tree, roots, species, clades, clusters_counter, outputDir)
        print("")
      
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gene_trees", help = "Directory conaining gene trees (with -d argument) or filename of a single gene tree to analyse (no -d argument)")
    parser.add_argument("-s", "--separator", choices=("dot", "dash", "second_dash", "3rd_dash", "hyphen"), help="Separator been species name and gene name in gene tree taxa")
    parser.add_argument("-S", "--Species_tree", help="Unrooted species tree in newick format")
    parser.add_argument("-d", "--directory", action="store_true", help="Process all trees in input directory")
#    parser.add_argument("--debug", action="store_true", help="Run in serial to enable easier debugging")
#    parser.add_argument("-o", "--output", action="store_true", help="Write out gene trees rooted at duplications")
    parser.set_defaults(Func=Main_Full)   
    args = parser.parse_args()
#    if args.output:
#        x = (args.gene_trees if args.directory else os.path.split(args.gene_trees)[0]) 
#        while x[-1] == "/":
#            x = x[:-1]
#        outputDir = x + "_rooted_duplications/" 
#        if not os.path.exists(outputDir): os.mkdir(outputDir)
    args.Func(args)
    
 
