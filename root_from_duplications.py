# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 14:34:46 2016

@author: david
"""

"""
Performance improvements:
- Use an unrooted tree representation so that the trees don't need to be rerooted at all
- convert all the gene names to species names once per tree

"""
import os
import glob
#import cProfile as profile
#import pstats
import tree
import argparse
import itertools
import multiprocessing as mp
from collections import Counter

def cmp_equal(exp, act):
    """exp - expected set of species
    act - actual set of species
    """
    return exp == act
    
def cmp_noExtra(exp, act):
    """exp - expected set of species
    act - actual set of species
    """
    return act.issubset(exp)

# criteria with which to compare sets of species
#compare = cmp_equal
compare = cmp_noExtra

# which clades must satisfy criteria
criteria = "crit_all"
#criteria = "crit_one"

spTreeFormat = 3
qSerial = False
nProcs = 64
#nProcs = 16

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
        nodes = [c] if c.is_leaf() else c.get_children() 
        return [ch.sp_down for ch in nodes]
        
    def get_grandchild_species_clades(self, c):
        return [self.get_child_species_clades(ch) for ch in c.get_children()]

    def get_up_grand_species_clades(self):
        """ 
        Exceptions:
        - if tree is non-binary
        """
        if self.node.is_root(): 
            children = self.node.get_children()
            if len(children) != 3: raise Exception("Expected binary tree")
            return [self.get_grandchild_species_clades(c) for c in children]            
        ancestors = self.node.get_ancestors()
        c = ancestors[0]
        # must be at least one ancestor
        ch_temp = c.get_children()
        if len(ch_temp) == 3: 
            # both easy
            c1, c2 = ch_temp[:2] if ch_temp[2] != self.node else ch_temp[1:3] if ch_temp[0] != self.node else (ch_temp[0], ch_temp[2])
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
            if len(ch_temp) == 3:
                # easy - both are children
                c21, c22 = ch_temp[:2] if ch_temp[2] != c else ch_temp[1:3] if ch_temp[0] != c else (ch_temp[0], ch_temp[2])
                c21clade = c21.sp_down
                c22clade = c22.sp_down
            else:
                c21 = ch_temp[0] if ch_temp[0] != c else ch_temp[1]
                c21clade = c21.sp_down
                c22clade = c2.sp_up     # c21clade = everything not in the others, work out at the end
            return [c1clades, [c21clade, c22clade]]
        else:
            raise Exception("Expected binary tree")
      
    def get_grandrelative_clades_stored(self, allTaxa):
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

def GetStoredSpeciesSets(node):
    children = node.get_children()
    if node.is_root():
        if len(children) != 3: return None
        return [ch.sp_down for ch in children]
    else:
        if len(children) != 2: return None
        return [ch.sp_down for ch in children] + [node.sp_up]

def GeneToSpecies_dash(g):
  return g.split("_", 1)[0]
  
def GeneToSpecies_secondDash(g):
  return "_".join(g.split("_", 2)[:2])
  
def GeneToSpecies_3rdDash(g):
  return "_".join(g.split("_", 3)[:3])
  
def GeneToSpecies_dot(g):
  return g.split(".", 1)[0]
                   
def LocalCheck_clades(clade1, clade2, expClades, GeneToSpecies):
    """Expected clades are now in tree structure going down two levels: [[A,B], [C,D]]
    Observed clades should match up with expected clades in terms of containing no unexpected species
    """   
    X = set.union(*expClades[0])
    Y = set.union(*expClades[1])
    x0, x1 = expClades[0] if len(expClades[0]) == 2 else (None, None) if len(expClades[0]) == 1 else (Exception, Exception)
    y0, y1 = expClades[1] if len(expClades[1]) == 2 else (None, None) if len(expClades[1]) == 1 else (Exception, Exception)
    if x0 == Exception or y0 == Exception: raise Exception
    matches = 0
    for actClades in [clade1, clade2]:
        if criteria == "crit_one":
            for clade in actClades:
                # left is (A u B) or (C u D) if length 1 else (A, B) or (AuB, AuB) or (C, D) or (CuD, CuD) if length 2 
                if len(clade) == 1:
                    c = clade[0]
                    if compare(X, c):
                        matches += 1
                        break
                    elif compare(Y, c):
                        matches += 1
                        break
                else:
                    c0 = clade[0] 
                    c1 = clade[1]
                    if x0 != None and ((compare(x0, c0) or compare(x1, c1)) or (compare(x1, c0) or compare(x0, c1))):
                        matches += 1
                        break
                    elif compare(X, c0) or compare(X, c1):
                        matches += 1
                        break
                    elif y0 != None and ((compare(y0, c0) or compare(y1, c1)) or (compare(y1, c0) or compare(y0, c1))):
                        matches += 1
                        break
                    elif compare(Y, c0) or compare(Y, c1):
                        matches += 1
                        break
            return (matches == 2)
        elif criteria == "crit_all":
            iUsed = None
            for clade in actClades:
                # left is (A u B) or (C u D) if length 1 else (A, B) or (AuB, AuB) or (C, D) or (CuD, CuD) if length 2 
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
        else:
            raise NotImplemented  
     
def SupportedHierachies(t, G, S, GeneToSpecies, species, dict_clades, clade_names, treeName):
    """
    Only get the species sets in the first instance as work out the clades as and when
    """
    supported = []
    # Pre-calcualte species sets on tree: traverse tree from leaves inwards
    StoreSpeciesSets(t, GeneToSpecies, G)
    for counter, n in enumerate(t.traverse()):
        if n.is_leaf(): continue
        # just get the list of putative descendant species at this point without worrying about grandchild clades
        spSets = GetStoredSpeciesSets(n)
        if spSets == None: continue # non-binary
        clades = None
        # check each of three directions to the root
        for i, j in itertools.combinations(range(3), 2):
            s1 = spSets[i]
            s2 = spSets[j]
            if min(len(s1), len(s2)) < 2: continue
            k1 = frozenset(species)
            k2 = frozenset(species)
            for kk in dict_clades: 
                if  s1.issubset(kk) and len(kk) < len(k1): k1 = kk
                if  s2.issubset(kk) and len(kk) < len(k2): k2 = kk
            if k1 != k2: continue
            if len(k1) == 1:
                # nothing to learn from a terminal species
                continue
            elif k1 == species:
                # putative ancient duplication
                continue
            elif all([not clade.isdisjoint(s1) for clade0 in dict_clades[k1] for clade in clade0]) and all([not clade.isdisjoint(s2) for clade0 in dict_clades[k1] for clade in clade0]):
                # Passed the check that the required species are present but at this point don't know where in the tree
                # Get grandchild clades as required (could still avoid the more costly up clade if it's not required)
                if clades == None:
                    N = Node(n)
                    clades = N.get_grandrelative_clades_stored(G)
                    if clades == None: break   # locally non-binary in vicinity of node, skip to next node
                if not LocalCheck_clades(clades[i], clades[j], dict_clades[k1], GeneToSpecies): continue
                supported.append(frozenset(k1))
    return supported   
    
"""
Parallelisation wrappers
================================================================================================================================
"""

def SupportedHierachies_wrapper(treeName, GeneToSpecies, species, dict_clades, clade_names):
    if not os.path.exists(treeName): return []
    t = tree.Tree(treeName, format=1)
    G = set(t.get_leaf_names())
    S = list(set(map(GeneToSpecies, G)))
    if len(S) < 4:
        return []
    result = SupportedHierachies(t, G, S, GeneToSpecies, species, dict_clades, clade_names, treeName)    
#    print(treeName)
    return result
    
def SupportedHierachies_wrapper2(args):
    return SupportedHierachies_wrapper(*args)
    
"""
End of Parallelisation wrappers
================================================================================================================================
"""
   
def AnalyseSpeciesTree(speciesTree):
    species = frozenset(speciesTree.get_leaf_names())
    parts = list(speciesTree.get_partitions())
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
        clade_names[p] = n.name + "_" + str(hash("".join(p)))[-4:]
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
    
def ParsimonyRoot(allSpecies, clades, supported_clusters):
    c = Counter(supported_clusters)
    contradictions = dict()
    for clade in clades:
        clade_p = allSpecies.difference(clade)
        against = 0
        for observed, n in c.items():
            if (not observed.issubset(clade)) and (not observed.issubset(clade_p)):
                against += n
        contradictions[clade] = against
    m = min(contradictions.values())
    n = len(supported_clusters)
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

def GetRoot(speciesTreeFN, treesDir, GeneToSpeciesMap, nProcessors, treeFmt=None):
    if treeFmt == None: treeFmt = spTreeFormat
    speciesTree = tree.Tree(speciesTreeFN, format=treeFmt)
    species, dict_clades, clade_names = AnalyseSpeciesTree(speciesTree)
    pool = mp.Pool(nProcessors, maxtasksperchild=1)       
    list_of_lists = pool.map(SupportedHierachies_wrapper2, [(fn, GeneToSpeciesMap, species, dict_clades, clade_names) for fn in glob.glob(treesDir + "/*")])
    clusters = []
    for l in list_of_lists:
        clusters.extend(l)
    roots, nSupport = ParsimonyRoot(species, dict_clades.keys(), clusters)
    roots = list(set(roots))
    speciesTrees_rootedFNs =[]
    for i, r in enumerate(roots):
        speciesTree = RootAtClade(speciesTree, r) 
        speciesTree_rootedFN = os.path.splitext(speciesTreeFN)[0] + "_%d_rooted.txt" % i 
#    speciesTree = LabelNodes()
        speciesTree.write(outfile=speciesTree_rootedFN, format=4)
        speciesTrees_rootedFNs.append(speciesTree_rootedFN)
    return roots, clusters, speciesTrees_rootedFNs, nSupport
      
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tree")
    parser.add_argument("-s", "--separator", choices=("dot", "dash", "second_dash", "3rd_dash"))
    parser.add_argument("-S", "--Species_tree")
    parser.add_argument("-d", "--directory", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
      
    GeneToSpecies = GeneToSpecies_dash
    if args.separator and args.separator == "dot":
        GeneToSpecies = GeneToSpecies_dot  
    elif args.separator and args.separator == "second_dash":
        GeneToSpecies = GeneToSpecies_secondDash  
    elif args.separator and args.separator == "3rd_dash":
        GeneToSpecies = GeneToSpecies_3rdDash  
    
    if not args.directory:
        speciesTree = tree.Tree(args.Species_tree, format=spTreeFormat)
        species, dict_clades, clade_names = AnalyseSpeciesTree(speciesTree)
        t = tree.Tree(args.input_tree, format=1)
        G = set(t.get_leaf_names())
        S = list(set(map(GeneToSpecies, G)))
        filename = 'profile_stats.stats'
#        profile.run('c = SupportedHierachies_wrapper(args.input_tree, GeneToSpecies, species, dict_clades, clade_names) ', filename)
        c = SupportedHierachies_wrapper(args.input_tree, GeneToSpecies, species, dict_clades, clade_names)      
        for cc in c: print(cc)
#        stats = pstats.Stats('profile_stats.stats')
#        stats.sort_stats('cumulative')
#        stats.print_stats(0.3)
    elif qSerial:
        speciesTree = tree.Tree(args.Species_tree, format=spTreeFormat)
        species, dict_clades, clade_names = AnalyseSpeciesTree(speciesTree)
        clusters = []
        for fn in glob.glob(args.input_tree + "/*"):
            c = SupportedHierachies_wrapper(fn, GeneToSpecies, species, dict_clades, clade_names)
            clusters.extend(c)
        roots, nSupport = ParsimonyRoot(species, dict_clades.keys(), clusters)
        if len(roots) > 1: print("Observed %d duplications. %d support the best roots and %d contradict them." % (len(clusters), nSupport, len(clusters) - nSupport))
        else: print("Observed %d duplications. %d support the best root and %d contradict it." % (len(clusters), nSupport, len(clusters) - nSupport))
        print("Best outgroup(s) for species tree:")
        for r in roots: 
            print(r)
        if args.verbose: PlotTree(speciesTree, args.input_tree, clusters)
    else:
        root, clusters, _, nSupport = GetRoot(args.Species_tree, args.input_tree, GeneToSpecies, nProcs, treeFmt = 1)
        for r in root: print(r)
        speciesTree = tree.Tree(args.Species_tree, format=1)
        if args.verbose: PlotTree(speciesTree, args.input_tree, clusters, qSimplePrint=False)
 
