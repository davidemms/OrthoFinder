# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:11:11 2017

@author: david

Perform directed 'reconciliation' first and then apply EggNOG method

1 - root gene trees on outgroup: unique one this time
2 - infer orthologues
"""
import os
import sys
import csv
import glob
import argparse
import operator
import itertools
import multiprocessing as mp
from collections import defaultdict, deque

from . import tree as tree_lib
from . import resolve, util, files, parallel_task_manager
try: 
    import queue
except ImportError:
    import Queue as queue   

PY2 = sys.version_info <= (3,)
csv_write_mode = 'wb' if PY2 else 'wt'
csv_append_mode = 'ab' if PY2 else 'at'
csv_read_mode = 'rb' if PY2 else 'rt'

debug = False   # HOGs

if not PY2:
    xrange = range

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


def SpeciesAndGene_dash(g):
  return g.split("_", 1)
    
def SpeciesAndGene_secondDash(g):
    a,b,c = g.split("_", 2)
    return (a+"_"+b, c)
  
def SpeciesAndGene_3rdDash(g):
    a,b,c,d = g.split("_", 3)
    return (a+"_"+b+"_"+c, d)
  
def SpeciesAndGene_dot(g):
  return g.split(".", 1)
  
def SpeciesAndGene_hyphen(g):
  return g.split("-", 1)
  
SpeciesAndGene_lookup = {GeneToSpecies_dash:SpeciesAndGene_dash, 
                        GeneToSpecies_secondDash:SpeciesAndGene_secondDash,
                        GeneToSpecies_3rdDash:SpeciesAndGene_3rdDash,
                        GeneToSpecies_dot:SpeciesAndGene_dot,
                        GeneToSpecies_hyphen:SpeciesAndGene_hyphen}
    
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

"""
HOGs
-------------------------------------------------------------------------------
""" 

def MRCA_node(t_rooted, taxa):
    return (t_rooted & next(taxon for taxon in taxa)) if len(taxa) == 1 else t_rooted.get_common_ancestor(taxa)

class HogWriter(object):
    def __init__(self, species_tree, species_tree_node_names, seq_ids, sp_ids, species_to_use):
        """
        Prepare files, get ready to write.
        species_tree_node_names - list of species tree nodes
        seq_ids - dict of sequence ids
        sp_ids - dict of species ids
        species_to_use - list of ints
        """
        self.seq_ids = seq_ids
        d = os.path.dirname(files.FileHandler.GetHierarchicalOrthogroupsFN("N0"))
        if not os.path.exists(d):
            os.mkdir(d)
        self.fhs = dict()
        self.iSps = list(map(str, sorted(species_to_use)))   # list of strings
        self.i_sp_to_index = {isp:i_col for i_col, isp in enumerate(self.iSps)}
        self.iHOG = defaultdict(int)
        self.species_tree = species_tree
        species_names = [sp_ids[i] for i in self.iSps]
        for name in species_tree_node_names:
            fn = files.FileHandler.GetHierarchicalOrthogroupsFN(name)
            self.fhs[name] = open(fn, csv_write_mode)
            util.writerow(self.fhs[name], ["HOG", "OG", "Gene Tree Parent Clade"] + species_names)
            self.fhs[name].flush()
        # Map from HOGs to genes that must be contained in them
        self.hog_contents = dict()  # sp_node_name = hog_name-> list of contents fo hog (internal nodes and leaves)
        for n in species_tree.traverse():
            desc = n.get_descendants()
            self.hog_contents[n.name] = set([int(nn.name) if nn.is_leaf() else nn.name for nn in desc])
        self.comp_nodes = get_comparable_nodes(self.species_tree)

    def get_hog_index(self, hog_name):
        i = self.iHOG[hog_name]
        self.iHOG[hog_name] += 1
        return i

    def write_hog_genes(self, genes, sp_node_name_list, og_name):
        if len(sp_node_name_list) == 0: return
        genes_per_species = defaultdict(list)
        for g in genes:
            isp, _ = g.split("_")
            genes_per_species[isp].append(self.seq_ids[g])
        i_hogs = [self.get_hog_index(sp_node_name) for sp_node_name in sp_node_name_list]
        row_genes = [", ".join(genes_per_species[isp]) for isp in self.iSps] 
        for i_hog, sp_node_name in zip(i_hogs, sp_node_name_list):
            util.writerow(self.fhs[sp_node_name], ["%s.HOG%07d" % (sp_node_name, i_hog),  og_name, "-"] + row_genes)


    def write_clade_v2(self, n, og_name, split_paralogous_clades_from_same_hog = False):
        """
        Look at parent node to know when to start, look at dups below to know when 
        to stop.
        - Current MRCA could be excluded either because it's already been done or
          because of a duplication below

        - Species-specific clades could still be HOGs if they are all that remain
          of that clade 
          Args:
            n - gene tree node
            og_name - name to use in file output
            split_paralogous_clades_from_same_hog - should clades which are within 
                the same HOG but are paralogous be split up
        """
        if n.is_leaf():
            return []
        if debug: print("\nTree node: %s" % n.name)
        if debug: print(n.sp_node)
        if (n.dup and n.sp_node == "N0"): 
            n.add_feature("done", set())
        # self.comp_nodes[n.sp_node] is the set of HOGs relevant to this node

        # Only skip doing HOGs for above if it is a dup and want to split paralogous clades from same HOG
        ch = n.get_children()
        if (split_paralogous_clades_from_same_hog and n.dup and (ch[0].sp_node == ch[1].sp_node)):
            # continue to record single-species orthogroups
            hogs_to_write = set() if n.sp_node.startswith("N") else self.comp_nodes[n.sp_node][0].copy()
        else:
            hogs_to_write = self.comp_nodes[n.sp_node][0].copy()
        
        # get scl & remove HOGs that can't be written yet due to duplications
        # 0. Get the scl units below this node in the gene tree
        # I.e. get genes (indexed by species) below each scl (the relevant gene tree nodes)
        genes_per_species_index = self.get_descendant_genes(n)

        # scl_mrca = {nn.sp_node for nn in scl if not nn.is_leaf()}
        if debug: print("Dups below: " + str(n.dups_below))
        stop_at_dups = lambda nn : nn.name in n.dups_below
        sp_node = self.species_tree & (n.sp_node)
        # don't need skip for dups, that's recorded in dups_below
        # traverse the species tree from the current node and record all nodes before hitting a duplication node from the gene tree
        hogs_to_write.update({nn.name for nn in sp_node.traverse('preorder', is_leaf_fn = stop_at_dups) if (not nn.is_leaf()) and (not nn.name in n.dups_below)})
        
        if not n.is_root():
            hogs_to_write.difference_update(n.up.done)

        # 2. Write HOGs
        n.add_feature("done", hogs_to_write if n.is_root() else n.up.done.union(hogs_to_write))
        if len(hogs_to_write) == 0:
            return []

        if debug: print(hogs_to_write)
        return self.get_hog_file_entries(hogs_to_write, genes_per_species_index, og_name, n.name)

    def get_descendant_genes(self, n):
        """
        Attempt at a simplified replacement to get_scl_units as shouldn't need to 
        care about which scl a gene belongs to.
        Args:
            n - node under consideration
        Returns:
            dict:sp_index->string of genes, comma separated
        """
        genes_per_species = defaultdict(list) # iCol (before 'name' columns) -> text string of genes
        genes = n.get_leaves()
        q_have_legitimate_gene = False # may be misplaced genes
        for g in genes:
            if "X" in g.features: continue
            q_have_legitimate_gene = True
            isp = g.name.split("_")[0]
            genes_per_species[self.i_sp_to_index[isp]].append(self.seq_ids[g.name])
        for k, v in genes_per_species.items():
            genes_per_species[k] = ", ".join(v)
        return genes_per_species

    def get_hog_file_entries(self, hogs_to_write, genes_per_species_index, og_name, gt_node_name):
        """
        Write the HOGs that can be determined from this gene tree node.
        Args:
            hogs_to_write - list of HOG names
            scl_units - dict:st_node_name->(dict:sp_index->genes string for HOGs file)   
                        e.g. st_node_name is N5 for internal node, 11 for leaf
            og_name - OG name
            gt_node_name - gene tree node name
        Implementation:
            - We have the HOGs that need writing plus knowledge of what scl units 
              each HOG should contain. For each hog take the intersection of what 
              we have with what the hog should contain.
        """
        ret = []
        for h in hogs_to_write:
            # print("HOG: " + h)
            q_empty = True
            # 2. We know the scl, these are the 'taxonomic units' available (clades or individual species in species tree for this node of the gene tree)
            # Note there can be at most one of each. Only a subset of these will fall under this HOG.
            units = self.hog_contents[h].intersection(genes_per_species_index.keys())
            # print("Units: " + str(units))
            genes_row = ["" for _ in self.iSps]
            # put the units into the row
            for isp in units:
                genes_row[isp] = genes_per_species_index[isp]
                q_empty = False
            if not q_empty: 
                # print((h, genes_row))
                ret.append((h, [og_name, gt_node_name] + genes_row))
                # self.writers[h].writerow(["%s.HOG%07d" % (h, self.get_hog_index(h)),  og_name, gt_node_name] + genes_row)
        return ret

    def close_files(self):
        for fh in self.fhs.values():
            fh.close()

    @staticmethod
    def get_skipped_nodes(n_sp_this, n_above_name, n_stop = None, n_gene=None):
        """
        Get the HOGs for the series of skipped species tree nodes above the current 
        node. 
        Args:
            n_sp_this - ete3 node from species tree 
            n_above_name - MRCA species tree node name for the gene tree node above
            n_stop - a HOG name above the MRCA that should be the last one added
            n_gene - ete3 node from the gene tree
        Implementation/Questions:
            - This is only used by the OGs with fewer than 4 taxa (and therefore 
              no tree) now
        """
        n = n_sp_this
        missed_sp_node_names = []
        if n.name == n_above_name:
            return missed_sp_node_names
        n = n.up
        while n is not None and n.name != n_above_name:
            missed_sp_node_names.append(n.name)
            if n.name == n_stop:
                break
            n = n.up
        # if above node is a duplication then we won't have written out a HOG for that, pass the next node 
        if n_stop is None:
            if n_gene is not None and n_gene.up is not None and n_gene.up.dup and n is not None:
                # then we also need to write the HOG above
                missed_sp_node_names.append(n.name)
        return missed_sp_node_names

    def mark_dups_below(self, tree):
        """
        Marks duplications below and at each node in feature 'dups_below'. This 
        determines the HOGs.
        Args:
            tree - ete3 gene tree with attributes 'dup' 'sp_node' for all non-terminals
        Returns
            tree - with attribute 'dups_below' on each node that is not a leaf or
                   a species-specific node and 'dup_level' on each node where n.dup=True
        """
        for n in tree.traverse('postorder'):
            if n.is_leaf():
                n.sp_node = n.name.split("_")[0]
                continue
            if n.dup:
                # get the MRCA level at which there is evidence of a duplication
                mrcas = [ch.name.split("_")[0] if ch.is_leaf() else ch.sp_node for ch in n.get_children()]
                if len(set(mrcas)) == 1 and len(mrcas) > 1:
                    n.add_feature('dup_level', mrcas[0])
                else:
                    # need two child nodes attesting to that level
                    l = self.get_evidenced_dup_level(mrcas)
                    if l is None:
                        n.dup = False
                        continue
                    n.add_feature('dup_level', l)
                # print((n.name, n.dup_level, mrcas))
            dups_below = set()
            for ch in n.get_children():
                if ch.is_leaf():
                    continue
                dups_below.update(ch.dups_below)
                # # don't care about terminal duplications
                # if ch.dup and ch.sp_node.startswith('N'):
                #     dups_below.add(ch.dup_level)
                if n.dup and n.dup_level.startswith('N'):
                    dups_below.add(n.dup_level)   # dups_below now includes current node too
            n.add_feature('dups_below', dups_below)
        return tree

    
    def get_evidenced_dup_level(self, mrcas):
        """
        Implementation
        - V3.1: Currently, we need a representative from X and Y in both ch1 & ch2.
        - What if we asked for evidence from each descendant clade of evidence of 
          a duplication?
            - Version that would be too stringent for duplication identification: 
              genetree =(ch1, ch2)n, sptree = (X,Y) and ask for a
              single representative of X that is in both ch1 & ch2 and similarly 
              for Y in ch1 & ch2.
            - Better version: V3.1 criterion (X & Y seen in ch1 & 2) plus two copies 
              of a gene from X and two copies of a gene from Y in clade n (don't 
              have to be correct topology such that they fall correctly in ch1 & 
              ch2, just evidence of duplicated genes).
                - Note, these two copies could arise from a separate & well-evidenced
                  lower duplication. The idea was that this would strike the right
                  balance, but can we do better? We'd have to identify the genes
                  below that go into each component of the interpretation of the 
                  tree. This is for another time, too much for now.

        - a. get MRCA for ch1 & ch2 and then the lower of the two
        - b. get all multi-copy species, get MRCA
        - Get lower of a & b
        """
        # try each in turn and see if it is supported
        attested = set()
        for l1, l2 in itertools.combinations(mrcas, 2):
            if l1 == l2:
                attested.add(l1)
            elif l2 in self.comp_nodes[l1][1]:
                attested.add(l2)
            elif l1 in self.comp_nodes[l2][1]:
                attested.add(l1)
        if len(attested) == 1:
            return attested.pop()
        elif len(attested) == 0:
            print("WARNING: Unexpected gene tree topology")
            print(mrcas)
            # raise Exception()
            return None
        else:
            # get the highest in the tree
            attested = list(attested)
            ancestor_lists = [(self.species_tree & a).get_ancestors() for a in attested]
            x = len(attested)
            for i in range(x):
                if all(attested[i] in ancestor_lists[j] for j in range(x) if j!=i):
                    return attested[i]
        print("WARNING: Unexpected gene tree topology 2")
        print(mrcas)
        # raise Exception()
        return None

    def WriteCachedHOGs(self, cached_hogs, lock_hogs):
        d = defaultdict(list)
        for h, row in cached_hogs:
            d[h].append(row)
        lock_hogs.acquire()
        try:
            for h, hog_rows in d.items():
                fh = self.fhs[h]
                for r in hog_rows:
                    hog_id = "%s.HOG%07d" % (h, self.get_hog_index(h))
                    util.writerow(fh, [hog_id, ] + r)
                fh.flush()
        finally:
            lock_hogs.release()

    @staticmethod
    def scl_fn(n):
        return n.is_leaf() or n.dup


def GetHOGs_from_tree(iog, tree, hog_writer, lock_hogs, q_split_paralogous_clades):
    og_name = "OG%07d" % iog
    if debug: print("\n===== %s =====" % og_name)
    try:
        tree = hog_writer.mark_dups_below(tree)
        cached_hogs = []
        for n in tree.traverse("preorder"):
            cached_hogs.extend(hog_writer.write_clade_v2(n, og_name, q_split_paralogous_clades))
        hog_writer.WriteCachedHOGs(cached_hogs, lock_hogs)
    except:
        print("WARNING: HOG analysis for %s failed" % og_name)
        print("Please report to https://github.com/davidemms/OrthoFinder/issues including \
SpeciesTree_rooted_ids.txt and Trees_ids/%s_tree_id.txt from WorkingDirectory/" % og_name)
        print(cached_hogs)
        raise


def get_highest_nodes(nodes, comp_nodes):
    """
    Returns the nodes closest to the root
    Args:
        nodes - the list of nodes to examine
        comp_nodes - dict:NX -> ( {closer to root}, {further from root} )
    """
    return  {n for n in nodes if not any(n in comp_nodes[n2][1] for n2 in nodes)}


def get_comparable_nodes(sp_tree):
    """
    Return a dictionary of comaprable nodes
    Node NX < NY if NX is on the path between NY and the root.
    If a node is not <, =, > another then they are incomparable
    Args:
        sp_tree - sp_tree with labelled nodes
    Returns:
        comp_nodes - dict:NX -> ( {n|n<NX}, {n|n>NX} ) i.e. (higher_nodes, lower_nodes)
    """
    comp_nodes = dict()
    for n in sp_tree.traverse('postorder'):
        nodes_below = set()
        if not n.is_leaf():
            for ch in n.get_children():
                if not ch.is_leaf():
                    nodes_below.update(ch.nodes_below)
                nodes_below.add(ch.name)
        above = set([nn.name for nn in n.get_ancestors()])
        n.add_feature('nodes_below', nodes_below)
        comp_nodes[n.name] = (above, nodes_below, above.union(nodes_below.union(set(n.name))))
    return comp_nodes

"""
Orthologs
-------------------------------------------------------------------------------
""" 

def OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, n1, n2):
    f_dup = len(sp_up.intersection(sett1)) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * len(sp_down.intersection(sett2)) * N_recip
    f_a = len(sp_up.intersection(sett1)) * (n2-len(sp_up.intersection(sett2))) * (n1-len(sp_down.intersection(sett1))) * len(sp_down.intersection(sett2)) * N_recip
    f_b = (n1-len(sp_up.intersection(sett1))) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * (n2-len(sp_down.intersection(sett2))) * N_recip
    choice = (f_dup, f_a, f_b)
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
        return [next(n for n in tree)] # arbitrary root if all genes are from the same species
    
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
                clades = [ch.inout_down for ch in nodes] if m.is_root() else ([m.inout_up] + [ch.inout_down for ch in m.get_children()])
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
    return [sorted(zip(scores_list, roots_list), key=lambda x: x[0], reverse=True)[0][1]]
                
def WriteQfO2(orthologues_list_pairs_list, outfilename, qAppend = True):
    """ takes a list where each entry is a pair, (genes1, genes2), which are orthologues of one another
    """
    with open(outfilename, 'a' if qAppend else 'w') as outfile:
        for gs1, gs2, _, _ in orthologues_list_pairs_list:
            for sp1, genes1 in gs1.items():
                for sp2, genes2 in gs2.items():
                    for g1 in genes1:
                        for g2 in genes2:
                            outfile.write("%s_%s\t%s_%s\n" % (sp1, g1, sp2, g2))
    
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
  
def OverlapSize(node, GeneToSpecies, suspect_genes):  
    descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()}.difference(suspect_genes) for n in node.get_children()]
    intersection = descendents[0].intersection(descendents[1])
    return len(intersection), intersection, descendents[0], descendents[1]

def ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours, GeneToSpecies, relOverlapCutoff=4):
    """
    Is an overlap suspicious and if so can it be resolved by identifying genes that are out of place?
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
        - The number of species in the overlap must be a 5th or less of the number of species in each clade - What if it's a single gene that's out of place? Won't make a difference then to the orthologs!
        - for each species with genes in both clades: the genes in one clade must all be more out of place (according to the 
          species tree) than all the gene from that species in the other tree
    """
    oSize = len(overlap)
    lsp0 = len(sp0)
    lsp1 = len(sp1)
    if (oSize == lsp0 or oSize == lsp1) or (relOverlapCutoff*oSize >= lsp0 and relOverlapCutoff*oSize >= lsp1): 
        return False, []
    # The overlap looks suspect, misplaced genes?
    # for each species, we'd need to be able to determine that all genes from A or all genes from B are misplaced
    genes_removed = []
    nA_removed = 0
    nB_removed = 0
    qResolved = True
    for sp in overlap:
        A = [g for g in ch[0].get_leaf_names() if GeneToSpecies(g) == sp]
        B = [g for g in ch[1].get_leaf_names() if GeneToSpecies(g) == sp]
        A_levels = []
        B_levels = []
        for X, level in zip((A,B),(A_levels, B_levels)):
            for g in X:
                gene_node = tree & g
                r = gene_node.up
                nextSpecies = set([GeneToSpecies(gg) for gg in r.get_leaf_names()])
                # having a gene from the same species isn't enough?? No, but we add to the count I think.
                while len(nextSpecies) == 1:
                    r = r.up
                    nextSpecies = set([GeneToSpecies(gg) for gg in r.get_leaf_names()])
                nextSpecies.remove(sp)
                # get the level
                # the sum of the closest and furthest expected distance topological distance for the closest genes in the gene tree (based on species tree topology)
                neigh = neighbours[sp]
                observed = [neigh[nSp] for nSp in nextSpecies]
                level.append(min(observed) + max(observed))
        qRemoveA = max(B_levels) + 2 < min(A_levels)   # if the clade is one step up the tree further way (min=max) then this gives +2. There's no way this is a problem                        
        qRemoveB = max(A_levels) + 2 < min(B_levels)                           
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
        return True, set(genes_removed)
    else:
        return False, set()
          
def GetRoot(tree, species_tree_rooted, GeneToSpecies):
        roots = GetRoots(tree, species_tree_rooted, GeneToSpecies)
        if len(roots) > 0:
            root_dists = [r.get_closest_leaf()[1] for r in roots]
            i, _ = max(enumerate(root_dists), key=operator.itemgetter(1))
            return roots[i]
        else:
            return None # single species tree

def CheckAndRootTree(treeFN, species_tree_rooted, GeneToSpecies):
    """
    Check that the tree can be analysed and rooted
    Root tree
    Returns None if this fails, i.e. checks: exists, has more than one gene, can be rooted
    """
    if (not os.path.exists(treeFN)) or os.stat(treeFN).st_size == 0: return None, False
    qHaveSupport = False
    try:
        tree = tree_lib.Tree(treeFN, format=2)
        qHaveSupport = True
    except:
        try:
            tree = tree_lib.Tree(treeFN)
        except:
            tree = tree_lib.Tree(treeFN, format=3)
    if len(tree) == 1: return None, False
    root = GetRoot(tree, species_tree_rooted, GeneToSpecies)
    if root == None: return None, False
    # Pick the first root for now
    if root != tree:
        tree.set_outgroup(root)
    return tree, qHaveSupport

def Orthologs_and_Suspect(ch, suspect_genes, misplaced_genes, SpeciesAndGene):
    """
    ch - the two child nodes that are orthologous
    suspect_genes - genes already identified as misplaced at lower levels
    misplaced_genes - genes identified as misplaced at this level

    Returns the tuple (o_0, o_1, os_0, os_1) where each element is a dictionary from species to genes from that species,
    the o are orthologs, the os are 'suspect' orthologs because the gene was previously identified as suspect
    """
    d = [defaultdict(list) for _ in range(2)]
    d_sus = [defaultdict(list) for _ in range(2)] 
    for node, di, d_susi in zip(ch, d, d_sus):
        for g in [g for g in node.get_leaf_names() if g not in misplaced_genes]:
            sp, seq = SpeciesAndGene(g)
            if g in suspect_genes:
                d_susi[sp].append(seq)
            else:
                di[sp].append(seq)
    return d[0], d[1], d_sus[0], d_sus[1]


def GetOrthologues_from_tree(iog, tree, species_tree_rooted, GeneToSpecies, neighbours,
                             q_get_dups=False, qNoRecon=False):
    """ 
    Returns:
        orthologues 
        tree - Each node of the tree has two features added: dup (bool) and sp_node (str)
        suspect_genes - set
        duplications - list of (sp_node_name, genes0, genes1)
    """
    og_name = "OG%07d" % iog
    n_species = len(species_tree_rooted)
    # max_genes_dups = 5*n_species
    # max_genes_text = (">%d genes" % max_genes_dups,)
    qPrune=False
    SpeciesAndGene = SpeciesAndGene_lookup[GeneToSpecies]
    orthologues = []
    duplications = []
    if not qNoRecon: tree = Resolve(tree, GeneToSpecies)
    if qPrune: tree.prune(tree.get_leaf_names())
    if len(tree) == 1: return set(orthologues), tree, set(), duplications
    """ At this point need to label the tree nodes """
    iNode = 1
    tree.name = "n0"
    suspect_genes = set()
    empty_set = set()
    # preorder traverse so that suspect genes can be identified first, before their closer orthologues are proposed.
    # Tree resolution has already been performed using a postorder traversal
    for n in tree.traverse('preorder'):
        if n.is_leaf(): continue
        if not n.is_root():
            n.name = "n%d" % iNode
            iNode += 1
        sp_present = None
        ch = n.get_children()
        if len(ch) == 2: 
            oSize, overlap, sp0, sp1 = OverlapSize(n, GeneToSpecies, suspect_genes)
            sp_present = sp0.union(sp1)
            stNode = MRCA_node(species_tree_rooted, sp_present)
            n.add_feature("sp_node", stNode.name)
            if oSize != 0:
                # this should be moved to the tree resolution step. Except that doesn't use the species tree, so can't
                qResolved, misplaced_genes = ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours, GeneToSpecies) 
                # label the removed genes
                for g in misplaced_genes:
                    nn = tree & g
                    nn.add_feature("X", True)
            else:
                misplaced_genes = empty_set
            dup = oSize != 0 and not qResolved
            n.add_feature("dup", dup)
            if dup:
                if q_get_dups:
                    # genes0 = ch[0].get_leaf_names() if len(ch[0]) <= max_genes_dups else max_genes_text
                    # genes1 = ch[1].get_leaf_names() if len(ch[1]) <= max_genes_dups else max_genes_text
                    genes0 = ch[0].get_leaf_names()
                    genes1 = ch[1].get_leaf_names()
                    duplications.append((stNode.name, n.name, float(oSize)/(len(stNode)), genes0, genes1))
            else:
                # sort out bad genes - no orthology for all the misplaced genes at this level (misplaced_genes). 
                # For previous levels, (suspect_genes) have their orthologues written to suspect orthologues file
                orthologues.append(Orthologs_and_Suspect(ch, suspect_genes, misplaced_genes, SpeciesAndGene))
                suspect_genes.update(misplaced_genes)
        elif len(ch) > 2:
            species = [{GeneToSpecies(l) for l in n_.get_leaf_names()} for n_ in ch]
            all_species = set.union(*species)
            stNode = MRCA_node(species_tree_rooted, all_species)
            n.add_feature("sp_node", stNode.name)
            # should skip if everything below this is a single species, but should write out the duplications
            if len(all_species) == 1:
                # genes = n.get_leaf_names() if len(n) <= max_genes_dups else max_genes_text
                genes = n.get_leaf_names()
                duplications.append((stNode.name, n.name, 1., genes, []))
                n.add_feature("dup", True)  
            else:
                dups = []
                for (n0, s0), (n1, s1) in itertools.combinations(zip(ch, species), 2):
                    if len(s0.intersection(s1)) == 0:
                        orthologues.append(Orthologs_and_Suspect((n0, n1), suspect_genes, empty_set, SpeciesAndGene))
                        dups.append(False)
                    else:
                        dups.append(True)
                if all(dups):
                    # genes = n.get_leaf_names() if len(n) <= max_genes_dups else max_genes_text
                    genes = n.get_leaf_names()
                    duplications.append((stNode.name, n.name, 1., genes, []))
                n.add_feature("dup", all(dups))
                # # if there are nodes below with same MRCA then dup (no HOGs) otherwise not dup (=> HOGS at this level)
                # descendant_nodes = [MRCA_node(species_tree_rooted, sp) for sp in species]
                # dup = any(down_sp_node == stNode for down_sp_node in descendant_nodes)
                # n.add_feature("dup", dup)  
                # print(n.name + ": dup3")
    return orthologues, tree, suspect_genes, duplications

def GetLinesForOlogFiles(orthologues_alltrees, speciesDict, iSpeciesToUse, sequenceDict, 
                            qContainsSuspectOlogs, olog_lines, olog_sus_lines):
    """
    Prepare teh lines of text for the pairwise ortholog files and the species-wise
    putative xenolog files
    Args:
        orthologues_alltrees - list of tuples (iog, (leavesL, leavesR, sus_leavesL, sus_leavesR))
        where each of the leavesL etc. is a dictionary, l : isp -> (iseq0, iseq0, ...)
          
    Implementation:
    Look at the genes and organise them per species. This is in contrast to the 
    first version which starts from the species pairs and looks for genes for each
    species pair.

    """
    nOrtho = util.nOrtho_sp(len(iSpeciesToUse))   
    sp_to_index = {str(sp):i for i, sp in enumerate(iSpeciesToUse)}
    for iog, orthologues_onetree in orthologues_alltrees:                   
        og = "OG%07d" % iog
        for leavesL, leavesR, leavesL_sus, leavesR_sus in orthologues_onetree:
            # suspect_genes are the genes which, for this level, the orthologues should be considered suspect as the gene appears misplaced (at this level)
            for spL, genesL in leavesL.items():
                spL_ = spL + "_"
                iL = sp_to_index[spL]
                nL = len(genesL)
                for spR, genesR in leavesR.items():
                    if spL == spR:
                        continue
                    spR_ = spR + "_"
                    iR = sp_to_index[spR]
                    nR = len(genesR)
                    textL = ", ".join([sequenceDict[spL_ + g] for g in genesL])
                    textR = ", ".join([sequenceDict[spR_ + g] for g in genesR])
                    olog_lines[iL][iR] += util.getrow((og, textL, textR))
                    olog_lines[iR][iL] += util.getrow((og, textR, textL))
                    nOrtho.n[iL, iR] += nL
                    nOrtho.n[iR, iL] += nR
                    if nL == 1 and nR == 1:
                        nOrtho.n_121[iL, iR] += 1
                        nOrtho.n_121[iR, iL] += 1
                    elif nL == 1:
                        nOrtho.n_12m[iL, iR] += 1
                        nOrtho.n_m21[iR, iL] += nR
                    elif nR == 1:
                        nOrtho.n_m21[iL, iR] += nL
                        nOrtho.n_12m[iR, iL] += 1
                    else:
                        nOrtho.n_m2m[iL, iR] += nL
                        nOrtho.n_m2m[iR, iL] += nR
            if not qContainsSuspectOlogs: continue
            # For each gene in suspect genes, must write out for any gene on other side
            # set up the two cases and regard the suspect genes as on '0' branch.
            # Iterate through the two assignments of which is branch 0
            # print((len(leavesL_sus), len(leavesR_sus)))
            leaves_sus = (leavesL_sus, leavesR_sus)
            leaves_norm = (leavesL, leavesR)
            for iPair in range(2):
                jPair = 1-iPair
                leaves_sus0 = leaves_sus[iPair]
                leaves_sus1 = leaves_sus[jPair]
                leaves_norm1 = leaves_norm[jPair]
                # Now deal with suspect genes on branch 0, other genes on 1 (no L & R)
                all_species1 = list(leaves_norm1.keys()) + list(leaves_sus1.keys())
                for sp0, genes0 in leaves_sus0.items():
                    if len(genes0) == 0:
                        continue    # should not occur any more
                    sp0_ = sp0 + "_"
                    i0 = sp_to_index[sp0]
                    for sp1 in all_species1:
                        if sp0 == sp1:
                            continue
                        sp1_ = sp1 + "_"
                        i1 = sp_to_index[sp1]
                        # don't query for elements that might not be there. That inserts and empty list
                        # which will cause problems with the algorithm
                        genes1 = leaves_norm1[sp1] if sp1 in leaves_norm1 else []
                        if sp1 in leaves_sus1: 
                            genes1 += leaves_sus1[sp1] 
                        if len(genes1) == 0:
                            continue    # should not occur any more
                        text0 = ", ".join([sequenceDict[sp0_ + g] for g in genes0])
                        text1 = ", ".join([sequenceDict[sp1_ + g] for g in genes1])
                        olog_sus_lines[i0] += util.getrow((og, text0, text1))
                        olog_sus_lines[i1] += util.getrow((og, text1, text0))
    return nOrtho
                                      
def Resolve(tree, GeneToSpecies):
    StoreSpeciesSets(tree, GeneToSpecies)
    for n in tree.traverse("postorder"):
        tree = resolve.resolve(n, GeneToSpecies)
    return tree

def GetSpeciesNeighbours(t):
    """
    Args: t = rooted species tree
    
    Returns:
    dict: species -> species_dict, such that species_dict: other_species -> toplogical_dist 
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
    neighbours = {sp:{other:n for n,others in enumerate(lev) for other in others} for sp, lev in levels.items()}
    return neighbours

def RootAndGetOrthologues_from_tree(iog, tree_fn, species_tree_rooted, GeneToSpecies, neighbours, qWrite=False, qNoRecon=False):
    rooted_tree_ids, qHaveSupport = CheckAndRootTree(tree_fn, species_tree_rooted, GeneToSpecies) # this can be parallelised easily
    if rooted_tree_ids is None: return
    orthologues, recon_tree, suspect_genes, dups = GetOrthologues_from_tree(iog, rooted_tree_ids, species_tree_rooted, GeneToSpecies, neighbours, qNoRecon=qNoRecon)
    if qWrite:
        directory = os.path.split(tree_fn)[0]
        WriteQfO2(orthologues, directory + "_Orthologues_M3/" + os.path.split(tree_fn)[1], qAppend=False)

def GetOrthologuesStandalone_Parallel(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree):
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    neighbours = GetSpeciesNeighbours(species_tree_rooted)
    args_queue = mp.Queue()
    for treeFn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")): args_queue.put((0, treeFn, species_tree_rooted, GeneToSpecies, neighbours))
    # Now need to root the tree first
    parallel_task_manager.RunMethodParallel(RootAndGetOrthologues_from_tree, args_queue, 16)

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
        # Now need to root the tree first
        RootAndGetOrthologues_from_tree(0, treeFn, species_tree_rooted, GeneToSpecies, neighbours)        

class OrthologsFiles(object):
    """wrapper to open all the orthologs files as once"""
    def __init__(self, directory, speciesDict, iSpeciesToUse, nSpecies, sp_to_index):
        self.d = directory
        self.speciesDict = speciesDict
        self.iSpeciesToUse = iSpeciesToUse
        self.nSpecies = nSpecies
        self.sp_to_index = sp_to_index
        self.dPutativeXenologs = files.FileHandler.GetPutativeXenelogsDir()
        self.ortholog_file_handles = [[None for _ in self.iSpeciesToUse] for _ in self.iSpeciesToUse]
        self.xenolog_file_handles = [None for _ in self.iSpeciesToUse]

    def __enter__(self):
        for i in xrange(self.nSpecies):
            sp0 = str(self.iSpeciesToUse[i])
            self.xenolog_file_handles[i] = open(self.dPutativeXenologs + "%s.tsv" % self.speciesDict[sp0], csv_append_mode)
            strsp0 = sp0 + "_"
            isp0 = self.sp_to_index[sp0]
            d0 = self.d + "Orthologues_" + self.speciesDict[sp0] + "/"
            for j in xrange(i, self.nSpecies):
                sp1 = str(self.iSpeciesToUse[j])
                if sp1 == sp0: continue
                strsp1 = sp1 + "_"
                isp1 = self.sp_to_index[sp1]
                d1 = self.d + "Orthologues_" + self.speciesDict[sp1] + "/"
                self.ortholog_file_handles[i][j] = open(d0 + '%s__v__%s.tsv' % (self.speciesDict[sp0], self.speciesDict[sp1]), csv_append_mode)
                self.ortholog_file_handles[j][i] = open(d1 + '%s__v__%s.tsv' % (self.speciesDict[sp1], self.speciesDict[sp0]), csv_append_mode)
        return self.ortholog_file_handles, self.xenolog_file_handles

    def __exit__(self, type, value, traceback):
        for fh in self.xenolog_file_handles:
            fh.close()
        for fh_list in self.ortholog_file_handles:
            for fh in fh_list:
                if fh is not None:
                    fh.close()

    @staticmethod
    def flush_olog_files(ortholog_file_handles):
        for i, handles in enumerate(ortholog_file_handles):
            for j, h in enumerate(handles):
                if i != j:
                    h.flush()

    @staticmethod
    def flush_xenolog_files(files_list):
        for h in files_list:
            h.flush()
                    

def InitialiseSuspectGenesDirs(nspecies, speciesIDs, speciesDict):
    files.FileHandler.GetSuspectGenesDir()  # creates the directory
    dSuspectOrthologues = files.FileHandler.GetPutativeXenelogsDir()
    for index1 in xrange(nspecies):
        with open(dSuspectOrthologues + '%s.tsv' % speciesDict[str(speciesIDs[index1])], csv_write_mode) as outfile:
            writer1 = csv.writer(outfile, delimiter="\t")
            writer1.writerow(("Orthogroup", speciesDict[str(speciesIDs[index1])], "Other"))

def WriteSuspectGenes(nspecies, speciesToUse, suspect_genes, speciesDict, SequenceDict):
    species = list(map(str, speciesToUse))
    dSuspectGenes = files.FileHandler.GetSuspectGenesDir()
    for index0 in xrange(nspecies):
        strsp0 = species[index0]
        strsp0_ = strsp0+"_"
        these_genes = [g for g in suspect_genes if g.startswith(strsp0_)]
        if len(these_genes) > 0:
            with open(dSuspectGenes + speciesDict[strsp0] + ".txt", csv_append_mode) as outfile:
                # not a CSV file so \n line endings are fine
                outfile.write("\n".join([SequenceDict[g] for g in these_genes]) + "\n")

def WriteDuplications(dups_file_handle, og_name, duplications, spIDs, seqIDs, stride_dups):
    """
    Args:
        duplications - list of (sp_node_id, gene_node_name, fraction, genes0, genes1)
    """
    for sp_node_id, gene_node_name, frac, genes0, genes1 in duplications:
        q_terminal = not sp_node_id.startswith("N")
        if stride_dups is None:
            isSTRIDE = "Terminal" if q_terminal else "Non-Terminal"
        else:
            isSTRIDE = "Terminal" if q_terminal else "Non-Terminal: STRIDE" if frozenset(genes0 + genes1) in stride_dups else "Non-Terminal"
        gene_list0 = ", ".join([seqIDs[g] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
        gene_list1 = ", ".join([seqIDs[g] for g in genes1])
        util.writerow(dups_file_handle, [og_name, spIDs[sp_node_id] if q_terminal else sp_node_id, gene_node_name, frac, isSTRIDE, gene_list0, gene_list1]) 

def DoOrthologuesForOrthoFinder(ogSet, species_tree_rooted_labelled, GeneToSpecies, stride_dups, qNoRecon, hog_writer, q_split_paralogous_clades, n_parallel):   
    """
    """
    try:
        # Create directory structure
        speciesDict = ogSet.SpeciesDict()
        SequenceDict = ogSet.SequenceDict()
        # Write directory and file structure
        nspecies = len(ogSet.speciesToUse)      
        dResultsOrthologues = files.FileHandler.GetOrthologuesDirectory()
        for index1 in xrange(nspecies):
            d = dResultsOrthologues + "Orthologues_" + speciesDict[str(ogSet.speciesToUse[index1])] + "/"
            if not os.path.exists(d): os.mkdir(d)     
            for index2 in xrange(nspecies):
                if index2 == index1: continue
                with open(d + '%s__v__%s.tsv' % (speciesDict[str(ogSet.speciesToUse[index1])], speciesDict[str(ogSet.speciesToUse[index2])]), csv_write_mode) as outfile:
                    writer1 = csv.writer(outfile, delimiter="\t")
                    writer1.writerow(("Orthogroup", speciesDict[str(ogSet.speciesToUse[index1])], speciesDict[str(ogSet.speciesToUse[index2])]))
        InitialiseSuspectGenesDirs(nspecies, ogSet.speciesToUse, speciesDict)
        neighbours = GetSpeciesNeighbours(species_tree_rooted_labelled)
        nOgs = len(ogSet.OGs()) 
        reconTreesRenamedDir = files.FileHandler.GetOGsReconTreeDir(True)
        spec_seq_dict = ogSet.Spec_SeqDict()
        sp_to_index = {str(sp):i for i, sp in enumerate(ogSet.speciesToUse)}

        # Infer orthologues and write them to file           
        with open(files.FileHandler.GetDuplicationsFN(), csv_write_mode) as outfile_dups, OrthologsFiles(dResultsOrthologues, speciesDict, ogSet.speciesToUse, nspecies, sp_to_index) as (ologs_file_handles, putative_xenolog_file_handles):
            util.writerow(outfile_dups, ["Orthogroup", "Species Tree Node", "Gene Tree Node", "Support", "Type",	"Genes 1", "Genes 2"])
            outfile_dups.flush()
            OrthologsFiles.flush_olog_files(ologs_file_handles)
            ta = TreeAnalyser(nOgs, dResultsOrthologues, reconTreesRenamedDir, species_tree_rooted_labelled, 
                              ogSet.speciesToUse, GeneToSpecies, SequenceDict, speciesDict, spec_seq_dict, 
                              neighbours, qNoRecon, outfile_dups, stride_dups, ologs_file_handles, 
                              putative_xenolog_file_handles, hog_writer, q_split_paralogous_clades)
            
            if n_parallel == 1:
                nOrthologues_SpPair = util.nOrtho_sp(nspecies)
                dummy_lock = mp.Lock()
                for iog in range(nOgs):
                    results = ta.AnalyseTree(iog) 
                    if results is None:
                        continue
                    nOrthologues_this, olog_lines, olog_sus_lines = results
                    if nOrthologues_this.n.sum() == 0 and sum(map(len, olog_sus_lines)) == 0: 
                        continue
                    nOrthologues_SpPair += nOrthologues_this
                    for i in range(nspecies):
                        for j in range(i+1, nspecies):
                            if len(olog_lines[i][j]) > 0:
                                # j is the largest (and intentionally changing quickest, which I think is best for the lock)
                                WriteOlogLinesToFile(ta.ologs_files_handles[i][j], olog_lines[i][j], dummy_lock)
                                WriteOlogLinesToFile(ta.ologs_files_handles[j][i], olog_lines[j][i], dummy_lock)
                        WriteOlogLinesToFile(ta.putative_xenolog_file_handles[i], olog_sus_lines[i], dummy_lock)
                # util.PrintTime("Done writing orthologs")
            else:
                args_queue = mp.Queue()
                for iog in range(nOgs):
                    args_queue.put(iog)
                nOrthologues_SpPair = RunOrthologsParallel(ta, len(ogSet.speciesToUse), args_queue, n_parallel)
    except IOError as e:
        if str(e).startswith("[Errno 24] Too many open files"):
            util.number_open_files_exception_advice(len(ogSet.speciesToUse), True)
            util.Fail()
        else:
            raise
    return nOrthologues_SpPair


class TreeAnalyser(object):
    def __init__(self, nOgs, dResultsOrthologues, reconTreesRenamedDir, species_tree_rooted_labelled, 
                speciesToUse, GeneToSpecies, SequenceDict, speciesDict, spec_seq_dict, 
                neighbours, qNoRecon, dups_file_handle, stride_dups, ologs_files_handles, 
                putative_xenolog_file_handles, hog_writer, q_split_paralogous_clades):
        self.nOgs = nOgs
        self.dResultsOrthologues = dResultsOrthologues
        self.reconTreesRenamedDir = reconTreesRenamedDir
        self.species_tree_rooted_labelled = species_tree_rooted_labelled
        self.speciesToUse = speciesToUse
        self.nspecies = len(self.speciesToUse)
        self.GeneToSpecies = GeneToSpecies
        self.SequenceDict = SequenceDict
        self.speciesDict = speciesDict
        self.spec_seq_dict = spec_seq_dict
        self.spec_seq_dict[">%d genes" % (5*self.nspecies)] = ">%d genes" % (5*self.nspecies)   # used in Duplications
        self.neighbours = neighbours
        self.qNoRecon = qNoRecon
        self.dups_file_handle = dups_file_handle
        self.stride_dups = stride_dups
        self.ologs_files_handles = ologs_files_handles
        self.putative_xenolog_file_handles = putative_xenolog_file_handles
        self.hog_writer = hog_writer
        self.q_split_paralogous_clades = q_split_paralogous_clades
        self.lock_ologs = [mp.Lock() for i in range(self.nspecies)]   # lock the larger of the two species index
        self.lock_dups = mp.Lock()
        self.lock_suspect = mp.Lock()
        self.lock_hogs = mp.Lock()

    def AnalyseTree(self, iog):
        try:
            if not os.path.exists(files.FileHandler.GetOGsTreeFN(iog)):
                return None
            og_name = "OG%07d" % iog
            n_species = len(self.speciesToUse)
            rooted_tree_ids, qHaveSupport = CheckAndRootTree(files.FileHandler.GetOGsTreeFN(iog), self.species_tree_rooted_labelled, self.GeneToSpecies) # this can be parallelised easily
            if rooted_tree_ids is None: 
                return None

            # Write rooted tree with accessions
            util.RenameTreeTaxa(rooted_tree_ids, files.FileHandler.GetOGsTreeFN(iog, True), 
                                self.spec_seq_dict, qSupport=qHaveSupport, qFixNegatives=True, qViaCopy=True)
            ologs, recon_tree, suspect_genes, dups = GetOrthologues_from_tree(iog, rooted_tree_ids, 
                                                        self.species_tree_rooted_labelled, self.GeneToSpecies, 
                                                        self.neighbours, q_get_dups=True, qNoRecon=self.qNoRecon)
            # Write Duplications
            self.lock_dups.acquire()
            try:                                                        
                WriteDuplications(self.dups_file_handle, og_name, dups, self.speciesDict, self.spec_seq_dict, self.stride_dups)
                self.dups_file_handle.flush()
            finally:
                self.lock_dups.release()
            # Write Suspect Genes
            if len(suspect_genes) > 0:
                self.lock_suspect.acquire()
                try:
                    WriteSuspectGenes(n_species, self.speciesToUse, suspect_genes, self.speciesDict, self.SequenceDict)
                finally:
                    self.lock_suspect.release()

            # Get Orthologues
            olog_lines = [["" for j in xrange(self.nspecies)] for i in xrange(self.nspecies)]
            olog_sus_lines = ["" for i in xrange(self.nspecies)]
            nOrthologues_SpPair = GetLinesForOlogFiles([(iog, ologs)], self.speciesDict, self.speciesToUse,
                                                    self.SequenceDict, len(suspect_genes) > 0, olog_lines, olog_sus_lines)
            GetHOGs_from_tree(iog, recon_tree, self.hog_writer, self.lock_hogs, self.q_split_paralogous_clades) 
            # don't relabel nodes, they've already been done
            util.RenameTreeTaxa(recon_tree, self.reconTreesRenamedDir + "OG%07d_tree.txt" % iog, self.spec_seq_dict, qSupport=False, qFixNegatives=True)
            if iog >= 0 and divmod(iog, 10 if self.nOgs <= 200 else 100 if self.nOgs <= 2000 else 1000)[1] == 0:
                util.PrintTime("Done %d of %d" % (iog, self.nOgs))
            return nOrthologues_SpPair, olog_lines, olog_sus_lines
        except:
            print("WARNING: Unknown error analysing tree %s" % og_name)
            olog_lines = [["" for j in xrange(self.nspecies)] for i in xrange(self.nspecies)]
            olog_sus_lines = ["" for i in xrange(self.nspecies)]
            return util.nOrtho_sp(n_species), olog_lines, olog_sus_lines

def Worker_RunOrthologsMethod(tree_analyser, nspecies, args_queue, results_queue, n_ologs_cache=100):
    try:
        nOrthologues_SpPair = util.nOrtho_sp(nspecies) 
        nCache = util.nOrtho_cache(nspecies) 
        olog_lines_tot = [["" for j in range(nspecies)] for i in range(nspecies)]
        olog_sus_lines_tot = ["" for i in range(nspecies)]
        while True:
            try:
                iog = args_queue.get(True, 1.)
                results = tree_analyser.AnalyseTree(iog)
                if results is None:
                    continue
                nOrtho, olog_lines, olog_sus_lines = results
                nOrthologues_SpPair += nOrtho
                nCache += nOrtho
                for i in range(nspecies):
                    olog_sus_lines_tot[i] += olog_sus_lines[i]
                    for j in range(nspecies):
                        olog_lines_tot[i][j] += olog_lines[i][j]
                # Now write those that we've collected enough lines for
                I,J = nCache.get_i_j_to_write(n_ologs_cache)
                for i, j in zip(I,J):
                    k_lock = max(i,j)
                    WriteOlogLinesToFile(tree_analyser.ologs_files_handles[i][j], olog_lines_tot[i][j], tree_analyser.lock_ologs[k_lock])
                    olog_lines_tot[i][j] = ""
            except parallel_task_manager.queue.Empty:
                for i in range(nspecies):
                    for j in range(i+1, nspecies):
                        # j is the largest (and intentionally changing quickest, which I think is best for the lock)
                        WriteOlogLinesToFile(tree_analyser.ologs_files_handles[i][j], olog_lines_tot[i][j], tree_analyser.lock_ologs[j])
                        WriteOlogLinesToFile(tree_analyser.ologs_files_handles[j][i], olog_lines_tot[j][i], tree_analyser.lock_ologs[j])
                    WriteOlogLinesToFile(tree_analyser.putative_xenolog_file_handles[i], olog_sus_lines_tot[i], tree_analyser.lock_suspect)
                results_queue.put(nOrthologues_SpPair)
                return
            except:
                print("WARNING: Unknown error, current orthogroup OG%07d" % iog)
                for i in range(nspecies):
                    for j in range(i+1, nspecies):
                        # j is the largest (and intentionally changing quickest, which I think is best for the lock)
                        WriteOlogLinesToFile(tree_analyser.ologs_files_handles[i][j], olog_lines_tot[i][j], tree_analyser.lock_ologs[j])
                        WriteOlogLinesToFile(tree_analyser.ologs_files_handles[j][i], olog_lines_tot[j][i], tree_analyser.lock_ologs[j])
                    WriteOlogLinesToFile(tree_analyser.putative_xenolog_file_handles[i], olog_sus_lines_tot[i], tree_analyser.lock_suspect)
                results_queue.put(nOrthologues_SpPair)
                return
    except Exception as e:
        print(e)
        print("WARNING: Unexpected error")
        return util.nOrtho_sp(nspecies) 
    return

def RunOrthologsParallel(tree_analyser, nspecies, args_queue, nProcesses):
    results_queue = mp.Queue()
    runningProcesses = [mp.Process(target=Worker_RunOrthologsMethod, args=(tree_analyser, nspecies, args_queue, results_queue)) for i_ in range(nProcesses)]
    for proc in runningProcesses:
        proc.start()
    nOrthologues_SpPair = util.nOrtho_sp(nspecies)
    n_remain = nProcesses
    while True:
        try:
            nOrtho = results_queue.get()  # block until an item is available
            nOrthologues_SpPair += nOrtho
            n_remain -= 1
            if n_remain == 0:
                break
        except parallel_task_manager.queue.Empty: 
            break
    return nOrthologues_SpPair


def WriteOlogLinesToFile(fh, text, lock):
    if len(text) == 0:
        return
    if debug: util.PrintTime("Waiting: %d" % os.getpid())
    lock.acquire()
    try:
        if debug: util.PrintTime("Acquired lock: %d" % os.getpid())
        fh.write(text)
        fh.flush()
    finally:
        lock.release()
        if debug: util.PrintTime("Released lock: %d" %  os.getpid())

def SortParallelFiles(n_parallel, speciesToUse, speciesDict):
    """
    Sort the lines by orthogroup for the files that were written in parallel:
    - Orthologs
    - Xenologs
    - Duplications
    - HOGs
    """
    species = [speciesDict[str(sp1)] for sp1 in speciesToUse]
    # Orthologs
    dResultsOrthologues = files.FileHandler.GetOrthologuesDirectory()
    ds = [dResultsOrthologues + "Orthologues_" + sp1 + "/"  for sp1 in species]
    fns_type = [(d + '%s__v__%s.tsv' % (sp1, sp2), "o") for sp1, d in zip(species, ds) for sp2 in species if sp1 != sp2]
    # Xenologs 
    dXenologs = files.FileHandler.GetPutativeXenelogsDir()
    fns_type.extend([(dXenologs + '%s.tsv' % sp1, "x") for sp1 in species])
    # HOGs
    fns_type.extend([(fn, "h") for fn in glob.glob(os.path.dirname(files.FileHandler.GetHierarchicalOrthogroupsFN("N0.tsv")) + "/*")])
    # Duplications
    fns_type.append((files.FileHandler.GetDuplicationsFN(), "d"))
    args_queue = mp.Queue()
    for x in fns_type:
        args_queue.put(x)
    parallel_task_manager.RunMethodParallel(SortFile, args_queue, n_parallel)

def SortFile(fn, f_type):
    """
    Sort the contents of the file by the first column and save back to the same filename.
    Args:
        fn - filename
        f_type - o, x, h or d for orthologs, xenologs, hogs or duplications
    """
    first_column_sort = lambda s : s.split("\t", 1)[0]
    if f_type == "h":
        # Need to renumber the hogs as the parallel numbering is incorrect
        with open(fn, csv_read_mode) as infile:
            try:
                header = next(infile)
            except StopIteration:
                return
            lines = []
            # remove incorrect HOG numbering
            for line in infile:
                lines.append(line.split("\t", 1)[-1])
            if len(lines) == 0:
                return
            hog_base = line.split(".", 1)[0]
        lines.sort(key = first_column_sort)
        with open(fn, csv_write_mode) as outfile:
            outfile.write(header)
            for ihog, l in enumerate(lines):
                outfile.write(hog_base + (".HOG%07d" % ihog) + "\t" + l)
    else:
        with open(fn, csv_read_mode) as infile:
            try:
                header = next(infile)
            except StopIteration:
                return
            lines = list(infile)
        lines.sort(key=first_column_sort)
        with open(fn, csv_write_mode) as outfile:
            outfile.write(header)
            outfile.write("".join(lines))

def GetOrthologues_from_phyldog_tree(iog, treeFN, GeneToSpecies, qWrite=False, dupsWriter=None, seqIDs=None, spIDs=None):
    """ if dupsWriter != None then seqIDs and spIDs must also be provided"""
    empty = set()
    orthologues = []
    if (not os.path.exists(treeFN)) or os.stat(treeFN).st_size == 0: return set(orthologues)
    tree = tree_lib.Tree(treeFN)
    if len(tree) == 1: return set(orthologues)
    """ At this point need to label the tree nodes """
    leaf_labels = dict()
    empty_dict = dict()
    for n in tree.traverse('preorder'):
        if n.is_leaf(): 
            leaf_labels[n.name] = ("n" + n.ND)
            continue
        else:
            n.name = "n" + n.ND
        ch = n.get_children()
        if len(ch) == 2:       
            oSize, overlap, sp0, sp1 = OverlapSize(n, GeneToSpecies, empty)
            if n.Ev == "D":
                if dupsWriter != None:
                    sp_present = sp0.union(sp1)
                    stNode = "N" + n.S
                    if len(sp_present) == 1:
                        isSTRIDE = "Terminal"
                    else:
                        isSTRIDE = "Non-Terminal"
                    dupsWriter.writerow(["OG%07d" % iog, spIDs[stNode] if len(stNode) == 1 else stNode, n.name, "-", isSTRIDE, ", ".join([seqIDs[g] for g in ch[0].get_leaf_names()]), ", ".join([seqIDs[g] for g in ch[1].get_leaf_names()])]) 
            else:
                d0 = defaultdict(list)
                for g in ch[0].get_leaf_names():
                    sp, seq = g.split("_")
                    d0[sp].append(seq)
                d1 = defaultdict(list)
                for g in ch[1].get_leaf_names():
                    sp, seq = g.split("_")
                    d1[sp].append(seq)
                orthologues.append((d0, d1, empty_dict, empty_dict))
        elif len(ch) > 2:
            print("Non-binary node")
            print((n.get_leaf_names()))
    if qWrite:
        directory = os.path.split(treeFN)[0]
        WriteQfO2(orthologues, directory + "/../Orthologues_M3/" + os.path.split(treeFN)[1], qAppend=False)
    return orthologues
    
def DoOrthologuesForOrthoFinder_Phyldog(ogSet, workingDirectory, GeneToSpecies, output_dir, reconTreesRenamedDir):    # Create directory structure
    speciesDict = ogSet.SpeciesDict()
    SequenceDict = ogSet.SequenceDict()
    # Write directory and file structure
    speciesIDs = ogSet.speciesToUse
    nspecies = len(speciesIDs)      
    for index1 in range(nspecies):
        d = output_dir + "Orthologues_" + speciesDict[str(speciesIDs[index1])] + "/"
        if not os.path.exists(d): os.mkdir(d)     
        for index2 in range(nspecies):
            if index2 == index1: continue
            with open(d + '%s__v__%s.tsv' % (speciesDict[str(speciesIDs[index1])], speciesDict[str(speciesIDs[index2])]), csv_write_mode) as outfile:
                writer1 = csv.writer(outfile, delimiter="\t")
                writer1.writerow(("Orthogroup", speciesDict[str(speciesIDs[index1])], speciesDict[str(speciesIDs[index2])]))
    with OrthologsFiles(dResultsOrthologues, speciesDict, ogSet.speciesToUse, nSpecies, sp_to_index) as ologs_files_handles, putative_xenolog_file_handles:
        nOgs = len(ogSet.OGs())
        nOrthologues_SpPair = util.nOrtho_sp(nspecies) 
        with open(files.FileHandler.GetDuplicationsFN(), csv_write_mode) as outfile:
            dupWriter = csv.writer(outfile, delimiter="\t")
            dupWriter.writerow(["Orthogroup", "Species Tree Node", "Gene Tree Node", "Support", "Type",	"Genes 1", "Genes 2"])
            for iog in range(nOgs):
                recon_tree = files.FileHandler.GetPhyldogOGResultsTreeFN(iog)
                orthologues = GetOrthologues_from_phyldog_tree(iog, recon_tree, GeneToSpecies, dupsWriter=dupWriter, seqIDs=ogSet.Spec_SeqDict(), spIDs=ogSet.SpeciesDict())
                allOrthologues = [(iog, orthologues)]
                util.RenameTreeTaxa(recon_tree, reconTreesRenamedDir + "OG%07d_tree.txt" % iog, ogSet.Spec_SeqDict(), qSupport=False, qFixNegatives=True, label='n') 
                if iog >= 0 and divmod(iog, 10 if nOgs <= 200 else 100 if nOgs <= 2000 else 1000)[1] == 0:
                    util.PrintTime("Done %d of %d" % (iog, nOgs))
                nOrthologues_SpPair += GetLinesForOlogFiles(allOrthologues, speciesDict, ogSet.speciesToUse, SequenceDict, False)
    return nOrthologues_SpPair
    
def RootAllTrees():
    speciesIDs = util.FirstWordExtractor("SpeciesIDs.txt").GetIDToNameDict()
    species_tree_rooted = tree.Tree("SpeciesTree_ids_0_rooted_unresolved.txt")
    GeneToSpecies = GeneToSpecies_dash
    for fn in glob.glob("Trees_ids/OG*txt"):
        print(("*** " + fn + " ***"))
        t = tree.Tree(fn)
        root = GetRoot(t, species_tree_rooted, GeneToSpecies)
        if root == None: 
            print(("Fail: " + fn))
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
    output_dir = output_dir + "_Orthologues_M3/"
    print(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    GetOrthologuesStandalone_Parallel(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree)
    # GetOrthologuesStandalone_Serial(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree)
#    RootTreeStandalone_Serial(args.trees_dir, args.rooted_species_tree, GeneToSpecies, output_dir, qSingleTree)
    util.Success()

