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
from collections import defaultdict

from . import tree as tree_lib
from . import resolve, util, files, parallel_task_manager

PY2 = sys.version_info <= (3,)
csv_write_mode = 'wb' if PY2 else 'wt'
csv_write_mode = 'ab' if PY2 else 'at'

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

def MRCA_node(t_rooted, taxa):
    return (t_rooted & next(taxon for taxon in taxa)) if len(taxa) == 1 else t_rooted.get_common_ancestor(taxa)

class HogWriter(object):
    def __init__(self, species_tree, species_tree_node_names, seq_ids, sp_ids):
        """
        Prepare files, get ready to write.
        species_tree_node_names - list of species tree nodes
        seq_ids - dict of sequence ids
        sp_ids - dict of species ids
        """
        self.seq_ids = seq_ids
        d = os.path.dirname(files.FileHandler.GetHierarchicalOrthogroupsFN("N0"))
        if not os.path.exists(d):
            os.mkdir(d)
        self.fhs = dict()
        self.writers = dict()
        self.iSps = list(sp_ids.keys())
        self.i_sp_to_index = {isp:i for i, isp in enumerate(self.iSps)}
        self.iHOG = defaultdict(int)
        self.species_tree = species_tree
        species_names = [sp_ids[i] for i in self.iSps]
        for name in species_tree_node_names:
            fn = files.FileHandler.GetHierarchicalOrthogroupsFN(name)
            self.fhs[name] = open(fn, csv_write_mode)
            self.writers[name] = csv.writer(self.fhs[name], delimiter="\t")
            self.writers[name].writerow(["HOG", "OG", "Gene Tree Parent Clade"] + species_names)
        # Map from HOGs to genes that must be contained in them
        self.hog_contents = dict()  # sp_node_name = hog_name-> list of contents fo hog (internal nodes and leaves)
        for n in species_tree.traverse():
            desc = n.get_descendants()
            self.hog_contents[n.name] = set([nn.name for nn in desc])

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
            self.writers[sp_node_name].writerow(["%s.HOG%07d" % (sp_node_name, i_hog),  og_name, "-"] + row_genes)

    def write_clade(self, n, og_name):
        """
        Write all the relevant HOGs for all the genes in the clade down to, but not beyond, the duplication nodes
        sp_node - the MRCA species tree node of the genes in the clade
        scl - List of 'Single-copy leaves', i.e. leaf or duplication nodes. They will each have the features 'dup' and 'sp_node'
        """
        # 0. Get the scl units below this node in the gene tree
        if n.is_leaf(): return
        scl = n.get_leaves(is_leaf_fn=self.scl_fn) # single-copy 'leaf'
        scl_units = dict()
        for unit in scl:
            # get all the genes from under here and sort them into a text string per species
            genes_per_species = defaultdict(list) # iCol (before 'name' columns) -> text string of genes
            genes = unit.get_leaves()
            q_have_legitimate_gene = False # may be misplaced genes
            for g in genes:
                if "X" in g.features: continue
                q_have_legitimate_gene = True
                isp = g.name.split("_")[0]
                genes_per_species[self.i_sp_to_index[isp]].append(self.seq_ids[g.name])
            if q_have_legitimate_gene:
                scl_units[unit.name.split("_")[0] if unit.is_leaf() else unit.sp_node] = {k:", ".join(v) for k,v in genes_per_species.items()}
                r = unit.name.split("_")[0] if unit.is_leaf() else unit.sp_node
        # 1. Get HOGs to write
        scl_mrca = {nn.sp_node for nn in scl if not nn.is_leaf()}
        sp_node_name = n.sp_node
        sp_node = self.species_tree & sp_node_name
        stop_at_dups = lambda nn : nn.name in scl_mrca
        hogs_to_write = [nn.name for nn in sp_node.traverse('preorder', is_leaf_fn = stop_at_dups) if not nn.is_leaf()] 
        hogs_to_write += self.get_skipped_nodes(sp_node, n.up.sp_node if n.up is not None else None, n)
        for h in hogs_to_write:
            q_empty = True
            # 2. We know the scl, these are the 'taxonomic units' available (clades or individual species in species tree for this node of the gene tree)
            # Note there can be at most one of each. Only a subset of these will fall under this HOG.
            units = self.hog_contents[h].intersection(scl_units.keys())
            genes_row = ["" for _ in self.iSps]
            # put the units into the row
            for u in units:
                for isp, genes_text in scl_units[u].items():
                    genes_row[isp] = genes_text
                    q_empty = False
            if not q_empty: 
                self.writers[h].writerow(["%s.HOG%07d" % (h, self.get_hog_index(h)),  og_name, n.name] + genes_row)

    def close_files(self):
        for fh in self.fhs.values():
            fh.close()

    @staticmethod
    def get_skipped_nodes(n_this, n_above_name, n_gene=None):
        """
        Write HOGs for the series of skipped species tree nodes
        Args:
            n_this - ete3 node from species tree 
            n_above - MRCA species tree node for the gene tree node above
            n_gene - ete3 node from the species tree
        """
        n = n_this
        missed_sp_node_names = []
        if n.name == n_above_name:
            return missed_sp_node_names
        n = n.up
        while n is not None and n.name != n_above_name:
            missed_sp_node_names.append(n.name)
            n = n.up
        # if above node is a duplication then we won't have written out a HOG for that, pass the next node 
        if n_gene is not None and n_gene.up is not None and n_gene.up.dup and n is not None:
            # then we also need to write the HOG above
            missed_sp_node_names.append(n.name)
        return missed_sp_node_names

    @staticmethod
    def scl_fn(n):
        return n.is_leaf() or n.dup

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

def GetHOGs_from_tree(iog, tree, hog_writer):
    """
    Implementation:
    - In orthologs pass, mark all duplication events (this could probably be done at the resolve step saving some time, but the tree nodes are being shifted around at this point so potentially hazardous)
    - In the HOGs pass, each duplication node is a stopping point, it is treated as an indivisible unit
    - Traverse from the root, to all duplication nodes- we look at the cldaes for the children of duplication nodes. Each duplication node is mapped to it's place on the species tree
        *** What if this conflicts with the topology of the species tree. I think it can't come strictly above a sister node 
            without a duplication node being involved, which stops us from doing anything incorrect.
    - At each iteration of the traverse, use the species tree to write out all the applicable HOGs, treating each leaf and
      each duplication event as an indivisible unit

      - See latex doc.
    """
    og_name = "OG%07d" % iog
    # First process the root (it will be processed below if it's a dup)
    if not (tree.dup and tree.sp_node == "N0"): 
        hog_writer.write_clade(tree, og_name)
    for n in tree.iter_search_nodes(dup=True):
        nodes_to_process = n.get_children()
        for n in nodes_to_process:
            hog_writer.write_clade(n, og_name)

def GetOrthologues_from_tree(iog, tree, species_tree_rooted, GeneToSpecies, neighbours, dupsWriter=None, seqIDs=None, spIDs=None, all_stride_dup_genes=None, qNoRecon=False):
    """ if dupsWriter != None then seqIDs and spIDs must also be provided

    Each node of the tree has two features added: dup (bool) and sp_node (str)
    """
    og_name = "OG%07d" % iog
    qPrune=True
    SpeciesAndGene = SpeciesAndGene_lookup[GeneToSpecies]
    orthologues = []
    if not qNoRecon: tree = Resolve(tree, GeneToSpecies)
    if qPrune: tree.prune(tree.get_leaf_names())
    if len(tree) == 1: return set(orthologues), tree, set()
    """ At this point need to label the tree nodes """
    iNode = 1
    tree.name = "n0"
    suspect_genes = set()
    empty_set = set()
    ilab = 0
    # preorder traverse so that suspect genes can be identified first, before their closer ortholgoues are proposed.
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
            ilab += 1
            if dup:
                if dupsWriter != None:
                    if len(sp_present) == 1:
                        isSTRIDE = "Terminal"
                    else:
                        isSTRIDE = "Non-Terminal" if all_stride_dup_genes == None else "Non-Terminal: STRIDE" if frozenset(n.get_leaf_names()) in all_stride_dup_genes else "Non-Terminal"
                    dupsWriter.writerow([og_name, spIDs[stNode.name] if len(stNode) == 1 else stNode.name, n.name, float(oSize)/(len(stNode)), isSTRIDE, ", ".join([seqIDs[g] for g in ch[0].get_leaf_names()]), ", ".join([seqIDs[g] for g in ch[1].get_leaf_names()])]) 
            else:
                # sort out bad genes - no orthology for all the misplaced genes at this level (misplaced_genes). 
                # For previous levels, (suspect_genes) have their orthologues written to suspect orthologues file
                orthologues.append(Orthologs_and_Suspect(ch, suspect_genes, misplaced_genes, SpeciesAndGene))
                suspect_genes.update(misplaced_genes)
        elif len(ch) > 2:
            species = [{GeneToSpecies(l) for l in n_.get_leaf_names()} for n_ in ch]
            stNode = MRCA_node(species_tree_rooted, set.union(species))
            n.add_feature("sp_node", stNode.name)
            n.add_feature("dup", False)
            for (n0, s0), (n1, s1) in itertools.combinations(zip(ch, species), 2):
                if len(s0.intersection(s1)) == 0:
                    orthologues.append(Orthologs_and_Suspect((n0, n1), suspect_genes, empty_set, SpeciesAndGene))
    return orthologues, tree, suspect_genes

def AppendOrthologuesToFiles(orthologues_alltrees, speciesDict, iSpeciesToUse, sequenceDict, resultsDir, ortholog_file_writers, suspect_genes_file_writers, qContainsSuspectOlogs):
    # Sort the orthologues according to species pairs
    sp_to_index = {str(sp):i for i, sp in enumerate(iSpeciesToUse)}
    nOrtho = util.nOrtho_sp(len(iSpeciesToUse))   
#    left = [[] for sp in species]  
#    right = [[] for sp in species]
    # reorder orthologues on a per-species basis
    if qContainsSuspectOlogs: dSuspect = files.FileHandler.GetPutativeXenelogsDir()
    nSpecies = len(iSpeciesToUse)
    for i in xrange(nSpecies):
        sp0 = str(iSpeciesToUse[i])
        if qContainsSuspectOlogs: 
            writer1_sus = suspect_genes_file_writers[i]
        strsp0 = sp0 + "_"
        isp0 = sp_to_index[sp0]
        for j in xrange(i, nSpecies):
            sp1 = str(iSpeciesToUse[j])
            if sp1 == sp0: continue
            strsp1 = sp1 + "_"
            isp1 = sp_to_index[sp1]
            if qContainsSuspectOlogs:
                writer2_sus = suspect_genes_file_writers[j]
            writer1 = ortholog_file_writers[i][j] 
            writer2 = ortholog_file_writers[j][i] 
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

def RootAndGetOrthologues_from_tree(iog, tree_fn, species_tree_rooted, GeneToSpecies, neighbours, qWrite=False, dupsWriter=None, seqIDs=None, spIDs=None, all_stride_dup_genes=None, qNoRecon=False):
    rooted_tree_ids, qHaveSupport = CheckAndRootTree(tree_fn, species_tree_rooted, GeneToSpecies) # this can be parallelised easily
    if rooted_tree_ids is None: return
    orthologues, recon_tree, suspect_genes = GetOrthologues_from_tree(iog, rooted_tree_ids, species_tree_rooted, GeneToSpecies, neighbours, dupsWriter=dupsWriter, seqIDs=seqIDs, spIDs=spIDs, all_stride_dup_genes=all_stride_dup_genes, qNoRecon=qNoRecon)
    if qWrite:
        directory = os.path.split(tree_fn)[0]
        WriteQfO2(orthologues, directory + "_Orthologues_M3/" + os.path.split(tree_fn)[1], qAppend=False)

def GetOrthologuesStandalone_Parallel(trees_dir, species_tree_rooted_fn, GeneToSpecies, output_dir, qSingleTree):
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    neighbours = GetSpeciesNeighbours(species_tree_rooted)
    args_queue = mp.Queue()
    for treeFn in glob.glob(trees_dir + ("*" if qSingleTree else "/*")): args_queue.put((0, treeFn, species_tree_rooted, GeneToSpecies, neighbours, True))
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
        RootAndGetOrthologues_from_tree(0, treeFn, species_tree_rooted, GeneToSpecies, neighbours, True)        

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
        self.suspect_genes_file_handles = [None for _ in self.iSpeciesToUse]

    def __enter__(self):
        ortholog_file_writers = [[None for _ in self.iSpeciesToUse] for _ in self.iSpeciesToUse]
        suspect_genes_file_writers = [None for _ in self.iSpeciesToUse]
        for i in xrange(self.nSpecies):
            sp0 = str(self.iSpeciesToUse[i])
            self.suspect_genes_file_handles[i] = open(self.dPutativeXenologs + "%s.tsv" % self.speciesDict[sp0], csv_write_mode)
            suspect_genes_file_writers[i] = csv.writer(self.suspect_genes_file_handles[i], delimiter="\t")
            strsp0 = sp0 + "_"
            isp0 = self.sp_to_index[sp0]
            d0 = self.d + "Orthologues_" + self.speciesDict[sp0] + "/"
            for j in xrange(i, self.nSpecies):
                sp1 = str(self.iSpeciesToUse[j])
                if sp1 == sp0: continue
                strsp1 = sp1 + "_"
                isp1 = self.sp_to_index[sp1]
                d1 = self.d + "Orthologues_" + self.speciesDict[sp1] + "/"
                self.ortholog_file_handles[i][j] = open(d0 + '%s__v__%s.tsv' % (self.speciesDict[sp0], self.speciesDict[sp1]), csv_write_mode)
                self.ortholog_file_handles[j][i] = open(d1 + '%s__v__%s.tsv' % (self.speciesDict[sp1], self.speciesDict[sp0]), csv_write_mode)
                ortholog_file_writers[i][j] = csv.writer(self.ortholog_file_handles[i][j], delimiter="\t")
                ortholog_file_writers[j][i] = csv.writer(self.ortholog_file_handles[j][i], delimiter="\t")
        return ortholog_file_writers, suspect_genes_file_writers

    def __exit__(self, type, value, traceback):
        for fh in self.suspect_genes_file_handles:
            fh.close()
        for fh_list in self.ortholog_file_handles:
            for fh in fh_list:
                if fh is not None:
                    fh.close()

def DoOrthologuesForOrthoFinder(ogSet, species_tree_rooted_labelled, GeneToSpecies, all_stride_dup_genes, qNoRecon, hog_writer):   
    """
    """
     # Create directory structure
    speciesDict = ogSet.SpeciesDict()
    SequenceDict = ogSet.SequenceDict()
    # Write directory and file structure
    qInitialisedSuspectGenesDirs = False
    speciesIDs = ogSet.speciesToUse
    nspecies = len(speciesIDs)      
    dResultsOrthologues = files.FileHandler.GetOrthologuesDirectory()
    for index1 in xrange(nspecies):
        d = dResultsOrthologues + "Orthologues_" + speciesDict[str(speciesIDs[index1])] + "/"
        if not os.path.exists(d): os.mkdir(d)     
        for index2 in xrange(nspecies):
            if index2 == index1: continue
            with open(d + '%s__v__%s.tsv' % (speciesDict[str(speciesIDs[index1])], speciesDict[str(speciesIDs[index2])]), csv_write_mode) as outfile:
                writer1 = csv.writer(outfile, delimiter="\t")
                writer1.writerow(("Orthogroup", speciesDict[str(speciesIDs[index1])], speciesDict[str(speciesIDs[index2])]))
    neighbours = GetSpeciesNeighbours(species_tree_rooted_labelled)
    # Infer orthologues and write them to file           
    nOgs = len(ogSet.OGs())
    nOrthologues_SpPair = util.nOrtho_sp(nspecies) 
    species = list(speciesDict.keys())
    reconTreesRenamedDir = files.FileHandler.GetOGsReconTreeDir(True)
    spec_seq_dict = ogSet.Spec_SeqDict()
    sp_to_index = {str(sp):i for i, sp in enumerate(ogSet.speciesToUse)}
    with open(files.FileHandler.GetDuplicationsFN(), csv_write_mode) as outfile, OrthologsFiles(dResultsOrthologues, speciesDict, ogSet.speciesToUse, nspecies, sp_to_index) as (ortholog_file_writers, suspect_genes_file_writers):
        dupWriter = csv.writer(outfile, delimiter="\t")
        dupWriter.writerow(["Orthogroup", "Species Tree Node", "Gene Tree Node", "Support", "Type",	"Genes 1", "Genes 2"])
        for iog in range(nOgs):
            # if iog != 3672: continue
            rooted_tree_ids, qHaveSupport = CheckAndRootTree(files.FileHandler.GetOGsTreeFN(iog), species_tree_rooted_labelled, GeneToSpecies) # this can be parallelised easily
            if rooted_tree_ids is None: continue
            # Write rooted tree with accessions
            util.RenameTreeTaxa(rooted_tree_ids, files.FileHandler.GetOGsTreeFN(iog, True), spec_seq_dict, qSupport=qHaveSupport, qFixNegatives=True, qViaCopy=True)
            orthologues, recon_tree, suspect_genes = GetOrthologues_from_tree(iog, rooted_tree_ids, species_tree_rooted_labelled, GeneToSpecies, neighbours, dupsWriter=dupWriter, seqIDs=spec_seq_dict, spIDs=ogSet.SpeciesDict(), all_stride_dup_genes=all_stride_dup_genes, qNoRecon=qNoRecon)
            GetHOGs_from_tree(iog, recon_tree, hog_writer)
            qContainsSuspectGenes = len(suspect_genes) > 0
            if (not qInitialisedSuspectGenesDirs) and qContainsSuspectGenes:
                qInitialisedSuspectGenesDirs = True
                dSuspectGenes = files.FileHandler.GetSuspectGenesDir()
                dSuspectOrthologues = files.FileHandler.GetPutativeXenelogsDir()
                for index1 in xrange(nspecies):
                    with open(dSuspectOrthologues + '%s.tsv' % speciesDict[str(speciesIDs[index1])], csv_write_mode) as outfile:
                        writer1 = csv.writer(outfile, delimiter="\t")
                        writer1.writerow(("Orthogroup", speciesDict[str(speciesIDs[index1])], "Other"))
            for index0 in xrange(nspecies):
                strsp0 = species[index0]
                strsp0_ = strsp0+"_"
                these_genes = [g for g in suspect_genes if g.startswith(strsp0_)]
                if len(these_genes) > 0:
                    with open(dSuspectGenes + speciesDict[strsp0] + ".txt", 'a') as outfile:
                        outfile.write("\n".join([SequenceDict[g] for g in these_genes]) + "\n")
            allOrthologues = [(iog, orthologues)]
            # don't relabel nodes, they've already been done
            util.RenameTreeTaxa(recon_tree, reconTreesRenamedDir + "OG%07d_tree.txt" % iog, spec_seq_dict, qSupport=False, qFixNegatives=True)
            if iog >= 0 and divmod(iog, 10 if nOgs <= 200 else 100 if nOgs <= 2000 else 1000)[1] == 0:
                util.PrintTime("Done %d of %d" % (iog, nOgs))
            nOrthologues_SpPair += AppendOrthologuesToFiles(allOrthologues, speciesDict, ogSet.speciesToUse, SequenceDict, dResultsOrthologues, ortholog_file_writers, suspect_genes_file_writers, qContainsSuspectGenes)
    return nOrthologues_SpPair


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
    with OrthologsFiles(dResultsOrthologues, speciesDict, ogSet.speciesToUse, nSpecies, sp_to_index) as ortholog_file_writers, suspect_genes_file_writers:
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
                nOrthologues_SpPair += AppendOrthologuesToFiles(allOrthologues, speciesDict, ogSet.speciesToUse, SequenceDict, output_dir, ortholog_file_writers, suspect_genes_file_writers, False)
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

