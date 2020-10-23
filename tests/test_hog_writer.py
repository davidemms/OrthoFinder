import os
import sys
import unittest
import argparse
from collections import defaultdict
import ete3

baseDir = os.path.dirname(os.path.realpath(__file__)) + os.sep

sys.path.append(baseDir + "../orthofinder")

import scripts_of.trees2ologs_of as t2o
t2o.files.FileHandler.rd1 = "/tmp/"

def species_tree():
    t = ete3.Tree("(((((1,5)N4,(((15,4)N7,(6,(9,((2,(12,7)N13)N11,(10,13)N12)N10)N9)N8)N6,3)N5)N3,11)N2,16)N1,8)N0;", format=1)
    names = ["N%d" % i for i in range(14)]
    return t, names

def gene_tree():
    t = ete3.Tree("((1_0, 5_1),4_3);")  # (C. elegans, Drosophila),Danio
    for n in t.traverse():
        n.add_feature('dup', False)
    n = t & '4_3'
    n.up.add_feature('sp_node', 'N3')
    return t

def gene_tree_dup_at_root_N12():
    gt = ete3.Tree("((10_0, 13_0),(10_1,13_1));")   # Mouse & Rat
    for n in gt.traverse():
        n.add_feature('dup', False)
    gt.dup = True
    gt.add_feature('sp_node', 'N12')
    for ch in gt.get_children():
        ch.add_feature('sp_node', 'N12')
    return gt

class TestHOGWriter(unittest.TestCase):
    def test_get_hogs_to_write(self):
        seq_ids = defaultdict(str)
        sp_ids = defaultdict(str)
        sp_tree, sp_tree_node_names = species_tree()
        hogw = t2o.HogWriter(sp_tree, sp_tree_node_names, seq_ids, sp_ids)
        gt = gene_tree()
        
        n = gt
        scl = n.get_leaves()    # they're just the leaves as no duplications
        hogs = hogw.get_hogs_to_write(n, scl)

        # all HOGs up and down from root=N3, when writing we'll find most are empty
        self.assertEqual(14, len(hogs))
        self.assertTrue('N0' in hogs)
        self.assertTrue('N7' in hogs)
        self.assertTrue('N12' in hogs)

    def test_get_skipped_nodes_root_duplication_not_N0(self):
        """
        What if the node above the same MRCA, the root, a duplication and 
        not N0? Then we don't write the higher HOGs? Yes, that's right. They're
        taken care of at the duplication node.
        """
        seq_ids = defaultdict(str)
        sp_ids = defaultdict(str)
        sp_tree, sp_tree_node_names = species_tree()
        hogw = t2o.HogWriter(sp_tree, sp_tree_node_names, seq_ids, sp_ids)
        gt = gene_tree_dup_at_root_N12()
        
        n_gene = gt.get_children()[0]
        n_sp_this = sp_tree & 'N12'
        sp_node_list = hogw.get_skipped_nodes(n_sp_this, 'N12', n_gene)

        self.assertEqual(0, len(sp_node_list))

    def test_get_hogs_to_write_root_duplication_not_N0(self):
        """
        How do we handle processing the root of a gene tree that is a duplication?
        """
        seq_ids = defaultdict(str)
        sp_ids = defaultdict(str)
        sp_tree, sp_tree_node_names = species_tree()
        hogw = t2o.HogWriter(sp_tree, sp_tree_node_names, seq_ids, sp_ids)
        gt = gene_tree_dup_at_root_N12()

        n = gt
        scl = n.get_leaves(is_leaf_fn=hogw.scl_fn)
        scl = [gt, ]     # just this node itself since it is a dup
        hogs = hogw.get_hogs_to_write(n, scl)

        # The original way that 'not writing N12' HOG was achieved was poor. It 
        # would be listed by 'get_hogs_to_write' as being required but then 'write_hogs'
        # would fail to find any genes for it (as the only scl_unit was N12, which
        # is not an expected child of N12) and would skip it as empty. It's far clearer
        # for it to be excluded here.
        self.assertEqual(9, len(hogs))  # should be all the ones up from N12 to root, not including N12!
        self.assertTrue('N10' in hogs) # first of these 
        self.assertTrue('N0' in hogs)  # last of these
        self.assertTrue('N12' not in hogs)  # not N12
        self.assertTrue('N11' not in hogs)  # N11 is a different part of the tree


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--test", help="Individual test to run")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print stdout from failing orthofinder run")
    
    args = parser.parse_args()
        
    qVerbose = args.verbose
    
    if args.test != None:
        suite = unittest.TestSuite()
        print("Running single test: %s" % args.test)
        suite.addTest(TestProgramCaller(args.test))
        runner = unittest.TextTestRunner()
        runner.run(suite)
    else:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestHOGWriter)
        unittest.TextTestRunner(verbosity=2).run(suite)  