import os
import sys
import unittest
import argparse
from collections import defaultdict
import ete3

baseDir = os.path.dirname(os.path.realpath(__file__)) + os.sep

sys.path.append(baseDir + "../orthofinder")


class TestGatehring(unittest.TestCase):
    def test_sum_cluster_hits(self):
        seq_ids = defaultdict(str)
        sp_ids = defaultdict(str)
        sp_tree, sp_tree_node_names = species_tree()
        hogw = t2o.HogWriter(sp_tree, sp_tree_node_names, seq_ids, sp_ids)
        empty_set = set()

        # act-assert
        gt = gene_tree()
        gt = hogw.mark_dups_below_v2(gt)
        self.assertTrue(all(n.dups_below == empty_set for n in gt.traverse() if not n.is_leaf()))

        # act-assert
        gt = gene_tree_dup_at_root_N12()
        gt = hogw.mark_dups_below_v2(gt)
        # root has a duplication
        self.assertEqual({'N12'}, gt.dups_below) 
        # none of the other nodes do
        self.assertTrue(all(n.dups_below == empty_set for n in gt.traverse() if not (n.is_leaf() or n.is_root())))

        # act-assert 
        gt = ete3.Tree("((1_0, (12_0, 7_7),(5_0,(12_1, 6_0)));")
        gt = hogw.mark_dups_below_v2(gt)
        self.assertEqual({'N3'}, gt.dups_below)

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