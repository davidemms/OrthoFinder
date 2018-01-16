# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 18:01:40 2018

@author: david
"""

import unittest

import tree
import consensus_tree as ct

import sys
sys.path.append("/home/david/workspace/p4/OrthoFinder/orthofinder/scripts/")
import tree_isomorphism as ti


taxa = "a b c d e".split()
taxa_index = {t:i for i, t in enumerate(taxa)}
a = ct.BitVector(taxa_index, "a")
b = ct.BitVector(taxa_index, "b")
c = ct.BitVector(taxa_index, "c")
d = ct.BitVector(taxa_index, "d")
e = ct.BitVector(taxa_index, "e")

class TestConsensusTree(unittest.TestCase):
    def test_BitVector(self):
        taxa = "a b c d e".split()
        taxa_index = {t:i for i, t in enumerate(taxa)}
        x = ct.BitVector(taxa_index, "a")
        self.assertTrue(x.Is("a"))
        self.assertFalse(x.Is("b"))
        self.assertFalse(x.Is("c"))
        self.assertFalse(x.Is("d"))
        self.assertFalse(x.Is("e"))
        
        y = ct.BitVector(taxa_index, "c")
        self.assertFalse(y.Is("a"))
        self.assertFalse(y.Is("b"))
        self.assertTrue(y.Is("c"))
        self.assertFalse(y.Is("d"))
        self.assertFalse(y.Is("e"))
        
        x.Add(y)
        self.assertTrue(x.Is("a"))
        self.assertFalse(x.Is("b"))
        self.assertTrue(x.Is("c"))
        self.assertFalse(x.Is("d"))
        self.assertFalse(x.Is("e"))
        
        z = ct.BitVector(taxa_index, "c")
        z2 = ct.BitVector(taxa_index, "d")
        z.Add(z2)
        z.Add(y)
        self.assertFalse(z.Is("a"))
        self.assertFalse(z.Is("b"))
        self.assertTrue(z.Is("c"))
        self.assertTrue(z.Is("d"))
        self.assertFalse(z.Is("e"))
        
    def test_UpdateSplits(self):
        all_splits = []
        t = tree.Tree("((a,b),(c,d));")
        taxa = "abcd"
        taxa_index = {t:i for i, t in enumerate(taxa)}
        ct.UpdateSplits(all_splits, t, taxa_index)
        self.assertEqual(len(all_splits), 5)
        s = map(bin, all_splits)
        self.assertTrue("0b1110" in s) # a
        self.assertTrue("0b10" in s)   # b
        self.assertTrue("0b1100" in s) #ab
        self.assertTrue("0b1000" in s) #abc
        self.assertTrue("0b100" in s) #abd
        
        
    def test_GetCompatibleSplits(self):
        x = ct.BitVector(taxa_index, "a")
        x.Add(ct.BitVector(taxa_index, "b"))
        x.Add(ct.BitVector(taxa_index, "c"))
        y = ct.BitVector(taxa_index, "a")
        y.Add(ct.BitVector(taxa_index, "b"))
        z = ct.BitVector(taxa_index, "b")
        z.Add(ct.BitVector(taxa_index, "e"))
        all_splits = [x.Canonical(), x.Canonical(), z.Canonical(), y.Canonical(), y.Canonical(), y.Canonical()]
        # x 0b00111 -> 0b11000
        # y 0b00011 -> 0b11100 
        com_sp = ct.GetCompatibleSplits(all_splits)
        self.assertEqual(len(com_sp),  2)
        self.assertEqual(bin(com_sp[0]), "0b11100")
        self.assertEqual(bin(com_sp[1]), "0b11000")
        
    def test_GetConstructTree(self):
        x = ct.BitVector(taxa_index, "a")
        x.Add(ct.BitVector(taxa_index, "b"))
        x.Add(ct.BitVector(taxa_index, "c"))
        y = ct.BitVector(taxa_index, "a")
        y.Add(ct.BitVector(taxa_index, "b"))
        z = ct.BitVector(taxa_index, "b")
        z.Add(ct.BitVector(taxa_index, "e"))
        all_splits = [x.Canonical(), x.Canonical(), z.Canonical(), y.Canonical(), y.Canonical(), y.Canonical()]
        com_sp = ct.GetCompatibleSplits(all_splits)
        t = ct.ConstructTree(com_sp, taxa_index, taxa)
        self.assertTrue(ti.IsIso_labelled_ete_nonbinary(t, tree.Tree("(((d,e),c),a,b);"), ti.Identity))
        
    def test_GetConstructTree2(self):
        x = ct.BitVector(taxa_index)
        x.Add(a)
        x.Add(b)
        y = ct.BitVector(taxa_index)
        y.Add(c)
        y.Add(e)
        com_sp = [x.Canonical(), y.Canonical()]
        t = ct.ConstructTree(com_sp, taxa_index, taxa)
        self.assertTrue(ti.IsIso_labelled_ete_nonbinary(t, tree.Tree("(a,b,((c,e),d));"), ti.Identity))
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestConsensusTree)
    unittest.TextTestRunner(verbosity=2).run(suite)   