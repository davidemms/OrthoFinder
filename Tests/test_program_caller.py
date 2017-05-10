#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:44:11 2017

@author: david
"""


import os
import sys
import ete2
import shutil
import filecmp
import argparse
import unittest

baseDir = os.path.dirname(os.path.realpath(__file__)) + os.sep

sys.path.append(baseDir + "../orthofinder/scripts")
import program_caller as pc
        
config_dir = baseDir + "Input/ConfigFiles/"
config_std = config_dir + "configure_cmds.txt"
config_alt = config_dir + "configure_cmds_alt.txt"

output_dir = baseDir + "Output_ProgramCaller/"

class TestProgramCaller(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)
        
    def test_normal_construct(self):
        c = pc.ProgramCaller(config_std)
        methods_msa = "mafft mergealign muscle".split()
        methods_trees = "fasttree raxml iqtree".split()
        msa = c.ListMSAMethods()
        self.assertEqual(len(msa), len(methods_msa))
        for m in methods_msa:
            self.assertTrue(m in msa, m)
        trees = c.ListTreeMethods()
        self.assertEqual(len(trees), len(methods_trees))
        for m in methods_trees:
            self.assertTrue(m in trees, m)
            
        c = pc.ProgramCaller(config_alt)
        methods_msa = "mafft mergealign bad_method".split()
        methods_trees = "fasttree iqtree false_method raxml_no_path".split()
        msa = c.ListMSAMethods()
        self.assertEqual(len(msa), len(methods_msa))
        for m in methods_msa:
            self.assertTrue(m in msa, m)
        trees = c.ListTreeMethods()
        self.assertEqual(len(trees), len(methods_trees))
        for m in methods_trees:
            self.assertTrue(m in trees, m)
        
    def test_GetMSAMethodCommand(self):
        workingDir =  "/home/david/WorkingDir/"
        c = pc.ProgramCaller(config_std)
        # test a method that uses the proposed output filename
        cmd, outfn = c.GetMSAMethodCommand("mafft", workingDir + "OG234.fa", workingDir + "Alignments/OG234.fa", "OG234_id")
        self.assertEqual(cmd, "mafft --localpair --maxiterate 1000 --anysymbol /home/david/WorkingDir/OG234.fa > /home/david/WorkingDir/Alignments/OG234.fa 2> /dev/null") 
        cmd, outfn = c.GetMSAMethodCommand("mafft", workingDir + "OGXYZ.fa", workingDir + "Alignments/OGXXX.fa", "OGYYY_id")
        self.assertEqual(cmd, "mafft --localpair --maxiterate 1000 --anysymbol /home/david/WorkingDir/OGXYZ.fa > /home/david/WorkingDir/Alignments/OGXXX.fa 2> /dev/null") 
            
        # test a method that generates its own output filename based on unique ID - don't actually test what the real output filename is here - need to run it for that
        cmd, outfn = c.GetMSAMethodCommand("mergealign", workingDir + "OG234.fa", workingDir + "Alignments/OG234.fa", "OG234_id")
        self.assertEqual(cmd, "python /home/david/software/MergeAlign/run_mergealign.py /home/david/WorkingDir/OG234.fa") 
        
#    def test_TestMSAMethod(self):
#        c = pc.ProgramCaller(config_alt)
#        # test a method that won't work
#        self.assertFalse(c.TestMSAMethod(output_dir, "bad_method"))
#        # test a method that ses the proposed output filename
#        self.assertTrue(c.TestMSAMethod(output_dir, "mafft"))
#        # test a method that generates its own output filename based on unique ID 
#        self.assertTrue(c.TestMSAMethod(output_dir, "mergealign"))
#        self.assertRaises(Exception, c.TestMSAMethod, output_dir, "made_up")
            
    def test_TestTreeMethod(self):
        c = pc.ProgramCaller(config_alt)
        # test a method that won't work
        self.assertFalse(c.TestTreeMethod(output_dir, "false_method"))
        # test a method that ses the proposed output filename
        self.assertTrue(c.TestTreeMethod(output_dir, "fasttree"))
        # test a method that generates its own output filename based on unique ID 
        self.assertTrue(c.TestTreeMethod(output_dir, "iqtree"))
        self.assertRaises(Exception, c.TestTreeMethod, output_dir, "made_up")
        # test a program not in path
        self.assertFalse(c.TestTreeMethod(output_dir, "raxml_no_path"))
    
        
    def test_CallMSAMethod(self):
        c = pc.ProgramCaller(config_alt)
        # Write an input file
        infn = output_dir + "Test.fa"
        with open(infn, 'wb') as outfile:
            outfile.write(">a1\nST\n>b2\nKL\n>c3\nSL\n>d4\nKT\n")
        
        # test a method that uses the proposed output filename
        exp_outfn = output_dir + "OG234.fa"
        if os.path.exists(exp_outfn):
            os.remove(exp_outfn)
        outfn = c.CallMSAMethod("mafft", infn, exp_outfn, "OGabc_id")
        self.assertTrue(outfn != None)
        self.assertEqual(outfn, exp_outfn)
        self.assertTrue(filecmp.cmp(infn, outfn))       # Alignment should be same as input
        
        # test a method that generates its own output filename based on input filename 
        exp_outfn = output_dir + "MergeAlign_alignment_Test.fa"
        if os.path.exists(exp_outfn):
            os.remove(outfn)
        outfn = c.CallMSAMethod("mergealign", infn, output_dir + "Alignments/OG234.fa", "OG234_id")
        self.assertEqual(outfn, exp_outfn)
        self.assertTrue(filecmp.cmp(baseDir + "ExpectedOutput/ProgramCaller/MergeAlign_alignment_Test.fa", outfn)) 
     
    def test_CallTreeMethod(self): 
        c = pc.ProgramCaller(config_alt)
        # Write an input file
        infn = output_dir + "Test.fa"
        with open(infn, 'wb') as outfile:
            outfile.write(">a1\nST\n>b2\nKL\n>c3\nSL\n>d4\nKT")
        
        # test a method that uses the proposed output filename
        exp_outfn = output_dir + "OG234.tre"
        if os.path.exists(exp_outfn):
            os.remove(exp_outfn)
        outfn = c.CallTreeMethod("fasttree", infn, exp_outfn, "OGabc_id")
        self.assertTrue(outfn != None)
        self.assertEqual(outfn, exp_outfn)
        expectedTree = "(a1:1.10312,c3:0.00055,(b2:0.00055,d4:1.46340)0.761:1.41871);"
        with open(outfn, 'rb') as infile: tree = infile.read().rstrip()
        self.assertEqual(expectedTree, tree)
            
        # test a method that generates its own output filename based on unique ID 
        exp_outfn = output_dir + "OG234_id.treefile"
        if os.path.exists(exp_outfn): os.remove(exp_outfn)
        outfn = c.CallTreeMethod("iqtree", infn, output_dir + "OG234.tre", "OG234_id")
        self.assertEqual(outfn, exp_outfn)
#        with open(outfn, 'rb') as infile: tree = infile.read().rstrip()
        expectedTree = ete2.Tree("(a1:0.7475350209,(b2:0.0000026604,c3:1.2318876921):1.279964,d4:0.0000022111);")
        actualTree = ete2.Tree(outfn)
        for n in expectedTree:
            x = (actualTree & n.name).dist
            self.assertLess(abs(x-n.dist)/n.dist, 0.3, (n.dist, (actualTree & n.name).dist))
        
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
        suite = unittest.TestLoader().loadTestsFromTestCase(TestProgramCaller)
        unittest.TextTestRunner(verbosity=2).run(suite)    