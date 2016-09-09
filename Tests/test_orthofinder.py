#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 11:25:10 2015

@author: david
"""

import os 
import time
import types
import datetime
import unittest
import subprocess
import shutil
import filecmp
import glob
import random
import csv
import argparse

__skipLongTests__ = False
qVerbose = False

baseDir = os.path.dirname(os.path.realpath(__file__)) + os.sep
qBinary = False
orthofinder = baseDir + "../orthofinder/orthofinder.py"
orthofinder_bin = baseDir + "../orthofinder/bin/orthofinder"
trees_for_orthogroups = baseDir + "../orthofinder/trees_from_MSA.py"
trees_for_orthogroups_bin = baseDir + "../orthofinder/bin/trees_from_MSA"
exampleFastaDir = baseDir + "Input/SmallExampleDataset/"
exampleBlastDir = baseDir + "Input/SmallExampleDataset_ExampleBlastDir/"

goldResultsDir_smallExample = baseDir + "ExpectedOutput/SmallExampleDataset/"
goldPrepareBlastDir = baseDir + "ExpectedOutput/SmallExampleDataset_PreparedForBlast/"

version = "1.0.4"
requiredBlastVersion = "2.2.28+"

citation = """When publishing work that uses OrthoFinder please cite:
    D.M. Emms & S. Kelly (2015), OrthoFinder: solving fundamental biases in whole genome comparisons
    dramatically improves orthogroup inference accuracy, Genome Biology 16:157.
""" 

expectedHelp="""OrthoFinder version %s Copyright (C) 2014 David Emms

    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software, and you are welcome to redistribute it under certain conditions.
    For details please see the License.md that came with this software.

Simple Usage
------------
python orthofinder.py -f fasta_directory [-t number_of_blast_threads] [-a number_of_orthofinder_threads]
    Infers orthogroups for the proteomes contained in fasta_directory running
    number_of_blast_threads in parallel for the BLAST searches and subsequently running
    number_of_orthofinder_threads in parallel for the OrthoFinder algorithm.

Advanced Usage
--------------
python orthofinder.py -f fasta_directory -p
    1. Prepares files for BLAST and prints the BLAST commands. Does not perform BLAST searches
    or infer orthogroups. Useful if you want to prepare the files in the form required by
    OrthoFinder but want to perform the BLAST searches using a job scheduler/on a cluster and
    then infer orthogroups using option 2.

python orthofinder.py -b precalculated_blast_results_directory [-a number_of_orthofinder_threads]
    2. Infers orthogroups using pre-calculated BLAST results. These can be after BLAST
    searches have been completed following the use of option 1 or using the WorkingDirectory
    from a previous OrthoFinder run. Species can be commented out with a '#' in the SpeciesIDs.txt
    file to exclude them from the analysis. See README file for details.

python orthofinder.py -b precalculated_blast_results_directory -f fasta_directory [-t number_of_blast_threads] [-a number_of_orthofinder_threads]
    3. Add species from fasta_directory to a previous OrthoFinder run where precalculated_blast_results_directory
    is the directory containing the BLAST results files etc. from the previous run.

Arguments
---------
-f fasta_directory, --fasta fasta_directory
    Predict orthogroups for the proteins in the fasta files in the fasta_directory

-b precalculated_blast_results_directory, --blast precalculated_blast_results_directory
    Predict orthogroups using the pre-calcualted BLAST results in precalculated_blast_results_directory.

-t number_of_blast_threads, --threads number_of_blast_threads
    The number of BLAST processes to be run simultaneously. This should be increased by the user to at least 
    the number of cores on the computer so as to minimise the time taken to perform the BLAST all-versus-all 
    queries. [Default is 16]

-a number_of_orthofinder_threads, --algthreads number_of_orthofinder_threads
    The number of threads to use for the OrthoFinder algorithm and MCL after BLAST searches have been completed. 
    Running the OrthoFinder algorithm with a number of threads simultaneously increases the RAM 
    requirements proportionally so be aware of the amount of RAM you have available (and see README file). 
    Additionally, as the algorithm implementation is very fast, file reading is likely to be the 
    limiting factor above about 5-10 threads and additional threads may have little effect other than 
    increase RAM requirements. [Default is 1]

-g, --groups
    Only infer orthogroups, do not infer gene trees of orthologues.

-I inflation_parameter, --inflation inflation_parameter
    Specify a non-default inflation parameter for MCL. [Default is 1.5]

-x speciesInfoFilename, --orthoxml speciesInfoFilename
    Output the orthogroups in the orthoxml format using the information in speciesInfoFilename.

-p , --prepare
    Only prepare the files in the format required by OrthoFinder and print out the BLAST searches that
    need to be performed but don't run BLAST or infer orthogroups

-h, --help
    Print this help text

When publishing work that uses OrthoFinder please cite:
    D.M. Emms & S. Kelly (2015), OrthoFinder: solving fundamental biases in whole genome comparisons
    dramatically improves orthogroup inference accuracy, Genome Biology 16:157.
""" % version

class CleanUp(object):
    """Cleans up after arbitrary code that could create any/all of the 'newFiles'
    and modify the 'modifiedFiles'
    
    Implementation:
        Makes copies of files in modifiedFiles
        when the context handler exits:
        - deletes any files from newFiles
        - uses the copies of modifiedFiles to revert them to their previous state
    """
    def __init__(self, newFiles, modifiedFiles, newDirs = [], qSaveFiles=False):
        """qSaveFiles is useful for debuging purposes
        """
        self.newFiles = newFiles
        self.modifiedFiles = modifiedFiles
        self.copies = []
        assert(types.ListType == type(newDirs)) # if it were a string the code could attempt to delete every file on computer
        self.newDirs = newDirs
        self.qSaveFiles = qSaveFiles
    def __enter__(self):
        for fn in self.modifiedFiles:
            copy = fn + "_bak%d" % random.randint(0, 999999)
            shutil.copy(fn, copy)
            self.copies.append(copy)
    def __exit__(self, type, value, traceback):
        if self.qSaveFiles and len(self.newFiles) != 0:
            saveDir = os.path.split(self.newFiles[0])[0] + "/SavedFiles/"
            os.mkdir(saveDir)
            for fn in self.modifiedFiles + self.newFiles:
                if os.path.exists(fn):
                    shutil.copy(fn, saveDir + os.path.split(fn)[1])
        for fn, copy in zip(self.modifiedFiles, self.copies):
            shutil.move(copy, fn)
        for fn in self.newFiles:
            if os.path.exists(fn): os.remove(fn)
        for d in self.newDirs:
            if self.qSaveFiles:
                while d[-1] == "/": d = d[:-1]
                shutil.move(d, d + "_bak/")
            else:
                shutil.rmtree(d)

class TestCommandLine(unittest.TestCase):
    """General expected (tested) behaviour 
    stdout:
        - Citation should be printed
        - Can check any other details later
        
    Results:
        - Prints name of results files
            - These are were they should be
        - Claimed results files should've only just been created
        - Results files should all contain the correct orthogroups
    """
    @classmethod
    def setUpClass(cls):
        capture = subprocess.Popen("blastp -version", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        stdout = "".join([x for x in capture.stdout])
        if requiredBlastVersion not in stdout:
            raise RuntimeError("Tests require BLAST version %s" % requiredBlastVersion)       
            
    def test_fromfasta_threads(self):
        currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        expectedCSVFile = currentResultsDir + "Orthogroups.csv"
        with CleanUp([], [], [currentResultsDir, ]):
            self.stdout, self.stderr = self.RunOrthoFinder("-f %s -t 4 -a 3 -g" % exampleFastaDir)
            self.CheckStandardRun(self.stdout, self.stderr, goldResultsDir_smallExample, expectedCSVFile)  
        self.test_passed = True         
           
    @unittest.skipIf(__skipLongTests__, "Only performing quick tests")     
    def test_fromfasta_full(self):
        currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        expectedCSVFile = currentResultsDir + "Orthogroups.csv"     
        with CleanUp([], [], [currentResultsDir, ]):
            self.stdout, self.stderr = self.RunOrthoFinder("--fasta %s -g" % exampleFastaDir)
            self.CheckStandardRun(self.stdout, self.stderr, goldResultsDir_smallExample, expectedCSVFile)  
        self.test_passed = True         
        
    def test_justPrepare(self):    
        self.currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        self.stdout, self.stderr = self.RunOrthoFinder("-f %s -p" % exampleFastaDir)
        expectedFiles = ['BlastDBSpecies0.phr', 'BlastDBSpecies0.pin', 'BlastDBSpecies0.psq', 
                         'BlastDBSpecies1.phr', 'BlastDBSpecies1.pin', 'BlastDBSpecies1.psq', 
                         'BlastDBSpecies2.phr', 'BlastDBSpecies2.pin', 'BlastDBSpecies2.psq', 
                         'Species0.fa', 'Species1.fa', 'Species2.fa',
                         'SequenceIDs.txt', 'SpeciesIDs.txt']
        # all files are present and have only just been created
        for fn in expectedFiles:
            fullFN = self.currentResultsDir + "WorkingDirectory/" + fn
            self.assertTrue(os.path.exists(fullFN))
            self.assertLess(time.time()-os.stat(fullFN).st_ctime, 20)
            if not fullFN.endswith(".pin"):
                self.assertTrue(filecmp.cmp(goldPrepareBlastDir + fn, fullFN), msg=fullFN)
            
        # no other files in directory or intermediate directory
        self.assertTrue(len(list(glob.glob(self.currentResultsDir + "WorkingDirectory/*"))), len(expectedFiles)) 
        self.assertTrue(len(list(glob.glob(self.currentResultsDir + "*"))), 1) # Only 'WorkingDirectory/"
        self.test_passed = True         
        
    def test_justPrepareFull(self):    
        self.currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        self.stdout, self.stderr = self.RunOrthoFinder("-f %s --prepare" % exampleFastaDir)
        expectedFiles = ['BlastDBSpecies0.phr', 'BlastDBSpecies0.pin', 'BlastDBSpecies0.psq', 
                         'BlastDBSpecies1.phr', 'BlastDBSpecies1.pin', 'BlastDBSpecies1.psq', 
                         'BlastDBSpecies2.phr', 'BlastDBSpecies2.pin', 'BlastDBSpecies2.psq', 
                         'Species0.fa', 'Species1.fa', 'Species2.fa',
                         'SequenceIDs.txt', 'SpeciesIDs.txt']
        # all files are present and have only just been created
        for fn in expectedFiles:
            fullFN = self.currentResultsDir + "WorkingDirectory/" + fn
            self.assertTrue(os.path.exists(fullFN))
            self.assertLess(time.time()-os.stat(fullFN).st_ctime, 20)
            if not fullFN.endswith(".pin"):
                self.assertTrue(filecmp.cmp(goldPrepareBlastDir + fn, fullFN), msg=fullFN)
            
        # no other files in directory or intermediate directory
        self.assertTrue(len(list(glob.glob(self.currentResultsDir + "WorkingDirectory/*"))), len(expectedFiles)) 
        self.assertTrue(len(list(glob.glob(self.currentResultsDir + "*"))), 1) # Only 'WorkingDirectory/"  
        self.test_passed = True         
        
    def test_fromblast(self):
        expectedCSVFile = exampleBlastDir + "Orthogroups.csv"
        newFiles = ("Orthogroups.csv Orthogroups_UnassignedGenes.csv Orthogroups.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt clusters_OrthoFinder_v%s_I1.5.txt OrthoFinder_v%s_graph.txt Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv" % (version, version, version)).split()
        newFiles = [exampleBlastDir + fn for fn in newFiles]
        with CleanUp(newFiles, []):
            self.stdout, self.stderr = self.RunOrthoFinder("-b %s -g" % exampleBlastDir)
            self.CheckStandardRun(self.stdout, self.stderr, goldResultsDir_smallExample, expectedCSVFile)  
        self.test_passed = True         
        
    def test_fromblast_full(self):
        expectedCSVFile = exampleBlastDir + "Orthogroups.csv"
        newFiles = ("Orthogroups.csv Orthogroups_UnassignedGenes.csv Orthogroups.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt clusters_OrthoFinder_v%s_I1.5.txt OrthoFinder_v%s_graph.txt Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv" % (version, version, version)).split()
        newFiles = [exampleBlastDir + fn for fn in newFiles]
        with CleanUp(newFiles, []):        
            self.stdout, self.stderr = self.RunOrthoFinder("--blast %s -g" % exampleBlastDir)
            self.CheckStandardRun(self.stdout, self.stderr, goldResultsDir_smallExample, expectedCSVFile)  
        self.test_passed = True         
        
    def test_fromblast_algthreads(self):
        expectedCSVFile = exampleBlastDir + "Orthogroups.csv"
        newFiles = ("Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv Orthogroups.csv Orthogroups_UnassignedGenes.csv Orthogroups.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt clusters_OrthoFinder_v%s_I1.5.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()
        newFiles = [exampleBlastDir + fn for fn in newFiles]
        with CleanUp(newFiles, []):
            self.stdout, self.stderr = self.RunOrthoFinder("-b %s -a 3 -g" % exampleBlastDir)
            self.CheckStandardRun(self.stdout, self.stderr, goldResultsDir_smallExample, expectedCSVFile)  
        self.test_passed = True         
        
    def test_blast_results_error(self):
        d = baseDir + "Input/SmallExampleDataset_ExampleBadBlast/"
        newFiles = [d + "%s%d_%d.pic" % (s, i,j) for i in xrange(1, 3) for j in xrange(3) for s in ["B", "BH"]]
        with CleanUp(newFiles, []):
            self.stdout, self.stderr = self.RunOrthoFinder("-a 2 -g -b " + d)
            self.assertTrue("Traceback" not in self.stderr)
            self.assertTrue("Offending line was:" in self.stderr)
            self.assertTrue("0_0	0_0	100.00	466" in self.stderr)
            self.assertTrue("Connected putatitive homologs" not in self.stdout)
            self.assertTrue("ERROR: An error occurred, please review previous error messages for more information." in self.stdout)  
        self.test_passed = True         
        
    def test_inflation(self):
        expectedCSVFile = exampleBlastDir + "Orthogroups.csv"
        newFiles = ("Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv Orthogroups.csv Orthogroups_UnassignedGenes.csv Orthogroups.txt clusters_OrthoFinder_v%s_I1.8.txt_id_pairs.txt clusters_OrthoFinder_v%s_I1.8.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()
        newFiles = [exampleBlastDir + fn for fn in newFiles]
        with CleanUp(newFiles, []):
            self.stdout, self.stderr = self.RunOrthoFinder("-I 1.8 -b %s -g" % exampleBlastDir)
            self.CheckStandardRun(self.stdout, self.stderr, baseDir + "ExpectedOutput/SmallExampleDataset_I1.8/", expectedCSVFile)  
        self.test_passed = True         
    
#    @unittest.skipIf(__skipLongTests__, "Only performing quick tests")       
#    def test_fromblastOrthobench(self):
#        goldResultsDir_orthobench = baseDir + "ExpectedOutput/Orthobench_blast/"
#        expectedCSVFileLocation = baseDir + "Input/Orthobench_blast/Orthogroups.csv"
#        self.currentResultsDir = None
#        expectedNewFiles = [baseDir + "Input/Orthobench_blast/" + x for x in "OrthoFinder_v0.4.0_graph.txt clusters_OrthoFinder_v0.4.0_I1.5.txt clusters_OrthoFinder_v0.4.0_I1.5.txt_id_pairs.txt".split()]
#        with CleanUp(expectedNewFiles, []):
#            self.stdout, self.stderr = self.RunOrthoFinder("-b %sInput/Orthobench_blast" % baseDir)
#            self.CheckStandardRun(stdout, stderr, goldResultsDir_orthobench, expectedCSVFileLocation, qDelete = True)  
#        self.test_passed = True         

    def test_help(self):
        self.stdout, self.stderr = self.RunOrthoFinder("")
        self.assertTrue(expectedHelp in self.stdout)
        self.assertEqual(len(self.stderr), 0)
        
        self.stdout, self.stderr = self.RunOrthoFinder("-h")
        self.assertTrue(expectedHelp in self.stdout)
        self.assertEqual(len(self.stderr), 0)
         
        self.stdout, self.stderr = self.RunOrthoFinder("--help")
        self.assertTrue(expectedHelp in self.stdout)
        self.assertEqual(len(self.stderr), 0)        
        self.test_passed = True     
         
#    def test_numberOfThreads(self):
#        pass  
#         
#    def test_numberOfThreads_noBlast(self):
#        # deal with it gracefully
#        pass
#    
#    def test_xmlOutput(self):
#        pass
         
    def test_addOneSpecies(self):
        expectedExtraFiles = [exampleBlastDir + fn for fn in ("Blast0_3.txt Blast3_0.txt Blast1_3.txt Blast3_1.txt Blast2_3.txt Blast3_2.txt Blast3_3.txt Species3.fa \
        Orthogroups.csv Orthogroups.txt Orthogroups_UnassignedGenes.csv \
        Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv \
        clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()]
        expectedChangedFiles = [exampleBlastDir + fn for fn in "SpeciesIDs.txt SequenceIDs.txt".split()]
        # cleanup afterwards including failed test
        goldDir = baseDir + "ExpectedOutput/AddOneSpecies/"
        expectedExtraDir = exampleBlastDir + "Orthologues_%s/" % datetime.date.today().strftime("%b%d") 
        with CleanUp(expectedExtraFiles, expectedChangedFiles, [expectedExtraDir]):        
            self.stdout, self.stderr = self.RunOrthoFinder("-b %s -f %s" % (exampleBlastDir, baseDir + "Input/ExtraFasta"))
            # check extra blast files
            # check extra fasta file: simple couple of checks to ensure all ok
            for fn in expectedExtraFiles:
                os.path.split(fn)[1]
                self.assertTrue(os.path.exists(fn), msg=fn)
                # mcl output files contain a variable header, these files are an implementation detail that I don't want to test (I want the final orthogroups to be correct)
                if "clusters" in os.path.split(fn)[1]: continue
                self.CompareFile(goldDir + os.path.split(fn)[1], fn)      
        self.test_passed = True
    
    def test_addTwoSpecies(self):
        expectedExtraFiles = [exampleBlastDir + fn for fn in ("Blast0_3.txt Blast3_0.txt Blast1_3.txt Blast3_1.txt Blast2_3.txt Blast3_2.txt Blast3_3.txt Species3.fa \
        Blast0_4.txt Blast4_0.txt Blast1_4.txt Blast4_1.txt Blast2_4.txt Blast4_2.txt Blast3_4.txt Blast4_3.txt Blast4_4.txt Species4.fa \
        Orthogroups.csv Orthogroups.txt Orthogroups_UnassignedGenes.csv \
        Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv \
        clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()]
        expectedChangedFiles = [exampleBlastDir + fn for fn in "SpeciesIDs.txt SequenceIDs.txt".split()]
        goldDir = baseDir + "ExpectedOutput/AddTwoSpecies/"
        with CleanUp(expectedExtraFiles, expectedChangedFiles):        
            self.stdout, self.stderr = self.RunOrthoFinder("-b %s -g -f %s" % (exampleBlastDir, baseDir + "Input/ExtraFasta2"))
            for fn in expectedExtraFiles:
                os.path.split(fn)[1]
                self.assertTrue(os.path.exists(fn), msg=fn)
                if "clusters" in os.path.split(fn)[1]: continue
                self.CompareFile(goldDir + os.path.split(fn)[1], fn)  
        self.test_passed = True
                    
    def test_addTwoSpecies_blastsRequired(self):  
        expectedExtraFiles = [exampleBlastDir + fn for fn in "Species3.fa Species4.fa".split()]
        expectedExtraFiles = expectedExtraFiles + [exampleBlastDir +  "BlastDBSpecies%d.%s" % (i, ext) for i in xrange(5) for ext in "phr pin psq".split()]
        expectedChangedFiles = [exampleBlastDir + fn for fn in "SpeciesIDs.txt SequenceIDs.txt".split()]
        # cleanup afterwards including failed test
        with CleanUp(expectedExtraFiles, expectedChangedFiles):        
            self.stdout, self.stderr = self.RunOrthoFinder("-b %s -f %s -p" % (exampleBlastDir, baseDir + "Input/ExtraFasta2"))
            original = [0, 1, 2]
            new = [3, 4]
            for i in new:
                for j in original:
                    assert("Blast%d_%d.txt" % (i,j) in self.stdout)
                    assert("Blast%d_%d.txt" % (j,i) in self.stdout)
                for j in new:
                    assert("Blast%d_%d.txt" % (i,j) in self.stdout)
                    assert("Blast%d_%d.txt" % (j,i) in self.stdout)
                
            assert(self.stdout.count("blastp") == 2*len(new)*len(original) + len(new)**2)     
        self.test_passed = True         
                
    def test_removeFirstSpecies(self):
        self.RemoveSpeciesTest(baseDir + "Input/ExampleDataset_removeFirst/",
                               baseDir + "ExpectedOutput/RemoveFirstSpecies/") 
                
    def test_removeMiddleSpecies(self):
        self.RemoveSpeciesTest(baseDir + "Input/ExampleDataset_removeMiddle/",
                               baseDir + "ExpectedOutput/RemoveMiddleSpecies/") 
    
    def test_removeLastSpecies(self):
        self.RemoveSpeciesTest(baseDir + "Input/ExampleDataset_removeLast/",
                               baseDir + "ExpectedOutput/RemoveLastSpecies/") 
    
    def RemoveSpeciesTest(self, inputDir, goldDir):
        """Working directory and results directory with correct files in"""
        requiredResults = [inputDir + fn for fn in "Orthogroups.csv Orthogroups_UnassignedGenes.csv Orthogroups.txt".split()]
        expectedExtraFiles = [inputDir + fn for fn in ("Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()]
        with CleanUp(expectedExtraFiles + requiredResults, []):
            self.stdout, self.stderr = self.RunOrthoFinder("-b %s" % inputDir)
            for fn in requiredResults:
                self.assertTrue(os.path.exists(fn), msg=fn)
                self.CompareFile(goldDir + os.path.split(fn)[1], fn)    
        self.test_passed = True         
    
#    def test_removeMultipleSpecies(self):
#        pass
    
    def test_removeOneAddOne(self):
        inputDir = baseDir + "Input/ExampleDataset_addOneRemoveOne/Results_Jan28/WorkingDirectory/"
        expectedExtraFiles = [inputDir + fn for fn in ("Blast0_3.txt Blast3_0.txt Blast1_3.txt Blast3_1.txt Blast2_3.txt Blast3_2.txt Blast3_3.txt Species3.fa \
        Orthogroups.csv Orthogroups.txt Orthogroups_UnassignedGenes.csv \
        Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv \
        clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()] 
        expectedChangedFiles = [inputDir + fn for fn in "SpeciesIDs.txt SequenceIDs.txt".split()]
        goldDir = baseDir + "ExpectedOutput/AddOneRemoveOne/"
        with CleanUp(expectedExtraFiles, expectedChangedFiles):        
            self.stdout, self.stderr = self.RunOrthoFinder("-b %s -f %s" % (inputDir, baseDir + "Input/ExampleDataset_addOneRemoveOne/ExtraFasta/"))
#            print(stdout)
#            print(stderr)
            for fn in expectedExtraFiles:
                os.path.split(fn)[1]
                self.assertTrue(os.path.exists(fn), msg=fn)
                if "Orthogroups" in os.path.split(fn)[1]:
                    self.CompareFile(goldDir + os.path.split(fn)[1], fn)  
            self.CompareFile(goldDir + "SpeciesIDs.txt", inputDir + "SpeciesIDs.txt") 
            self.CompareFile(goldDir + "SequenceIDs.txt", inputDir + "SequenceIDs.txt")  
        self.test_passed = True         
    
    
    def test_removeOneAddOne_fullAnalysis(self):
        inputDir = baseDir + "Input/AddOneRemoveOne_FullAnalysis/Results_Sep09/WorkingDirectory/"
        extraBlast = [inputDir + "Blast%d_4.txt" % i for i in xrange(5)] + [inputDir + "Blast4_%d.txt" % i for i in xrange(4)]
        expectedExtraFiles = extraBlast + [inputDir + fn for fn in ("Orthogroups.csv Orthogroups.txt Orthogroups_UnassignedGenes.csv \
        Statistics_PerSpecies.csv Statistics_Overall.csv Orthogroups_SpeciesOverlaps.csv \
        clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()] 
        expectedChangedFiles = [inputDir + fn for fn in "SpeciesIDs.txt SequenceIDs.txt".split()]
        goldDir = baseDir + "ExpectedOutput/AddOneRemoveOne_FullAnalysis/"
        expExtraDir = [inputDir + "Orthologues_%s/" % datetime.date.today().strftime("%b%d") ]
        with CleanUp(expectedExtraFiles, expectedChangedFiles, expExtraDir):        
            self.stdout, self.stderr = self.RunOrthoFinder("-b %s -f %s" % (inputDir, baseDir + "Input/AddOneRemoveOne_FullAnalysis/ExtraFasta/"))
            for fn in expectedExtraFiles:
                os.path.split(fn)[1]
                self.assertTrue(os.path.exists(fn), msg=fn) 
            self.assertTrue(citation in self.stdout)
            egOrthologuesFN = expExtraDir[0] + "Orthologues/Orthologues_Mycoplasma_alkalescens/Mycoplasma_alkalescens__v__Mycoplasma_hyopneumoniae.csv"
            self.assertTrue(os.path.exists(egOrthologuesFN))
            self.CompareFile_95Percent(goldDir + "Mycoplasma_alkalescens__v__Mycoplasma_hyopneumoniae.csv", egOrthologuesFN) 
            nTrees = 366
            for i in xrange(nTrees):
                self.assertTrue(os.path.exists(expExtraDir[0] + "Gene_Trees/OG%07d_tree.txt" % i))
                self.assertGreater(os.stat(expExtraDir[0] + "Gene_Trees/OG%07d_tree.txt" % i).st_size, 200)
            self.assertTrue(os.path.exists(expExtraDir[0] + "SpeciesTree_rooted.txt"))
            self.assertGreater(os.stat(expExtraDir[0] + "SpeciesTree_rooted.txt").st_size, 150)
        self.test_passed = True      
        
#    def test_addMultipleSpecies_prepare(selg):
#        pass
#    
#    def test_incorrectlyPreparedBlast(self):
#        pass
#    
#    def test_errorsWithFastaFiles(self):
#        pass
#
#    def test_differentCommandOrder(self):
#        pass
#
#    def test_withWithoutSeparator(self):
#        # with and without trailing directory separator
#        pass
#    
#    def test_timeForBenchmark(self):
#        pass
    
    """ Test that all Sequence files are identical to expected and that the expected Alignment & Tree files all exist.
    Don't require that the Alignment & Tree files are identical, they will change with different versions of programs used"""
    def test_trees(self):
        expectedDirs = "Alignments Trees Sequences".split()
        newDirs = [baseDir + "Input/SmallExampleDataset_forTrees/Results_Jan28/" + d +"/" for d in expectedDirs]
        goldDirs = [baseDir + "ExpectedOutput/SmallExampleDataset_trees/" + d + "/" for d in expectedDirs]
        nTrees = 427
        with CleanUp([], [], newDirs):   
            self.stdout, self.stderr = self.RunTrees(baseDir + "Input/SmallExampleDataset_forTrees/Results_Jan28/")
            for i in xrange(nTrees):
                self.assertTrue(os.path.exists(newDirs[0] + "/OG%07d.fa" % i), msg=str(i))
                self.assertTrue(os.path.exists(newDirs[1] + "/OG%07d_tree.txt" % i))
            goldD = goldDirs[2]
            newD = newDirs[2]
            for goldFN in glob.glob(goldD + "*"):
                newFN = newD + os.path.split(goldFN)[1]
                self.assertTrue(os.path.exists(newFN), msg=goldFN)
                self.CompareFile(goldFN, newFN)     
        self.test_passed = True         
    
    def test_treesAfterOption_b(self):
        expectedDirs = "Alignments Trees Sequences".split()
        newDirs = [baseDir + "Input/SmallExampleDataset_forTreesFromBlast/" + d +"/" for d in expectedDirs]
        goldDirs = [baseDir + "ExpectedOutput/SmallExampleDataset_trees/" + d +"/" for d in expectedDirs]
        with CleanUp([], [], newDirs):        
            self.stdout, self.stderr = self.RunTrees(baseDir + "Input/SmallExampleDataset_forTreesFromBlast/")
            for goldD, newD in zip(goldDirs, newDirs):
                for goldFN in glob.glob(goldD + "*"):
                    newFN = newD + os.path.split(goldFN)[1]
                    self.assertTrue(os.path.exists(newFN), msg=goldFN)
                    if newD[:-1].rsplit("/", 1)[-1] == "Sequences":
                        self.CompareFile(goldFN, newFN)       
        self.test_passed = True         
    
    def test_treesMultipleResults(self):
        expectedDirs = "Alignments Trees Sequences".split()
        newDirs = [baseDir + "Input/MultipleResults/Results_Jan28/WorkingDirectory/" + d +"/" for d in expectedDirs]
        goldDirs = [baseDir + "ExpectedOutput/MultipleResultsInDirectory/" + d +"/" for d in expectedDirs]
        with CleanUp([], [], newDirs):        
            self.stdout, self.stderr = self.RunTrees(baseDir + "Input/MultipleResults/Results_Jan28/WorkingDirectory/clusters_OrthoFinder_v0.4.0_I1.5.txt_id_pairs.txt")
            for goldD, newD in zip(goldDirs, newDirs):
                for goldFN in glob.glob(goldD + "*"):
                    newFN = newD + os.path.split(goldFN)[1]
                    self.assertTrue(os.path.exists(newFN), msg=goldFN)
                    if newD[:-1].rsplit("/", 1)[-1] == "Sequences":
                        self.CompareFile(goldFN, newFN)      
        self.test_passed = True         
        
    def test_treesSubset_unspecified(self):
        self.stdout, self.stderr = self.RunTrees(baseDir + "/Input/Trees_OneSpeciesRemoved")
        self.assertTrue("ERROR: Results from multiple OrthoFinder runs found" in self.stdout)
        self.assertTrue("Please run with only one set of results in directories or specifiy the specific clusters_OrthoFinder_*.txt_id_pairs.txt file on the command line" in self.stdout)  
        self.test_passed = True         
    
    def test_treesResultsChoice_subset(self):   
        expectedDirs = "Alignments Trees Sequences".split()
        newDirs = [baseDir + "/Input/Trees_OneSpeciesRemoved/" + d +"/" for d in expectedDirs]
        goldDirs = [baseDir + "ExpectedOutput/SmallExampleDataset_trees/" + d + "/" for d in expectedDirs]
        nTrees = 427
        with CleanUp([], [], newDirs):   
            self.stdout, self.stderr = self.RunTrees(baseDir + "/Input/Trees_OneSpeciesRemoved/clusters_OrthoFinder_v0.6.1_I1.5_1.txt_id_pairs.txt")
            for i in xrange(nTrees):
                self.assertTrue(os.path.exists(newDirs[0] + "/OG%07d.fa" % i), msg=str(i))
                self.assertTrue(os.path.exists(newDirs[1] + "/OG%07d_tree.txt" % i))
            goldD = goldDirs[2]
            newD = newDirs[2]
            for goldFN in glob.glob(goldD + "*"):
                newFN = newD + os.path.split(goldFN)[1]
                self.assertTrue(os.path.exists(newFN), msg=goldFN)  
        self.test_passed = True         
#    
    def test_treesResultsChoice_full(self):
        dirs = ["Trees/", "Alignments/", "Sequences/"]
        goldDirs = [baseDir + "ExpectedOutput/Trees_OneSpeciesRemoved/" + d for d in dirs]
        expectedDirs = [baseDir + "/Input/Trees_OneSpeciesRemoved/" + d for d in dirs]
        nExpected = [536, 536, 1330]
        fnPattern = [baseDir + "/Input/Trees_OneSpeciesRemoved/" + p for p in ["Trees/OG%07d_tree.txt", "Alignments/OG%07d.fa", "Sequences/OG%07d.fa"]]
        with CleanUp([], [],  expectedDirs):
            self.stdout, self.stderr = self.RunTrees(baseDir + "/Input/Trees_OneSpeciesRemoved/clusters_OrthoFinder_v0.6.1_I1.5.txt_id_pairs.txt")
            for goldD, expD, n, fnPat in zip(goldDirs, expectedDirs, nExpected, fnPattern):
                for i in xrange(n):
                    self.assertTrue(os.path.exists(fnPat % i))
                # Only test the contents of a limited number
                for goldFN in glob.glob(goldD + "*"):
                    newFN = expD + os.path.split(goldFN)[1]
                    self.CompareFile(goldFN, newFN)   
        self.test_passed = True         
                    
    def test_ProblemCharacters_issue28(self):
        # Can't have ( ) : ,
        inputDir = baseDir + "Input/DisallowedCharacters/"
        with CleanUp([], [], [inputDir + x for x in "Trees/ Sequences/ Alignments/".split()]):        
            self.stdout, self.stderr = self.RunTrees(inputDir)
            fn = inputDir + "Trees/OG0000000_tree.txt"
            self.assertTrue(os.path.exists(fn))
            with open(fn, 'rb') as infile:
                l = infile.readline()
            # Tree exists: open braces, comma, close, end of tree
            self.assertTrue("(" in l)
            self.assertTrue(")" in l)
            self.assertTrue("," in l)
            self.assertTrue(";" in l)  
        self.test_passed = True      
        
    def test_DuplicateAccession(self):
        currentResultsDir = baseDir + "Input/SmallExampleDataset_dup_acc/" + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        with CleanUp([], [], [currentResultsDir]):        
            self.stdout, self.stderr = self.RunOrthoFinder("-f %s" % baseDir + "Input/SmallExampleDataset_dup_acc/")
            # Tree exists
            expectedFN = baseDir + "Input/SmallExampleDataset_dup_acc/" + ("Results_%s/" % datetime.date.today().strftime("%b%d")) + ("Orthologues_%s/" % datetime.date.today().strftime("%b%d")) + "Gene_Trees/OG0000000_tree.txt"
            self.assertTrue(os.path.exists(expectedFN))
            # No error
            self.assertTrue("ERROR" not in self.stdout)
            self.assertTrue(len(self.stderr) == 0)
            # run completes
            self.assertTrue("Orthogroup statistics:"  in self.stdout)
        
        
#    def test_treesExtraSpecies(self):
#        pass
        
    def setUp(self):
        self.currentResultsDir = None
        self.test_passed = False
        self.stdout = None
        self.stderr = None
        
    def tearDown(self):
        self.CleanCurrentResultsDir()
        if not self.test_passed:
            if qVerbose: print(self.stdout)
            print(self.stderr)
        
    def RunOrthoFinder(self, commands):
        if qBinary:
            capture = subprocess.Popen("%s %s" % (orthofinder_bin, commands), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        else:
            capture = subprocess.Popen("python %s %s" % (orthofinder, commands), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        stdout = "".join([x for x in capture.stdout])
        stderr = "".join([x for x in capture.stderr])
        return stdout, stderr
        
    def RunTrees(self, commands):
        if qBinary:
            capture = subprocess.Popen("%s -t 8 %s" % (trees_for_orthogroups_bin, commands), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        else:
            capture = subprocess.Popen("python %s -t 8 %s" % (trees_for_orthogroups, commands), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        stdout = "".join([x for x in capture.stdout])
        stderr = "".join([x for x in capture.stderr])
        return stdout, stderr
        
    def CheckStandardRun(self, stdout, stderr, goldResultsDir, expectedCSVFileLocation):
        # Citation
        self.assertTrue(citation in stdout)
        self.assertTrue("Connected putatitive homologs" in stdout) # see checks for errors, in that case this shouldn't be here
        
        # Results
        resultsLine = "Fifty percent of all genes were in orthogroups"
        self.assertTrue(resultsLine in stdout)
        self.assertLess(time.time()-os.stat(expectedCSVFileLocation).st_ctime, 60)
        
        # Results - orthogroups correct
        resultsDir = os.path.split(expectedCSVFileLocation)[0] + "/"
        self.CompareFile(goldResultsDir + "Orthogroups.csv", resultsDir + "Orthogroups.csv")
        self.CompareFile(goldResultsDir + "Orthogroups_UnassignedGenes.csv", resultsDir + "Orthogroups_UnassignedGenes.csv")
        self.CompareFile(goldResultsDir + "Orthogroups.txt", resultsDir + "Orthogroups.txt")
         
    def CleanCurrentResultsDir(self):
        if self.currentResultsDir == None: return
        # safety check, make sure directory was created in last 15 mins
        if time.time()-os.stat(self.currentResultsDir).st_ctime > 15*60:
            print("WARNING: Directory seems to predate test, it won't be deleted: %s" % self.currentResultsDir)
        else:
            shutil.rmtree(self.currentResultsDir)
            
    def CompareBlast(self, fn_gold, fn_actual, nMaxCompareLines = 10):
        with open(fn_gold, 'rb') as in1, open(fn_actual, 'rb') as in2:
            r1 = csv.reader(in1, delimiter="\t")
            r2 = csv.reader(in2, delimiter="\t")
            for ln1, ln2, _ in zip(r1, r2, xrange(nMaxCompareLines)):
                if ln1[0] != ln2[0]: 
                    shutil.copy(fn_actual, baseDir + "FailedOutput/" + os.path.split(fn_actual)[1]) 
                    self.assertTrue(False, msg=fn_gold) 
                if ln1[1] != ln2[1]: 
                    shutil.copy(fn_actual, baseDir + "FailedOutput/" + os.path.split(fn_actual)[1]) 
                    self.assertTrue(False, msg=fn_gold) 
                if abs(float(ln1[11]) - float(ln2[11]))  > 2.0: 
                    shutil.copy(fn_actual, baseDir + "FailedOutput/" + os.path.split(fn_actual)[1]) 
                    self.assertTrue(False, msg=fn_gold) 
        
    def CompareFile(self, fn_gold, fn_actual):
        fn_gold = os.path.split(fn_gold)[0] + "/" + os.path.split(fn_gold)[1].replace(version, "0.4.0")
        f = os.path.split(fn_actual)[1]                
        if "Blast" in f:
            self.CompareBlast(fn_gold, fn_actual)
        elif "Statistics" in f:
            self.CompareStatsFile(fn_gold, fn_actual)
        else: 
            if not filecmp.cmp(fn_gold, fn_actual):
                shutil.copy(fn_actual, baseDir + "FailedOutput/" + os.path.split(fn_actual)[1]) 
                self.assertTrue(False, msg=fn_gold) 
        
    def CompareFile_95Percent(self, fn_gold, fn_actual):
        with open(fn_gold, 'rb') as g, open(fn_actual, 'rb') as a:
            G = set(g.readlines())
            A = set(a.readlines())
        result = float(len(G.symmetric_difference(A))) / float((min(len(A), len(G))))
        self.assertLess(result, 0.05)
        
    def CompareStatsFile(self, fn_gold, fn_actual):
        with open(fn_gold, 'rb') as f_gold, open(fn_actual, 'rb') as f_actual:
            for gold, actual in zip(f_gold, f_actual):
                if "Date" in gold: continue
                if gold != actual:
                    shutil.copy(fn_actual, baseDir + "FailedOutput/" + os.path.split(fn_actual)[1]) 
                    self.assertTrue(False, msg=fn_gold) 
            
""" 
Test to add:
    - Extra results files when already files in directory (and clusters and graph file)
    - multiple -s arguments passed 
    - multiple -f arguments passed

"""  
      
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--binaries", action="store_true", help="Run tests on binary files")
    parser.add_argument("-t", "--test", help="Individual test to run")
    parser.add_argument("-d", "--dir", help="Test program in specified directory")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print stdout from failing orthofinder run")
    
    args = parser.parse_args()
    if args.dir:
        orthofinder = args.dir + "/orthofinder.py"
        orthofinder_bin = os.path.splitext(orthofinder)[0]
        trees_for_orthogroups = args.dir + "/trees_from_MSA.py"
        trees_for_orthogroups_bin = os.path.splitext(trees_for_orthogroups)[0]
        
    qBinary = args.binaries
    qVerbose = args.verbose
    print("Testing:")
    print("  " + orthofinder_bin if qBinary else orthofinder)
    print("  " + trees_for_orthogroups_bin if qBinary else trees_for_orthogroups)
    print("")
    
    if args.test != None:
        suite = unittest.TestSuite()
        print("Running single test: %s" % args.test)
        suite.addTest(TestCommandLine(args.test))
        runner = unittest.TextTestRunner()
        runner.run(suite)
    else:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestCommandLine)
        unittest.TextTestRunner(verbosity=2).run(suite)        
    
