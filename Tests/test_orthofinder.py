#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 11:25:10 2015

@author: david
"""

import os 
import sys
import time
import datetime
import unittest
import subprocess
import shutil
import filecmp
import glob
import random
import csv

__skipLongTests__ = False

baseDir = os.path.dirname(os.path.realpath(__file__)) + os.sep
orthofinder = baseDir + "../orthofinder.py"
trees_for_orthogroups = baseDir + "../trees_for_orthogroups.py"
exampleFastaDir = baseDir + "Input/SmallExampleDataset/"
exampleBlastDir = baseDir + "Input/SmallExampleDataset_ExampleBlastDir/"

goldResultsDir_smallExample = baseDir + "ExpectedOutput/SmallExampleDataset/"
goldPrepareBlastDir = baseDir + "ExpectedOutput/SmallExampleDataset_PreparedForBlast/"

version = "0.6.1"
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
    Infers orthologous groups for the proteomes contained in fasta_directory running
    number_of_blast_threads in parallel for the BLAST searches and subsequently running
    number_of_orthofinder_threads in parallel for the OrthoFinder algorithm.

Advanced Usage
--------------
python orthofinder.py -f fasta_directory -p
    1. Prepares files for BLAST and prints the BLAST commands. Does not perform BLAST searches
    or infer orthologous groups. Useful if you want to prepare the files in the form required by
    OrthoFinder but want to perform the BLAST searches using a job scheduler/on a cluster and
    then infer orthologous groups using option 2.

python orthofinder.py -b precalculated_blast_results_directory [-a number_of_orthofinder_threads]
    2. Infers orthologous groups using pre-calculated BLAST results. These can be after BLAST
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

-I inflation_parameter, --inflation inflation_parameter
    Specify a non-default inflation parameter for MCL. [Default is 1.5]

-x speciesInfoFilename, --orthoxml speciesInfoFilename
    Output the orthogroups in the orthoxml format using the information in speciesInfoFilename.

-p , --prepare
    Only prepare the files in the format required by OrthoFinder and print out the BLAST searches that
    need to be performed but don't run BLAST or infer orthologous groups

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
        - Results files should all contain the correct orthologous groups
    """
    @classmethod
    def setUpClass(cls):
        capture = subprocess.Popen("blastp -version", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        stdout = "".join([x for x in capture.stdout])
        if requiredBlastVersion not in stdout:
            raise RuntimeError("Tests require BLAST version %s" % requiredBlastVersion)
        
    
#    @unittest.skipIf(__skipLongTests__, "Only performing quick tests")     
    def test_fromfasta(self):
        currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        expectedCSVFile = currentResultsDir + "OrthologousGroups.csv"
        with CleanUp([], [], [currentResultsDir, ]):
            stdout, stderr = self.RunOrthoFinder("-f %s" % exampleFastaDir)
            self.CheckStandardRun(stdout, stderr, goldResultsDir_smallExample, expectedCSVFile)
            
    def test_fromfasta_threads(self):
        currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        expectedCSVFile = currentResultsDir + "OrthologousGroups.csv"
        with CleanUp([], [], [currentResultsDir, ]):
            stdout, stderr = self.RunOrthoFinder("-f %s -t 4 -a 3" % exampleFastaDir)
            self.CheckStandardRun(stdout, stderr, goldResultsDir_smallExample, expectedCSVFile)
           
    @unittest.skipIf(__skipLongTests__, "Only performing quick tests")     
    def test_fromfasta_full(self):
        currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        expectedCSVFile = currentResultsDir + "OrthologousGroups.csv"     
        with CleanUp([], [], [currentResultsDir, ]):
            stdout, stderr = self.RunOrthoFinder("--fasta %s" % exampleFastaDir)
            self.CheckStandardRun(stdout, stderr, goldResultsDir_smallExample, expectedCSVFile)
        
    def test_justPrepare(self):    
        self.currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        stdout, stderr = self.RunOrthoFinder("-f %s -p" % exampleFastaDir)
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
        
    def test_justPrepareFull(self):    
        self.currentResultsDir = exampleFastaDir + "Results_%s/" % datetime.date.today().strftime("%b%d") 
        stdout, stderr = self.RunOrthoFinder("-f %s --prepare" % exampleFastaDir)
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
        
    def test_fromblast(self):
        expectedCSVFile = exampleBlastDir + "OrthologousGroups.csv"
        newFiles = ("OrthologousGroups.csv OrthologousGroups_UnassignedGenes.csv OrthologousGroups.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt clusters_OrthoFinder_v%s_I1.5.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()
        newFiles = [exampleBlastDir + fn for fn in newFiles]
        with CleanUp(newFiles, []):
            stdout, stderr = self.RunOrthoFinder("-b %s" % exampleBlastDir)
            self.CheckStandardRun(stdout, stderr, goldResultsDir_smallExample, expectedCSVFile)
        
    def test_fromblast_full(self):
        expectedCSVFile = exampleBlastDir + "OrthologousGroups.csv"
        newFiles = ("OrthologousGroups.csv OrthologousGroups_UnassignedGenes.csv OrthologousGroups.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt clusters_OrthoFinder_v%s_I1.5.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()
        newFiles = [exampleBlastDir + fn for fn in newFiles]
        with CleanUp(newFiles, []):        
            stdout, stderr = self.RunOrthoFinder("--blast %s" % exampleBlastDir)
            self.CheckStandardRun(stdout, stderr, goldResultsDir_smallExample, expectedCSVFile)
        
    def test_fromblast_algthreads(self):
        expectedCSVFile = exampleBlastDir + "OrthologousGroups.csv"
        newFiles = ("OrthologousGroups.csv OrthologousGroups_UnassignedGenes.csv OrthologousGroups.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt clusters_OrthoFinder_v%s_I1.5.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()
        newFiles = [exampleBlastDir + fn for fn in newFiles]
        with CleanUp(newFiles, []):
            stdout, stderr = self.RunOrthoFinder("-b %s -a 3" % exampleBlastDir)
            self.CheckStandardRun(stdout, stderr, goldResultsDir_smallExample, expectedCSVFile)
        
    def test_blast_results_error(self):
        with CleanUp([], []):
            stdout, stderr = self.RunOrthoFinder("-a 2 -b " + baseDir + "Input/SmallExampleDataset_ExampleBadBlast/")
            self.assertTrue("Traceback" not in stderr)
            self.assertTrue("Offending line was:" in stderr)
            self.assertTrue("0_0	0_0	100.00	466" in stderr)
            self.assertTrue("Connected putatitive homologs" not in stdout)
            self.assertTrue("ERROR: An error occurred, please review previous error messages for more information." in stdout)
        
    def test_inflation(self):
        expectedCSVFile = exampleBlastDir + "OrthologousGroups.csv"
        newFiles = ("OrthologousGroups.csv OrthologousGroups_UnassignedGenes.csv OrthologousGroups.txt clusters_OrthoFinder_v%s_I1.8.txt_id_pairs.txt clusters_OrthoFinder_v%s_I1.8.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()
        newFiles = [exampleBlastDir + fn for fn in newFiles]
        with CleanUp(newFiles, []):
            stdout, stderr = self.RunOrthoFinder("-I 1.8 -b %s" % exampleBlastDir)
            self.CheckStandardRun(stdout, stderr, baseDir + "ExpectedOutput/SmallExampleDataset_I1.8/", expectedCSVFile)
    
#    @unittest.skipIf(__skipLongTests__, "Only performing quick tests")       
#    def test_fromblastOrthobench(self):
#        goldResultsDir_orthobench = baseDir + "ExpectedOutput/Orthobench_blast/"
#        expectedCSVFileLocation = baseDir + "Input/Orthobench_blast/OrthologousGroups.csv"
#        self.currentResultsDir = None
#        expectedNewFiles = [baseDir + "Input/Orthobench_blast/" + x for x in "OrthoFinder_v0.4.0_graph.txt clusters_OrthoFinder_v0.4.0_I1.5.txt clusters_OrthoFinder_v0.4.0_I1.5.txt_id_pairs.txt".split()]
#        with CleanUp(expectedNewFiles, []):
#            stdout, stderr = self.RunOrthoFinder("-b %sInput/Orthobench_blast" % baseDir)
#            self.CheckStandardRun(stdout, stderr, goldResultsDir_orthobench, expectedCSVFileLocation, qDelete = True)

    def test_help(self):
         stdout, stderr = self.RunOrthoFinder("")
         self.assertTrue(expectedHelp in stdout)
         self.assertEqual(len(stderr), 0)
         
         stdout, stderr = self.RunOrthoFinder("-h")
         self.assertTrue(expectedHelp in stdout)
         self.assertEqual(len(stderr), 0)
         
         stdout, stderr = self.RunOrthoFinder("--help")
         self.assertTrue(expectedHelp in stdout)
         self.assertEqual(len(stderr), 0)        
         
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
        OrthologousGroups.csv OrthologousGroups.txt OrthologousGroups_UnassignedGenes.csv \
        clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()]
        expectedChangedFiles = [exampleBlastDir + fn for fn in "SpeciesIDs.txt SequenceIDs.txt".split()]
        # cleanup afterwards including failed test
        goldDir = baseDir + "ExpectedOutput/AddOneSpecies/"
        with CleanUp(expectedExtraFiles, expectedChangedFiles):        
            stdout, stderr = self.RunOrthoFinder("-b %s -f %s" % (exampleBlastDir, baseDir + "Input/ExtraFasta"))
            # check extra blast files
            # check extra fasta file: simple couple of checks to ensure all ok
            for fn in expectedExtraFiles:
                os.path.split(fn)[1]
                self.assertTrue(os.path.exists(fn), msg=fn)
                # mcl output files contain a variable header, these files are an implementation detail that I don't want to test (I want the final orthogroups to be correct)
                if "clusters" in os.path.split(fn)[1]: continue
                if "Blast" in os.path.split(fn)[1]:
                    self.CompareBlast(goldDir + os.path.split(fn)[1], fn)
                else:
                    self.CompareFile(goldDir + os.path.split(fn)[1].replace(version, "0.4.0"), fn)               
    
    def test_addTwoSpecies(self):
        expectedExtraFiles = [exampleBlastDir + fn for fn in ("Blast0_3.txt Blast3_0.txt Blast1_3.txt Blast3_1.txt Blast2_3.txt Blast3_2.txt Blast3_3.txt Species3.fa \
        Blast0_4.txt Blast4_0.txt Blast1_4.txt Blast4_1.txt Blast2_4.txt Blast4_2.txt Blast3_4.txt Blast4_3.txt Blast4_4.txt Species4.fa \
        OrthologousGroups.csv OrthologousGroups.txt OrthologousGroups_UnassignedGenes.csv \
        clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()]
        expectedChangedFiles = [exampleBlastDir + fn for fn in "SpeciesIDs.txt SequenceIDs.txt".split()]
        # cleanup afterwards including failed test
        goldDir = baseDir + "ExpectedOutput/AddTwoSpecies/"
        with CleanUp(expectedExtraFiles, expectedChangedFiles):        
            stdout, stderr = self.RunOrthoFinder("-b %s -f %s" % (exampleBlastDir, baseDir + "Input/ExtraFasta2"))
            for fn in expectedExtraFiles:
                os.path.split(fn)[1]
                self.assertTrue(os.path.exists(fn), msg=fn)
                if "clusters" in os.path.split(fn)[1]: continue
                if "Blast" in os.path.split(fn)[1]:
                    self.CompareBlast(goldDir + os.path.split(fn)[1], fn)
                else:
                    self.CompareFile(goldDir + os.path.split(fn)[1].replace(version, "0.4.0"), fn)           
                
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
        requiredResults = [inputDir + fn for fn in "OrthologousGroups.csv OrthologousGroups_UnassignedGenes.csv OrthologousGroups.txt".split()]
        expectedExtraFiles = [inputDir + fn for fn in ("clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()]
        with CleanUp(expectedExtraFiles + requiredResults, []):
            stdout, stderr = self.RunOrthoFinder("-b %s" % inputDir)
            for fn in requiredResults:
                self.assertTrue(os.path.exists(fn), msg=fn)
                self.CompareFile(goldDir + os.path.split(fn)[1], fn)  
    
#    def test_removeMultipleSpecies(self):
#        pass
    
    def test_removeOneAddOne(self):
        inputDir = baseDir + "Input/ExampleDataset_addOneRemoveOne/Results_Jan28/WorkingDirectory/"
        expectedExtraFiles = [inputDir + fn for fn in ("Blast0_3.txt Blast3_0.txt Blast1_3.txt Blast3_1.txt Blast2_3.txt Blast3_2.txt Blast3_3.txt Species3.fa \
        OrthologousGroups.csv OrthologousGroups.txt OrthologousGroups_UnassignedGenes.csv \
        clusters_OrthoFinder_v%s_I1.5.txt clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt OrthoFinder_v%s_graph.txt" % (version, version, version)).split()] 
        expectedChangedFiles = [inputDir + fn for fn in "SpeciesIDs.txt SequenceIDs.txt".split()]
        goldDir = baseDir + "ExpectedOutput/AddOneRemoveOne/"
        with CleanUp(expectedExtraFiles, expectedChangedFiles):        
            stdout, stderr = self.RunOrthoFinder("-b %s -f %s" % (inputDir, baseDir + "Input/ExampleDataset_addOneRemoveOne/ExtraFasta/"))
#            print(stdout)
#            print(stderr)
            for fn in expectedExtraFiles:
                os.path.split(fn)[1]
                self.assertTrue(os.path.exists(fn), msg=fn)
                if "OrthologousGroups" in os.path.split(fn)[1]:
                    self.CompareFile(goldDir + os.path.split(fn)[1], fn)  
            self.CompareFile(goldDir + "SpeciesIDs.txt", inputDir + "SpeciesIDs.txt") 
            self.CompareFile(goldDir + "SequenceIDs.txt", inputDir + "SequenceIDs.txt")
    
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
            stdout, stderr = self.RunTrees(baseDir + "Input/SmallExampleDataset_forTrees/Results_Jan28/")
            for i in xrange(nTrees):
                self.assertTrue(os.path.exists(newDirs[0] + "/OG%07d.fa" % i), msg=str(i))
                self.assertTrue(os.path.exists(newDirs[1] + "/OG%07d_tree.txt" % i))
            goldD = goldDirs[2]
            newD = newDirs[2]
            for goldFN in glob.glob(goldD + "*"):
                newFN = newD + os.path.split(goldFN)[1]
                self.assertTrue(os.path.exists(newFN), msg=goldFN)
                self.CompareFile(goldFN, newFN)   
    
    def test_treesAfterOption_b(self):
        expectedDirs = "Alignments Trees Sequences".split()
        newDirs = [baseDir + "Input/SmallExampleDataset_forTreesFromBlast/" + d +"/" for d in expectedDirs]
        goldDirs = [baseDir + "ExpectedOutput/SmallExampleDataset_trees/" + d +"/" for d in expectedDirs]
        with CleanUp([], [], newDirs):        
            stdout, stderr = self.RunTrees(baseDir + "Input/SmallExampleDataset_forTreesFromBlast/")
            for goldD, newD in zip(goldDirs, newDirs):
                for goldFN in glob.glob(goldD + "*"):
                    newFN = newD + os.path.split(goldFN)[1]
                    self.assertTrue(os.path.exists(newFN), msg=goldFN)
                    if newD[:-1].rsplit("/", 1)[-1] == "Sequences":
                        self.CompareFile(goldFN, newFN)     
    
    def test_treesMultipleResults(self):
        expectedDirs = "Alignments Trees Sequences".split()
        newDirs = [baseDir + "Input/MultipleResults/Results_Jan28/WorkingDirectory/" + d +"/" for d in expectedDirs]
        goldDirs = [baseDir + "ExpectedOutput/MultipleResultsInDirectory/" + d +"/" for d in expectedDirs]
        with CleanUp([], [], newDirs):        
            stdout, stderr = self.RunTrees(baseDir + "Input/MultipleResults/Results_Jan28/WorkingDirectory/clusters_OrthoFinder_v0.4.0_I1.5.txt_id_pairs.txt")
            for goldD, newD in zip(goldDirs, newDirs):
                for goldFN in glob.glob(goldD + "*"):
                    newFN = newD + os.path.split(goldFN)[1]
                    self.assertTrue(os.path.exists(newFN), msg=goldFN)
                    if newD[:-1].rsplit("/", 1)[-1] == "Sequences":
                        self.CompareFile(goldFN, newFN)    
        
    def test_treesSubset_unspecified(self):
        stdout, stderr = self.RunTrees(baseDir + "/Input/Trees_OneSpeciesRemoved")
        self.assertTrue("ERROR: Results from multiple OrthoFinder runs found" in stdout)
        self.assertTrue("Please run with only one set of results in directories or specifiy the specific clusters_OrthoFinder_*.txt_id_pairs.txt file on the command line" in stdout)
    
    def test_treesResultsChoice_subset(self):   
        expectedDirs = "Alignments Trees Sequences".split()
        newDirs = [baseDir + "/Input/Trees_OneSpeciesRemoved/" + d +"/" for d in expectedDirs]
        goldDirs = [baseDir + "ExpectedOutput/SmallExampleDataset_trees/" + d + "/" for d in expectedDirs]
        nTrees = 427
        with CleanUp([], [], newDirs):   
            stdout, stderr = self.RunTrees(baseDir + "/Input/Trees_OneSpeciesRemoved/clusters_OrthoFinder_v0.6.1_I1.5_1.txt_id_pairs.txt")
            for i in xrange(nTrees):
                self.assertTrue(os.path.exists(newDirs[0] + "/OG%07d.fa" % i), msg=str(i))
                self.assertTrue(os.path.exists(newDirs[1] + "/OG%07d_tree.txt" % i))
            goldD = goldDirs[2]
            newD = newDirs[2]
            for goldFN in glob.glob(goldD + "*"):
                newFN = newD + os.path.split(goldFN)[1]
                self.assertTrue(os.path.exists(newFN), msg=goldFN)
#    
    def test_treesResultsChoice_full(self):
        dirs = ["Trees/", "Alignments/", "Sequences/"]
        goldDirs = [baseDir + "ExpectedOutput/Trees_OneSpeciesRemoved/" + d for d in dirs]
        expectedDirs = [baseDir + "/Input/Trees_OneSpeciesRemoved/" + d for d in dirs]
        nExpected = [536, 536, 1330]
        fnPattern = [baseDir + "/Input/Trees_OneSpeciesRemoved/" + p for p in ["Trees/OG%07d_tree.txt", "Alignments/OG%07d.fa", "Sequences/OG%07d.fa"]]
        with CleanUp([], [],  expectedDirs):
            stdout, stderr = self.RunTrees(baseDir + "/Input/Trees_OneSpeciesRemoved/clusters_OrthoFinder_v0.6.1_I1.5.txt_id_pairs.txt")
            for goldD, expD, n, fnPat in zip(goldDirs, expectedDirs, nExpected, fnPattern):
                for i in xrange(n):
                    self.assertTrue(os.path.exists(fnPat % i))
                # Only test the contents of a limited number
                for goldFN in glob.glob(goldD + "*"):
                    newFN = expD + os.path.split(goldFN)[1]
                    self.CompareFile(goldFN, newFN) 
        
#    def test_treesExtraSpecies(self):
#        pass
        
    def setUp(self):
        self.currentResultsDir = None
        
    def tearDown(self):
        self.CleanCurrentResultsDir()
        
    def RunOrthoFinder(self, commands):
        capture = subprocess.Popen("python %s %s" % (orthofinder, commands), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        stdout = "".join([x for x in capture.stdout])
        stderr = "".join([x for x in capture.stderr])
        return stdout, stderr
        
    def RunTrees(self, commands):
        capture = subprocess.Popen("python %s -t 8 %s" % (trees_for_orthogroups, commands), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        stdout = "".join([x for x in capture.stdout])
        stderr = "".join([x for x in capture.stderr])
        return stdout, stderr
        
    def CheckStandardRun(self, stdout, stderr, goldResultsDir, expectedCSVFileLocation):
        # Citation
        self.assertTrue(citation in stdout)
        self.assertTrue("Connected putatitive homologs" in stdout) # see checks for errors, in that case this shouldn't be here
        
        # Results
        resultsLine = """Orthologous groups have been written to tab-delimited files:
   %s""" % expectedCSVFileLocation
        self.assertTrue(resultsLine in stdout)
        self.assertLess(time.time()-os.stat(expectedCSVFileLocation).st_ctime, 60)
        
        # Results - orthogroups correct
        resultsDir = os.path.split(expectedCSVFileLocation)[0] + "/"
        self.CompareFile(goldResultsDir + "OrthologousGroups.csv", resultsDir + "OrthologousGroups.csv")
        self.CompareFile(goldResultsDir + "OrthologousGroups_UnassignedGenes.csv", resultsDir + "OrthologousGroups_UnassignedGenes.csv")
        self.CompareFile(goldResultsDir + "OrthologousGroups.txt", resultsDir + "OrthologousGroups.txt")
         
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
        if not filecmp.cmp(fn_gold, fn_actual):
            shutil.copy(fn_actual, baseDir + "FailedOutput/" + os.path.split(fn_actual)[1]) 
            self.assertTrue(False, msg=fn_gold) 
            
""" 
Test:
    - Extra results files when already files in directory (and clusters and graph file)
    - multiple -s arguments passed 
    - multiple -f arguments passed

"""        
if __name__ == "__main__":
    if len(sys.argv) == 2: 
        # run single test
        suite = unittest.TestSuite()
        suite.addTest(TestCommandLine(sys.argv[1]))
        runner = unittest.TextTestRunner()
        runner.run(suite)
    else:
        # run suite
        suite = unittest.TestLoader().loadTestsFromTestCase(TestCommandLine)
        unittest.TextTestRunner(verbosity=2).run(suite)
    
