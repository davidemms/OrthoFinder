#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2016 David Emms
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
# david_emms@hotmail.comhor: david


"""
Handles location of all input and output files
Users:
1. Call InitialiseFileHandler
2. Interact with FileHandler object (which is initialised by the above method)
    - it is an instance of __Files_new_dont_manually_create__

The code also contains the help class: PreviousFilesLocator (and child classes of it)
"""

import os
import sys
import glob
import time
import shutil
import datetime

from . import util

class SpeciesInfo(object):
    def __init__(self):
        self.speciesToUse = []           #       seqsInfo.iSpeciesToUse   - which to include for this analysis 
        self.nSpAll = None               #       seqsInfo.nSpAll => 0, 1, ..., nSpAll - 1 are valid species indices
        self.iFirstNewSpecies = None     #       iFirstNew   => (0, 1, ..., iFirstNew-1) are from previous and (iFirstNew, iFirstNew+1, ..., nSpecies-1) are the new species indices
    def __str__(self):
        return str((self.speciesToUse, self.nSpAll, self.iFirstNewSpecies))
   
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """
                    
class __Files_new_dont_manually_create__(object):    
    def __init__(self):
        self.baseOgFormat = "OG%07d"
        self.wd_base = []               # Base: blast, species & sequence IDs, species fasta files - should not request this and then write here
        self.wd_current = None          # Location to write out any new files
        self.wd_trees = None            # Location of working dir containing tree files 
        self.rd1 = None
        self.fileIdentifierString = "OrthoFinder"
        self.clustersFilename = None
        self.iResultsVersion = None
        self.nondefaultPickleDir = None
        self.speciesTreeRootedIDsFN = None
        self.multipleRootedSpeciesTreesDir = None
        self.species_ids_corrected = None
        # to be modified as appropriate
     
    """ ========================================================================================== """
    # RefactorDS - FileHandler
    def CreateOutputDirFromStart_new(self, fasta_dir, base, user_name = None, old_wd_base_list=None):
        """
        The initial difference will be that results will go in OrthoFinder/Results_DATE or USER_SPECIFIED/RESULTS_DATE
        whereas before they went in Results_DATE or USER_SPECIFIED.
        
        If this is a composite analysis (-f + -b) then old_wd_base_list != None 
        
        old_wd_base_list - first item is the WD from a previous analysis to be extended. If this extended other
          ones itself then there will be other items in the list.
        """
        if user_name == None:
            self.rd1 = util.CreateNewWorkingDirectory(base + "Results_")
        else:
            self.rd1 = util.CreateNewWorkingDirectory(base + "Results_" + user_name, qDate=False)
        self.wd_current = self.rd1 + "WorkingDirectory/"
        os.mkdir(self.wd_current)
        self.wd_base = [self.wd_current]        
        if old_wd_base_list != None:
            shutil.copy(old_wd_base_list[0] + "SpeciesIDs.txt", self.wd_current + "SpeciesIDs.txt")
            shutil.copy(old_wd_base_list[0] + "SequenceIDs.txt", self.wd_current + "SequenceIDs.txt")
            # Log the first wd in list, this can then be followed back to previous ones
            # Log file - point to WD at start of chain which contains the new species
            # wd_base_list - should contain current directory and then previous linked directories
            with open(self.wd_current + "previous_wd.txt", 'w') as outfile: outfile.write(old_wd_base_list[0] + "\n")
            self.wd_base.extend(old_wd_base_list)
        self.wd_trees = self.wd_current
        self.StartLog()
        
    # RefactorDS - PreviousFilesLocator
    def StartFromOrthogroupsOrSequenceSearch(self, wd_base_list, base, clustersFilename_pairs=None, user_name = None, userSpeciesTree=None):
        """
        NEed to initialise:
        wd_base
        wd_trees
        wd_current
        """
        if len(self.wd_base) != 0: raise Exception("Changing WorkingDirectory1")
        self.wd_base = wd_base_list
        if clustersFilename_pairs != None: self.clustersFilename = clustersFilename_pairs[:-len("_id_pairs.txt")]
        if user_name == None:
            self.rd1 = util.CreateNewWorkingDirectory(base + "Results_")
        else:
            self.rd1 = util.CreateNewWorkingDirectory(base + "Results_" + user_name, qDate=False)
        self.wd_current = self.rd1 + "WorkingDirectory/"
        os.mkdir(self.wd_current)
        with open(self.rd1 + "Log.txt", 'w'):
            pass
        self.wd_trees = self.wd_current
        self.StartLog()
    
    
    def StartFromTrees(self, 
                       wd1_list, 
                       wd2,
                       base, 
                       clustersFilename_pairs,
                       speciesTreeFN, 
                       qIsUSerSpeciesTree,
                       user_name=None):
        """
        Convert user species tree here if necessary
        For OF species tree copy it to location given by FileHandler
        For user species tree, this must be done immediately by OF code
        """
        self.wd_base = wd1_list
        self.wd_trees = wd2
        if user_name == None:
            self.rd1 = util.CreateNewWorkingDirectory(base + "Results_")
        else:
            self.rd1 = util.CreateNewWorkingDirectory(base + "Results_" + user_name, qDate=False)
        self.wd_current = self.rd1 + "WorkingDirectory/"
        os.mkdir(self.wd_current)
        self.clustersFilename = clustersFilename_pairs[:-len("_id_pairs.txt")]
        self.StartLog()
        if not qIsUSerSpeciesTree:
            shutil.copy(speciesTreeFN, self.GetSpeciesTreeIDsRootedFN())
        self.WriteToLog("Species Tree: %s\n" % speciesTreeFN)
        self.LogWorkingDirectoryTrees()
                                         
    def CreateOutputDirectories(self, options, previous_files_locator, base_dir, fastaDir=None):
        if options.qStartFromFasta and options.qStartFromBlast:
            wd1 = previous_files_locator.GetStartFromBlast()
            self.CreateOutputDirFromStart_new(fastaDir, base_dir, user_name=options.name, old_wd_base_list = wd1)
        
        elif options.qStartFromFasta:
            self.CreateOutputDirFromStart_new(fastaDir, base_dir, user_name=options.name)
            
        elif options.qStartFromBlast: 
            wd1 = previous_files_locator.GetStartFromBlast()
            self.StartFromOrthogroupsOrSequenceSearch(wd1, 
                                                      base_dir,
                                                      user_name=options.name)  
            
        elif options.qStartFromGroups:
            wd1, clustersFilename_pairs = previous_files_locator.GetStartFromOGs()
            self.StartFromOrthogroupsOrSequenceSearch(wd1, 
                                                      base_dir,
                                                      clustersFilename_pairs, 
                                                      user_name=options.name)
                                                      
        elif options.qStartFromTrees:
            wd1, clustersFilename_pairs, wd_trees, speciesTreeFN = previous_files_locator.GetStartFromTrees()
            if options.speciesTreeFN != None:
                qIsUserSpeciesTree = True
                speciesTreeFN = options.speciesTreeFN
            elif speciesTreeFN != None:
                qIsUserSpeciesTree = False
            else:
                print("ERROR: Could not find species tree")
                util.Fail()
            self.StartFromTrees(wd1, 
                                wd_trees,
                                base_dir, 
                                clustersFilename_pairs,
                                speciesTreeFN, 
                                qIsUserSpeciesTree,
                                user_name=options.name)
        if (options.qStartFromGroups or options.qStartFromTrees) and previous_files_locator.species_ids_lines != None:
            # In only these cases, it's possible that the SpeciesIDs.txt file is out of sync and the version in the previous log should be used instead
            self.CreateCorrectedSpeciesIDsFile(previous_files_locator.species_ids_lines)

    def CreateCorrectedSpeciesIDsFile(self, species_ids_lines):
        self.species_ids_corrected = self.wd_current + "SpeciesIDs.txt"
        with open(self.species_ids_corrected, 'w') as outfile:
            outfile.write(species_ids_lines)
                                                      
    """ ========================================================================================== """
       
       
    # RefactorDS - FileHandler
    def SetNondefaultPickleDir(self, d):
        self.pickleDir = d
   
    def GetPickleDir(self):
        if self.nondefaultPickleDir != None: 
            d = self.pickleDir
        else:
            d = self.wd_current + "pickle/"
        if not os.path.exists(d): os.mkdir(d)
        return d
            
    """ Standard Methods
        ========================================================================================== """                       
    def LogSpecies(self):
        text = "\nSpecies used: \n"
        fn = self.GetSpeciesIDsFN()
        with open(fn, 'r') as infile:
            text += "".join(infile.readlines())
        self.WriteToLog(text + "\n")
        
    """ Standard Directories
        ========================================================================================== """
    
    def GetWorkingDirectory1_Read(self):
        if len(self.wd_base) == 0: raise Exception("No wd1")
        return self.wd_base 
        
    def GetWorkingDirectory_Write(self):
        if self.wd_current == None: raise Exception("No wd_current")
        return self.wd_current 
        
    def GetResultsDirectory1(self):
        if self.rd1 == None: raise Exception("No rd1")
        return self.rd1 
        
    def GetResultsDirectory2(self):
        if self.rd1 == None: raise Exception("No rd1")
        return self.rd1 
        
    def GetOrthologuesDirectory(self):
        """"Where the directories of species orthologues are"""
        if self.rd1 == None: raise Exception("No rd1")
        d = self.rd1 + "Orthologues/"
        if not os.path.exists(d): os.mkdir(d)
        return d
        
    """ Orthogroups files 
        ========================================================================================== """
    
    def GetSpeciesIDsFN(self):
        if self.species_ids_corrected != None:
            return self.species_ids_corrected
        return self.wd_base[0] + "SpeciesIDs.txt"
        
    def GetSequenceIDsFN(self):
        # It is always in the first of the 'extension' directories (as this is the relevant one)
        return self.wd_base[0] + "SequenceIDs.txt"
    
#    def GetSpeciesIDsFN(self):
#        if self.species_ids_corrected != None:
#            return self.species_ids_corrected
#        if len(self.wd_base) == 0: raise Exception("No wd1")
#        for d in self.wd_base:
#            fn = d + "SpeciesIDs.txt"
#            if os.path.exists(fn): return fn
#        raise Exception(fn + " not found")
#        
#    def GetSequenceIDsFN(self):
#        if len(self.wd_base) == 0: raise Exception("No wd1")
#        for d in self.wd_base:
#            fn = d + "SequenceIDs.txt"
#            if os.path.exists(fn): return fn
#        raise Exception(fn + " not found")
        
    def GetSpeciesSeqsDir(self):
        if len(self.wd_base) == 0: raise Exception("No wd1")
        return self.wd_base 
        
    def GetSpeciesFastaFN(self, iSpecies, qForCreation=False):
        """
        qForCreation: A path is required at which the file should be created (don't search for it)
        """
        if len(self.wd_base) == 0: raise Exception("No wd1")
        if qForCreation:
            return "%sSpecies%d.fa" % (self.wd_base[0], iSpecies)
        for d in self.wd_base:
            fn = "%sSpecies%d.fa" % (d, iSpecies)
            if os.path.exists(fn): return fn
        raise Exception(fn + " not found")
        
    def GetSortedSpeciesFastaFiles(self):
        if len(self.wd_base) == 0: raise Exception("No wd1")
        fastaFilenames = []
        for d in self.wd_base:
            fastaFilenames.extend(glob.glob(d + "Species*.fa"))
        speciesIndices = []
        for f in fastaFilenames:
            start = f.rfind("Species")
            speciesIndices.append(int(f[start+7:-3]))
        indices, sortedFasta = util.SortArrayPairByFirst(speciesIndices, fastaFilenames)
        return sortedFasta  
        
    def GetSpeciesDatabaseN(self, iSpecies, program="Blast"):
        return "%s%sDBSpecies%d" % (self.wd_current, program, iSpecies)
    
    def GetBlastResultsDir(self):
        return self.wd_base
        
    def GetBlastResultsFN(self, iSpeciesSearch, jSpeciesDB, qForCreation=False):
        if len(self.wd_base) == 0: raise Exception("No wd1")
        if qForCreation: return "%sBlast%d_%d.txt" % (self.wd_base[0], iSpeciesSearch, jSpeciesDB)     
        for d in self.wd_base:
            fn = "%sBlast%d_%d.txt" % (d, iSpeciesSearch, jSpeciesDB)
            if os.path.exists(fn) or os.path.exists(fn + ".gz"): return fn
        raise Exception(fn + " not found")
        
    def GetGraphFilename(self):
        if self.wd_current == None: raise Exception("No wd_current")
        return self.wd_current + "%s_graph.txt" % self.fileIdentifierString
        
    def CreateUnusedClustersFN(self, mclInflation):
        if self.wd_current == None: raise Exception("No wd_current")
        self.clustersFilename, self.iResultsVersion = util.GetUnusedFilename(self.wd_current  + "clusters_%s_I%0.1f" % (self.fileIdentifierString, mclInflation), ".txt")
        return self.clustersFilename, self.clustersFilename + "_id_pairs.txt"
        
    def SetClustersFN(self, pairsFN):
        self.clustersFilename = pairsFN[:-len("_id_pairs.txt")]
        log = "Orthogroups used: %s\n\n" % self.clustersFilename
        self.WriteToLog(log)
        
    def GetClustersFN(self):
        return self.clustersFilename + "_id_pairs.txt"
        
    """ Orthologues files
        ========================================================================================== """
    
    def GetResultsSeqsDir_SingleCopy(self):
        d = self.rd1 + "Single_Copy_Orthologue_Sequences/"
        if not os.path.exists(d): os.mkdir(d)
        return d        
        
    def GetResultsSeqsDir(self):
        return self.rd1 + "Orthogroup_Sequences/"
        
    def GetResultsAlignDir(self):
        return self.rd1 + "MultipleSequenceAlignments/"
        
    def GetResultsTreesDir(self):
        return self.rd1 + "Gene_Trees/"
    
    def GetOGsSeqFN(self, iOG, qResults=False):
        if qResults:
            return self.rd1 + "Orthogroup_Sequences/" + (self.baseOgFormat % iOG) + ".fa"
        else:
            return self.wd_current + "Sequences_ids/" + (self.baseOgFormat % iOG) + ".fa"
            
    def GetOGsAlignFN(self, iOG, qResults=False):
        if qResults:
            return self.rd1 + "MultipleSequenceAlignments/" + (self.baseOgFormat % iOG) + ".fa"
        else:
            return self.wd_current + "Alignments_ids/" + (self.baseOgFormat % iOG) + ".fa"
            
    def GetOGsTreeFN(self, iOG, qResults=False):
        if qResults:
            return self.rd1 + "Gene_Trees/" + (self.baseOgFormat % iOG) + "_tree.txt"
        else:
            return self.wd_trees + "Trees_ids/" + (self.baseOgFormat % iOG) + "_tree_id.txt"   
        
    def GetSpeciesTreeConcatAlignFN(self, qResults=False):
        if qResults:
            return self.rd1 + "MultipleSequenceAlignments/" + "SpeciesTreeAlignment.fa"
        else:
            return self.wd_current + "Alignments_ids/SpeciesTreeAlignment.fa"  
        
    def GetSpeciesTreeMatrixFN(self, qPutInWorkingDir = False):
        if qPutInWorkingDir:
            return self.wd_current + "SpeciesMatrix.phy"
        else:
            return self.wd_current + "Distances/SpeciesMatrix.phy"
            
    def GetSpeciesTreeUnrootedFN(self, qAccessions=False):
        if qAccessions:
            return self.wd_trees + "SpeciesTree_unrooted.txt"
        else: 
            return self.wd_trees + "SpeciesTree_unrooted_ids.txt"  
                        
    def GetSpeciesTreeIDsRootedFN(self):
        return self.wd_current + "SpeciesTree_rooted_ids.txt"
            
    def GetSpeciesTreeResultsFN(self, i, qUnique):
        """
        The results species tree (rooted, accessions, support values)
        i: index for species tree, starting at 0
        qUnique: bool, has a unique root been identified (as it may not be known exatly which branch the root belongs on)
        E.g. if there were just one species tree, the correct call would be GetSpeciesTreeResultsFN(0,True)
        """
        d = self.rd1 + "Species_Tree/"
        if not os.path.exists(d): os.mkdir(d)
        if qUnique:
            return d + "SpeciesTree_rooted.txt"
        else:
            if not self.multipleRootedSpeciesTreesDir:
                self.multipleRootedSpeciesTreesDir = d + "Potential_Rooted_Species_Trees/"
                if not os.path.exists(self.multipleRootedSpeciesTreesDir): os.mkdir(self.multipleRootedSpeciesTreesDir)
            return self.multipleRootedSpeciesTreesDir + "SpeciesTree_rooted_at_outgroup_%d.txt" % i    
            
    def GetSpeciesTreeResultsNodeLabelsFN(self):
        return self.GetSpeciesTreeResultsFN(0, True)[:-4] + "_node_labels.txt"
        
    def GetOGsDistMatFN(self, iOG):
        return self.wd_current + "Distances/OG%07d.phy" % iOG
        
    def GetSpeciesDict(self):
        d = util.FullAccession(self.GetSpeciesIDsFN()).GetIDToNameDict()
        return {k:v.rsplit(".",1)[0] for k,v in d.items()}

    def GetHierarchicalOrthogroupsFN(self, sp_node_name):
        return self.rd1 + "Phylogenetic_Hierarchical_Orthogroups/%s.tsv" % sp_node_name
        
    """ ========================================================================================== """
            
    def GetOGsTreeDir(self, qResults=False):
        if qResults:
            return self.rd1 + "Gene_Trees/" 
        else:
            return self.wd_trees + "Trees_ids/" 
            
    def GetOGsReconTreeDir(self, qResults=False):
        if qResults:
            d = self.rd1 + "Resolved_Gene_Trees/" 
            if not os.path.exists(d): os.mkdir(d)
            return d
        else:
            raise NotImplemented() 
            
    def GetOGsReconTreeFN(self, iOG):
        return self.rd1 + "Resolved_Gene_Trees/OG%07d_tree.txt" % iOG
            
    def GetPhyldogWorkingDirectory(self):
        d = self.wd_current + "phyldog/"
        if not os.path.exists(d): os.mkdir(d)
        return d
            
    def GetPhyldogOGResultsTreeFN(self, i):
        return self.wd_current + "phyldog/OG%07d.ReconciledTree.txt" % i
    
    """ ========================================================================================== """
         
    def CleanWorkingDir2(self):
        dirs = ['Distances/']
        for d in dirs:
            dFull = self.wd_current + d
            if os.path.exists(dFull): 
                try:
                    shutil.rmtree(dFull)
                except OSError:
                    time.sleep(1)
                    shutil.rmtree(dFull, True)  # shutil / NFS bug - ignore errors, it's less crucial that the files are deleted
                    
    """ ************************************************************************************************************************* """
            
# RefactorDS - FileHandler 
    """ Standard Methods ========================================================================================== """  
    def LogFailAndExit(self, text=""):
        if text != "": print(text)
        self.WriteToLog("\nERROR: An error occurred\n" + text)
        util.Fail()
             
    def WriteToLog(self, text, qWithTime=False):
        prepend = ""
        if qWithTime:
            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + " : "
        with open(self.rd1 + "Log.txt", 'a') as outfile:
            outfile.write(prepend + text)
    
    def StartLog(self):
        self.WriteToLog("Started OrthoFinder version " + util.version + "\n", True)
        text = "Command Line: " + " ".join(sys.argv) + "\n\n"
        text += "WorkingDirectory_Base: %s\n" % self.wd_base[0]
        self.WriteToLog(text)
        if self.clustersFilename != None:self.LogOGs()
    
    def LogOGs(self):
        self.WriteToLog("FN_Orthogroups: %s\n" % (self.clustersFilename + "_id_pairs.txt"))
    
    def LogWorkingDirectoryTrees(self):
        self.WriteToLog("WorkingDirectory_Trees: %s\n" % self.wd_trees)
        
    def MakeResultsDirectory2(self, tree_generation_method, stop_after="", append_name=""):
        """
        Args
        tree_method: msa, dendroblast, phyldog (determines the directory structure that will be created)
        stop_after: seqs, align
        """
        # RefactorDS - need to change where it puts things
        if self.rd1 == None: raise Exception("No rd1") 
        self.wd_trees = self.wd_current
        os.mkdir(self.rd1 + "Orthologues/")
        if tree_generation_method == "msa":
            for i, d in enumerate([self.GetResultsSeqsDir(), self.wd_current + "Sequences_ids/", self.GetResultsAlignDir(), self.wd_current + "Alignments_ids/", self.GetResultsTreesDir(), self.wd_current + "Trees_ids/"]):
                if stop_after == "seqs" and i == 2: break 
                if stop_after == "align" and i == 4: break 
                if not os.path.exists(d): os.mkdir(d)
        elif tree_generation_method == "dendroblast":
            for i, d in enumerate([self.wd_current + "Distances/", self.GetResultsTreesDir(), self.wd_current + "Trees_ids/"]):
                if not os.path.exists(d): os.mkdir(d)
    
    def GetOrthogroupResultsFNBase(self):
        if self.rd1 == None: 
            raise Exception("No rd1")
        if self.iResultsVersion == None:
            raise Exception("Base results identifier has not been created")
        d = self.rd1 + "Orthogroups/"
        if not os.path.exists(d): os.mkdir(d)
        return d + "Orthogroups" + ("" if self.iResultsVersion == 0 else "_%d" % self.iResultsVersion)
        
    def GetOGsStatsResultsDirectory(self):
        d = self.rd1 + "Comparative_Genomics_Statistics/"
        if not os.path.exists(d): os.mkdir(d)
        return d
        
    def GetDuplicationsFN(self):
        d = self.rd1 + "Gene_Duplication_Events/"
        if not os.path.exists(d): os.mkdir(d)
        return d + "Duplications.tsv"
    
    def GetSuspectGenesDir(self):
        d = self.rd1 + "Phylogenetically_Misplaced_Genes/"
        if not os.path.exists(d): os.mkdir(d)
        return d
        
    def GetPutativeXenelogsDir(self):
        d = self.rd1 + "Putative_Xenologs/"
        if not os.path.exists(d): os.mkdir(d)
        return d    

FileHandler = __Files_new_dont_manually_create__()
                    
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """

class Unprocessable(Exception):
    pass

class PreviousFilesLocator(object):
    def __init__(self):
        self.wd_base_prev = []
        self.clustersFilename_pairs = None
        self.wd_trees = None
        self.home_for_results = None
        self.speciesTreeRootedIDsFN = None
        self.species_ids_lines = None
                
    def GetHomeForResults(self):
        return self.home_for_results

    def GetStartFromBlast(self):
        return self.wd_base_prev

    def GetStartFromOGs(self):
        return self.wd_base_prev, self.clustersFilename_pairs

    def GetStartFromTrees(self):
        return self.wd_base_prev, self.clustersFilename_pairs, self.wd_trees, self.speciesTreeRootedIDsFN
        
""" ************************************************************************************************************************* """

class PreviousFilesLocator_new(PreviousFilesLocator):
    def __init__(self, options, continuationDir):
        PreviousFilesLocator.__init__(self)
        if not continuationDir.endswith("/"): continuationDir += "/"
        self.home_for_results = continuationDir + "../"
        if (options.qStartFromFasta and not options.qStartFromBlast):
            # there are no files to find
            return
        if not self._IsNewDirStructure(continuationDir): raise Unprocessable("Input directory structure is not processable as new structure")
        self._ProcessLog(continuationDir + "/Log.txt")
    
    def _IsNewDirStructure(self, inputDir):
        return os.path.exists(inputDir + "/Log.txt") 
    
    def _ProcessLog(self, logFN):
        """
        Get all relevant data from log file. 
        Checks the paths ssaved do exist still
        Should work with relevant paths to allow directory to move
        Other methods can then check that the data required for a particualr run is available
        """
        with open(logFN, 'r') as infile:
            for line in infile:
                if line.startswith("Species used:"):
                    self.species_ids_lines = ""
                    line = next(infile)
                    while line.rstrip() != "":
                        self.species_ids_lines += line
                        line = next(infile)
                wd_base_str = "WorkingDirectory_Base: "
                wd_trees_str = "WorkingDirectory_Trees: "
                clusters_str = "FN_Orthogroups: "
                if line.startswith(wd_base_str): 
                    wd_base_anchor = line.rstrip()[len(wd_base_str):]
                    if not os.path.exists(wd_base_anchor):
                        # try to see if it's a relative directory to current one
                        path, d_wd = os.path.split(wd_base_anchor[:-1])
                        path, d_res = os.path.split(path)
                        wd_base_anchor = os.path.split(logFN)[0] + ("/../%s/%s/" % (d_res, d_wd))
                        if not os.path.exists(wd_base_anchor):
                            print("ERROR: Missing directory: %s" % wd_base_anchor)
                            util.Fail()
                    self.wd_base_prev = self.GetWDBaseChain(wd_base_anchor)
                if line.startswith(clusters_str): 
                    clusters_fn_full_path = line.rstrip()[len(clusters_str):]
                    self.clustersFilename_pairs = clusters_fn_full_path 
                    if not os.path.exists(self.clustersFilename_pairs):
                        # try to see if it's a relative directory to current one
                        path, clusters_fn = os.path.split(self.clustersFilename_pairs)
                        path, d_wd = os.path.split(path)
                        path, d_res = os.path.split(path)
                        self.clustersFilename_pairs = os.path.split(logFN)[0] + ("/../%s/%s/%s" %(d_res, d_wd, clusters_fn))
                        if not os.path.exists(self.clustersFilename_pairs):
                            print("ERROR: Missing orthogroups file: %s or %s" % (self.clustersFilename_pairs, clusters_fn_full_path))
                            util.Fail()
#                    self._GetOGsFile(wd_ogs_path)
                if line.startswith(wd_trees_str): 
                    self.wd_trees = line.rstrip()[len(wd_trees_str):]
                    if not os.path.exists(self.wd_trees):
                        # try to see if it's a relative directory to current one
                        path, d_wd = os.path.split(self.wd_trees[:-1])
                        path, d_res = os.path.split(path)
                        self.wd_trees = os.path.split(logFN)[0] + ("/../%s/%s/" % (d_res, d_wd))
                        if not os.path.exists(self.wd_trees):
                            print("ERROR: Missing directory: %s" % self.wd_trees)
                            util.Fail()
                    self.speciesTreeRootedIDsFN = self.wd_trees + "SpeciesTree_rooted_ids.txt" 
                        
                            
    def GetWDBaseChain(self, wd_base_anchor):
        chain = [wd_base_anchor]
        while os.path.exists(chain[-1] + "previous_wd.txt"):
            with open(chain[-1] + "previous_wd.txt", 'r') as infile:
                wd = infile.readline().rstrip()
                if not os.path.exists(wd):
                    # try to see if it's a relative directory to current one
                    path, d_wd = os.path.split(wd[:-1])
                    path, d_res = os.path.split(path)
                    wd = wd_base_anchor + ("/../../%s/%s/" % (d_res, d_wd))
                chain.append(wd)
        return chain
                
            
""" ************************************************************************************************************************* """

class PreviousFilesLocator_old(PreviousFilesLocator):
    def __init__(self, options, continuationDir):
        PreviousFilesLocator.__init__(self)
        if not continuationDir.endswith("/"): continuationDir += "/"
        self.home_for_results = continuationDir + "OrthoFinder/"
        if options.qStartFromGroups or options.qStartFromTrees:
            # User can specify it using clusters_id_pairs file, process this first to get the workingDirectory
            ogs_dir = continuationDir + "../" if options.qStartFromTrees else continuationDir
            self.wd_base_prev, self.orthofinderResultsDir, self.clustersFilename_pairs = self._GetOGsFile(ogs_dir)
            if options.qStartFromTrees:
                self._FindFromTrees(continuationDir, options.speciesTreeFN)
        elif options.qStartFromBlast:
            if self._IsWorkingDirectory(continuationDir): 
                self.wd_base_prev = continuationDir
            elif self._IsWorkingDirectory(continuationDir + "WorkingDirectory/"):
                self.wd_base_prev = continuationDir + "WorkingDirectory/"
            else:
                self.wd_base_prev = continuationDir   # nothing much to do, set this as the one to try and fail later
        self.wd_base_prev = [self.wd_base_prev]
                
    def _GetOGsFile(self, userArg):
        """returns the WorkingDirectory, ResultsDirectory and clusters_id_pairs filename"""
        qSpecifiedResultsFile = False
        if userArg == None:
            print("ERROR: orthofinder_results_directory has not been specified")
            util.Fail()
        if os.path.isfile(userArg):
            fn = os.path.split(userArg)[1]
            if ("clusters_OrthoFinder_" not in fn) or ("txt_id_pairs.txt" not in fn):
                print("ERROR:\n    %s\nis neither a directory or a clusters_OrthoFinder_*.txt_id_pairs.txt file." % userArg)
                util.Fail()
            qSpecifiedResultsFile = True
            # user has specified specific results file
        elif userArg[-1] != os.path.sep: 
            userArg += os.path.sep
        
        # find required files
        if qSpecifiedResultsFile:
            orthofinderWorkingDir = os.path.split(userArg)[0] + os.sep
            if not self._IsWorkingDirectory(orthofinderWorkingDir):
                print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s" % orthofinderWorkingDir)
                util.Fail()
        else:
            orthofinderWorkingDir = os.path.split(userArg)[0] if qSpecifiedResultsFile else userArg
            if not self._IsWorkingDirectory(orthofinderWorkingDir):
                orthofinderWorkingDir = userArg + "WorkingDirectory" + os.sep   
                if not self._IsWorkingDirectory(orthofinderWorkingDir):
                    print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s\nor\n   %s\n" % (userArg, orthofinderWorkingDir))
                    util.Fail()
                
        if qSpecifiedResultsFile:
            print("\nUsing orthogroups in file:\n    %s" % userArg)
            return orthofinderWorkingDir, orthofinderWorkingDir, userArg
        else:     
            # identify orthogroups file
            clustersFiles = glob.glob(orthofinderWorkingDir + "clusters_OrthoFinder_*.txt_id_pairs.txt")
            orthogroupFiles = glob.glob(orthofinderWorkingDir + "OrthologousGroups*.txt") + glob.glob(orthofinderWorkingDir + "Orthogroups*.txt")
            if orthofinderWorkingDir != userArg:
                orthogroupFiles += glob.glob(userArg + "OrthologousGroups*.txt")
                orthogroupFiles += glob.glob(userArg + "Orthogroups*.txt")
            # User may have specified a WorkingDirectory and results could be in directory above
            if len(orthogroupFiles) < len(clustersFiles):
                orthogroupFiles += glob.glob(userArg + ".." + os.sep + "OrthologousGroups*.txt")
                orthogroupFiles += glob.glob(userArg + ".." + os.sep + "Orthogroups*.txt")
            clustersFiles = sorted(clustersFiles)
            orthogroupFiles = sorted(orthogroupFiles)
            if len(clustersFiles) > 1 or len(orthogroupFiles) > 1:
                print("ERROR: Results from multiple OrthoFinder runs found\n")
                print("Tab-delimiter Orthogroups*.txt/OrthologousGroups*.txt files:")
                for fn in orthogroupFiles:
                    print("    " + fn)
                print("With corresponding cluster files:")
                for fn in clustersFiles:
                    print("    " + fn)
                print("\nPlease run with only one set of results in directories or specifiy the specific clusters_OrthoFinder_*.txt_id_pairs.txt file on the command line")
                util.Fail()        
                
            if len(clustersFiles) != 1 or len(orthogroupFiles) != 1:
                print("ERROR: Results not found in <orthofinder_results_directory> or <orthofinder_results_directory>/WorkingDirectory")
                print("\nCould not find:\n    Orthogroups*.txt/OrthologousGroups*.txt\nor\n    clusters_OrthoFinder_*.txt_id_pairs.txt")
                util.Fail()
                
            print("\nUsing orthogroups in file:\n    %s" % orthogroupFiles[0])
            print("and corresponding clusters file:\n    %s" % clustersFiles[0])
            return orthofinderWorkingDir, userArg, clustersFiles[0]
            
    def _IsWorkingDirectory(self, orthofinderWorkingDir):
        ok = True
        ok = ok and len(glob.glob(orthofinderWorkingDir + "clusters_OrthoFinder_*.txt_id_pairs.txt")) > 0
        ok = ok and len(glob.glob(orthofinderWorkingDir + "Species*.fa")) > 0
        return ok
        
    def _FindFromTrees(self, orthologuesDir, userSpeciesTree):
        """
        if userSpeciesTree == None: Use existing tree
        """
        print("\nFind from trees:")
        print((orthologuesDir, userSpeciesTree))
        self.wd_trees = orthologuesDir + "WorkingDirectory/"
        # Find species tree
        if userSpeciesTree == None:
            possibilities = ["SpeciesTree_ids_0_rooted.txt", "SpeciesTree_ids_1_rooted.txt", "SpeciesTree_user_ids.txt", "SpeciesTree_unrooted_0_rooted.txt", "STAG_SpeciesTree_ids_0_rooted.txt"] # etc (only need to determine if unique)
            nTrees = 0
            for p in possibilities:
                for d in [self.wd_trees, self.wd_trees + "Trees_ids/"]:
                    fn = d + p
                    if os.path.exists(fn): 
                        nTrees += 1
                        speciesTree_fn = fn
            if nTrees == 0:
                print("\nERROR: There is a problem with the specified directory. The rooted species tree %s or %s is not present." % (possibilities[0], possibilities[2]))
                print("Please rectify the problem or alternatively use the -s option to specify the species tree to use.\n")
                util.Fail()
            if nTrees > 1:
                print("\nERROR: There is more than one rooted species tree in the specified directory structure. Please use the -s option to specify which species tree should be used\n")
                util.Fail()
            self.speciesTreeRootedIDsFN = speciesTree_fn
        else:
            if not os.path.exists(userSpeciesTree):
                print("\nERROR: %s does not exist\n" % userSpeciesTree)
                util.Fail()
            self.speciesTreeRootedIDsFN = userSpeciesTree
            
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """

def InitialiseFileHandler(options, fastaDir=None, continuationDir=None, resultsDir_nonDefault=None, pickleDir_nonDefault=None):
    """
    Creates a file handler object which will determine the location of all the files:
    Results will be under the user specified directory of the default results location. Defaults:
        - New, from start: 
                FastaDir/OrthoFinder/Results_Date
            or
                resultsDir_nonDefault/Results_Date
        - New, continuation: Existing_OrthoFinder_Dir/Results_Date
        - Old, continuation: 
                ContinuationDir/OrthoFinder/Results_Date
            or 
                resultsDir_nonDefault/Results_Date
        
    
    Implementation
    1. Working out if an old directory structure is being used
    2. Construct and appropriate PreviousFilesLocator if necessary - this locates all required files
    3. Pass this to FileHandler - this creates the directory structure required for this run
    4. if error: print and exit
    5. Return FileHandler
    
    Tasks:
    - Switch this round, I can tell if it's and old or new directory right from the start - read log and check info present,
    perhaps just psss it to the new file handler and let it decide if everything is there
    """
    # 1 & 2
    # If starting from scratch, no need for a PreviousFileLocator
    if options.qStartFromFasta and not options.qStartFromBlast:
        pfl = None
        base_dir = resultsDir_nonDefault if resultsDir_nonDefault != None else fastaDir + "OrthoFinder/"
    else:
        try:
            # Try to process these as the new directory structure
            pfl = PreviousFilesLocator_new(options, continuationDir)
            # don't create any new directory, it already exists
            base_dir = pfl.GetHomeForResults()
        except Unprocessable:
            pfl = PreviousFilesLocator_old(options, continuationDir)
            base_dir = resultsDir_nonDefault if resultsDir_nonDefault != None else pfl.GetHomeForResults()
    if not os.path.exists(base_dir): os.mkdir(base_dir)
    # 3 
    # RefactorDS - this might be suitable as a constructor now
    # base_dir - should now exist
    
    """The previous file locator should decide where the output directory should be rooted? 
    Rules:
    - If starting from Fasta then Fasta/OrthoFinder/Results_Date
          - Or, SpecifiedDirectory/Results_Date if user specified
    - If starting from a previous new-structure directory (TopLevel/Results_X) then TopLevel/Results_Date
    - If starting from a previous old-structure directory then, as high up as we can go and still be in the directory structure:
        - Fasta/Results_OldDate/OrthoFinder/Results_Date
    """
    FileHandler.CreateOutputDirectories(options, pfl, base_dir, fastaDir)    
        