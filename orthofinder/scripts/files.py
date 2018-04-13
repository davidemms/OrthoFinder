#!/usr/bin/env python
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
"""
import os
import glob
import time
import shutil

import util


class Directories(object):
    def __init__(self):
        self.resultsDir = None           # directory for orthogroup results files
        self.workingDir = None           # Orthogroup inference workingDir
                                         # Will need to store 3 bits of information in total
    
    def IDsFilename(self):
        return self.workingDir + "SequenceIDs.txt"
    def SpeciesIdsFilename(self):
        return self.workingDir + "SpeciesIDs.txt"

class SpeciesInfo(object):
    def __init__(self):
        self.speciesToUse = []           #       seqsInfo.iSpeciesToUse   - which to include for this analysis 
        self.nSpAll = None               #       seqsInfo.nSpAll => 0, 1, ..., nSpAll - 1 are valid species indices
        self.iFirstNewSpecies = None     #       iFirstNew   => (0, 1, ..., iFirstNew-1) are from previous and (iFirstNew, iFirstNew+1, ..., nSpecies-1) are the new species indices
    def __str__(self):
        return str((self.speciesToUse, self.nSpAll, self.iFirstNewSpecies))
        

#class Files_singleton(object):
#    class __Singleton(object):
#        pass
#    instance = None
#    
#    def __init__(self):
#        if not Files_singleton.instance:
#            Files_singleton.instance = Files_singleton.__Singleton()
            
class Files_singleton(object):    
    def __init__(self):
        self.baseOgFormat = "OG%07d"
        self.wd1 = None
        self.rd1 = None
        self.ologdir = None
        self.fastaDir = None
        self.fileIdentifierString = "OrthoFinder"
        self.clustersFilename = None
        self.iResultsVersion = None
        self.resultsBaseFilename = None
        self.nondefaultPickleDir = None
        self.speciesTreeRootedIDsFN = None
        self.multipleRootedSpeciesTreesDir = None
    
    """ ========================================================================================== """
    """ Three options (currently) for initialisation:
    1. CreateOutputDirFromWorkingDirectory1   - for working from a directory of BLAST results
    2. CreateOutputDirFromInputFastaDir       - for working from a new analysis of a directory of fasta files
    3. CreateNonDefaultResultsDirectory1      - for a user specified results directory
    """
    def CreateOutputDirFromExistingDirs(self, wd1, rd1=None):
        if self.wd1 != None: raise Exception("Changing WorkingDirectory1")
        self.wd1 = wd1
        if rd1 == None: 
            self.rd1 = self.wd1
        else:
            self.rd1 = rd1

    def CreateOutputDirFromInputFastaDir(self, d, append_name = ""):
        if self.wd1 != None: raise Exception("Changing WorkingDirectory1")
        if not d.endswith("/"): d+="/"        
        self.rd1 = util.CreateNewWorkingDirectory(d + "Results_" + ("" if append_name == "" else append_name + "_"))    
        self.wd1 = self.rd1 + "WorkingDirectory/" 
        os.mkdir(self.wd1)
    
    def CreateNonDefaultResultsDirectory1(self, d):
        if self.wd1 != None: raise Exception("Changing WorkingDirectory1")
        if d[-1] != "/":
            d += "/"
        self.rd1 = d 
        self.wd1 = self.rd1 + "WorkingDirectory/"
        os.mkdir(self.rd1)
        os.mkdir(self.wd1)

    def CreateOutputDirForAnalysisFromTrees(self, orthologuesDir, userSpeciesTree):
        """
        if userSpeciesTree == None: Use existing tree
        """
        self.rd1 = orthologuesDir + "../"
        self.wd2 = orthologuesDir + "WorkingDirectory/"
        # Find species tree
        if userSpeciesTree == None:
            possibilities = ["SpeciesTree_ids_0_rooted.txt", "SpeciesTree_ids_1_rooted.txt", "SpeciesTree_user_ids.txt", "SpeciesTree_unrooted_0_rooted.txt", "STAG_SpeciesTree_ids_0_rooted.txt"] # etc (only need to determine if unique)
            nTrees = 0
            for p in possibilities:
                for d in [self.wd2, self.wd2 + "Trees_ids/"]:
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
        resultsDir_new = orthologuesDir + "New_Analysis_From_Trees"      # for the Orthologues_Species/ directories
        self.rd2 = util.CreateNewWorkingDirectory(resultsDir_new + "_")
        self.ologdir = self.rd2  + "Orthologues/"
        os.mkdir(self.ologdir)
        self.wd1, self.rd1, self.clustersFilename = util.GetOGsFile(self.rd1)
        if self.clustersFilename.endswith("_id_pairs.txt"):
            self.clustersFilename = self.clustersFilename[:-len("_id_pairs.txt")]
        
    """ ========================================================================================== """
   
    def SetNondefaultPickleDir(self, d):
        self.pickleDir = d
   
    def GetPickleDir(self):
        if self.nondefaultPickleDir != None: 
            d = self.pickleDir
        else:
            d = self.wd1 + "pickle/"
        if not os.path.exists(d): os.mkdir(d)
        return d
   
    def MakeResultsDirectory2(self, tree_generation_method, stop_after="", append_name=""):
        """
        Args
        tree_method: msa, dendroblast, phyldog (determines the directory structure that will be created)
        stop_after: seqs, align
        """
        if self.rd1 == None: raise Exception("No rd1") 
        self.rd2 = util.CreateNewWorkingDirectory(self.GetResultsDirectory1() + "Orthologues_" + ("" if append_name == "" else append_name + "_"))   
        self.wd2 = self.rd2 + "WorkingDirectory/"
        os.mkdir(self.wd2)
        os.mkdir(self.rd2 + "Orthologues/")
        if tree_generation_method == "msa":
            for i, d in enumerate([self.rd2 + "Sequences/", self.wd2 + "Sequences_ids/", self.rd2 + "Alignments/", self.wd2 + "Alignments_ids/", self.rd2 + "Gene_Trees/", self.wd2 + "Trees_ids/"]):
                if stop_after == "seqs" and i == 2: break 
                if stop_after == "align" and i == 4: break 
                if not os.path.exists(d): os.mkdir(d)
        elif tree_generation_method == "dendroblast":
            for i, d in enumerate([self.wd2 + "Distances/", self.rd2 + "Gene_Trees/", self.wd2 + "Trees_ids/"]):
                if not os.path.exists(d): os.mkdir(d)
        
    """ Standard Dirctories
        ========================================================================================== """
    
    def GetWorkingDirectory1(self):
        if self.wd1 == None: raise Exception("No wd1")
        return self.wd1 
        
    def GetWorkingDirectory2(self):
        if self.wd2 == None: raise Exception("No wd2")
        return self.wd2 
        
    def GetResultsDirectory1(self):
        if self.rd1 == None: raise Exception("No rd1")
        return self.rd1 
        
    def GetResultsDirectory2(self):
        if self.rd2 == None: raise Exception("No rd2")
        return self.rd2 
        
    def GetOrthologuesDirectory(self):
        """"Where the directories of species orthologues are"""
        if self.rd2 == None: raise Exception("No rd2")
        return self.rd2 + "Orthologues/"
        
    """ Orthogroups files 
        ========================================================================================== """
        
    def GetSpeciesIDsFN(self):
        if self.wd1 == None: raise Exception("No wd1")
        return self.wd1 + "SpeciesIDs.txt"
        
    def GetSequenceIDsFN(self):
        if self.wd1 == None: raise Exception("No wd1")
        return self.wd1 + "SequenceIDs.txt"
        
    def GetSpeciesSeqsDir(self):
        if self.wd1 == None: raise Exception("No wd1")
        return self.wd1 
        
    def GetSpeciesFastaFN(self, iSpecies):
        if self.wd1 == None: raise Exception("No wd1")
        return "%sSpecies%d.fa" % (self.wd1, iSpecies)
        
    def GetSortedSpeciesFastaFiles(self):
        if self.wd1 == None: raise Exception("No wd1")
        fastaFilenames = glob.glob(self.wd1 + "Species*.fa")
        speciesIndices = []
        for f in fastaFilenames:
            start = f.rfind("Species")
            speciesIndices.append(int(f[start+7:-3]))
        indices, sortedFasta = util.SortArrayPairByFirst(speciesIndices, fastaFilenames)
        return sortedFasta  
        
    def GetSpeciesDatabaseN(self, iSpecies, program="Blast"):
        if self.wd1 == None: raise Exception("No wd1")
        return "%s%sDBSpecies%d" % (self.wd1, program, iSpecies)
        
    def GetBlastResultsDir(self):
        return self.wd1
        
    def GetBlastResultsFN(self, iSpeciesSearch, jSpeciesDB):
        if self.wd1 == None: raise Exception("No wd1")
        return "%sBlast%d_%d.txt" % (self.wd1, iSpeciesSearch, jSpeciesDB)
        
    def GetGraphFilename(self):
        if self.wd1 == None: raise Exception("No wd1")
        return self.wd1 + "%s_graph.txt" % self.fileIdentifierString
        
    def CreateUnusedClustersFN(self, mclInflation):
        if self.wd1 == None: raise Exception("No wd1")
        self.clustersFilename, self.iResultsVersion = util.GetUnusedFilename(self.wd1  + "clusters_%s_I%0.1f" % (self.fileIdentifierString, mclInflation), ".txt")
        return self.clustersFilename, self.clustersFilename + "_id_pairs.txt"
        
    def SetClustersFN(self, pairsFN):
        self.clustersFilename = pairsFN[:-len("_id_pairs.txt")]
        
    def GetClustersFN(self):
        return self.clustersFilename + "_id_pairs.txt"
        
    def GetResultsFNBase(self):
        if self.wd1 == None: 
            raise Exception("No wd1")
        if self.iResultsVersion == None:
            raise Exception("Base results identifier has not been created")
        return self.rd1 + "Orthogroups" + ("" if self.iResultsVersion == 0 else "_%d" % self.iResultsVersion)
        
    """ Orthologues files
        ========================================================================================== """
        
    def GetResultsSeqsDir(self):
        return self.rd2 + "Sequences/"
        
    def GetResultsAlignDir(self):
        return self.rd2 + "Alignments/"
        
    def GetResultsTreesDir(self):
        return self.rd2 + "Gene_Trees/"
        
    def GetOlogStatsDir(self):
        return self.rd2
        
    def GetPutativeXenelogsDir(self):
        d = self.GetOrthologuesDirectory() + "Putative_Xenologues/"
        if not os.path.exists(d): os.mkdir(d)
        return d
    
    def GetOGsSeqFN(self, iOG, qResults=False):
        if qResults:
            return self.rd2 + "Sequences/" + (self.baseOgFormat % iOG) + ".fa"
        else:
            return self.wd2 + "Sequences_ids/" + (self.baseOgFormat % iOG) + ".fa"
            
    def GetOGsAlignFN(self, iOG, qResults=False):
        if qResults:
            return self.rd2 + "Alignments/" + (self.baseOgFormat % iOG) + ".fa"
        else:
            return self.wd2 + "Alignments_ids/" + (self.baseOgFormat % iOG) + ".fa"
            
    def GetOGsTreeFN(self, iOG, qResults=False):
        if qResults:
            return self.rd2 + "Gene_Trees/" + (self.baseOgFormat % iOG) + "_tree.txt"
        else:
            return self.wd2 + "Trees_ids/" + (self.baseOgFormat % iOG) + "_tree_id.txt"   
        
    def GetSpeciesTreeConcatAlignFN(self, qResults=False):
        if qResults:
            return self.rd2 + "Alignments/SpeciesTreeAlignment.fa"
        else:
            return self.wd2 + "Alignments_ids/SpeciesTreeAlignment.fa"  
        
    def GetSpeciesTreeMatrixFN(self, qPutInWorkingDir = False):
        if qPutInWorkingDir:
            return self.wd2 + "SpeciesMatrix.phy"
        else:
            return self.wd2 + "Distances/SpeciesMatrix.phy"
            
    def GetSpeciesTreeUnrootedFN(self, qAccessions=False):
        if qAccessions:
            return self.wd2 + "SpeciesTree_unrooted.txt"
        else: 
            return self.wd2 + "Trees_ids/SpeciesTree_unrooted_id.txt"  
            
    def SetSpeciesTreeRootedFN(self, fn):
        self.speciesTreeRootedIDsFN = fn
            
    def GetSpeciesTreeRootedFN(self):
        return self.speciesTreeRootedIDsFN
        
    def GetSpeciesTreeResultsFN(self, i, qUnique):
        """
        The results species tree (rooted, accessions, support values)
        i: index for species tree, starting at 0
        qUnique: bool, has a unique root been identified (as it may not be known exatly which branch the root belongs on)
        E.g. if there were just one species tree, the correct call would be GetSpeciesTreeResultsFN(0,True)
        """
        if qUnique:
            return self.rd2 + "SpeciesTree_rooted.txt"
        else:
            if not self.multipleRootedSpeciesTreesDir:
                self.multipleRootedSpeciesTreesDir = self.rd2 + "Potential_Rooted_Species_Trees/"
                if not os.path.exists(self.multipleRootedSpeciesTreesDir): os.mkdir(self.multipleRootedSpeciesTreesDir)
            return self.multipleRootedSpeciesTreesDir + "SpeciesTree_rooted_at_outgroup_%d.txt" % i
        
    def GetSpeciesTreeUserSupplied_idsFN(self):
        return self.wd2 + "SpeciesTree_UserSupplied_Rooted_IDs.txt"
        
    def GetOGsDistMatFN(self, iOG):
        return self.wd2 + "Distances/OG%07d.phy" % iOG
        
    def GetSpeciesDict(self):
        d = util.FullAccession(self.GetSpeciesIDsFN()).GetIDToNameDict()
        return {k:v.rsplit(".",1)[0] for k,v in d.items()}
        
    """ ========================================================================================== """
            
    def GetOGsTreeDir(self, qResults=False):
        if qResults:
            return self.rd2 + "Gene_Trees/" 
        else:
            return self.wd2 + "Trees_ids/" 
            
    def GetOGsReconTreeDir(self, qResults=False):
        if qResults:
            d = self.rd2 + "Recon_Gene_Trees/" 
            if not os.path.exists(d): os.mkdir(d)
            return d
        else:
            raise NotImplemented() 
            
    def GetPhyldogWorkingDirectory(self):
        d = self.wd2 + "phyldog/"
        if not os.path.exists(d): os.mkdir(d)
        return d
            
    def GetPhyldogOGResultsTreeFN(self, i):
        return self.wd2 + "phyldog/OG%07d.ReconciledTree.txt" % i
            
    def GetDuplicationsFN(self):
        return self.rd2 + "Duplications.csv"
        
    
    """ ========================================================================================== """
         
    def CleanWorkingDir2(self):
        dirs = ['Distances/']
        for d in dirs:
            dFull = self.wd2 + d
            if os.path.exists(dFull): 
                try:
                    shutil.rmtree(dFull)
                except OSError:
                    time.sleep(1)
                    shutil.rmtree(dFull, True)  # shutil / NFS bug - ignore errors, it's less crucial that the files are deleted
                    
FileHandler = Files_singleton()
        
    

        