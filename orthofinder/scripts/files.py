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

import util


class Directories(object):
    def __init__(self):
        self.resultsDir = None           # directory for orthogroup results files
        self.workingDir = None           # Orthogroup inference workingDir
        self.separatePickleDir = None
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
        self.wd1 = None
        self.rd1 = None
        self.fastaDir = None
        self.fileIdentifierString = "OrthoFinder"
        self.clustersFilename = None
        self.iResultsVersion = None
        self.resultsBaseFilename = None
    
    """ Two options (currently) for creating:
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

    def CreateOutputDirFromInputFastaDir(self, d, name = ""):
        if self.wd1 != None: raise Exception("Changing WorkingDirectory1")
        if not d.endswith("/"): d+="/"        
        self.rd1 = util.CreateNewWorkingDirectory(d + "Results_" + ("" if name == "" else name + "_"))    
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

    
    """ ========================================================================================== """
        
    def GetWorkingDirectory1(self):
        if self.wd1 == None: raise Exception("No wd1")
        return self.wd1 
        
    def GetResultsDirectory1(self):
        if self.rd1 == None: raise Exception("No rd1")
        return self.rd1 
        
    def GetSpeciesIDsFN(self):
        if self.wd1 == None: raise Exception("No wd1")
        return self.wd1 + "SpeciesIDs.txt"
        
    def GetSequenceIDsFN(self):
        if self.wd1 == None: raise Exception("No wd1")
        return self.wd1 + "SequenceIDs.txt"
        
    def GetSortedSpeciesFastaFiles(self):
        if self.wd1 == None: raise Exception("No wd1")
        fastaFilenames = glob.glob(self.wd1 + "Species*.fa")
        speciesIndices = []
        for f in fastaFilenames:
            start = f.rfind("Species")
            speciesIndices.append(int(f[start+7:-3]))
        indices, sortedFasta = util.SortArrayPairByFirst(speciesIndices, fastaFilenames)
        return sortedFasta  
        
    def GetSpeciesFastaFN(self, iSpecies):
        if self.wd1 == None: raise Exception("No wd1")
        return "%sSpecies%d.fa" % (self.wd1, iSpecies)
        
    def GetSpeciesDatabaseN(self, iSpecies, program="Blast"):
        if self.wd1 == None: raise Exception("No wd1")
        return "%s%sDBSpecies%d" % (self.wd1, program, iSpecies)
        
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
        
    def GetResultsFNBase(self):
        if self.wd1 == None: 
            raise Exception("No wd1")
        if self.iResultsVersion == None:
            raise Exception("Base results identifier has not been created")
        return self.rd1 + "Orthogroups" + ("" if self.iResultsVersion == 0 else "_%d" % self.iResultsVersion)
             
FileHandler = Files_singleton()
        
    

        