#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014 David Emms
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
# david_emms@hotmail.com

import sys                                      # Y
import subprocess                               # Y
import os                                       # Y
import glob                                     # Y
import multiprocessing                          # optional  (problems on OpenBSD)
import itertools                                # Y
import datetime                                 # Y
from scipy.optimize import curve_fit            # install
import numpy as np                              # install
import csv                                      # Y
import scipy.sparse as sparse                   # install
import os.path                                  # Y
import numpy.core.numeric as numeric            # install
import cPickle as pic                           # Y
import time                                     # Y
from collections import defaultdict             # Y
import xml.etree.ElementTree as ET              # Y
from xml.etree.ElementTree import SubElement    # Y
from xml.dom import minidom                     # Y

version = "0.4.0"
fastaExtensions = {"fa", "faa", "fasta", "fas"}
if sys.platform.startswith("linux"):
    with open(os.devnull, "w") as f:
        subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f) # get round problem with python multiprocessing library that can set all cpu affinities to a single cpu

"""
Utilities
-------------------------------------------------------------------------------
"""
def RunCommand(command):
    subprocess.call(command)
        
def RunCommandReport(command):
    util.PrintTime("Running command: %s" % " ".join(command))
    RunCommand(command)
    util.PrintTime("Finished command: %s" % " ".join(command))
   
class util:
    @staticmethod
    def GetDirectoryName(baseDirName, dateString, i):
        if i == 0:
            return baseDirName + dateString + os.sep
        else:
            return baseDirName + dateString + ("_%d" % i) + os.sep
    
    """Call GetNameForNewWorkingDirectory before a call to CreateNewWorkingDirectory to find out what directory will be created"""
    @staticmethod
    def CreateNewWorkingDirectory(baseDirectoryName):
        dateStr = datetime.date.today().strftime("%b%d") 
        iAppend = 0
        newDirectoryName = util.GetDirectoryName(baseDirectoryName, dateStr, iAppend)
        while os.path.exists(newDirectoryName):
            iAppend += 1
            newDirectoryName = util.GetDirectoryName(baseDirectoryName, dateStr, iAppend)
        os.mkdir(newDirectoryName)
        return newDirectoryName
    
    @staticmethod
    def GetUnusedFilename(baseFilename, ext):
        iAppend = 0
        newFilename = baseFilename + ext
        while os.path.exists(newFilename):
            iAppend += 1
            newFilename = baseFilename + ("_%d" % iAppend) + ext
        return newFilename, iAppend
    
    @staticmethod
    def PrintTime(message):
        print(str(datetime.datetime.now()) + " : " + message)  
           
    @staticmethod
    def SortArrayPairByFirst(useForSortAr, keepAlignedAr, qLargestFirst=False):
        sortedTuples = sorted(zip(useForSortAr, keepAlignedAr), reverse=qLargestFirst)
        useForSortAr = [i for i, j in sortedTuples]
        keepAlignedAr = [j for i, j in sortedTuples]
        return useForSortAr, keepAlignedAr
           
    @staticmethod
    def PrintNoNewLine(text):
        sys.stdout.write(text)
     
    @staticmethod
    def SortFastaFilenames(fastaFilenames):
        speciesIndices = []
        for f in fastaFilenames:
            start = f.rfind("Species")
            speciesIndices.append(int(f[start+7:-3]))
        indices, sortedFasta = util.SortArrayPairByFirst(speciesIndices, fastaFilenames)
        return sortedFasta    
    
"""
IDExtractor
-------------------------------------------------------------------------------
"""
class IDExtractor(object):
    """IDExtractor deals with the fact that for different datasets a user will
    want to extract a unique sequence ID from the fasta file accessions uin different 
    ways."""
    def GetIDToNameDict(self):
        raise NotImplementedError("Should not be implemented")
    def GetNameToIDDict(self):
        raise NotImplementedError("Should not be implemented")

class FullAccession(IDExtractor):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        self.nameToIDDict = dict()
        with open(idsFilename, 'rb') as idsFile:
            for line in idsFile:
                id, accession = line.rstrip().split(": ", 1)
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession                
                self.nameToIDDict[accession] = id 
                
    def GetIDToNameDict(self):
        return self.idToNameDict
        
    def GetNameToIDDict(self):
        return self.nameToIDDict
                
class FirstWordExtractor(IDExtractor):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        self.nameToIDDict = dict()
        with open(idsFilename, 'rb') as idsFile:
            for line in idsFile:
                id, rest = line.split(": ", 1)
                accession = rest.split(None, 1)[0]
                if accession in self.nameToIDDict:
                    raise RuntimeError("A duplicate accession was found using just first part: % s" % accession)
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession                
                self.nameToIDDict[accession] = id   
                
    def GetIDToNameDict(self):
        return self.idToNameDict
        
    def GetNameToIDDict(self):
        return self.nameToIDDict

"""
MCL
-------------------------------------------------------------------------------
"""

class MCL:
    @staticmethod
    def GetPredictedOGs(clustersFilename):
        predictedOGs = []
        nOGsString = ""
        qContainsProfiles = False
        with open(clustersFilename, 'rb') as clusterFile:
            header = True
            og = set()
            for line in clusterFile:
                if header:
                    if line.count("begin"):
                        header = False
                else:
                    if line.find(")") != -1:
                        break
                    if line[-2] == "$":
                        line = line[:-3]
                    if line[0] == " ":
                        # continuation of group
                        x = line.split()
                        y = [x_ for x_ in x if not x_.startswith('Prof')]
                        og = og.union(y)
                    else:
                        # new OG
                        if len(og) != 0:
                            predictedOGs.append(og)
                        nOGsString, line = line.split(" ", 1)
                        x = line.split()
                        y = [x_ for x_ in x if not x_.startswith('Prof')]
                        if len(x) != len(y):
                            qContainsProfiles = True
                        og = set(y)
            if len(og) > 0:
                predictedOGs.append(og)
        if not qContainsProfiles:
            assert(len(predictedOGs) == int(nOGsString) + 1)
        return predictedOGs
        
    @staticmethod
    def GetSingleID(speciesStartingIndices, seq, speciesToUse):    
        iSpecies, iSeq = map(int, seq.split("_"))
        offset = speciesStartingIndices[speciesToUse.index(iSpecies)]
        return iSeq + offset
        
    @staticmethod
    def GetIDPair(speciesStartingIndices, singleID, speciesToUse):   
        for i, startingIndex in enumerate(speciesStartingIndices):
            if startingIndex > singleID:
                return "%d_%d" % (speciesToUse[i-1], singleID - speciesStartingIndices[i-1])
        return "%d_%d" % (speciesToUse[-1], singleID - speciesStartingIndices[len(speciesStartingIndices)-1]) 
    
    @staticmethod
    def ConvertSingleIDsToIDPair(speciesStartingIndices, clustersFilename, newFilename, speciesToUse):
        with open(clustersFilename, 'rb') as clusterFile, open(newFilename, "wb") as output:
            header = True
            for line in clusterFile:
                appendDollar = False
                initialText = ""
                idsString = ""
                ids = []
                if header:
                    output.write(line)
                    if line.count("begin"):
                        header = False
                else:
                    if line.find(")") != -1:
                        output.write(line)
                        break
                    if line[-2] == "$":
                        line = line[:-3]
                        appendDollar = True
                    if line[0] != " ":
                        initialText, line = line.split(None, 1)
                    # continuation of group
                    ids = line.split()
                    for id in ids:
                        idsString += MCL.GetIDPair(speciesStartingIndices, int(id), speciesToUse) + " "
                    output.write(initialText + "      " + idsString)
                    if appendDollar:
                        output.write("$\n")
                    else:
                        output.write("\n")
                        
    @staticmethod
    def CreateOGs(predictedOGs, outputFilename, idDict):
        with open(outputFilename, 'wb') as outputFile:
            for iOg, og in enumerate(predictedOGs):
                outputFile.write("OG%07d:" % iOg)
                for seq in og:
                    outputFile.write(" " + idDict[seq])
                outputFile.write("\n")
      
    @staticmethod            
    def prettify(elem):
        """Return a pretty-printed XML string for the Element.
        """
        rough_string = ET.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")
    
    @staticmethod            
    def WriteOrthoXML(speciesInfo, predictedOGs, numbersOfSequences, idDict, orthoxmlFilename, speciesToUse):
        """ speciesInfo: ordered array for which each element has
            fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
        """
                
        
        # Write OrthoXML file
        root = ET.Element("orthoXML")
        root.set('xsi:schemaLocation', "http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd")
        root.set('originVersion', version)
        root.set('origin', 'OrthoFinder')
        root.set('version', "0.3")
        root.set('xmlns:xsi', "http://www.w3.org/2001/XMLSchema-instance")
        #notes = SubElement(root, 'notes')

        # Species: details of source of genomes and sequences they contain
        speciesStartingIndices = []
        iGene_all = 0
        for iSpecies, (species, nSeqs, thisSpeciesInfo) in enumerate(zip(speciesInfo, numbersOfSequences, speciesInfo)):
            speciesNode = SubElement(root, 'species')
            speciesNode.set('NCBITaxId', thisSpeciesInfo[2])           # required
            speciesNode.set('name', thisSpeciesInfo[1])                # required
            speciesDatabaseNode = SubElement(speciesNode, "database")
            speciesDatabaseNode.set('name', thisSpeciesInfo[3])            # required
            speciesDatabaseNode.set('version', thisSpeciesInfo[4])         # required
#            speciesDatabaseNode.set('geneLink', "")        # skip
#            speciesDatabaseNode.set('protLink', "")        # skip
#            speciesDatabaseNode.set('transcriptLink', "")  # skip
            allGenesNode = SubElement(speciesDatabaseNode, "genes")
            speciesStartingIndices.append(iGene_all)
            for iGene_species in xrange(nSeqs):
                geneNode = SubElement(allGenesNode, 'gene')
                geneNode.set("geneId", idDict["%d_%d" % (iSpecies , iGene_species)])  
                geneNode.set('id', str(iGene_all))       # required
#                geneNode.set("protID", "")  # skip
                iGene_all += 1
                
        # Scores tag - unused
#            scoresNode = SubElement(root, 'scores')        # skip

        # Orthogroups
        allGroupsNode = SubElement(root, 'groups')
        for iOg, og in enumerate(predictedOGs):
            groupNode = SubElement(allGroupsNode, 'orthologGroup')
            groupNode.set('id', str(iOg))
#                groupScoreNode = SubElement(groupNode, 'score')    # skip
#                groupScoreNode.set('id', "")                       # skip
#                groupScoreNode.set('value', "")                    # skip
#                SubElement(groupNode, 'property')                  # skip
            for seq in og:
                geneNode = SubElement(groupNode, 'geneRef')
                geneNode.set('id', str(MCL.GetSingleID(speciesStartingIndices, seq, speciesToUse)))
#                    SubElement(geneNode, 'score')                  # skip
        with open(orthoxmlFilename, 'wb') as orthoxmlFile:
#            ET.ElementTree(root).write(orthoxmlFile)
            orthoxmlFile.write(MCL.prettify(root))
        print("Orthologous groups have been written to orthoxml file:\n   %s" % orthoxmlFilename)
                        
    @staticmethod                       
    def RunMCL(graphFilename, clustersFilename, inflation = 1.5):
        command = ["mcl", graphFilename, "-I", "1.5", "-o", clustersFilename]
        RunCommand(command)
        util.PrintTime("Ran MCL")  

"""
scnorm
-------------------------------------------------------------------------------
"""
class scnorm:
    @staticmethod
    def loglinear(x, a, b):
        return a*np.log10(x)+b     
    
    @staticmethod
    def GetLengthArraysForMatrix(m, len_i, len_j):
        I, J = m.nonzero()
        scores = [v for row in m.data for v in row]     # use fact that it's lil
        Li = np.array(len_i[I])
        Lj = np.array(len_j[J])
        return Li, Lj, scores
        
    @staticmethod
    def GetTopPercentileOfScores(L, S, percentileToKeep):
        # Get the top x% of hits at each length
        nScores = len(S)
        t_sort = sorted(zip(L, range(nScores)))
        indices = [j for i, j in t_sort]
        s_sorted = [S[i] for i in indices]
        l_sorted = [L[i] for i in indices]
        if nScores < 100:
            # then we can't split them into bins, return all for fitting
            return l_sorted, s_sorted
        nInBins = 1000 if nScores > 5000 else (200 if nScores > 1000 else 20)
        nBins, remainder = divmod(nScores, nInBins)
        topScores = []
        topLengths = []
        for i in xrange(nBins):
            first = i*nInBins
            last = min((i+1)*nInBins-1, nScores - 1)
            theseLengths = l_sorted[first:last+1]
            theseScores = s_sorted[first:last+1]
            cutOff = np.percentile(theseScores, percentileToKeep)
            lengthsToKeep = [thisL for thisL, thisScore in zip(theseLengths, theseScores) if thisScore >= cutOff]
            topLengths.extend(lengthsToKeep)
            topScores.extend([thisScore for thisL, thisScore in zip(theseLengths, theseScores) if thisScore >= cutOff])
        return topLengths, topScores
        
    @staticmethod
    def CalculateFittingParameters(Lf, S):
        pars,covar =  curve_fit(scnorm.loglinear, Lf, np.log10(S))
        return pars
           
    @staticmethod   
    def NormaliseScoresByLogLengthProduct(b, Lq, Lh, params): 
        rangeq = range(len(Lq))
        rangeh = range(len(Lh))
        li_vals = Lq**(-params[0])
        lj_vals = Lh**(-params[0])
        li_matrix = sparse.csr_matrix((li_vals, (rangeq, rangeq)))
        lj_matrix = sparse.csr_matrix((lj_vals, (rangeh, rangeh)))
        return sparse.lil_matrix(10**(-params[1]) * li_matrix * b * lj_matrix)

 
"""
BlastFileProcessor
-------------------------------------------------------------------------------
"""   
class BlastFileProcessor(object): 
    def __init__(self, filesDirectory, speciesToUse, nSeqs, nSpecies, speciesStartingIndices):     
        self.filesDirectory = filesDirectory
        self.speciesToUse = speciesToUse
        self.nSeqs = nSeqs
        self.nSpecies = nSpecies
        self.speciesStartingIndices = speciesStartingIndices
        self.sep = "_"
        self.tol = 1e-3
        
    def GetBH_s(self, pairwiseScoresMatrices, iSpecies):
        nSeqs_i = self.NumberOfSequences(iSpecies)
        bestHitForSequence = -1*np.ones(nSeqs_i)
        H = [None for i_ in xrange(self.nSpecies)] # create array of Nones to be replace by matrices
        for j in xrange(self.nSpecies):
            if iSpecies == j:
                # identify orthologs then come back to paralogs
                continue
            W = pairwiseScoresMatrices[j]
            I = []
            J = []
            for kRow in xrange(nSeqs_i):
                values=W.getrowview(kRow)
                if values.nnz == 0:
                    continue
                m = max(values.data[0])
                bestHitForSequence[kRow] = m if m > bestHitForSequence[kRow] else bestHitForSequence[kRow]
                # get all above this value with tolerance
                temp = [index for index, value in zip(values.rows[0], values.data[0]) if value > m - self.tol]
                J.extend(temp)
                I.extend(kRow * np.ones(len(temp), dtype=np.dtype(int)))
            H[j] = sparse.csr_matrix((np.ones(len(I)), (I, J)), shape=W.get_shape())
        # now look for paralogs
        I = []
        J = []
        W = pairwiseScoresMatrices[iSpecies]
        for kRow in xrange(nSeqs_i):
            values=W.getrowview(kRow)
            if values.nnz == 0:
                continue
            temp = [index for index, value in zip(values.rows[0], values.data[0]) if value > bestHitForSequence[kRow] - self.tol]
            J.extend(temp)
            I.extend(kRow * np.ones(len(temp), dtype=np.dtype(int)))
        H[iSpecies] = sparse.csr_matrix((np.ones(len(I)), (I, J)), shape=W.get_shape())
        return H
       
    def MatrixAnd(self, H):       
        for i in xrange(self.nSpecies):
            for j in xrange(i + 1):
                H[i][j] = H[i][j].multiply(H[j][i].transpose())
                if i != j:
                    H[j][i] = H[i][j].transpose() 
        return H
       
    @staticmethod             
    def MatricesAnd_s(Xarr, Yarr):
        Zarr = []
        for x, y in zip(Xarr, Yarr):
            Zarr.append(x.multiply(y))
        return Zarr
        
    @staticmethod             
    def MatricesAndTr_s(Xarr, Yarr):
        Zarr = []
        for x, y in zip(Xarr, Yarr):
            Zarr.append(x.multiply(y.transpose()))
        return Zarr
    
    @staticmethod        
    def GetNumberOfSequencesInFile(filename):
        lastIDLine = ""
        sequenceStartingIndices = []
        currentSpecies = 0
        sequenceStartingIndices.append(0)
        with open(filename) as file:
            count = 0
            for line in file:
                if len(line) > 1 and line[0] == ">":
                    count+=1
                    lastIDLine = line
                    thisSpecies = int(line[1:].split("_", 1)[0])
                    if thisSpecies != currentSpecies:
                        sequenceStartingIndices.append(count-1)
                        currentSpecies = thisSpecies 
        nSpecies = int(lastIDLine[1:].split("_")[0]) + 1
        return count, nSpecies, sequenceStartingIndices
        
    def GetIDPairFromString(self, line):
        return map(int, line.split(self.sep))
        
    def NumberOfSequences(self, iSpecies):
        return (self.speciesStartingIndices[iSpecies+1] if iSpecies != self.nSpecies-1 else self.nSeqs) - self.speciesStartingIndices[iSpecies] 
  
    @staticmethod   
    def GetNumberOfSequencesInFileFromDir(inputDirectory, speciesToUse):
        sequenceStartingIndices = [0]
        count = 0
        for i, iFasta in enumerate(speciesToUse):
            fastaFilename = inputDirectory + "Species%d.fa" % iFasta
            with open(fastaFilename) as file:
                for line in file:
                    if len(line) > 1 and line[0] == ">":
                        count+=1
            sequenceStartingIndices.append(count)
        sequenceStartingIndices = sequenceStartingIndices[:-1]
        nSpecies = len(speciesToUse)
        return count, nSpecies, sequenceStartingIndices
        
    def GetSequenceLengths(self):                
        sequenceLengths = []
        for iSpecies, iFasta in enumerate(self.speciesToUse):
            sequenceLengths.append(np.zeros(self.NumberOfSequences(iSpecies)))
            fastaFilename = self.filesDirectory + "Species%d.fa" % iFasta
            currentSequenceLength = 0
            iCurrentSequence = -1
            qFirstLine = True
            with open(fastaFilename) as file:
                for row in file:
                    if len(row) > 1 and row[0] == ">":    
                        if qFirstLine:
                            qFirstLine = False
                        else:
                            sequenceLengths[iSpecies][iCurrentSequence] = currentSequenceLength
                            currentSequenceLength = 0
                        _, iCurrentSequence = self.GetIDPairFromString(row[1:])
                    else:
                        currentSequenceLength += len(row.rstrip())
            sequenceLengths[iSpecies][iCurrentSequence] = currentSequenceLength
        return sequenceLengths
                   
    def GetBLAST6Scores(self, iSpecies, jSpecies): 
        nSeqs_i = self.NumberOfSequences(iSpecies)
        nSeqs_j = self.NumberOfSequences(jSpecies)
        B = sparse.lil_matrix((nSeqs_i, nSeqs_j))
        row = ""
        try:
            with open(self.filesDirectory + "Blast%d_%d.txt" % (self.speciesToUse[iSpecies], self.speciesToUse[jSpecies]), 'rb') as blastfile:
                blastreader = csv.reader(blastfile, delimiter='\t')
                for row in blastreader:    
                    # Get hit and query IDs
                    try:
                        species1ID, sequence1ID = map(int, row[0].split(self.sep, 1)) 
                        species2ID, sequence2ID = map(int, row[1].split(self.sep, 1))     
                    except (IndexError, ValueError):
                        sys.stderr.write("\nERROR: Query or hit sequence ID in BLAST results file was missing or incorrectly formatted.\n")
                        raise
                    # Get bit score for pair
                    try:
                        score = float(row[11])   
                    except (IndexError, ValueError):
                        sys.stderr.write("\nERROR: 12th field in BLAST results file line should be the bit-score for the hit\n")
                        raise
                    qSameSequence = (species1ID == species2ID and sequence1ID == sequence2ID)
                    if qSameSequence:
                        continue
                    # store bit score
                    try:
                        if score > B[sequence1ID, sequence2ID]: 
                            B[sequence1ID, sequence2ID] = score   
                    except IndexError:
                        def ord(n):
                            return str(n)+("th" if 4<=n%100<=20 else {1:"st",2:"nd",3:"rd"}.get(n%10, "th"))
#                        sys.stderr.write("\nError in input files, expected only %d sequences in species %d and %d sequences in species %d but found a hit in the Blast%d_%d.txt between sequence %d_%d (i.e. %s sequence in species) and sequence %d_%d (i.e. %s sequence in species)\n" %  (nSeqs_i, iSpecies, nSeqs_j, jSpecies, iSpecies, jSpecies, iSpecies, sequence1ID, ord(sequence1ID+1), jSpecies, sequence2ID, ord(sequence2ID+1)))
                        sys.stderr.write("\nERROR: Inconsistent input files.\n")
                        kSpecies, nSeqs_k, sequencekID = (iSpecies,  nSeqs_i, sequence1ID) if sequence1ID >= nSeqs_i else (jSpecies,  nSeqs_j, sequence2ID)
                        sys.stderr.write("Species%d.fa contains only %d sequences " % (kSpecies,  nSeqs_k)) 
                        sys.stderr.write("but found a query/hit in the Blast%d_%d.txt for sequence %d_%d (i.e. %s sequence in species %d).\n" %  (iSpecies, jSpecies, kSpecies, sequencekID, ord(sequencekID+1), kSpecies))
                        Fail()
        except Exception:
            sys.stderr.write("Malformatted line in %sBlast%d_%d.txt\nOffending line was:\n" % (self.filesDirectory, iSpecies, jSpecies))
            sys.stderr.write("\t".join(row) + "\n")
            Fail()
        return B       

"""
WaterfallMethod
-------------------------------------------------------------------------------
"""   
class WaterfallMethod:
    def __init__(self, inputDirectory, outputDirectory, speciesToUse, nSeqs, nSpecies, speciesStartingIndices):
        self.thisBfp = BlastFileProcessor(inputDirectory, speciesToUse, nSeqs, nSpecies, speciesStartingIndices)    
        self.outputDir = outputDirectory 
        self.speciesToUse = speciesToUse
        if not os.path.exists(self.outputDir):
           os.mkdir(self.outputDir)  
        self.picProtocol = 1
        self.totalDump = 0.
        self.totalLoad = 0.
        
    def RunWaterfallMethod(self, graphFilename):
        util.PrintTime("Started")   
        Lengths = self.thisBfp.GetSequenceLengths()
        util.PrintTime("Got sequence lengths")
        util.PrintTime("Initial processing of each species")
        # process up to the best hits for each species
        for iSpecies in xrange(self.thisBfp.nSpecies):
            Bi = []
            for jSpecies in xrange(self.thisBfp.nSpecies):
                Bij = self.thisBfp.GetBLAST6Scores(iSpecies, jSpecies)  
                Bij = self.NormaliseScores(Bij, Lengths, iSpecies, jSpecies)
                Bi.append(Bij)
            self.DumpMatrixArray("B", Bi, iSpecies)
            BH = self.thisBfp.GetBH_s(Bi, iSpecies)
            self.DumpMatrixArray("BH", BH, iSpecies)
            util.PrintTime("Initial processing of species %d complete" % iSpecies)
        
        for iSpecies in xrange(self.thisBfp.nSpecies):
            # calculate RBH for species i
            BHix = self.LoadMatrixArray("BH", iSpecies)
            BHxi = self.LoadMatrixArray("BH", iSpecies, row=False)
            RBHi = self.thisBfp.MatricesAndTr_s(BHix, BHxi)   # twice as much work as before (only did upper triangular before)
            B = self.LoadMatrixArray("B", iSpecies)
            connect = self.ConnectAllBetterThanAnOrtholog_s(RBHi, B, iSpecies) 
            self.DumpMatrixArray("connect", connect, iSpecies) 
        util.PrintTime("Connected putatitive homologs") 
        
        with open(graphFilename, 'wb') as graphFile:
            graphFile.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n" % (self.thisBfp.nSeqs, self.thisBfp.nSeqs)) 
            graphFile.write("\n(mclmatrix\nbegin\n\n")  
            for iSpec in xrange(self.thisBfp.nSpecies):
                # calculate the 2-way connections for one query species
                connect2 = []
                for jSpec in xrange(self.thisBfp.nSpecies):
                    m1 = self.LoadMatrix("connect", iSpec, jSpec)
                    m2tr = numeric.transpose(self.LoadMatrix("connect", jSpec, iSpec))
                    connect2.append(m1 + m2tr)
                B = self.LoadMatrixArray("B", iSpec)
                B_connect = self.thisBfp.MatricesAnd_s(connect2, B)
                util.PrintTime("Writen final scores for species %d to graph file" % iSpec)
                
                W = [b.sorted_indices().tolil() for b in B_connect]
                for query in xrange(self.thisBfp.NumberOfSequences(iSpec)):
                    offset = self.thisBfp.speciesStartingIndices[iSpec]
                    graphFile.write("%d    " % (offset + query))
                    for jSpec in xrange(self.thisBfp.nSpecies):
                        row = W[jSpec].getrowview(query)
                        jOffset = self.thisBfp.speciesStartingIndices[jSpec]
                        for j, value in zip(row.rows[0], row.data[0]):
                            graphFile.write("%d:%.3f " % (j + jOffset, value))
                    graphFile.write("$\n")
                        
            graphFile.write(")\n")

        # delete pic files
        self.DeleteMatrices()
        
    def DeleteMatrices(self):
        for f in glob.glob(self.outputDir + "B*_*.pic"):
            os.remove(f)
        for f in glob.glob(self.outputDir + "connect*_*.pic"):
            os.remove(f)

    def DumpMatrix(self, name, m, iSpecies, jSpecies):
        start = time.time()
        with open(self.outputDir + "%s%d_%d.pic" % (name, iSpecies, jSpecies), 'wb') as picFile:
            pic.dump(m, picFile, protocol=self.picProtocol)
        self.totalDump += (time.time() - start)
        
    def DumpMatrixArray(self, name, matrixArray, iSpecies):
        for jSpecies, m in enumerate(matrixArray):
            self.DumpMatrix(name, m, iSpecies, jSpecies)

    def LoadMatrix(self, name, iSpecies, jSpecies):   
        start = time.time()
        with open(self.outputDir + "%s%d_%d.pic" % (name, iSpecies, jSpecies), 'rb') as picFile:  
            M = pic.load(picFile)
        self.totalLoad += (time.time() - start)
        return M
            
    def LoadMatrixArray(self, name, iSpecies, row=True):
        matrixArray = []
        for jSpecies in xrange(self.thisBfp.nSpecies):
            if row == True:
                matrixArray.append(self.LoadMatrix(name, iSpecies, jSpecies))
            else:
                matrixArray.append(self.LoadMatrix(name, jSpecies, iSpecies))
        return matrixArray
        
    def NormaliseScores(self, B, Lengths, iSpecies, jSpecies):              
        Li, Lj, scores = scnorm.GetLengthArraysForMatrix(B, Lengths[iSpecies], Lengths[jSpecies])
        Lf = Li * Lj     
        topLf, topScores = scnorm.GetTopPercentileOfScores(Lf, scores, 95)   
        if len(topScores) > 1:
            fittingParameters = scnorm.CalculateFittingParameters(topLf, topScores)  
            return scnorm.NormaliseScoresByLogLengthProduct(B, Lengths[iSpecies], Lengths[jSpecies], fittingParameters)
        else:
            print("WARNING: Too few hits between species %d and species %d to normalise the scores, these hits will be ignored" % (iSpecies, jSpecies))
            return sparse.lil_matrix(B.get_shape())

    def GetMostDistant_s(self, RBH, B, iSpec):
        mostDistant = numeric.transpose(np.ones(self.thisBfp.NumberOfSequences(iSpec))*1e9)
        for kSpec in xrange(self.thisBfp.nSpecies):
            B[kSpec] = B[kSpec].tocsr()
            if iSpec == kSpec:
                continue
            I, J = RBH[kSpec].nonzero()
            if len(I) > 0:
                mostDistant[I] = np.minimum(B[kSpec][I, J], mostDistant[I])
        return mostDistant
    
    def ConnectAllBetterThanCutoff_s(self, B, mostDistant, iSpec):
        connect = []
        nSeqs_i = self.thisBfp.NumberOfSequences(iSpec)
        for jSpec in xrange(self.thisBfp.nSpecies):
            M=B[jSpec].tolil()
            if iSpec != jSpec:
                IIJJ = [(i,j) for i, (valueRow, indexRow) in enumerate(zip(M.data, M.rows)) for j, v in zip(indexRow, valueRow) if v >= mostDistant[i]]
            else:
                IIJJ = [(i,j) for i, (valueRow, indexRow) in enumerate(zip(M.data, M.rows)) for j, v in zip(indexRow, valueRow) if (i != j) and v >= mostDistant[i]]
            II = [i for (i, j) in IIJJ]
            JJ = [j for (i, j) in IIJJ]
            onesArray = np.ones(len(IIJJ))
            mat = sparse.csr_matrix( (onesArray,  (II, JJ)), shape=(nSeqs_i,  self.thisBfp.NumberOfSequences(jSpec)))
            connect.append(mat)
        return connect
      
    def ConnectAllBetterThanAnOrtholog_s(self, RBH, B, iSpec):        
        mostDistant = self.GetMostDistant_s(RBH, B, iSpec) 
        connect = self.ConnectAllBetterThanCutoff_s(B, mostDistant, iSpec)
        return connect

"""
OrthoFinder
-------------------------------------------------------------------------------
"""   
nBlastDefault = 16

def CanRunCommand(command, qAllowStderr = False):
    util.PrintNoNewLine("Test can run \"%s\"" % command)       # print without newline
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]
    if len(stdout) > 0 and (qAllowStderr or len(stderr) == 0):
        print(" - ok")
        return True
    else:
        print(" - failed")
        return False

def CanRunBLAST():
    if CanRunCommand("makeblastdb -help") and CanRunCommand("blastp -help"):
        return True
    else:
        print("ERROR: Cannot run BLAST+")
        print("Please check BLAST+ is installed and that the executables are in the system path\n")
        return False

def CanRunMCL():
    command = "mcl -h"
    if CanRunCommand(command):
        return True
    else:
        print("ERROR: Cannot run MCL with the command \"%s\"" % command)
        print("Please check MCL is installed and in the system path\n")
        return False
        
def PrintCitation():
    print("""\nWhen publishing work that uses OrthoFinder please cite:
    D.M. Emms & S. Kelly (2015), OrthoFinder: solving fundamental biases in whole genome comparisons
    dramatically improves orthogroup inference accuracy, Genome Biology 16:157.\n""")   
    
def PrintHelp():  
    print("Simple Usage") 
    print("------------")
    print("python orthofinder.py -f fasta_directory [-t max_number_of_threads]")
    print("    Infers orthologous groups for the proteomes contained in fasta_directory running")
    print("    max_number_of_threads in parallel for the BLAST searches.")
    print("")    
    print("Advanced Usage")
    print("--------------")
    print("python orthofinder.py -f fasta_directory -p")
    print("    1. Prepares files for BLAST and prints the BLAST commands. Does not perform BLAST searches")
    print("    or infer othologous groups. Useful if you want to prepare the files in the form required by")
    print("    OrthoFinder but want to perform the BLAST searches using a job scheduler/on a cluster and")
    print("    then infer orthologous groups using option 2.")
    print("")
    print("python orthofinder.py -b precalculated_blast_results_directory ")
    print("    2. Infers orthologous groups using pre-calculated BLAST results. These can be after BLAST")
    print("    searches have been completed following the use of option 1 or using the WorkingDirectory")
    print("    from a previous OrthoFinder run. Species can be commented out with a '#' in the SpeciesIDs.txt")
    print("    file to exclude them from the analysis. See README file for details.")
    print("")
    print("python orthofinder.py -b precalculated_blast_results_directory -f fasta_directory")
    print("    3. Add species from fasta_directory to a previous OrthoFinder run where precalculated_blast_results_directory")
    print("    is the directory containing the BLAST results files etc. from the previous run.")
    print("")
#    print("python orthofinder.py -s precalculated_blast_results_directory")
#    print("    4. Use a subset of species from a previous OrthoFinder run. The precalculated_blast_results_directory")
#    print("    should contain the BLAST results files etc. plus an extra, user-supplied file 'SpeciesSubset.txt'")
#    print("    specifying which species should be included. See README file for details.")
#    print("")
    
#    print("Arguments:\n")
#    print("""-f fasta_directory, --fasta fasta_directory
#   Predict orthogroups for the genes in the fasta files in the fasta_directory\n""")
#    
#    print("""-b precalculated_blast_results_directory, --blast precalculated_blast_results_directory
#   Predict orthogroups using the pre-calcualted BLAST results in precalculated_blast_results_directory.
#   The directory must contain the BLAST results files, fasta files with IDs for the accessions, 
#   SequenceIDs.txt and SpeciesIDs.txt in the formats described in the README file.\n""")
#    
#    print("""-t max_number_of_threads, --threads max_number_of_threads
#   The maximum number of BLAST processes to be run simultaneously. The deafult is %d but this 
#   should be increased by the user to at least the number of cores on the computer so as to 
#   minimise the time taken to perform the BLAST all-versus-all queries.\n""" % nBlastDefault)
#    
#    print("""-x speciesInfoFilename, --orthoxml speciesInfoFilename
#   Output the orthogroups in the orthoxml format using the information in speciesInfoFilename.\n""")
        
    print("Arguments")
    print("---------")
    print("""-f fasta_directory, --fasta fasta_directory
    Predict orthogroups for the proteins in the fasta files in the fasta_directory\n""")
    
    print("""-b precalculated_blast_results_directory, --blast precalculated_blast_results_directory
    Predict orthogroups using the pre-calcualted BLAST results in precalculated_blast_results_directory.\n""")
    
    print("""-t max_number_of_threads, --threads max_number_of_threads
    The maximum number of BLAST processes to be run simultaneously. The deafult is %d but this 
    should be increased by the user to at least the number of cores on the computer so as to 
    minimise the time taken to perform the BLAST all-versus-all queries.\n""" % nBlastDefault)
    
    print("""-x speciesInfoFilename, --orthoxml speciesInfoFilename
    Output the orthogroups in the orthoxml format using the information in speciesInfoFilename.\n""")
    
    print("""-p , --prepare
    Only prepare the files in the format required by OrthoFinder and print out the BLAST searches that
    need to be performed but don't run BLAST or infer orthologous groups\n""" )
    
    print("""-h, --help
    Print this help text""")
    PrintCitation()

def GetDirectoryArgument(arg, args):
    if len(args) == 0:
        print("Missing option for command line argument %s" % arg)
        Fail()
    directory = os.path.abspath(args.pop(0))
    if directory[-1] != os.sep:
        directory += os.sep
    return directory

def Fail():
    sys.exit()
    
def AssignIDsToSequences(fastaDirectory, outputDirectory):
    idsFilename = outputDirectory + "SequenceIDs.txt"
    speciesFilename = outputDirectory + "SpeciesIDs.txt"
    iSeq = 0
    iSpecies = 0
    # check if SpeciesIDs.txt already exists
    if os.path.exists(speciesFilename):
        with open(speciesFilename, 'rb') as infile:
            for line in infile: pass
        if line.startswith("#"): line = line[1:]
        iSpecies = int(line.split(":")[0]) + 1
    originalFastaFilenames = sorted([f for f in os.listdir(fastaDirectory) if os.path.isfile(os.path.join(fastaDirectory,f))])
    originalFastaFilenames = [f for f in originalFastaFilenames if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in fastaExtensions]
    returnFilenames = []
    newSpeciesIDs = []
    with open(idsFilename, 'ab') as idsFile, open(speciesFilename, 'ab') as speciesFile:
        for fastaFilename in originalFastaFilenames:
            newSpeciesIDs.append(iSpecies)
            outputFastaFilename = outputDirectory + "Species%d.fa" % iSpecies
            outputFasta = open(outputFastaFilename, 'wb')
            returnFilenames.append(outputFastaFilename)            
            fastaFilename = fastaFilename.rstrip()
            speciesFile.write("%d: %s\n" % (iSpecies, fastaFilename))
            baseFilename, extension = os.path.splitext(fastaFilename)
            with open(fastaDirectory + os.sep + fastaFilename, 'rb') as fastaFile:
                for line in fastaFile:
                    if len(line) > 0 and line[0] == ">":
                        newID = "%d_%d" % (iSpecies, iSeq)
                        idsFile.write("%s: %s" % (newID, line[1:]))
                        outputFasta.write(">%s\n" % newID)    
                        iSeq += 1
                    else:
                        outputFasta.write(line)
                outputFasta.write("\n")
            iSpecies += 1
            iSeq = 0
            outputFasta.close()
    if len(originalFastaFilenames) > 0: outputFasta.close()
    return returnFilenames, originalFastaFilenames, idsFilename, speciesFilename, newSpeciesIDs

def AnalyseSequences(inputDir, outputDir, speciesToUse, nSeqs, nSpecies, speciesStartingIndices, graphFilename):
    wfAlg = WaterfallMethod(inputDir, outputDir, speciesToUse, nSeqs, nSpecies, speciesStartingIndices)
    wfAlg.RunWaterfallMethod(graphFilename)

def WriteOrthogroupFiles(ogs, idsFilenames, resultsBaseFilename, clustersFilename_pairs):
    outputFN = resultsBaseFilename + ".txt"
    try:
        fullDict = dict()
        for idsFilename in idsFilenames:
            idExtract = FirstWordExtractor(idsFilename)
            idDict = idExtract.GetIDToNameDict()
            fullDict.update(idDict)
        MCL.CreateOGs(ogs, outputFN, fullDict)
    except KeyError as e:
        sys.stderr.write("ERROR: Sequence ID not found in %s\n" % idsFilename)
        sys.stderr.write(str(e) + "\n")
        Fail()        
    except RuntimeError as error:
        print(error.message)
        if error.message.startswith("ERROR"):
            print("ERROR: %s contains a duplicate ID. The IDs for the orthologous groups in %s will not be replaced with the sequence accessions. If %s was prepared manually then please check the IDs are correct. " % (idsFilename, clustersFilename_pairs, idsFilename))
            Fail()
        else:
            print("Tried to use only the first part of the accession in order to list the sequences in each orthologous group more concisely but these were not unique. Will use the full accession line instead.")     
            try:
                fullDict = dict()
                for idsFilename in idsFilenames:
                    idExtract = FullAccession(idsFilename)
                    idDict = idExtract.GetIDToNameDict()
                    fullDict.update(idDict)
                MCL.CreateOGs(ogs, outputFN, fullDict)   
            except:
                print("ERROR: %s contains a duplicate ID. The IDs for the orthologous groups in %s will not be replaced with the sequence accessions. If %s was prepared manually then please check the IDs are correct. " % (idsFilename, clustersFilename_pairs, idsFilename))
                Fail()
    return fullDict

def GetSpeciesToUse(speciesIDsFN):
    speciesToUse = []
    with open(speciesIDsFN, 'rb') as speciesF:
        for line in speciesF:
            if len(line) == 0 or line[0] == "#": continue
            speciesToUse.append(int(line.split(":")[0]))
    return speciesToUse

def CreateOrthogroupTable(ogs, 
                          idToNameDict, 
                          speciesFilename, 
                          speciesToUse,
                          resultsBaseFilename):
    
    speciesNamesDict = dict()
    with open(speciesFilename, 'rb') as speciesNamesFile:
        for line in speciesNamesFile:
            if line.startswith("#"): continue
            short, full = line.rstrip().split(": ")
            speciesNamesDict[int(short)] = full    
    nSpecies = len(speciesNamesDict) 
    
    ogs_names = [[idToNameDict[seq] for seq in og] for og in ogs]
    ogs_ints = [[map(int, sequence.split("_")) for sequence in og] for og in ogs]

    # write out
    outputFilename = resultsBaseFilename + ".csv"
    singleGeneFilename = resultsBaseFilename + "_UnassignedGenes.csv"
    with open(outputFilename, 'wb') as outputFile, open(singleGeneFilename, 'wb') as singleGeneFile:
        fileWriter = csv.writer(outputFile, delimiter="\t")
        singleGeneWriter = csv.writer(singleGeneFile, delimiter="\t")
        for writer in [fileWriter, singleGeneWriter]:
            row = [""] + [speciesNamesDict[index] for index in speciesToUse]
            writer.writerow(row)
        
        for iOg, (og, og_names) in enumerate(zip(ogs_ints, ogs_names)):
            ogDict = defaultdict(list)
            rows = ["OG%07d" % iOg]
            thisOutputWriter = fileWriter
            # separate it into sequences from each species
            if len(og) == 1:
                rows.extend(['' for x in xrange(nSpecies)])
                rows[speciesToUse.index(og[0][0]) + 1] = og_names[0]
                thisOutputWriter = singleGeneWriter
            else:
                for (iSpecies, iSequence), name in zip(og, og_names):
                    ogDict[speciesToUse.index(iSpecies)].append(name)
                for iSpecies in xrange(nSpecies):
                    rows.append(", ".join(ogDict[iSpecies]))
            thisOutputWriter.writerow(rows)
    print("""Orthologous groups have been written to tab-delimited files:\n   %s\n   %s""" % (outputFilename, singleGeneFilename))
    print("""And in OrthoMCL format:\n   %s""" % (outputFilename[:-3] + "txt"))
            
#def GetOrderedBlastCommands(fastaFilenames, blastDBs, workingDirectory):
def GetOrderedBlastCommands(previousFastaFiles, newFastaFiles, workingDir):
    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands 
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing 
    the BLAST commands.
    """
    iSpeciesPrevious = [int(fn[fn.rfind("Species") + 7:].split(".")[0]) for fn in previousFastaFiles]
    iSpeciesNew = [int(fn[fn.rfind("Species") + 7:].split(".")[0]) for fn in newFastaFiles]
    nSeqs = {i:BlastFileProcessor.GetNumberOfSequencesInFile(workingDir + "Species%d.fa" % i)[0] for i in (iSpeciesPrevious+iSpeciesNew)}
    speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew)] + \
                   [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesPrevious)] + \
                   [(i, j) for i, j in itertools.product(iSpeciesPrevious, iSpeciesNew)] 
    taskSizes = [nSeqs[i]*nSeqs[j] for i,j in speciesPairs]
    taskSizes, speciesPairs = util.SortArrayPairByFirst(taskSizes, speciesPairs, True)
    commands = [["blastp", "-outfmt", "6", "-evalue", "0.001", "-query", workingDir + "Species%d.fa" % iFasta, "-db", workingDir + "BlastDBSpecies%d" % iDB, "-out", "%sBlast%d_%d.txt" % (workingDir, iFasta, iDB)]
                    for iFasta, iDB in speciesPairs]               

    return commands
"""
Main
-------------------------------------------------------------------------------
"""   

if __name__ == "__main__":
    print("\nOrthoFinder version %s Copyright (C) 2014 David Emms\n" % version)
    print("""    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software, and you are welcome to redistribute it under certain conditions.
    For details please see the License.md that came with this software.\n""")
    if len(sys.argv) == 1 or sys.argv[1] == "--help" or sys.argv[1] == "help" or sys.argv[1] == "-h":
        PrintHelp()
        sys.exit()
             
    # Control
    nBlast = nBlastDefault
    qUsePrecalculatedBlast = False  # remove, just store BLAST to do
    qUseFastaFiles = False  # local to argument checking
    qXML = False
    qOnlyPrepare = False
#    qUseSubset = False # forget once arguments parsed and have decided what the BLAST files are (Effectively part of BLASTFileProcessor)
                     
    # Files         
    fastaDir = None             # Where to find the original fasta files        # replace with list of files
    resultsDir = None           # Where to put the redults files
    workingDir_previous = None  # where to find original Blast, processed fasta etc
    workingDir = ""             # where to put all new files 
    previousFastaFiles = []
    newFastaFiles = []
    speciesToUse = []
    
    args = sys.argv[1:]
    while len(args) > 0:
        arg = args.pop(0)    
        if arg == "-f" or arg == "--fasta":
            if qUseFastaFiles:
                print("Repeated argument: -f/--fasta")
                Fail()
            qUseFastaFiles = True
            fastaDir = GetDirectoryArgument(arg, args)
        elif arg == "-b" or arg == "--blast":
            if qUsePrecalculatedBlast:
                print("Repeated argument: -b/--blast")
                Fail()
            qUsePrecalculatedBlast = True
            workingDir_previous = GetDirectoryArgument(arg, args)
        elif arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument -t")
                Fail()
            arg = args.pop(0)
            try:
                nBlast = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s" % arg)
                Fail()    
        elif arg == "-x" or arg == "--orthoxml":  
            if qXML:
                print("Repeated argument: -x/--orthoxml")
                Fail()
            qXML = True
            if len(args) == 0:
                print("Missing option for command line argument %s" % arg)
                Fail()
            speciesInfoFilename = args.pop(0)
#        elif arg == "-s" or arg == "--subset":  
#            if qUseSubset:
#                print("Repeated argument: -s/--subset")
#                Fail()
#            qUseSubset = True
#            qUsePrecalculatedBlast = True
#            workingDir_previous = GetDirectoryArgument(arg, args)
        elif arg == "-p" or arg == "--prepare":
            qOnlyPrepare = True
        elif arg == "-h" or arg == "--help":
            PrintHelp()
            sys.exit()
        else:
            print("Unrecognised argument: %s\n" % arg)
            Fail()            
    
    # check argument combinations   
    if qUseFastaFiles and qUsePrecalculatedBlast:
        print("Adding new species in %s to existing analysis in %s" % (fastaDir, workingDir_previous))
    
    # if using previous results, check everything is ok
    if qUsePrecalculatedBlast:
        speciesIdsFilename = workingDir_previous + "SpeciesIDs.txt"
        if not os.path.exists(speciesIdsFilename):
            print("%s file must be provided if using previously calculated BLAST results" % speciesIdsFilename)
            Fail()
        speciesToUse = GetSpeciesToUse(speciesIdsFilename)
        workingDir = os.path.abspath(workingDir) + os.sep
        if resultsDir == None: 
            workingDir = workingDir_previous
            resultsDir = workingDir_previous
        else:
            workingDir = os.path.abspath(workingDir) + os.sep + "WorkingDirectory" + os.sep
        
        # check BLAST results directory exists
        if not os.path.exists(workingDir_previous):
            print("Previous/Pre-calculated BLAST results directory does not exist: %s\n" % workingDir_previous)
            Fail()
     
        # check fasta files are present 
        previousFastaFiles = util.SortFastaFilenames(glob.glob(workingDir_previous + "Species*.fa"))
        if len(previousFastaFiles) == 0:
            print("No processed fasta files in the supplied previous working directory: %s\n" % workingDir_previous)
            Fail()
        tokens = previousFastaFiles[-1][:-3].split("Species")
        lastFastaNumberString = tokens[-1]
        iLastFasta = 0
        nFasta = len(previousFastaFiles)
        try:
            iLastFasta = int(lastFastaNumberString)
        except:
            print("Filenames for processed fasta files are incorrect: %s\n" % previousFastaFiles[-1])
            Fail()
        if nFasta != iLastFasta + 1:
            print("Not all expected fasta files are present. Index of last fasta file is %s but found %d fasta files.\n" % (lastFastaNumberString, len(previousFastaFiles)))
            Fail()
        
        # check BLAST files
        qHaveBlast = True
        for iSpecies in speciesToUse:
            for jSpecies in speciesToUse:
                filename = "%sBlast%d_%d.txt" % (workingDir_previous, iSpecies, jSpecies) 
                if not os.path.exists(filename):
                    print("BLAST results file is missing: %s" % filename)
                    qHaveBlast = False
        if not qHaveBlast: Fail()
                    
        # check SequenceIDs.txt and SpeciesIDs.txt files are present
        idsFilename = workingDir_previous + "SequenceIDs.txt"
        if not os.path.exists(idsFilename):
            print("%s file must be provided if using previous calculated BLAST results" % idsFilename)
            Fail()
        print("Using previously calculated BLAST results in %s" % workingDir_previous)       
    else:
        # - create working directory
        if resultsDir == None: resultsDir = util.CreateNewWorkingDirectory(fastaDir + "Results_")
        workingDir = resultsDir + "WorkingDirectory" + os.sep
        os.mkdir(workingDir)
     
     
    # check for BLAST+ and MCL - else instruct how to install and add to path
    print("\n1. Checking required programs are installed")
    print("-------------------------------------------")
    if (not qUsePrecalculatedBlast) and (not CanRunBLAST()):
        Fail()
    if not CanRunMCL():
        Fail()
        
    # - rename sequences with unique, simple identifiers
    print("\n2. Temporarily renaming sequences with unique, simple identifiers")
    print( "------------------------------------------------------------------")
    if not qUseFastaFiles:
        print("Skipping")
    else:
        newFastaFiles, userFastaFilenames, idsFilename, speciesIdsFilename, newSpeciesIDs = AssignIDsToSequences(fastaDir, workingDir)
        speciesToUse = speciesToUse + newSpeciesIDs
        print("Done")
     
    
    if qXML:   
        print("\n2b. Reading species information file")
        print( "-------------------------------------")        
        # do this now so that we can alert user to any errors prior to running the algorithm
        # speciesInfo:  name, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
        userFastaFilenames = previousFastaFiles + userFastaFilenames
        speciesInfo = [[] for i_ in speciesToUse]
        userFastaFilenames_justNames = [name for path, name in map(os.path.split, userFastaFilenames)]
        fastaFileIndices = {filename:iSpecies for iSpecies, filename in enumerate(userFastaFilenames_justNames)}
        with open(speciesInfoFilename, 'rb') as speciesInfoFile:
            reader = csv.reader(speciesInfoFile, delimiter = "\t")
            for iLine, line in enumerate(reader):
                if len(line) != 5:
                    # allow for an extra empty line at the end
                    if len(line) == 0 and iLine == len(userFastaFilenames_justNames):
                        continue
                    print("ERROR")
                    print("Species information file %s line %d is incorrectly formatted." % (speciesInfoFilename, iLine + 1))        
                    print("File should be contain one line per species")
                    print("Each line should contain 5 tab-delimited fields:")
                    print("  fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseFastaFilename")
                    print("See README file for more information.")
                    Fail() 
                fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile = line
                try:
                    iSpecies = fastaFileIndices[fastaFilename]
                except KeyError:
                    print("ERROR")
                    print("%s from line %d of the species information file was not one of the" % (fastaFilename, iLine+1))
                    print("input fasta files. The input fasta files were:")
                    for filename in userFastaFilenames_justNames:
                        print(filename)
                    print("Please provide information for each of these species in the species information file")
                    Fail() 
                speciesInfo[iSpecies] = line   
        # check information has been provided for all species
        speciesMissing = False        
        for fastaFilename, iSpecies in fastaFileIndices.items():
            if speciesInfo[iSpecies] == []:
                if not speciesMissing:
                    print("ERROR")
                    print("Species information file %s does not contain information for all species." % speciesInfoFilename) 
                    print("Information is missing for:") 
                    speciesMissing = True
                print(fastaFilename)
        if speciesMissing:
            Fail()
     
    print("\n3. Dividing up work for BLAST for parallel processing")
    print(  "-----------------------------------------------------")
    if not qUseFastaFiles:
        print("Skipping")
    else:
        print("\n3a. Creating BLAST databases")
        print(  "----------------------------")
        for iSp in speciesToUse:
            command = ["makeblastdb", "-dbtype", "prot", "-in", workingDir + "Species%d.fa" % iSp, "-out", workingDir + "BlastDBSpecies%d" % iSp]
            RunCommand(command)    
    
    if qOnlyPrepare:
        print("\n4. BLAST commands that must be run")
        print(  "----------------------------------")
    else:        
        print("\n4. Running BLAST all-versus-all")
        print(  "-------------------------------")
    if not qUseFastaFiles:
        print("Skipping")
    else:
        if qUsePrecalculatedBlast:
            print("Only running new BLAST searches")
        commands = GetOrderedBlastCommands(previousFastaFiles, newFastaFiles, workingDir)
        if qOnlyPrepare:
            for command in commands:
                print(" ".join(command))
            sys.exit()
        print("Maximum number of BLAST processer: %d" % nBlast)
        util.PrintTime("This may take some time....")  
        pool = multiprocessing.Pool(nBlast) 
        pool.map(RunCommandReport, commands)   
        print("Done!")  
        # remove BLAST databases
        for f in glob.glob(workingDir + "BlastDBSpecies*"):
            os.remove(f)


    # Run Algorithm, cluster and output cluster files with original accessions
    print("\n5. Running OrthoFinder algorithm")
    print(  "--------------------------------")
    fileIdentifierString = "OrthoFinder_v%s" % version
    graphFilename = workingDir + "%s_graph.txt" % fileIdentifierString
    nSeqs, nSpecies, speciesStartingIndices = BlastFileProcessor.GetNumberOfSequencesInFileFromDir(workingDir_previous if qUsePrecalculatedBlast else workingDir, speciesToUse)
    # it's important to free up the memory from python used for processing the genomes
    # before launching MCL becuase both use sizeable ammounts of memory. The only
    # way I can find to do this is to launch the memory intensive python code 
    # as separate process that exitsbefore MCL is launched.
#    AnalyseSequences(workingDir, nSeqs, nSpecies, speciesStartingIndices, graphFilename)  # without launching a new process
    p = multiprocessing.Process(target=AnalyseSequences, args=(workingDir_previous if qUsePrecalculatedBlast else workingDir, workingDir, speciesToUse, nSeqs, nSpecies, speciesStartingIndices, graphFilename))
    p.start()
    p.join()
    if p.exitcode != 0:
        Fail()
    
    
    
    # 5b. MCL     
    inflation = 1.5
    clustersFilename, iResultsVersion = util.GetUnusedFilename(workingDir + "clusters_%s_I%0.1f" % (fileIdentifierString, inflation), ".txt")
    MCL.RunMCL(graphFilename, clustersFilename, inflation)
    clustersFilename_pairs = clustersFilename + "_id_pairs.txt"
    MCL.ConvertSingleIDsToIDPair(speciesStartingIndices, clustersFilename, clustersFilename_pairs, speciesToUse)  
    # delete single ID filename
#        os.remove(clustersFilename)   
    
    print("\n6. Creating files for Orthologous Groups")
    print(  "----------------------------------------")
    PrintCitation()
    ogs = MCL.GetPredictedOGs(clustersFilename_pairs)
    resultsBaseFilename = util.GetUnusedFilename(resultsDir + "OrthologousGroups", ".csv")[:-4]         # remove .csv from base filename
    resultsBaseFilename = resultsDir + "OrthologousGroups" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion)
    idsDict = WriteOrthogroupFiles(ogs, [idsFilename], resultsBaseFilename, clustersFilename_pairs)
    CreateOrthogroupTable(ogs, idsDict, speciesIdsFilename, speciesToUse, resultsBaseFilename)
    if qXML:
        numbersOfSequences = list(np.diff(speciesStartingIndices))
        numbersOfSequences.append(nSeqs - speciesStartingIndices[-1])
        orthoxmlFilename = resultsBaseFilename + ".orthoxml"
        MCL.WriteOrthoXML(speciesInfo, ogs, numbersOfSequences, idsDict, orthoxmlFilename, speciesToUse)
    print("\n")
