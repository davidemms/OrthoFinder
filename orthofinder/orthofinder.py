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
import multiprocessing as mp                    # optional  (problems on OpenBSD)
import itertools                                # Y
import datetime                                 # Y
from collections import Counter                 # Y
from scipy.optimize import curve_fit            # install
import numpy as np                              # install
import csv                                      # Y
import scipy.sparse as sparse                   # install
import os.path                                  # Y
import numpy.core.numeric as numeric            # install
from collections import defaultdict             # Y
import xml.etree.ElementTree as ET              # Y
from xml.etree.ElementTree import SubElement    # Y
from xml.dom import minidom                     # Y
import Queue                                    # Y
import warnings                                 # Y

import scripts.mcl as MCLread
import scripts.blast_file_processor as BlastFileProcessor
from scripts import util, matrices, get_orthologues

fastaExtensions = {"fa", "faa", "fasta", "fas"}
if sys.platform.startswith("linux"):
    with open(os.devnull, "w") as f:
        subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f) # get round problem with python multiprocessing library that can set all cpu affinities to a single cpu

# Fix LD_LIBRARY_PATH when using pyinstaller 
my_env = os.environ.copy()
if getattr(sys, 'frozen', False):
    if 'LD_LIBRARY_PATH_ORIG' in my_env:
        my_env['LD_LIBRARY_PATH'] = my_env['LD_LIBRARY_PATH_ORIG']  
    else:
        my_env['LD_LIBRARY_PATH'] = ''  
    if 'DYLD_LIBRARY_PATH_ORIG' in my_env:
        my_env['DYLD_LIBRARY_PATH'] = my_env['DYLD_LIBRARY_PATH_ORIG']  
    else:
        my_env['DYLD_LIBRARY_PATH'] = ''      
         
def RunBlastDBCommand(command):
    capture = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]
    nLines_success= 10
    if len(stdout) > nLines_success or len(stderr) > 0:
        print("\nWarning:")
        print("".join(stdout[2:]))
        if len(stderr) > 0: print(stderr)
            
def SpeciesNameDict(speciesIDsFN):
    speciesNamesDict = dict()
    with open(speciesIDsFN, 'rb') as speciesNamesFile:
        for line in speciesNamesFile:
            if line.startswith("#"): continue
            short, full = line.rstrip().split(": ")
            speciesNamesDict[int(short)] = full.split(".", 1)[0]
    return speciesNamesDict
    
"""
MCL
-------------------------------------------------------------------------------
"""    
class MCL:
    @staticmethod
    def CreateOGs(predictedOGs, outputFilename, idDict):
        with open(outputFilename, 'wb') as outputFile:
            for iOg, og in enumerate(predictedOGs):
                outputFile.write("OG%07d: " % iOg)
                accessions = sorted([idDict[seq] for seq in og])
                outputFile.write(" ".join(accessions))
                outputFile.write("\n")
      
    @staticmethod
    def prettify(elem):
        """Return a pretty-printed XML string for the Element.
        """
        rough_string = ET.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")
       
    @staticmethod         
    def WriteOrthoXML(speciesInfo, predictedOGs, nSequencesDict, idDict, orthoxmlFilename, speciesToUse):
        """ speciesInfo: ordered array for which each element has
            fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
        """
        # Write OrthoXML file
        root = ET.Element("orthoXML")
        root.set('xsi:schemaLocation', "http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd")
        root.set('originVersion', util.version)
        root.set('origin', 'OrthoFinder')
        root.set('version', "0.3")
        root.set('xmlns:xsi', "http://www.w3.org/2001/XMLSchema-instance")
        #notes = SubElement(root, 'notes')
    
        # Species: details of source of genomes and sequences they contain
        speciesStartingIndices = []
        iGene_all = 0
        for iPos, thisSpeciesInfo in enumerate(speciesInfo):
            iSpecies = speciesToUse[iPos]
            nSeqs = nSequencesDict[iSpecies]
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
                geneNode.set('id', str(MCLread.GetSingleID(speciesStartingIndices, seq, speciesToUse)))
    #                    SubElement(geneNode, 'score')                  # skip
        with open(orthoxmlFilename, 'wb') as orthoxmlFile:
    #            ET.ElementTree(root).write(orthoxmlFile)
            orthoxmlFile.write(MCL.prettify(root))
        print("Orthogroups have been written to orthoxml file:\n   %s" % orthoxmlFilename)
            
    @staticmethod               
    def RunMCL(graphFilename, clustersFilename, nProcesses, inflation):
        nProcesses = 4 if nProcesses > 4 else nProcesses    # MCL appears to take *longer* as more than 4 processes are used
        command = ["mcl", graphFilename, "-I", str(inflation), "-o", clustersFilename, "-te", str(nProcesses), "-V", "all"]
        util.RunCommand(command)
        util.PrintTime("Ran MCL")  
    
    @staticmethod
    def WriteOrthogroupFiles(ogs, idsFilenames, resultsBaseFilename, clustersFilename_pairs):
        outputFN = resultsBaseFilename + ".txt"
        try:
            fullDict = dict()
            for idsFilename in idsFilenames:
                idExtract = util.FirstWordExtractor(idsFilename)
                idDict = idExtract.GetIDToNameDict()
                fullDict.update(idDict)
            MCL.CreateOGs(ogs, outputFN, fullDict)
        except KeyError as e:
            sys.stderr.write("ERROR: Sequence ID not found in %s\n" % idsFilename)
            sys.stderr.write(str(e) + "\n")
            util.Fail()        
        except RuntimeError as error:
            print(error.message)
            if error.message.startswith("ERROR"):
                print("ERROR: %s contains a duplicate ID. The IDs for the orthogroups in %s will not be replaced with the sequence accessions. If %s was prepared manually then please check the IDs are correct. " % (idsFilename, clustersFilename_pairs, idsFilename))
                util.Fail()
            else:
                print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup\nmore concisely but these were not unique. The full accession line will be used instead.\n")     
                try:
                    fullDict = dict()
                    for idsFilename in idsFilenames:
                        idExtract = util.FullAccession(idsFilename)
                        idDict = idExtract.GetIDToNameDict()
                        fullDict.update(idDict)
                    MCL.CreateOGs(ogs, outputFN, fullDict)   
                except:
                    print("ERROR: %s contains a duplicate ID. The IDs for the orthogroups in %s will not be replaced with the sequence accessions. If %s was prepared manually then please check the IDs are correct. " % (idsFilename, clustersFilename_pairs, idsFilename))
                    util.Fail()
        return fullDict
    
    @staticmethod
    def CreateOrthogroupTable(ogs, 
                              idToNameDict, 
                              speciesNamesDict, 
                              speciesToUse,
                              resultsBaseFilename):
        
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
                row = ["OG%07d" % iOg]
                thisOutputWriter = fileWriter
                # separate it into sequences from each species
                if len(og) == 1:
                    row.extend(['' for x in xrange(nSpecies)])
                    row[speciesToUse.index(og[0][0]) + 1] = og_names[0]
                    thisOutputWriter = singleGeneWriter
                else:
                    for (iSpecies, iSequence), name in zip(og, og_names):
                        ogDict[speciesToUse.index(iSpecies)].append(name)
                    for iSpecies in xrange(nSpecies):
                        row.append(", ".join(sorted(ogDict[iSpecies])))
                thisOutputWriter.writerow(row)
        resultsFilesString = "Orthogroups have been written to tab-delimited files:\n   %s\n   %s (OrthoMCL format)\n   %s" % (outputFilename, outputFilename[:-3] + "txt", singleGeneFilename)
        return resultsFilesString

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
RunInfo
-------------------------------------------------------------------------------
"""     

def GetSequenceLengths(seqsInfo, fileInfo):                
    sequenceLengths = []
    for iSpecies, iFasta in enumerate(seqsInfo.speciesToUse):
        sequenceLengths.append(np.zeros(seqsInfo.nSeqsPerSpecies[iFasta]))
        fastaFilename = fileInfo.workingDir + "Species%d.fa" % iFasta
        currentSequenceLength = 0
        iCurrentSequence = -1
        qFirstLine = True
        with open(fastaFilename) as infile:
            for row in infile:
                if len(row) > 1 and row[0] == ">":    
                    if qFirstLine:
                        qFirstLine = False
                    else:
                        sequenceLengths[iSpecies][iCurrentSequence] = currentSequenceLength
                        currentSequenceLength = 0
                    _, iCurrentSequence = util.GetIDPairFromString(row[1:])
                else:
                    currentSequenceLength += len(row.rstrip())
        sequenceLengths[iSpecies][iCurrentSequence] = currentSequenceLength
    return sequenceLengths
  
# Redundant?  
def GetNumberOfSequencesInFile(filename):
    count = 0
    with open(filename) as infile:
        for line in infile:
            if line.startswith(">"): count+=1
    return count

""" Question: Do I want to do all BLASTs or just the required ones? It's got to be all BLASTs I think. They could potentially be 
run after the clustering has finished."""
def GetOrderedBlastCommands(seqsInfo, dirs):
    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands 
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing 
    the BLAST commands.
    """
    iSpeciesPrevious = range(dirs.iFirstNewSpecies)
    iSpeciesNew = range(dirs.iFirstNewSpecies, dirs.nSpAll)
    speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew)] + \
                   [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesPrevious)] + \
                   [(i, j) for i, j in itertools.product(iSpeciesPrevious, iSpeciesNew)] 
    taskSizes = [seqsInfo.nSeqsPerSpecies[i]*seqsInfo.nSeqsPerSpecies[j] for i,j in speciesPairs]
    taskSizes, speciesPairs = util.SortArrayPairByFirst(taskSizes, speciesPairs, True)
    commands = [["blastp", "-outfmt", "6", "-evalue", "0.001", "-query", dirs.workingDir + "Species%d.fa" % iFasta, "-db", dirs.workingDir + "BlastDBSpecies%d" % iDB, "-out", "%sBlast%d_%d.txt" % (dirs.workingDir, iFasta, iDB)]
                    for iFasta, iDB in speciesPairs]               
    return commands     

"""
Matrices
-------------------------------------------------------------------------------
""" 

def DeleteMatrices(fileInfo):
    for f in glob.glob(fileInfo.workingDir + "B*_*.pic"):
        if os.path.exists(f): os.remove(f)
    for f in glob.glob(fileInfo.workingDir + "connect*_*.pic"):
        if os.path.exists(f): os.remove(f)
            
def GetBH_s(pairwiseScoresMatrices, seqsInfo, iSpecies, tol=1e-3):
    nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpecies]]
    bestHitForSequence = -1*np.ones(nSeqs_i)
    H = [None for i_ in xrange(seqsInfo.nSpecies)] # create array of Nones to be replace by matrices
    for j in xrange(seqsInfo.nSpecies):
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
            temp = [index for index, value in zip(values.rows[0], values.data[0]) if value > m - tol]
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
        temp = [index for index, value in zip(values.rows[0], values.data[0]) if value > bestHitForSequence[kRow] - tol]
        J.extend(temp)
        I.extend(kRow * np.ones(len(temp), dtype=np.dtype(int)))
    H[iSpecies] = sparse.csr_matrix((np.ones(len(I)), (I, J)), shape=W.get_shape())
    return H
    
    
"""
WaterfallMethod
-------------------------------------------------------------------------------
""" 

def WriteGraph_perSpecies(args):
    seqsInfo, fileInfo, iSpec = args            
    # calculate the 2-way connections for one query species
    with open(fileInfo.graphFilename + "_%d" % iSpec, 'wb') as graphFile:
        connect2 = []
        for jSpec in xrange(seqsInfo.nSpecies):
            m1 = matrices.LoadMatrix("connect", fileInfo, iSpec, jSpec)
            m2tr = numeric.transpose(matrices.LoadMatrix("connect", fileInfo, jSpec, iSpec))
            connect2.append(m1 + m2tr)
        B = matrices.LoadMatrixArray("B", fileInfo, seqsInfo, iSpec)
        B_connect = matrices.MatricesAnd_s(connect2, B)
        
        W = [b.sorted_indices().tolil() for b in B_connect]
        for query in xrange(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]):
            offset = seqsInfo.seqStartingIndices[iSpec]
            graphFile.write("%d    " % (offset + query))
            for jSpec in xrange(seqsInfo.nSpecies):
                row = W[jSpec].getrowview(query)
                jOffset = seqsInfo.seqStartingIndices[jSpec]
                for j, value in zip(row.rows[0], row.data[0]):
                    graphFile.write("%d:%.3f " % (j + jOffset, value))
            graphFile.write("$\n")
        if iSpec == (seqsInfo.nSpecies - 1): graphFile.write(")\n")
        util.PrintTime("Writen final scores for species %d to graph file" % iSpec)
            
            
class WaterfallMethod:    
    @staticmethod
    def NormaliseScores(B, Lengths, iSpecies, jSpecies):    
        Li, Lj, scores = scnorm.GetLengthArraysForMatrix(B, Lengths[iSpecies], Lengths[jSpecies])
        Lf = Li * Lj     
        topLf, topScores = scnorm.GetTopPercentileOfScores(Lf, scores, 95)   
        if len(topScores) > 1:
            fittingParameters = scnorm.CalculateFittingParameters(topLf, topScores)  
            return scnorm.NormaliseScoresByLogLengthProduct(B, Lengths[iSpecies], Lengths[jSpecies], fittingParameters)
        else:
            print("WARNING: Too few hits between species %d and species %d to normalise the scores, these hits will be ignored" % (iSpecies, jSpecies))
            return sparse.lil_matrix(B.get_shape())
            
    @staticmethod
    def ProcessBlastHits(seqsInfo, fileInfo, Lengths, iSpecies):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            # process up to the best hits for each species
            Bi = []
            for jSpecies in xrange(seqsInfo.nSpecies):
                Bij = BlastFileProcessor.GetBLAST6Scores(seqsInfo, fileInfo, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[jSpecies])  
                Bij = WaterfallMethod.NormaliseScores(Bij, Lengths, iSpecies, jSpecies)
                Bi.append(Bij)
            matrices.DumpMatrixArray("B", Bi, fileInfo, iSpecies)
            BH = GetBH_s(Bi, seqsInfo, iSpecies)
            matrices.DumpMatrixArray("BH", BH, fileInfo, iSpecies)
            util.PrintTime("Initial processing of species %d complete" % iSpecies)
        
    @staticmethod 
    def Worker_ProcessBlastHits(cmd_queue):
        while True:
            try:
                args = cmd_queue.get(True, 1)
                WaterfallMethod.ProcessBlastHits(*args)
            except Queue.Empty:
                return 

    @staticmethod
    def ConnectCognates(seqsInfo, fileInfo, iSpecies): 
        # calculate RBH for species i
        BHix = matrices.LoadMatrixArray("BH", fileInfo, seqsInfo, iSpecies)
        BHxi = matrices.LoadMatrixArray("BH", fileInfo, seqsInfo, iSpecies, row=False)
        RBHi = matrices.MatricesAndTr_s(BHix, BHxi)   # twice as much work as before (only did upper triangular before)
        B = matrices.LoadMatrixArray("B", fileInfo, seqsInfo, iSpecies)
        connect = WaterfallMethod.ConnectAllBetterThanAnOrtholog_s(RBHi, B, seqsInfo, iSpecies) 
        matrices.DumpMatrixArray("connect", connect, fileInfo, iSpecies)
            
    @staticmethod 
    def Worker_ConnectCognates(cmd_queue):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            while True:
                try:
                    args = cmd_queue.get(True, 1)
                    WaterfallMethod.ConnectCognates(*args)
                except Queue.Empty:
                    return  
                                   
    @staticmethod
    def WriteGraphParallel(seqsInfo, fileInfo, nProcess):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            with open(fileInfo.graphFilename + "_header", 'wb') as graphFile:
                graphFile.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n" % (seqsInfo.nSeqs, seqsInfo.nSeqs)) 
                graphFile.write("\n(mclmatrix\nbegin\n\n") 
            pool = mp.Pool(nProcess)
            pool.map(WriteGraph_perSpecies, [(seqsInfo, fileInfo, iSpec) for iSpec in xrange(seqsInfo.nSpecies)])
            subprocess.call("cat " + fileInfo.graphFilename + "_header " + " ".join([fileInfo.graphFilename + "_%d" % iSp for iSp in xrange(seqsInfo.nSpecies)]) + " > " + fileInfo.graphFilename, shell=True)
            # Cleanup
            os.remove(fileInfo.graphFilename + "_header")
            for iSp in xrange(seqsInfo.nSpecies): os.remove(fileInfo.graphFilename + "_%d" % iSp)
            DeleteMatrices(fileInfo) 
                
    @staticmethod
    def GetMostDistant_s(RBH, B, seqsInfo, iSpec):
        mostDistant = numeric.transpose(np.ones(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]])*1e9)
        for kSpec in xrange(seqsInfo.nSpecies):
            B[kSpec] = B[kSpec].tocsr()
            if iSpec == kSpec:
                continue
            I, J = RBH[kSpec].nonzero()
            if len(I) > 0:
                mostDistant[I] = np.minimum(B[kSpec][I, J], mostDistant[I])
        return mostDistant
    
    @staticmethod
    def ConnectAllBetterThanCutoff_s(B, mostDistant, seqsInfo, iSpec):
        connect = []
        nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]
        for jSpec in xrange(seqsInfo.nSpecies):
            M=B[jSpec].tolil()
            if iSpec != jSpec:
                IIJJ = [(i,j) for i, (valueRow, indexRow) in enumerate(zip(M.data, M.rows)) for j, v in zip(indexRow, valueRow) if v >= mostDistant[i]]
            else:
                IIJJ = [(i,j) for i, (valueRow, indexRow) in enumerate(zip(M.data, M.rows)) for j, v in zip(indexRow, valueRow) if (i != j) and v >= mostDistant[i]]
            II = [i for (i, j) in IIJJ]
            JJ = [j for (i, j) in IIJJ]
            onesArray = np.ones(len(IIJJ))
            mat = sparse.csr_matrix( (onesArray,  (II, JJ)), shape=(nSeqs_i,  seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[jSpec]]))
            connect.append(mat)
        return connect
    
    @staticmethod
    def ConnectAllBetterThanAnOrtholog_s(RBH, B, seqsInfo, iSpec):        
        mostDistant = WaterfallMethod.GetMostDistant_s(RBH, B, seqsInfo, iSpec) 
        connect = WaterfallMethod.ConnectAllBetterThanCutoff_s(B, mostDistant, seqsInfo, iSpec)
        return connect

"""
Stats
-------------------------------------------------------------------------------
"""
def OrthogroupsMatrix(iSpecies, properOGs):
    speciesIndexDict = {iSp:iCol for iCol, iSp in enumerate(iSpecies)}
    nSpecies = len(iSpecies)
    nGroups = len(properOGs)
    # (i, j)-th entry of ogMatrix gives the number of genes from i in orthogroup j
    ogMatrix = np.zeros((nGroups, nSpecies)) 
    for i_og, og in enumerate(properOGs):
        for species, _ in og:
            ogMatrix[i_og, speciesIndexDict[species]] += 1
    return ogMatrix
  
def Stats_SpeciesOverlaps(fn, speciesNamesDict, iSpecies, speciesPresence):
    """ Number of orthogroups in which each species-pair is present. Called by Stats"""
    with open(fn, 'wb') as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow([""] + [speciesNamesDict[index] for index in iSpecies])
        for iSp in iSpecies:
            overlap = [len([1 for og in speciesPresence if (iSp in og and jSp in og)]) for jSp in iSpecies]
            writer.writerow([speciesNamesDict[iSp]] + overlap)
 
def Stats_SizeTable(writer_sum, writer_sp, properOGs, allGenesCounter, iSpecies, speciesPresence):
    """ Overall and per-species histogram tables of orthogroup sizes. Called by Stats"""
    bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 21, 51, 101, 151, 201, 501, 1001, 9e99]
    writer_sum.writerow([])
    writer_sp.writerow([])
    nGenesPerOG = [len(og) for og in properOGs]
    nGenesPerOGPerSpecies = [[len([1 for g in og if g[0] == iSp]) for og in properOGs] for iSp in iSpecies]
    counters_GenesPerOGPerSpecies = [Counter(spCount) for spCount in nGenesPerOGPerSpecies]
    counters_GenesPerOG = Counter(nGenesPerOG)
    nSp = len(iSpecies)
    nOGs = len(properOGs)
    nGenesTotal = sum([len(og) for og in properOGs])
    percentFormat = "%0.1f"
    writer_sum.writerow(["Average number of genes per-species in orthogroup", "Number of orthogroups", "Percentage of orthogroups", "Number of genes", "Percentage of genes"])
    # Per-species tables
    table_NO = [["Number of genes per-species in orthogroup"] + ["Number of orthogroups" for _ in iSpecies]]      # number of orthogroups
    table_PO = [["Number of genes per-species in orthogroup"] + ["Percentage of orthogroups" for _ in iSpecies]]  # percentage of orthogroups
    table_NG = [["Number of genes per-species in orthogroup"] + ["Number of genes" for _ in iSpecies]]            # Number of genes
    table_PG = [["Number of genes per-species in orthogroup"] + ["Percentage of genes" for _ in iSpecies]]        # percentage of genes
    for start, end in zip(bins, bins[1:]):
        binName = "<1" if start == 0 else ("'%d" % start) if start+1 == end else "'%d+" % start if end == 9e99 else "%d-%d" % (start, end-1)
        nOrthogroups = sum([count for size, count in counters_GenesPerOG.items() if start*nSp<=size<end*nSp])
        nGenes = sum([size*count for size, count in counters_GenesPerOG.items() if start*nSp<=size<end*nSp])
        row_sum = [binName, nOrthogroups, percentFormat % (100.*nOrthogroups/nOGs), nGenes, percentFormat % (100.*nGenes/nGenesTotal)]
        writer_sum.writerow(row_sum)        
        # Per-species stats
        if binName == "<1": binName = "'0"
        nOrthogroups_ps = [sum([number for size, number in c.items() if start<=size<end]) for c in counters_GenesPerOGPerSpecies]
        nGenes_ps = [sum([size*number for size, number in c.items() if start<=size<end]) for c in counters_GenesPerOGPerSpecies]
        table_NO.append([binName] + nOrthogroups_ps)
        table_PO.append([binName] + [percentFormat % (100.*n/nOGs) for n in nOrthogroups_ps])
        table_NG.append([binName] + nGenes_ps)
        table_PG.append([binName] + [percentFormat % (100.*n/allGenesCounter[iSp]) for iSp, n in zip(iSpecies, nGenes_ps)])
    writer_sp.writerow([])
    for r in table_NO: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_PO: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_NG: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_PG: writer_sp.writerow(r)
        
    # Species presence
    n = map(len, speciesPresence)
    writer_sum.writerow([])
    writer_sum.writerow(["Number of species in orthogroup", "Number of orthogroups"])
    for i in xrange(1, nSp+1):
        writer_sum.writerow([i, n.count(i)])

def Stats(ogs, speciesNamesDict, iSpecies, resultsDir, iResultsVersion):
    """ Top-level method for calcualtion of stats for the orthogroups"""
    allOgs = [[map(int, g.split("_")) for g in og] for og in ogs]
    properOGs = [og for og in allOgs if len(og) > 1]
    allGenes = [g for og in allOgs for g in og]
    filename_sp = resultsDir +  "Statistics_PerSpecies" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".csv"
    filename_sum = resultsDir +  "Statistics_Overall" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".csv"
    filename_overlap = resultsDir +  "Orthogroups_SpeciesOverlaps" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".csv"
    percentFormat = "%0.1f"
    with open(filename_sp, 'wb') as outfile_species, open(filename_sum, 'wb') as outfile_sum:
        writer_sp = csv.writer(outfile_species, delimiter="\t")
        writer_sum = csv.writer(outfile_sum, delimiter="\t")
        # header
        writer_sp.writerow([""] + [speciesNamesDict[index] for index in iSpecies])
        
        # Number of genes
        allGenesCounter = Counter([g[0] for g in allGenes])
        nGenes = sum(allGenesCounter.values())
#        writer_sum.writerow(["Number of species", len(iSpecies)])
        writer_sp.writerow(["Number of genes"] + [allGenesCounter[iSp] for iSp in iSpecies])
        writer_sum.writerow(["Number of genes", nGenes])
        
        # Number of assigned/unassigned genes
        assignedGenesCounter = Counter([g[0] for og in properOGs for g in og])
        nAssigned = sum(assignedGenesCounter.values())
        writer_sp.writerow(["Number of genes in orthogroups"] + [assignedGenesCounter[iSp] for iSp in iSpecies])
        writer_sum.writerow(["Number of genes in orthogroups"] + [nAssigned])
        writer_sp.writerow(["Number of unassigned genes"] + [allGenesCounter[iSp] - assignedGenesCounter[iSp] for iSp in iSpecies])
        writer_sum.writerow(["Number of unassigned genes"] + [nGenes - nAssigned])     
        # Percentages
        pAssigned = 100.*nAssigned/nGenes
        writer_sp.writerow(["Percentage of genes in orthogroups"] + [percentFormat % (100.*assignedGenesCounter[iSp]/allGenesCounter[iSp]) for iSp in iSpecies])
        writer_sum.writerow(["Percentage of genes in orthogroups", percentFormat % pAssigned])   
        writer_sp.writerow(["Percentage of unassigned genes"] + [percentFormat % (100*(1.-(float(assignedGenesCounter[iSp])/allGenesCounter[iSp]))) for iSp in iSpecies])
        writer_sum.writerow(["Percentage of unassigned genes", percentFormat % (100*(1.-(float(nAssigned)/nGenes)))])
        
        # Number of Orthogroups
        speciesPresence = [set([g[0] for g in og]) for og in properOGs]
        nOgs = len(properOGs)
        writer_sum.writerow(["Number of orthogroups", nOgs])
        writer_sp.writerow(["Number of orthogroups containing species"] + [sum([iSp in og_sp for og_sp in speciesPresence]) for iSp in iSpecies])
        writer_sp.writerow(["Percentage of orthogroups containing species"] + [percentFormat % (100.*sum([iSp in og_sp for og_sp in speciesPresence])/len(properOGs)) for iSp in iSpecies])
        
        # Species specific orthogroups - orthogroups-based
        speciesSpecificOGsCounter = Counter([next(iter(og_sp)) for og_sp in speciesPresence if len(og_sp) == 1])
        writer_sp.writerow(["Number of species-specific orthogroups"] + [speciesSpecificOGsCounter[iSp] for iSp in iSpecies])
        writer_sum.writerow(["Number of species-specific orthogroups", sum(speciesSpecificOGsCounter.values())])
        
        # Species specific orthogroups - gene-based
        iSpeciesSpecificOGs = [i for i, og_sp in enumerate(speciesPresence) if len(og_sp) == 1] 
        iSpSpecificOGsGeneCounts = [sum([len(properOGs[iog]) for iog in iSpeciesSpecificOGs if properOGs[iog][0][0] == iSp]) for iSp in iSpecies]
        writer_sp.writerow(["Number of genes in species-specific orthogroups"] + iSpSpecificOGsGeneCounts)
        writer_sum.writerow(["Number of genes in species-specific orthogroups", sum(iSpSpecificOGsGeneCounts)])
        writer_sp.writerow(["Percentage of genes in species-specific orthogroups"] + [percentFormat % (100.*n_ss/allGenesCounter[iSp]) for n_ss, iSp in zip(iSpSpecificOGsGeneCounts, iSpecies)])
        writer_sum.writerow(["Percentage of genes in species-specific orthogroups", percentFormat % (100.*sum(iSpSpecificOGsGeneCounts)/nGenes)])
        
        # 'averages'
        l = list(reversed(map(len, properOGs)))
        writer_sum.writerow(["Mean orthogroup size", "%0.1f" % np.mean(l)])
        writer_sum.writerow(["Median orthogroup size", np.median(l)])
        L = np.cumsum(l)
        j, _ = next((i, x) for i, x in enumerate(L) if x > nAssigned/2)
        writer_sum.writerow(["G50 (assigned genes)",l[j]])
        l2 = list(reversed(map(len, ogs)))
        L2 = np.cumsum(l2)
        j2, _ = next((i, x) for i, x in enumerate(L2) if x > nGenes/2)
        G50 = l2[j2]
        writer_sum.writerow(["G50 (all genes)", G50])
        writer_sum.writerow(["O50 (assigned genes)", len(l) - j])
        O50 = len(l2) - j2
        writer_sum.writerow(["O50 (all genes)", O50])
        
        # Single-copy orthogroups
        ogMatrix = OrthogroupsMatrix(iSpecies, properOGs)
        nSpecies = len(iSpecies)
        nPresent = (ogMatrix > np.zeros((1, nSpecies))).sum(1)
        nCompleteOGs = list(nPresent).count(nSpecies)
        singleCopyOGs = (ogMatrix == np.ones((1, nSpecies))).all(1).nonzero()[0]   
        nSingleCopy = len(singleCopyOGs)
        writer_sum.writerow(["Number of orthogroups with all species present", nCompleteOGs])
        writer_sum.writerow(["Number of single-copy orthogroups", nSingleCopy])
        
        # Results filenames
        writer_sum.writerow(["Date", str(datetime.datetime.now()).split()[0]])
        writer_sum.writerow(["Orthogroups file", "Orthogroups" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".csv"])
        writer_sum.writerow(["Unassigned genes file", "Orthogroups" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + "_UnassignedGenes.csv"])
        writer_sum.writerow(["Per-species statistics", os.path.split(filename_sp)[1]])
        writer_sum.writerow(["Overall statistics", os.path.split(filename_sum)[1]])
        writer_sum.writerow(["Orthogroups shared between species", os.path.split(filename_overlap)[1]])
        
        # Sizes
        Stats_SizeTable(writer_sum, writer_sp, properOGs, allGenesCounter, iSpecies, speciesPresence)
        Stats_SpeciesOverlaps(filename_overlap, speciesNamesDict, iSpecies, speciesPresence)

    statsFiles = "Orthogroup statistics:\n"
    statsFiles += "   " + "   ".join([os.path.split(fn)[1] for fn in [filename_sp, filename_sum, filename_overlap]])
    summaryText = """OrthoFinder assigned %d genes (%0.1f%% of total) to %d orthogroups. Fifty percent of all genes were in orthogroups 
with %d or more genes (G50 was %d) and were contained in the largest %d orthogroups (O50 was %d). There were %d 
orthogroups with all species present and %d of these consisted entirely of single-copy genes.""" % (nAssigned, pAssigned, nOgs, G50, G50, O50, O50, nCompleteOGs, nSingleCopy)
    return summaryText, statsFiles
          

"""
OrthoFinder
-------------------------------------------------------------------------------
"""   
g_mclInflation = 1.5

def CanRunBLAST():
    if util.CanRunCommand("makeblastdb -help") and util.CanRunCommand("blastp -help"):
        return True
    else:
        print("ERROR: Cannot run BLAST+")
        print("Please check BLAST+ is installed and that the executables are in the system path\n")
        return False

def CanRunMCL():
    command = "mcl -h"
    if util.CanRunCommand(command):
        return True
    else:
        print("ERROR: Cannot run MCL with the command \"%s\"" % command)
        print("Please check MCL is installed and in the system path\n")
        return False
    
def PrintHelp():  
    print("=== Simple Usage ===\n") 
    print("orthofinder -f fasta_directory [-t n_blast_threads]")
    print("\n    Infers orthogroups for the proteomes contained in fasta_directory using")
    print("    n_blast_threads in parallel for the BLAST searches and tree inference.")
    print("")         
    print("=== Arguments ===")
#    print("")         
    print("Control where analysis starts/finishes:")
    print("")
    print("""-f fasta_directory, --fasta fasta_directory
    Predict orthogroups for the proteins in the fasta files in the fasta_directory\n""")
    
    print("""-b blast_results_directory, --blast blast_results_directory
    Predict orthogroups using the pre-calcualted BLAST results in blast_results_directory.\n""")
    
    print("""-og, --only-groups
    Only infer orthogroups, do not infer gene trees of orthologues.\n""")
    
    print("""-op, --only-prepare
    Only prepare the BLAST input files in the format required by OrthoFinder.\n""" )
    
    print("""-fg, --from-groups
    Infer gene trees and orthologues starting from pre-computed orthogroups.\n""")
    
    print("""-ft, --from-trees
    Infer orthologues starting from pre-computed gene trees.\n""")
    
    print("Additional arguments:")
    print("")
    
    print("""-t n_blast_threads, --threads n_blast_threads
    The number of BLAST processes to be run simultaneously. [Default is %d]\n""" % util.nThreadsDefault)
    
    print("""-a n_orthofinder_threads, --algthreads n_orthofinder_threads
    The number of threads to use for the less readily parallelised parts of the OrthoFinder algorithm.
    There are a number of trade-offs involved, see manual for details. [Default is %d]\n""" % util.nAlgDefault)
    
#    print("""-s, --species-tree
#    Use user suplied rooted species tree for gene-tree/species-tree reconciliation.\n""")
    
    print("""-I inflation_parameter, --inflation inflation_parameter
    Specify a non-default inflation parameter for MCL. Not recommended. [Default is %0.1f]\n""" % g_mclInflation)
    
    print("""-x speciesInfoFilename, --orthoxml speciesInfoFilename
    Output the orthogroups in the orthoxml format using the information in speciesInfoFilename.\n""")
    
    print("""-s rootedSpeciesTree, --speciestree rootedSpeciesTree
    Use rootedSpeciesTree for gene-tree/species-tree reconciliation (i.e. orthologue inference).\n""")
        
    print("""-h, --help
    Print this help text""")
    print("")
    util.PrintCitation() 
    
"""
Main
-------------------------------------------------------------------------------
"""   

def GetDirectoryArgument(arg, args):
    if len(args) == 0:
        print("Missing option for command line argument %s" % arg)
        util.Fail()
    directory = os.path.abspath(args.pop(0))
    if not os.path.isfile(directory) and directory[-1] != os.sep: 
        directory += os.sep
    return directory

#def GetOrthogroupsDirectory(suppliedDir, options):
#    """
#    Possible directory structures
#    1. Default: 
#        FastaFiles/Results_<date>/                          <- Orthogroups spreadsheets                 
#        FastaFiles/Results_<date>/WorkingDirectory/         <- Sequence and BLAST files
#        FastaFiles/Results_<date>/Orthologues_<date>/       <- Orthologues
#        FastaFiles/Results_<date>/Orthologues_<date>/WorkingDirectory/, Trees/, Orthologues 
#    2. From BLAST: 
#        <MainDirectory>/                                    <- Orthogroups spreadsheets / Sequence and BLAST files
#        FastaFiles/Results_<date>/WorkingDirectory/
#        FastaFiles/Results_<date>/Orthologues_<date>/
#        FastaFiles/Results_<date>/Orthologues_<date>/WorkingDirectory/, Trees/, Orthologues 
#    """
   
# Control
class Options(object):#
    def __init__(self):
        self.nBlast = util.nThreadsDefault
        self.nProcessAlg = util.nAlgDefault
        self.qStartFromBlast = False  # remove, just store BLAST to do
        self.qStartFromFasta = False  # local to argument checking
        self.qStartFromGroups = False
        self.qStartFromTrees = False
        self.qStopAfterPrepare = False
        self.qStopAfterGroups = False
        self.speciesXMLInfoFN = None
        self.speciesTreeFN = None
        self.mclInflation = g_mclInflation
    
    def what(self):
        for k, v in self.__dict__.items():
            if v == True:
                print(k)

class Directories(object):
    def __init__(self):
        self.resultsDir = None           # directory for orthogroup results files
        self.workingDir = None           # Orthogroup inference workingDir
                                         # Will need to store 3 bits of information in total    
        self.speciesToUse = []           #       seqsInfo.iSpeciesToUse   - which to include for this analysis 
        self.nSpAll = None               #       seqsInfo.nSpAll => 0, 1, ..., nSpAll - 1 are valid species indices
        self.iFirstNewSpecies = None     #       iFirstNew   => (0, 1, ..., iFirstNew-1) are from previous and (iFirstNew, iFirstNew+1, ..., nSpecies-1) are the new species indices
    
    def IDsFilename(self):
        return self.workingDir + "SequenceIDs.txt"
    def SpeciesIdsFilename(self):
        return self.workingDir + "SpeciesIDs.txt"
                                 
def ProcessArgs():
    """ 
    Workflow
    | 1. Fasta Files | 2.  Prepare files    | 3.   Blast    | 4. Orthogroups    | 5.   Gene Trees     | 6.   Reconciliations/Orthologues   |

    Options
    Start from:
    -f: 1,2,..,6    (start from fasta files, --fasta)
    -b: 4,5,6       (start from blast results, --blast)
    -fg: 5,6         (start from orthogroups/do orthologue workflow, --from-groups)
    -ft: 6           (start from gene tree/do reconciliation, --from-trees)
    Stop at:
    -op: 2           (only prepare, --only-prepare)
    -og: 4           (orthogroups, --only-groups)
    """
    if len(sys.argv) == 1 or sys.argv[1] == "--help" or sys.argv[1] == "help" or sys.argv[1] == "-h":
        PrintHelp()
        sys.exit() 

    options = Options()
    fastaDir = None
    workingDir = None
    orthologuesDir = None
    
    """
    -f: store fastaDir
    -b: store workingDir
    -fg: store orthologuesDir 
    -ft: store orthologuesDir 
    + xml: speciesXMLInfoFN
    """    
    
    args = sys.argv[1:]
    while len(args) > 0:
        arg = args.pop(0)    
        if arg == "-f" or arg == "--fasta":
            if options.qStartFromFasta:
                print("Repeated argument: -f/--fasta")
                util.Fail()
            options.qStartFromFasta = True
            fastaDir = GetDirectoryArgument(arg, args)
        elif arg == "-b" or arg == "--blast":
            if options.qStartFromBlast:
                print("Repeated argument: -b/--blast")
                util.Fail()
            options.qStartFromBlast = True
            workingDir = GetDirectoryArgument(arg, args)
        elif arg == "-fg" or arg == "--from-groups":
            if options.qStartFromGroups:
                print("Repeated argument: -fg/--from-groups")
                util.Fail()
            options.qStartFromGroups = True
            workingDir = GetDirectoryArgument(arg, args)
        elif arg == "-ft" or arg == "--from-trees":
            if options.qStartFromTrees:
                print("Repeated argument: -ft/--from-trees")
                util.Fail()
            options.qStartFromTrees = True
            orthologuesDir = GetDirectoryArgument(arg, args)
        elif arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument -t")
                util.Fail()
            arg = args.pop(0)
            try:
                options.nBlast = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s" % arg)
                util.Fail()    
        elif arg == "-a" or arg == "--algthreads":
            if len(args) == 0:
                print("Missing option for command line argument -a")
                util.Fail()
            arg = args.pop(0)
            try:
                options.nProcessAlg = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s" % arg)
                util.Fail()   
        elif arg == "-I" or arg == "--inflation":
            if len(args) == 0:
                print("Missing option for command line argument -I")
                util.Fail()
            arg = args.pop(0)
            try:
                options.mclInflation = float(arg)
            except:
                print("Incorrect argument for MCL inflation parameter: %s" % arg)
                util.Fail()    
        elif arg == "-x" or arg == "--orthoxml":  
            if options.speciesXMLInfoFN:
                print("Repeated argument: -x/--orthoxml")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s" % arg)
                util.Fail()
            options.speciesXMLInfoFN = args.pop(0)
        elif arg == "-s" or arg == "--speciestree":  
            if options.speciesXMLInfoFN:
                print("Repeated argument: -s/--speciestree")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s" % arg)
                util.Fail()
            options.speciesTreeFN = args.pop(0)
        elif arg == "-op" or arg == "--only-prepare":
            options.qStopAfterPrepare = True
        elif arg == "-og" or arg == "--only-groups":
            options.qStopAfterGroups = True
        elif arg == "-o" or arg == "--orthologues":
            options.qOrthogroups = False
        elif arg == "-h" or arg == "--help":
            PrintHelp()
            sys.exit()
        else:
            print("Unrecognised argument: %s\n" % arg)
            util.Fail()    
    
    # check argument combinations           
    if options.qStartFromTrees and options.qStartFromFasta:
        print("ERROR: Incompatible arguments, -f (start from fasta files) and -g (start from orthogroups)")
        util.Fail()
        
    if options.qStartFromTrees and options.qStartFromBlast:
        print("ERROR: Incompatible arguments, -b (start from pre-calcualted BLAST results) and -g (start from orthogroups)")
        util.Fail()      

    print("%d thread(s) for highly parallel tasks (BLAST searches etc.)" % options.nBlast)
    print("%d thread(s) for OrthoFinder algorithm" % options.nProcessAlg)
    return options, fastaDir, workingDir, orthologuesDir              

def GetXMLSpeciesInfo(dirs, options):
    # speciesInfo:  name, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
    util.PrintUnderline("Reading species information file")
    # do this now so that we can alert user to any errors prior to running the algorithm
    speciesInfo = [[] for i_ in dirs.speciesToUse]
    speciesNamesDict = SpeciesNameDict(dirs.SpeciesIdsFilename())
    speciesRevDict = {v:k for k,v in speciesNamesDict.items()}
    userFastaFilenames = [os.path.split(speciesNamesDict[i])[1] for i in dirs.speciesToUse]
    with open(options.speciesXMLInfoFN, 'rb') as speciesInfoFile:
        reader = csv.reader(speciesInfoFile, delimiter = "\t")
        for iLine, line in enumerate(reader):
            if len(line) != 5:
                # allow for an extra empty line at the end
                if len(line) == 0 and iLine == len(userFastaFilenames):
                    continue
                print("ERROR")
                print("Species information file %s line %d is incorrectly formatted." % (options.speciesXMLInfoFN, iLine + 1))        
                print("File should be contain one line per species")
                print("Each line should contain 5 tab-delimited fields:")
                print("  fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseFastaFilename")
                print("See README file for more information.")
                util.Fail() 
            fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile = line
            try:
                iSpecies = speciesRevDict[os.path.splitext(fastaFilename)[0]]
            except KeyError:
                print("Skipping %s from line %d as it is not being used in this analysis" % (fastaFilename, iLine+1))
                continue
            speciesInfo[dirs.speciesToUse.index(iSpecies)] = line   
    # check information has been provided for all species
    speciesMissing = False        
    for iPos, iSpecies in enumerate(dirs.speciesToUse):
        if speciesInfo[iPos] == []:
            if not speciesMissing:
                print("ERROR")
                print("Species information file %s does not contain information for all species." % options.speciesXMLInfoFN) 
                print("Information is missing for:") 
                speciesMissing = True
            print(speciesNamesDict[iSpecies])
    if speciesMissing:
        util.Fail()
    return speciesInfo

def CheckDependencies(options):
    util.PrintUnderline("1. Checking required programs are installed")
    if (not options.qStartFromBlast) and (not CanRunBLAST()):
        util.Fail()
    if not CanRunMCL():
        util.Fail()

def DoOrthogroups(options, dirs, seqsInfo):
    # Run Algorithm, cluster and output cluster files with original accessions
    util.PrintUnderline("Running OrthoFinder algorithm")
    fileIdentifierString = "OrthoFinder_v%s" % util.version
    graphFilename = dirs.workingDir + "%s_graph.txt" % fileIdentifierString
    # it's important to free up the memory from python used for processing the genomes
    # before launching MCL becuase both use sizeable ammounts of memory. The only
    # way I can find to do this is to launch the memory intensive python code 
    # as separate process that exitsbefore MCL is launched.
    fileInfo = util.FileInfo(workingDir = dirs.workingDir, graphFilename=graphFilename) 
    Lengths = GetSequenceLengths(seqsInfo, fileInfo)
    
    # Process BLAST hits
    util.PrintTime("Initial processing of each species")
    cmd_queue = mp.Queue()
    for iSpecies in xrange(seqsInfo.nSpecies):
        cmd_queue.put((seqsInfo, fileInfo, Lengths, iSpecies))
    runningProcesses = [mp.Process(target=WaterfallMethod.Worker_ProcessBlastHits, args=(cmd_queue, )) for i_ in xrange(options.nProcessAlg)]
    for proc in runningProcesses:
        proc.start()
    util.ManageQueue(runningProcesses, cmd_queue)
    
    cmd_queue = mp.Queue()
    for iSpecies in xrange(seqsInfo.nSpecies):
        cmd_queue.put((seqsInfo, fileInfo, iSpecies))
    runningProcesses = [mp.Process(target=WaterfallMethod.Worker_ConnectCognates, args=(cmd_queue, )) for i_ in xrange(options.nProcessAlg)]
    for proc in runningProcesses:
        proc.start()
    util.ManageQueue(runningProcesses, cmd_queue)
    
    util.PrintTime("Connected putatitive homologs") 
    WaterfallMethod.WriteGraphParallel(seqsInfo, fileInfo, options.nProcessAlg)
    
    # 5b. MCL     
    clustersFilename, iResultsVersion = util.GetUnusedFilename(dirs.workingDir  + "clusters_%s_I%0.1f" % (fileIdentifierString, options.mclInflation), ".txt")
    MCL.RunMCL(graphFilename, clustersFilename, options.nProcessAlg, options.mclInflation)
    clustersFilename_pairs = clustersFilename + "_id_pairs.txt"
    MCLread.ConvertSingleIDsToIDPair(seqsInfo, clustersFilename, clustersFilename_pairs)   
    
    util.PrintUnderline("Writing orthogroups to file")
    if options.qStopAfterGroups: util.PrintCitation()
    ogs = MCLread.GetPredictedOGs(clustersFilename_pairs)
    resultsBaseFilename = util.GetUnusedFilename(dirs.resultsDir + "Orthogroups", ".csv")[:-4]         # remove .csv from base filename
    resultsBaseFilename = dirs.resultsDir + "Orthogroups" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion)
    idsDict = MCL.WriteOrthogroupFiles(ogs, [dirs.IDsFilename()], resultsBaseFilename, clustersFilename_pairs)
    speciesNamesDict = SpeciesNameDict(dirs.SpeciesIdsFilename())
    orthogroupsResultsFilesString = MCL.CreateOrthogroupTable(ogs, idsDict, speciesNamesDict, dirs.speciesToUse, resultsBaseFilename)
    print(orthogroupsResultsFilesString)
    summaryText, statsFile = Stats(ogs, speciesNamesDict, dirs.speciesToUse, dirs.resultsDir, iResultsVersion)
    if options.speciesXMLInfoFN:
        MCL.WriteOrthoXML(speciesInfo, ogs, seqsInfo.nSeqsPerSpecies, idsDict, resultsBaseFilename + ".orthoxml", dirs.speciesToUse)
    return clustersFilename_pairs, statsFile, summaryText, orthogroupsResultsFilesString

# 0
def ProcessPreviousFiles(workingDir):
    """Checks for:
    SpeciesIDs.txt
    Species*.fa
    Blast*txt
    SequenceIDs.txt
    
    Checks which species should be included
    """
    dirs = Directories()
    dirs.workingDir = workingDir
    if not os.path.exists(dirs.SpeciesIdsFilename()):
        print("%s file must be provided if using previously calculated BLAST results" % dirs.SpeciesIdsFilename())
        util.Fail()
    dirs.speciesToUse, dirs.nSpAll = util.GetSpeciesToUse(dirs.SpeciesIdsFilename())
    dirs.resultsDir = dirs.workingDir
    
    # check BLAST results directory exists
    if not os.path.exists(dirs.workingDir):
        print("Previous/Pre-calculated BLAST results directory does not exist: %s\n" % dirs.workingDir)
        util.Fail()
 
    # check fasta files are present 
    previousFastaFiles = util.SortFastaFilenames(glob.glob(dirs.workingDir + "Species*.fa"))
    if len(previousFastaFiles) == 0:
        print("No processed fasta files in the supplied previous working directory: %s\n" % dirs.workingDir)
        util.Fail()
    tokens = previousFastaFiles[-1][:-3].split("Species")
    lastFastaNumberString = tokens[-1]
    iLastFasta = 0
    nFasta = len(previousFastaFiles)
    try:
        iLastFasta = int(lastFastaNumberString)
    except:
        print("Filenames for processed fasta files are incorrect: %s\n" % previousFastaFiles[-1])
        util.Fail()
    if nFasta != iLastFasta + 1:
        print("Not all expected fasta files are present. Index of last fasta file is %s but found %d fasta files.\n" % (lastFastaNumberString, len(previousFastaFiles)))
        util.Fail()
    
    # check BLAST files
    qHaveBlast = True
    for iSpecies in dirs.speciesToUse:
        for jSpecies in dirs.speciesToUse:
            filename = "%sBlast%d_%d.txt" % (dirs.workingDir, iSpecies, jSpecies) 
            if not os.path.exists(filename):
                print("BLAST results file is missing: %s" % filename)
                qHaveBlast = False
    if not qHaveBlast: util.Fail()
                
    # check SequenceIDs.txt and SpeciesIDs.txt files are present
    if not os.path.exists(dirs.IDsFilename()):
        print("%s file must be provided if using previous calculated BLAST results" % dirs.IDsFilename())
        util.Fail()
    return dirs

# 6
def CreateBlastDatabases(dirs):
    nDB = max(dirs.speciesToUse) + 1
    for iSp in xrange(nDB):
        command = ["makeblastdb", "-dbtype", "prot", "-in", dirs.workingDir + "Species%d.fa" % iSp, "-out", dirs.workingDir + "BlastDBSpecies%d" % iSp]
        util.PrintTime("Creating Blast database %d of %d" % (iSp + 1, nDB))
        RunBlastDBCommand(command) 

# 7
def RunBlast(options, dirs, seqsInfo):
    if options.qStopAfterPrepare:
        util.PrintUnderline("BLAST commands that must be run")
    else:        
        util.PrintUnderline("Running BLAST all-versus-all")
    commands = GetOrderedBlastCommands(seqsInfo, dirs)
    if options.qStopAfterPrepare:
        for command in commands:
            print(" ".join(command))
        sys.exit()
    print("Using %d thread(s)" % options.nBlast)
    util.PrintTime("This may take some time....")  
    cmd_queue = mp.Queue()
    for iCmd, cmd in enumerate(commands):
        cmd_queue.put((iCmd+1, cmd))           
    runningProcesses = [mp.Process(target=util.Worker_RunCommand, args=(cmd_queue, options.nBlast, len(commands))) for i_ in xrange(options.nBlast)]
    for proc in runningProcesses:
        proc.start()#
    for proc in runningProcesses:
        while proc.is_alive():
            proc.join()
    # remove BLAST databases
    for f in glob.glob(dirs.workingDir + "BlastDBSpecies*"):
        os.remove(f)

def GetOrthogroupsFile(userDir):
    return util.GetOGsFile(userDir)

# 9
def GetOrthologues(dirs, options, clustersFilename_pairs, orthogroupsResultsFilesString=None):
    util.PrintUnderline("Running Orthologue Prediction", True)
    orthologuesResultsFilesString = get_orthologues.GetOrthologues(dirs.workingDir, dirs.resultsDir, dirs.speciesToUse, dirs.nSpAll, clustersFilename_pairs, options.nBlast, options.speciesTreeFN)
    if None != orthogroupsResultsFilesString: print(orthogroupsResultsFilesString)
    print(orthologuesResultsFilesString.rstrip())    

def GetOrthologues_FromTrees(orthologuesDir, userSpeciesTreeFN = None):
    groupsDir = orthologuesDir + "../"
    workingDir = orthologuesDir + "WorkingDirectory/"
    speciesTree_fn = userSpeciesTreeFN if userSpeciesTreeFN else workingDir + "Trees_ids/SpeciesTree_ids_0_rooted.txt"
    qUserSpTree = userSpeciesTreeFN != None
    return get_orthologues.JustOrthologues(groupsDir, speciesTree_fn, workingDir, qUserSpTree)
 
def ProcessesNewFasta(fastaDir, existingDirs=None):
    """
    Process fasta files and return a Directory object with all paths completed.
    """
    # Check files present
    if not os.path.exists(fastaDir):
        print("\nDirectory does not exist: %s" % fastaDir)
        util.Fail()
    originalFastaFilenames = sorted([f for f in os.listdir(fastaDir) if os.path.isfile(os.path.join(fastaDir,f))])
    originalFastaFilenames = [f for f in originalFastaFilenames if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in fastaExtensions]
    if len(originalFastaFilenames) == 0:
        print("\nNo fasta files found in supplied directory: %s" % fastaDir)
        util.Fail()
    if None == existingDirs:
        dirs = Directories()
        dirs.resultsDir = util.CreateNewWorkingDirectory(fastaDir + "Results_")
        dirs.workingDir = dirs.resultsDir + "WorkingDirectory" + os.sep
        os.mkdir(dirs.workingDir)
    else:
        dirs = existingDirs
    iSeq = 0
    iSpecies = 0
    # check if SpeciesIDs.txt already exists
    if os.path.exists(dirs.SpeciesIdsFilename()):
        with open(dirs.SpeciesIdsFilename(), 'rb') as infile:
            for line in infile: pass
        if line.startswith("#"): line = line[1:]
        iSpecies = int(line.split(":")[0]) + 1
    dirs.iFirstNewSpecies = iSpecies
    newSpeciesIDs = []
    with open(dirs.IDsFilename(), 'ab') as idsFile, open(dirs.SpeciesIdsFilename(), 'ab') as speciesFile:
        for fastaFilename in originalFastaFilenames:
            newSpeciesIDs.append(iSpecies)
            outputFastaFilename = dirs.workingDir + "Species%d.fa" % iSpecies
            outputFasta = open(outputFastaFilename, 'wb')
            fastaFilename = fastaFilename.rstrip()
            speciesFile.write("%d: %s\n" % (iSpecies, fastaFilename))
            baseFilename, extension = os.path.splitext(fastaFilename)
            with open(fastaDir + os.sep + fastaFilename, 'rb') as fastaFile:
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
    dirs.speciesToUse = dirs.speciesToUse + newSpeciesIDs
    dirs.nSpAll = max(dirs.speciesToUse) + 1      # will be one of the new species
    return dirs
            
def CheckOptions(options, dirs):
    """Check any optional arguments are valid once we know what species are in the analysis
    - user supplied species tree
    """
    if options.speciesTreeFN:
        expSpecies = SpeciesNameDict(dirs.SpeciesIdsFilename()).values()
        get_orthologues.CheckUserSpeciesTree(options.speciesTreeFN, expSpecies)

if __name__ == "__main__":    
    print("OrthoFinder version %s Copyright (C) 2014 David Emms\n" % util.version)
    print("""    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software, and you are welcome to redistribute it under certain conditions.
    For details please see the License.md that came with this software.\n""")
    options, fastaDir, workingDir, orthologuesDir = ProcessArgs()  
    # 2.
    CheckDependencies(options) 
    
    # if using previous Trees etc., check these are all present - Job for get_orthologues
    if options.qStartFromBlast and options.qStartFromFasta:
        # 0. Check Files
        dirs = ProcessPreviousFiles(workingDir)
        print("Adding new species in %s to existing analysis in %s" % (fastaDir, dirs.workingDir))
        # 3. 
        dirs = ProcessesNewFasta(fastaDir, dirs)
        CheckOptions(options, dirs)
        # 4.
        seqsInfo = util.GetSeqsInfo(dirs.workingDir, dirs.speciesToUse, dirs.nSpAll)
        # 5.
        if options.speciesXMLInfoFN:   
            speciesInfo = GetXMLSpeciesInfo(dirs, options)
        # 6.    
        util.PrintUnderline("Dividing up work for BLAST for parallel processing")
        CreateBlastDatabases(dirs)
        # 7.  
        RunBlast(options, dirs, seqsInfo)
        # 8.
        clustersFilename_pairs, statsFile, summaryText, orthogroupsResultsFilesString = DoOrthogroups(options, dirs, seqsInfo)
        # 9.
        if not options.qStopAfterGroups:
            GetOrthologues(dirs, options, clustersFilename_pairs, orthogroupsResultsFilesString)
        # 10.
        print("\n" + statsFile + "\n\n" + summaryText) 
        util.PrintCitation()
            
    elif options.qStartFromFasta:
        # 3. 
        dirs = ProcessesNewFasta(fastaDir)
        CheckOptions(options, dirs)
        # 4
        seqsInfo = util.GetSeqsInfo(dirs.workingDir, dirs.speciesToUse, dirs.nSpAll)
        # 5.
        if options.speciesXMLInfoFN:   
            speciesInfo = GetXMLSpeciesInfo(dirs, options)
        # 6.    
        util.PrintUnderline("Dividing up work for BLAST for parallel processing")
        CreateBlastDatabases(dirs)
        # 7. 
        RunBlast(options, dirs, seqsInfo)
        # 8.  
        clustersFilename_pairs, statsFile, summaryText, orthogroupsResultsFilesString = DoOrthogroups(options, dirs, seqsInfo)    
        # 9. 
        if not options.qStopAfterGroups:
            GetOrthologues(dirs, options, clustersFilename_pairs, orthogroupsResultsFilesString)
        # 10.
        print("\n" + statsFile + "\n\n" + summaryText) 
        util.PrintCitation()
        
    elif options.qStartFromBlast:
        # 0.
        dirs = ProcessPreviousFiles(workingDir)
        print("Using previously calculated BLAST results in %s" % dirs.workingDir) 
        CheckOptions(options, dirs)
        # 4.
        seqsInfo = util.GetSeqsInfo(dirs.workingDir, dirs.speciesToUse, dirs.nSpAll)
        # 5.
        if options.speciesXMLInfoFN:   
            speciesInfo = GetXMLSpeciesInfo(dirs, options)
        # 8        
        clustersFilename_pairs, statsFile, summaryText, orthogroupsResultsFilesString = DoOrthogroups(options, dirs, seqsInfo)    
        # 9
        if not options.qStopAfterGroups:
            GetOrthologues(dirs, options, clustersFilename_pairs, orthogroupsResultsFilesString)
        # 10
        print("\n" + statsFile + "\n\n" + summaryText) 
        util.PrintCitation() 
    elif options.qStartFromGroups:
        # 0.
        orthofinderWorkingDir, orthofinderResultsDir, clustersFilename_pairs = GetOrthogroupsFile(workingDir)   
        dirs = ProcessPreviousFiles(orthofinderWorkingDir)
        dirs.resultsDir = orthofinderResultsDir       
        CheckOptions(options, dirs)
        # 9
        GetOrthologues(dirs, options, clustersFilename_pairs)
        # 10
        util.PrintCitation() 
    elif options.qStartFromTrees:
        # Don't call CheckOptions(options, dirs), Check species tree in GetOrthologues_FromTrees
        summaryText = GetOrthologues_FromTrees(orthologuesDir, options.speciesTreeFN)
        print(summaryText) 
        util.PrintCitation() 
    else:
        raise NotImplementedError
        
        
     

    

    

