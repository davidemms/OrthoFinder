#!/usr/bin/env python3
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

# first import parallel task manager to minimise RAM overhead for small processes
from __future__ import absolute_import
from . import parallel_task_manager

import os                                       # Y
os.environ["OPENBLAS_NUM_THREADS"] = "1"    # fix issue with numpy/openblas. Will mean that single threaded options aren't automatically parallelised 

import sys                                      # Y
import subprocess                               # Y
import glob                                     # Y
import shutil                                   # Y
import time                                     # Y
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
try: 
    import queue
except ImportError:
    import Queue as queue                       # Y
import warnings                                 # Y

PY2 = sys.version_info <= (3,)
csv_write_mode = 'wb' if PY2 else 'wt'

from . import blast_file_processor, files, mcl, util, matrices, orthologues, program_caller, trees_msa 

# Get directory containing script/bundle
if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
max_int = sys.maxsize
ok = False
while not ok:
    try:
        csv.field_size_limit(max_int)
        ok = True
    except OverflowError:
        max_int = int(max_int/10)
sys.setrecursionlimit(10**6)
    
fastaExtensions = {"fa", "faa", "fasta", "fas", "pep"}
# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f) 

my_env = os.environ.copy()
# use orthofinder supplied executables by preference
my_env['PATH'] = os.path.join(__location__, 'bin:') + my_env['PATH']
# Fix LD_LIBRARY_PATH when using pyinstaller 
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
    capture = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env, shell=True)
    stdout, stderr = capture.communicate()
    try:
        stdout = stdout.decode()
        stderr = stderr.decode()
    except (UnicodeDecodeError, AttributeError):
        stdout = stdout.encode()
        stderr = stderr.encode()
    n_stdout_lines = stdout.count("\n")
    n_stderr_lines = stderr.count("\n")
    nLines_success= 12
    if n_stdout_lines > nLines_success or n_stderr_lines > 0 or capture.returncode != 0:
        print("\nWARNING: Likely problem with input FASTA files")
        if capture.returncode != 0:
            print("makeblastdb returned an error code: %d" % capture.returncode)
        else:
            print("makeblastdb produced unexpected output")
        print("Command: %s" % " ".join(command))
        print("stdout:\n-------")
        print(stdout)
        if len(stderr) > 0:
            print("stderr:\n-------")
            print(stderr)
            
def SpeciesNameDict(speciesIDsFN):
    speciesNamesDict = dict()
    with open(speciesIDsFN, 'r') as speciesNamesFile:
        for line in speciesNamesFile:
            if line.startswith("#"): continue
            line = line.rstrip()
            if not line: continue
            short, full = line.split(": ")
            speciesNamesDict[int(short)] = full.rsplit(".", 1)[0]
    return speciesNamesDict
    
"""
MCL
-------------------------------------------------------------------------------
"""    
class MCL:
    @staticmethod
    def CreateOGs(predictedOGs, outputFilename, idDict):
        with open(outputFilename, 'w') as outputFile:
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
            for iGene_species in range(nSeqs):
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
                geneNode.set('id', str(mcl.GetSingleID(speciesStartingIndices, seq, speciesToUse)))
    #                    SubElement(geneNode, 'score')                  # skip
        with open(orthoxmlFilename, 'w') as orthoxmlFile:
    #            ET.ElementTree(root).write(orthoxmlFile)
            orthoxmlFile.write(MCL.prettify(root))
        print("Orthogroups have been written to orthoxml file:\n   %s" % orthoxmlFilename)
            
    @staticmethod               
    def RunMCL(graphFilename, clustersFilename, nProcesses, inflation):
        command = " ".join(["mcl", graphFilename, "-I", str(inflation), "-o", clustersFilename, "-te", str(nProcesses), "-V", "all"])
        parallel_task_manager.RunCommand(command, qPrintOnError=True)
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
            files.FileHandler.LogFailAndExit(("ERROR: Sequence ID not found in %s\n" % idsFilename) + str(e) + "\n")        
        except RuntimeError as error:
            print(str(error))
            if str(error).startswith("ERROR"):
                err_text = "ERROR: %s contains a duplicate ID. The IDs for the orthogroups in %s will not be replaced with the sequence accessions. If %s was prepared manually then please check the IDs are correct. " % (idsFilename, clustersFilename_pairs, idsFilename)
                files.FileHandler.LogFailAndExit(err_text)
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
                    err_text = "ERROR: %s contains a duplicate ID. The IDs for the orthogroups in %s will not be replaced with the sequence accessions. This is probably because the same accession was used more than once in your input FASTA files. However, if %s was prepared manually then you may need to check that file instead." % (idsFilename, clustersFilename_pairs, idsFilename) 
                    files.FileHandler.LogFailAndExit(err_text)
        return fullDict
    
    @staticmethod
    def CreateOrthogroupTable(ogs, 
                              idToNameDict, 
                              speciesNamesDict, 
                              speciesToUse,
                              resultsBaseFilename):
        
        nSpecies = len(speciesNamesDict) 
        
        ogs_names = [[idToNameDict[seq] for seq in og] for og in ogs]
        ogs_ints = [[list(map(int, sequence.split("_"))) for sequence in og] for og in ogs]
    
        # write out
        outputFilename = resultsBaseFilename + ".tsv"
        outputFilename_counts = resultsBaseFilename + ".GeneCount.tsv"
        singleGeneFilename = resultsBaseFilename + "_UnassignedGenes.tsv"
        with open(outputFilename, csv_write_mode) as outputFile, open(singleGeneFilename, csv_write_mode) as singleGeneFile, open(outputFilename_counts, csv_write_mode) as outFile_counts:
            fileWriter = csv.writer(outputFile, delimiter="\t")
            fileWriter_counts = csv.writer(outFile_counts, delimiter="\t")
            singleGeneWriter = csv.writer(singleGeneFile, delimiter="\t")
            for writer in [fileWriter, singleGeneWriter]:
                row = ["Orthogroup"] + [speciesNamesDict[index] for index in speciesToUse]
                writer.writerow(row)
            fileWriter_counts.writerow(row + ['Total'])
            
            for iOg, (og, og_names) in enumerate(zip(ogs_ints, ogs_names)):
                ogDict = defaultdict(list)
                row = ["OG%07d" % iOg]
                thisOutputWriter = fileWriter
                # separate it into sequences from each species
                if len(og) == 1:
                    row.extend(['' for x in range(nSpecies)])
                    row[speciesToUse.index(og[0][0]) + 1] = og_names[0]
                    thisOutputWriter = singleGeneWriter
                else:
                    for (iSpecies, iSequence), name in zip(og, og_names):
                        ogDict[speciesToUse.index(iSpecies)].append(name)
                    for iSpecies in range(nSpecies):
                        row.append(", ".join(sorted(ogDict[iSpecies])))
                    counts = Counter([iSpecies for iSpecies, _ in og])
                    counts_row = [counts[iSpecies] for iSpecies in speciesToUse]
                    fileWriter_counts.writerow(row[:1] + counts_row + [sum(counts_row)])
                thisOutputWriter.writerow(row)

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
        for i in range(nBins):
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
        rangeq = list(range(len(Lq)))
        rangeh = list(range(len(Lh)))
        li_vals = Lq**(-params[0])
        lj_vals = Lh**(-params[0])
        li_matrix = sparse.csr_matrix((li_vals, (rangeq, rangeq)))
        lj_matrix = sparse.csr_matrix((lj_vals, (rangeh, rangeh)))
        return sparse.lil_matrix(10**(-params[1]) * li_matrix * b * lj_matrix)

"""
RunInfo
-------------------------------------------------------------------------------
"""     

def GetSequenceLengths(seqsInfo):                
    sequenceLengths = []
    for iSpecies, iFasta in enumerate(seqsInfo.speciesToUse):
        sequenceLengths.append(np.zeros(seqsInfo.nSeqsPerSpecies[iFasta]))
        fastaFilename = files.FileHandler.GetSpeciesFastaFN(iFasta)
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
def GetOrderedSearchCommands(seqsInfo, speciesInfoObj, qDoubleBlast, search_program, prog_caller):
    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands 
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing 
    the BLAST commands.
    """
    iSpeciesPrevious = list(range(speciesInfoObj.iFirstNewSpecies))
    iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew) if (qDoubleBlast or i <=j)] + \
                   [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesPrevious) if (qDoubleBlast or i <=j)] + \
                   [(i, j) for i, j in itertools.product(iSpeciesPrevious, iSpeciesNew) if (qDoubleBlast or i <=j)] 
    taskSizes = [seqsInfo.nSeqsPerSpecies[i]*seqsInfo.nSeqsPerSpecies[j] for i,j in speciesPairs]
    taskSizes, speciesPairs = util.SortArrayPairByFirst(taskSizes, speciesPairs, True)
    if search_program == "blast":
        commands = [" ".join(["blastp", "-outfmt", "6", "-evalue", "0.001", "-query", files.FileHandler.GetSpeciesFastaFN(iFasta), "-db", files.FileHandler.GetSpeciesDatabaseN(iDB), "-out", files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)]) for iFasta, iDB in speciesPairs]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(search_program, files.FileHandler.GetSpeciesFastaFN(iFasta), files.FileHandler.GetSpeciesDatabaseN(iDB, search_program), files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)) for iFasta, iDB in speciesPairs]
    return commands     

"""
Matrices
-------------------------------------------------------------------------------
""" 
            
def GetBH_s(pairwiseScoresMatrices, seqsInfo, iSpecies, tol=1e-3):
    nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpecies]]
    bestHitForSequence = -1*np.ones(nSeqs_i)
    H = [None for i_ in range(seqsInfo.nSpecies)] # create array of Nones to be replace by matrices
    for j in range(seqsInfo.nSpecies):
        if iSpecies == j:
            # identify orthologs then come back to paralogs
            continue
        W = pairwiseScoresMatrices[j]
        I = []
        J = []
        for kRow in range(nSeqs_i):
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
    for kRow in range(nSeqs_i):
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
    seqsInfo, graphFN, iSpec, d_pickle = args            
    # calculate the 2-way connections for one query species
    with open(graphFN + "_%d" % iSpec, 'w') as graphFile:
        connect2 = []
        for jSpec in range(seqsInfo.nSpecies):
            m1 = matrices.LoadMatrix("connect", iSpec, jSpec, d_pickle)
            m2tr = numeric.transpose(matrices.LoadMatrix("connect", jSpec, iSpec, d_pickle))
            connect2.append(m1 + m2tr)
            del m1, m2tr
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpec, d_pickle)
        B_connect = matrices.MatricesAnd_s(connect2, B)
        del B, connect2
        
        W = [b.sorted_indices().tolil() for b in B_connect]
        del B_connect
        for query in range(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]):
            offset = seqsInfo.seqStartingIndices[iSpec]
            graphFile.write("%d    " % (offset + query))
            for jSpec in range(seqsInfo.nSpecies):
                row = W[jSpec].getrowview(query)
                jOffset = seqsInfo.seqStartingIndices[jSpec]
                for j, value in zip(row.rows[0], row.data[0]):
                    graphFile.write("%d:%.3f " % (j + jOffset, value))
            graphFile.write("$\n")
        if iSpec == (seqsInfo.nSpecies - 1): graphFile.write(")\n")
        util.PrintTime("Written final scores for species %d to graph file" % iSpec)
            
            
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
    def ProcessBlastHits(seqsInfo, blastDir_list, Lengths, iSpecies, d_pickle, qDoubleBlast):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            # process up to the best hits for each species
            Bi = []
            for jSpecies in range(seqsInfo.nSpecies):
                Bij = blast_file_processor.GetBLAST6Scores(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[jSpecies], qDoubleBlast=qDoubleBlast)  
                Bij = WaterfallMethod.NormaliseScores(Bij, Lengths, iSpecies, jSpecies)
                Bi.append(Bij)
            matrices.DumpMatrixArray("B", Bi, iSpecies, d_pickle)
            BH = GetBH_s(Bi, seqsInfo, iSpecies)
            matrices.DumpMatrixArray("BH", BH, iSpecies, d_pickle)
            util.PrintTime("Initial processing of species %d complete" % iSpecies)
        
    @staticmethod 
    def Worker_ProcessBlastHits(cmd_queue, d_pickle, qDoubleBlast):
        while True:
            try:
                args = cmd_queue.get(True, 1)
                WaterfallMethod.ProcessBlastHits(*args, d_pickle=d_pickle, qDoubleBlast=qDoubleBlast)
            except queue.Empty:
                return 
            except Exception:
                seqsInfo, _, _, iSpecies = args
                i = seqsInfo.speciesToUse[iSpecies]
                print("ERROR: Error processing files Blast%d_*" % i)
                raise

    @staticmethod
    def ConnectCognates(seqsInfo, iSpecies, d_pickle): 
        # calculate RBH for species i
        BHix = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle)
        BHxi = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle, row=False)
        RBHi = matrices.MatricesAndTr_s(BHix, BHxi)   # twice as much work as before (only did upper triangular before)
        del BHix, BHxi
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpecies, d_pickle)
        connect = WaterfallMethod.ConnectAllBetterThanAnOrtholog_s(RBHi, B, seqsInfo, iSpecies) 
        matrices.DumpMatrixArray("connect", connect, iSpecies, d_pickle)
            
    @staticmethod 
    def Worker_ConnectCognates(cmd_queue, d_pickle):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            while True:
                try:
                    args = cmd_queue.get(True, 1)
                    WaterfallMethod.ConnectCognates(*args, d_pickle=d_pickle)
                except queue.Empty:
                    return  
                                   
    @staticmethod
    def WriteGraphParallel(seqsInfo, nProcess):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            with open(files.FileHandler.GetGraphFilename(), 'w') as graphFile:
                graphFile.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n" % (seqsInfo.nSeqs, seqsInfo.nSeqs)) 
                graphFile.write("\n(mclmatrix\nbegin\n\n") 
            pool = mp.Pool(nProcess)
            graphFN = files.FileHandler.GetGraphFilename()
            pool.map(WriteGraph_perSpecies, [(seqsInfo, graphFN, iSpec, files.FileHandler.GetPickleDir()) for iSpec in range(seqsInfo.nSpecies)])
            for iSp in range(seqsInfo.nSpecies):
                subprocess.call("cat " + graphFN + "_%d" % iSp + " >> " + graphFN, shell=True)
                os.remove(graphFN + "_%d" % iSp)
            # Cleanup
            pool.close()
            matrices.DeleteMatrices("B", files.FileHandler.GetPickleDir()) 
            matrices.DeleteMatrices("connect", files.FileHandler.GetPickleDir()) 
    
    @staticmethod
    def GetMostDistant_s(RBH, B, seqsInfo, iSpec):
        # most distant RBB - as cut-off for connecting to other genes
        mostDistant = numeric.transpose(np.ones(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]])*1e9)
        # Best hit in another species: species-specific paralogues will now be connected - closer than any gene in any other species
        bestHit = numeric.transpose(np.zeros(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]))
        for kSpec in range(seqsInfo.nSpecies):
            B[kSpec] = B[kSpec].tocsr()
            if iSpec == kSpec:
                continue
            bestHit = np.maximum(bestHit, matrices.sparse_max_row(B[kSpec]))
            I, J = RBH[kSpec].nonzero()
            if len(I) > 0:
                mostDistant[I] = np.minimum(B[kSpec][I, J], mostDistant[I])
        # anything that doesn't have an RBB, set to distance to the closest gene in another species. I.e. it will hit just that gene and all genes closer to it in the its own species
        I = mostDistant > 1e8
        mostDistant[I] = bestHit[I] + 1e-6   # to connect to one in it's own species it must be closer than other species. We can deal with hits outside the species later 
        return mostDistant

    @staticmethod
    def ConnectAllBetterThanCutoff_s(B, mostDistant, seqsInfo, iSpec):
        connect = []
        nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]
        for jSpec in range(seqsInfo.nSpecies):
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
    with open(fn, csv_write_mode) as outfile:
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
    n = list(map(len, speciesPresence))
    writer_sum.writerow([])
    writer_sum.writerow(["Number of species in orthogroup", "Number of orthogroups"])
    for i in range(1, nSp+1):
        writer_sum.writerow([i, n.count(i)])

def Stats(ogs, speciesNamesDict, iSpecies, iResultsVersion):
    """ Top-level method for calculation of stats for the orthogroups"""
    allOgs = [[list(map(int, g.split("_"))) for g in og] for og in ogs]
    properOGs = [og for og in allOgs if len(og) > 1]
    allGenes = [g for og in allOgs for g in og]
    ogStatsResultsDir = files.FileHandler.GetOGsStatsResultsDirectory()
    filename_sp = ogStatsResultsDir +  "Statistics_PerSpecies" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".tsv"
    filename_sum = ogStatsResultsDir +  "Statistics_Overall" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".tsv"
    filename_overlap = ogStatsResultsDir +  "Orthogroups_SpeciesOverlaps" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".tsv"
    filename_single_copy = files.FileHandler.GetOrthogroupResultsFNBase() + "_SingleCopyOrthologues.txt"
    percentFormat = "%0.1f"
    with open(filename_sp, csv_write_mode) as outfile_species, open(filename_sum, csv_write_mode) as outfile_sum:
        writer_sp = csv.writer(outfile_species, delimiter="\t")
        writer_sum = csv.writer(outfile_sum, delimiter="\t")
        # header
        writer_sp.writerow([""] + [speciesNamesDict[index] for index in iSpecies])
        
        # Number of genes
        allGenesCounter = Counter([g[0] for g in allGenes])
        nGenes = sum(allGenesCounter.values())
        writer_sum.writerow(["Number of species", len(iSpecies)])
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
        writer_sp.writerow(["Percentage of orthogroups containing species"] + [percentFormat % ((100.*sum([iSp in og_sp for og_sp in speciesPresence])/len(properOGs)) if len(properOGs) > 0 else 0.) for iSp in iSpecies])
        
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
        l = list(reversed(list(map(len, properOGs))))
        writer_sum.writerow(["Mean orthogroup size", "%0.1f" % np.mean(l)])
        writer_sum.writerow(["Median orthogroup size", np.median(l)])
        L = np.cumsum(l)
        j, _ = next((i, x) for i, x in enumerate(L) if x > nAssigned/2)
        writer_sum.writerow(["G50 (assigned genes)",l[j]])
        l2 = list(reversed(list(map(len, ogs))))
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
        with open(filename_single_copy, 'w') as outfile_singlecopy:
            outfile_singlecopy.write("\n".join(["OG%07d" % i_ for i_ in singleCopyOGs]))
        # Link single-copy orthologues
        f =  files.FileHandler.GetOGsSeqFN
        in_fn = [f(i, True) for i in singleCopyOGs]
        g_fmt = files.FileHandler.GetResultsSeqsDir_SingleCopy() + files.FileHandler.baseOgFormat + ".fa"
        out_fn =[g_fmt % i for i in singleCopyOGs]
        for i, o in zip(in_fn, out_fn):
            shutil.copy(i, o)
            
        # Results filenames
        writer_sum.writerow(["Date", str(datetime.datetime.now()).split()[0]])
        writer_sum.writerow(["Orthogroups file", "Orthogroups" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".tsv"])
        writer_sum.writerow(["Unassigned genes file", "Orthogroups" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + "_UnassignedGenes.tsv"])
        writer_sum.writerow(["Per-species statistics", os.path.split(filename_sp)[1]])
        writer_sum.writerow(["Overall statistics", os.path.split(filename_sum)[1]])
        writer_sum.writerow(["Orthogroups shared between species", os.path.split(filename_overlap)[1]])
        
        # Sizes
        Stats_SizeTable(writer_sum, writer_sp, properOGs, allGenesCounter, iSpecies, speciesPresence)
        Stats_SpeciesOverlaps(filename_overlap, speciesNamesDict, iSpecies, speciesPresence)

    summaryText = """OrthoFinder assigned %d genes (%0.1f%% of total) to %d orthogroups. Fifty percent of all genes were in orthogroups with %d or more genes (G50 was %d) and were contained in the largest %d orthogroups (O50 was %d). There were %d orthogroups with all species present and %d of these consisted entirely of single-copy genes.""" % (nAssigned, pAssigned, nOgs, G50, G50, O50, O50, nCompleteOGs, nSingleCopy)
    print(summaryText)

"""
OrthoFinder
-------------------------------------------------------------------------------
"""   
g_mclInflation = 1.5

def CanRunBLAST():
    if parallel_task_manager.CanRunCommand("makeblastdb -help") and parallel_task_manager.CanRunCommand("blastp -help"):
        return True
    else:
        print("ERROR: Cannot run BLAST+")
        print("Please check BLAST+ is installed and that the executables are in the system path\n")
        return False

def CanRunMCL():
    command = "mcl -h"
    if parallel_task_manager.CanRunCommand(command):
        return True
    else:
        print("ERROR: Cannot run MCL with the command \"%s\"" % command)
        print("Please check MCL is installed and in the system path\n")
        return False
    
def GetProgramCaller():
    config_file = os.path.join(__location__, 'config.json') 
    pc = program_caller.ProgramCaller(config_file if os.path.exists(config_file) else None)
    config_file_user = os.path.expanduser("~/config_orthofinder_user.json")
    if os.path.exists(config_file_user):
        pc_user = program_caller.ProgramCaller(config_file_user)
        pc.Add(pc_user)
    return pc

def PrintHelp(prog_caller):  
    msa_ops = prog_caller.ListMSAMethods()
    tree_ops = prog_caller.ListTreeMethods()
    search_ops = prog_caller.ListSearchMethods()
    
    print("SIMPLE USAGE:") 
    print("Run full OrthoFinder analysis on FASTA format proteomes in <dir>")
    print("  orthofinder [options] -f <dir>")   
    print("")          
    print("Add new species in <dir1> to previous run in <dir2> and run new analysis")
    print("  orthofinder [options] -f <dir1> -b <dir2>")
    print("") 
      
    print("OPTIONS:")
    print(" -t <int>        Number of parallel sequence search threads [Default = %d]" % util.nThreadsDefault)
    print(" -a <int>        Number of parallel analysis threads")
    print(" -d              Input is DNA sequences")
    print(" -M <txt>        Method for gene tree inference. Options 'dendroblast' & 'msa'")
    print("                 [Default = dendroblast]")
    print(" -S <txt>        Sequence search program [Default = diamond]")
    print("                 Options: " + ", ".join(['blast'] + search_ops))
    print(" -A <txt>        MSA program, requires '-M msa' [Default = mafft]")
    print("                 Options: " + ", ".join(msa_ops))
    print(" -T <txt>        Tree inference method, requires '-M msa' [Default = fasttree]")
    print("                 Options: " + ", ".join(tree_ops)) 
#    print(" -R <txt>        Tree reconciliation method [Default = of_recon]")
#    print("                 Options: of_recon, dlcpar, dlcpar_convergedsearch")
    print(" -s <file>       User-specified rooted species tree")
    print(" -I <int>        MCL inflation parameter [Default = %0.1f]" % g_mclInflation)
    print(" -x <file>       Info for outputting results in OrthoXML format")
    print(" -p <dir>        Write the temporary pickle files to <dir>")
    print(" -1              Only perform one-way sequence search")
    print(" -X              Don't add species names to sequence IDs")
    print(" -y              Split paralogous clades below root of a HOG into separate HOGs")
    print(" -z              Don't trim MSAs (columns>=90% gap, min. alignment length 500)")
    print(" -n <txt>        Name to append to the results directory")  
    print(" -o <txt>        Non-default results directory")  
    print(" -h              Print this help text")

    print("")    
    print("WORKFLOW STOPPING OPTIONS:")   
    print(" -op             Stop after preparing input files for BLAST" )
    print(" -og             Stop after inferring orthogroups")
    print(" -os             Stop after writing sequence files for orthogroups")
    print("                 (requires '-M msa')")
    print(" -oa             Stop after inferring alignments for orthogroups")
    print("                 (requires '-M msa')")
    print(" -ot             Stop after inferring gene trees for orthogroups " )
   
    print("")   
    print("WORKFLOW RESTART COMMANDS:") 
    print(" -b  <dir>         Start OrthoFinder from pre-computed BLAST results in <dir>")   
    print(" -fg <dir>         Start OrthoFinder from pre-computed orthogroups in <dir>")
    print(" -ft <dir>         Start OrthoFinder from pre-computed gene trees in <dir>")
    
    print("")
    print("LICENSE:")
    print(" Distributed under the GNU General Public License (GPLv3). See License.md")
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
    if not os.path.exists(directory):
        print("Specified directory doesn't exist: %s" % directory)
        util.Fail()
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
        self.nProcessAlg = None
        self.qStartFromBlast = False  # remove, just store BLAST to do
        self.qStartFromFasta = False  # local to argument checking
        self.qStartFromGroups = False
        self.qStartFromTrees = False
        self.qStopAfterPrepare = False
        self.qStopAfterGroups = False
        self.qStopAfterSeqs = False
        self.qStopAfterAlignments = False
        self.qStopAfterTrees = False
        self.qMSATrees = False
        self.qAddSpeciesToIDs = True
        self.qTrim = True
        self.search_program = "diamond"
        self.msa_program = "mafft"
        self.tree_program = "fasttree"
        self.recon_method = "of_recon"
        self.name = None   # name to identify this set of results
        self.qDoubleBlast = True
        self.qSplitParaClades = False
        self.qPhyldog = False
        self.speciesXMLInfoFN = None
        self.speciesTreeFN = None
        self.mclInflation = g_mclInflation
        self.dna = False
    
    def what(self):
        for k, v in self.__dict__.items():
            if v == True:
                print(k)
                                 
def ProcessArgs(prog_caller, args):
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
    if len(args) == 0 or args[0] == "--help" or args[0] == "help" or args[0] == "-h":
        PrintHelp(prog_caller)
        util.Success() 

    options = Options()
    fastaDir = None
    continuationDir = None
    resultsDir_nonDefault = None
    pickleDir_nonDefault = None
    q_selected_msa_options = False
    q_selected_search_option = False
    
    """
    -f: store fastaDir
    -b: store workingDir
    -fg: store orthologuesDir 
    -ft: store orthologuesDir 
    + xml: speciesXMLInfoFN
    """    
    
    while len(args) > 0:
        arg = args.pop(0)    
        if arg == "-f" or arg == "--fasta":
            if options.qStartFromFasta:
                print("Repeated argument: -f/--fasta\n")
                util.Fail()
            options.qStartFromFasta = True
            fastaDir = GetDirectoryArgument(arg, args)
        elif arg == "-b" or arg == "--blast":
            if options.qStartFromBlast:
                print("Repeated argument: -b/--blast\n")
                util.Fail()
            options.qStartFromBlast = True
            continuationDir = GetDirectoryArgument(arg, args)
        elif arg == "-fg" or arg == "--from-groups":
            if options.qStartFromGroups:
                print("Repeated argument: -fg/--from-groups\n")
                util.Fail()
            options.qStartFromGroups = True
            continuationDir = GetDirectoryArgument(arg, args)
        elif arg == "-ft" or arg == "--from-trees":
            if options.qStartFromTrees:
                print("Repeated argument: -ft/--from-trees\n")
                util.Fail()
            options.qStartFromTrees = True
            continuationDir = GetDirectoryArgument(arg, args)
        elif arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nBlast = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s\n" % arg)
                util.Fail()    
        elif arg == "-a" or arg == "--algthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nProcessAlg = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s\n" % arg)
                util.Fail()   
        elif arg == "-1":
            options.qDoubleBlast = False
        elif arg == "-d" or arg == "--dna":
            options.dna = True
            if not q_selected_search_option:
                options.search_program = "blast_nucl"
        elif arg == "-X":
            options.qAddSpeciesToIDs = False
        elif arg == "-y":
            options.qSplitParaClades = True
        elif arg == "-z":
            options.qTrim = False
        elif arg == "-I" or arg == "--inflation":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.mclInflation = float(arg)
            except:
                print("Incorrect argument for MCL inflation parameter: %s\n" % arg)
                util.Fail()    
        elif arg == "-x" or arg == "--orthoxml":  
            if options.speciesXMLInfoFN:
                print("Repeated argument: -x/--orthoxml")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.speciesXMLInfoFN = args.pop(0)
        elif arg == "-n" or arg == "--name":  
            if options.name:
                print("Repeated argument: -n/--name")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.name = args.pop(0)
            while options.name.endswith("/"): options.name = options.name[:-1]
            if any([symbol in options.name for symbol in [" ", "/"]]): 
                print("Invalid symbol for command line argument %s\n" % arg)
                util.Fail()
        elif arg == "-o" or arg == "--output":  
            if resultsDir_nonDefault != None:
                print("Repeated argument: -o/--output")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            resultsDir_nonDefault = args.pop(0)
            while resultsDir_nonDefault.endswith("/"): resultsDir_nonDefault = resultsDir_nonDefault[:-1]
            resultsDir_nonDefault += "/"
            if os.path.exists(resultsDir_nonDefault):
                print("ERROR: non-default output directory already exists: %s\n" % resultsDir_nonDefault)
                util.Fail()
            if " " in resultsDir_nonDefault:
                print("ERROR: non-default output directory cannot include spaces: %s\n" % resultsDir_nonDefault)
                util.Fail()
            checkDirName = resultsDir_nonDefault
            while checkDirName.endswith("/"):
                checkDirName = checkDirName[:-1]
            path, newDir = os.path.split(checkDirName)
            if path != "" and not os.path.exists(path):
                print("ERROR: location '%s' for results directory '%s' does not exist.\n" % (path, newDir))
                util.Fail()
        elif arg == "-s" or arg == "--speciestree":  
            if options.speciesXMLInfoFN:
                print("Repeated argument: -s/--speciestree")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.speciesTreeFN = args.pop(0)
        elif arg == "-S" or arg == "--search":
            choices = ['blast'] + prog_caller.ListSearchMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.search_program = arg
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()
        elif arg == "-M" or arg == "--method":
            arg_M_or_msa = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg == "msa": 
                options.qMSATrees = True
            elif arg == "phyldog": 
                options.qPhyldog = True
                options.recon_method = "phyldog"
                options.qMSATrees = False
            elif arg == "dendroblast": options.qMSATrees = False    
            else:
                print("Invalid argument for option %s: %s" % (arg_M_or_msa, arg))
                print("Valid options are 'dendroblast' and 'msa'\n")
                util.Fail()
        elif arg == "-A" or arg == "--msa_program":
            choices = ['mafft'] + prog_caller.ListMSAMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.msa_program = arg
                q_selected_msa_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()
        elif arg == "-T" or arg == "--tree_program":
            choices = ['fasttree'] + prog_caller.ListTreeMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.tree_program = arg
                q_selected_msa_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()
        elif arg == "-R" or arg == "--recon_method":
            choices = ['of_recon', 'dlcpar', 'dlcpar_convergedsearch', 'only_overlap']
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.recon_method = arg
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()
        elif arg == "-p":
            pickleDir_nonDefault = GetDirectoryArgument(arg, args)
        elif arg == "-op" or arg == "--only-prepare":
            options.qStopAfterPrepare = True
        elif arg == "-og" or arg == "--only-groups":
            options.qStopAfterGroups = True
        elif arg == "-os" or arg == "--only-seqs":
            options.qStopAfterSeqs = True
        elif arg == "-oa" or arg == "--only-alignments":
            options.qStopAfterAlignments = True
        elif arg == "-ot" or arg == "--only-trees":
            options.qStopAfterTrees = True
        elif arg == "-h" or arg == "--help":
            PrintHelp(prog_caller)
            util.Success()
        else:
            print("Unrecognised argument: %s\n" % arg)
            util.Fail()    
    
    # set a default for number of algorithm threads
    if options.nProcessAlg is None:
        options.nProcessAlg = min(16, max(1, int(options.nBlast/8)))

    # check argument combinations       
    if not (options.qStartFromFasta or options.qStartFromBlast or options.qStartFromGroups or options.qStartFromTrees):
        print("ERROR: Please specify the input directory for OrthoFinder using one of the options: '-f', '-b', '-fg' or '-ft'.")
        util.Fail()
    
    if options.qStartFromFasta and (options.qStartFromTrees or options.qStartFromGroups):
        print("ERROR: Incompatible arguments, -f (start from fasta files) and" + (" -fg (start from orthogroups)" if options.qStartFromGroups else " -ft (start from trees)"))
        util.Fail()
        
    if options.qStartFromBlast and (options.qStartFromTrees or options.qStartFromGroups):
        print("ERROR: Incompatible arguments, -b (start from pre-calcualted BLAST results) and" + (" -fg (start from orthogroups)" if options.qStartFromGroups else " -ft (start from trees)"))
        util.Fail()      

    if options.qStartFromTrees and options.qStartFromGroups:
        print("ERROR: Incompatible arguments, -fg (start from orthogroups) and -ft (start from trees)")
        util.Fail()    

    if options.qStopAfterSeqs and (not options.qMSATrees):
        print("ERROR: Argument '-os' (stop after sequences) also requires option '-M msa'")
        util.Fail()   

    if options.qStopAfterAlignments and (not options.qMSATrees):
        print("ERROR: Argument '-oa' (stop after alignments) also requires option '-M msa'")
        util.Fail()     

    if q_selected_msa_options and (not options.qMSATrees and not options.qPhyldog):
        print("ERROR: Argument '-A' or '-T' (multiple sequence alignment/tree inference program) also requires option '-M msa'")
        util.Fail()       
        
    if options.qPhyldog and (not options.speciesTreeFN):
        print("ERROR: Phyldog currently needs a species tree to be provided")
        util.Fail()          

    if resultsDir_nonDefault != None and ((not options.qStartFromFasta) or options.qStartFromBlast):
        print("ERROR: Incompatible arguments, -o (non-default output directory) can only be used with a new OrthoFinder run using option '-f'")
        util.Fail()       
        
    if options.search_program not in (prog_caller.ListSearchMethods() + ['blast']):
        print("ERROR: Search program (%s) not configured in config.json file" % options.search_program)
        util.Fail()
        
    util.PrintTime("Starting OrthoFinder %s" % util.version)    
    print("%d thread(s) for highly parallel tasks (BLAST searches etc.)" % options.nBlast)
    print("%d thread(s) for OrthoFinder algorithm" % options.nProcessAlg)
    return options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault            

def GetXMLSpeciesInfo(seqsInfoObj, options):
    # speciesInfo:  name, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
    util.PrintUnderline("Reading species information file")
    # do this now so that we can alert user to any errors prior to running the algorithm
    speciesXML = [[] for i_ in seqsInfoObj.speciesToUse]
    speciesNamesDict = SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
    speciesRevDict = {v:k for k,v in speciesNamesDict.items()}
    userFastaFilenames = [os.path.split(speciesNamesDict[i])[1] for i in seqsInfoObj.speciesToUse]
    with open(options.speciesXMLInfoFN, 'r') as speciesInfoFile:
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
            speciesXML[seqsInfoObj.speciesToUse.index(iSpecies)] = line   
    # check information has been provided for all species
    speciesMissing = False        
    for iPos, iSpecies in enumerate(seqsInfoObj.speciesToUse):
        if speciesXML[iPos] == []:
            if not speciesMissing:
                print("ERROR")
                print("Species information file %s does not contain information for all species." % options.speciesXMLInfoFN)
                print("Information is missing for:") 
                speciesMissing = True
            print(speciesNamesDict[iSpecies])
    if speciesMissing:
        util.Fail()
    return speciesXML

def IDsFileOK(filename):
    """
    It is best to detect any issues with input files at start, perform all required checks here
    """
    with open(filename, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if len(line) == 0: continue
            tokens = line.split(": ", 1)
            if len(tokens) !=2 or len(tokens[1]) == 0:
                return False, line
    return True, None

def CheckDependencies(options, prog_caller, dirForTempFiles):
    util.PrintUnderline("Checking required programs are installed")
    if (options.qStartFromFasta):
        if options.search_program == "blast":
            if not CanRunBLAST(): util.Fail()
        elif not prog_caller.TestSearchMethod(dirForTempFiles, options.search_program):
            print("\nERROR: Cannot run %s" % options.search_program)
            print("Format of make database command:")
            print("  " + prog_caller.GetSearchMethodCommand_DB(options.search_program, "INPUT", "OUTPUT"))
            print("ERROR: Cannot run %s" % options.search_program)
            print("Format of search database command:")
            print("  " + prog_caller.GetSearchMethodCommand_Search(options.search_program, "INPUT", "DATABASE", "OUTPUT"))
            print("Please check %s is installed and that the executables are in the system path\n" % options.search_program)
            util.Fail()
    if (options.qStartFromFasta or options.qStartFromBlast) and not CanRunMCL():
        util.Fail()
    if not (options.qStopAfterPrepare or options.qStopAfterSeqs or options.qStopAfterGroups):
        if not orthologues.CanRunOrthologueDependencies(dirForTempFiles, 
                                                            options.qMSATrees, 
                                                            options.qPhyldog, 
                                                            options.qStopAfterTrees, 
                                                            options.msa_program, 
                                                            options.tree_program, 
                                                            options.recon_method,
                                                            prog_caller, 
                                                            options.qStopAfterAlignments):
            print("Dependencies have been met for inference of orthogroups but not for the subsequent orthologue inference.")
            print("Either install the required dependencies or use the option '-og' to stop the analysis after the inference of orthogroups.\n")
            util.Fail()

def DoOrthogroups(options, speciesInfoObj, seqsInfo):
    # Run Algorithm, cluster and output cluster files with original accessions
    util.PrintUnderline("Running OrthoFinder algorithm")
    # it's important to free up the memory from python used for processing the genomes
    # before launching MCL becuase both use sizeable ammounts of memory. The only
    # way I can find to do this is to launch the memory intensive python code 
    # as separate process that exits before MCL is launched.
    Lengths = GetSequenceLengths(seqsInfo)
    
    # Process BLAST hits
    util.PrintTime("Initial processing of each species")
    cmd_queue = mp.Queue()
    blastDir_list = files.FileHandler.GetBlastResultsDir()
    for iSpecies in range(seqsInfo.nSpecies):
        cmd_queue.put((seqsInfo, blastDir_list, Lengths, iSpecies))
    files.FileHandler.GetPickleDir()     # create the pickle directory before the parallel processing to prevent a race condition
    runningProcesses = [mp.Process(target=WaterfallMethod.Worker_ProcessBlastHits, args=(cmd_queue, files.FileHandler.GetPickleDir(), options.qDoubleBlast)) for i_ in range(options.nProcessAlg)]
    for proc in runningProcesses:
        proc.start()
    parallel_task_manager.ManageQueue(runningProcesses, cmd_queue)
    
    cmd_queue = mp.Queue()
    for iSpecies in range(seqsInfo.nSpecies):
        cmd_queue.put((seqsInfo, iSpecies))
    runningProcesses = [mp.Process(target=WaterfallMethod.Worker_ConnectCognates, args=(cmd_queue, files.FileHandler.GetPickleDir())) for i_ in range(options.nProcessAlg)]
    for proc in runningProcesses:
        proc.start()
    parallel_task_manager.ManageQueue(runningProcesses, cmd_queue)
    
    util.PrintTime("Connected putative homologues") 
    WaterfallMethod.WriteGraphParallel(seqsInfo, options.nProcessAlg)
    
    # 5b. MCL     
    clustersFilename, clustersFilename_pairs = files.FileHandler.CreateUnusedClustersFN(options.mclInflation) 
    graphFilename = files.FileHandler.GetGraphFilename() 
    MCL.RunMCL(graphFilename, clustersFilename, options.nProcessAlg, options.mclInflation)
    mcl.ConvertSingleIDsToIDPair(seqsInfo, clustersFilename, clustersFilename_pairs)   
    
    util.PrintUnderline("Writing orthogroups to file")
    ogs = mcl.GetPredictedOGs(clustersFilename_pairs)
    
    resultsBaseFilename = files.FileHandler.GetOrthogroupResultsFNBase()
    idsDict = MCL.WriteOrthogroupFiles(ogs, [files.FileHandler.GetSequenceIDsFN()], resultsBaseFilename, clustersFilename_pairs)
    speciesNamesDict = SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
    MCL.CreateOrthogroupTable(ogs, idsDict, speciesNamesDict, speciesInfoObj.speciesToUse, resultsBaseFilename)
    
    # Write Orthogroup FASTA files    
    ogSet = orthologues.OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll, options.qAddSpeciesToIDs, idExtractor = util.FirstWordExtractor)
    treeGen = trees_msa.TreesForOrthogroups(None, None, None)
    fastaWriter = trees_msa.FastaWriter(files.FileHandler.GetSpeciesSeqsDir(), speciesInfoObj.speciesToUse)
    d_seqs = files.FileHandler.GetResultsSeqsDir()
    if not os.path.exists(d_seqs): os.mkdir(d_seqs)
    treeGen.WriteFastaFiles(fastaWriter, ogSet.OGs(qInclAll=True), idsDict, False)
    
    Stats(ogs, speciesNamesDict, speciesInfoObj.speciesToUse, files.FileHandler.iResultsVersion)
    if options.speciesXMLInfoFN:
        MCL.WriteOrthoXML(speciesXML, ogs, seqsInfo.nSeqsPerSpecies, idsDict, resultsBaseFilename + ".orthoxml", speciesInfoObj.speciesToUse)
    print("")
    util.PrintTime("Done orthogroups")
    files.FileHandler.LogOGs()

# 0
def ProcessPreviousFiles(workingDir_list, qDoubleBlast):
    """Checks for:
    workingDir should be the WorkingDirectory containing Blast*.txt files
    
    SpeciesIDs.txt
    Species*.fa
    Blast*txt
    SequenceIDs.txt
    
    Checks which species should be included
    
    """
    # check BLAST results directory exists
    if not os.path.exists(workingDir_list[0]):
        err_text = "ERROR: Previous/Pre-calculated BLAST results directory does not exist: %s\n" % workingDir_list[0]
        files.FileHandler.LogFailAndExit(err_text)
        
    speciesInfo = files.SpeciesInfo()
    if not os.path.exists(files.FileHandler.GetSpeciesIDsFN()):
        err_text = "ERROR: %s file must be provided if using previously calculated BLAST results" % files.FileHandler.GetSpeciesIDsFN()
        files.FileHandler.LogFailAndExit(err_text)
    file_ok, err_line = IDsFileOK(files.FileHandler.GetSpeciesIDsFN())
    if not file_ok: 
        files.FileHandler.LogFailAndExit("ERROR: %s file contains a blank accession. Line:\n %s" % (files.FileHandler.GetSpeciesIDsFN(), err_line))
    speciesInfo.speciesToUse, speciesInfo.nSpAll, speciesToUse_names = util.GetSpeciesToUse(files.FileHandler.GetSpeciesIDsFN())
 
    # check fasta files are present 
    previousFastaFiles = files.FileHandler.GetSortedSpeciesFastaFiles()
    if len(previousFastaFiles) == 0:
        err_text = "ERROR: No processed fasta files in the supplied previous working directories:\n" + "\n".join(workingDir_list) + "\n"
        files.FileHandler.LogFailAndExit(err_text)
    tokens = previousFastaFiles[-1][:-3].split("Species")
    lastFastaNumberString = tokens[-1]
    iLastFasta = 0
    nFasta = len(previousFastaFiles)
    try:
        iLastFasta = int(lastFastaNumberString)
    except:
        files.FileHandler.LogFailAndExit("ERROR: Filenames for processed fasta files are incorrect: %s\n" % previousFastaFiles[-1])
    if nFasta != iLastFasta + 1:
        files.FileHandler.LogFailAndExit("ERROR: Not all expected fasta files are present. Index of last fasta file is %s but found %d fasta files.\n" % (lastFastaNumberString, len(previousFastaFiles)))
    
    # check BLAST files
    blast_fns_triangular = [files.FileHandler.GetBlastResultsFN(iSpecies, jSpecies) for iSpecies in speciesInfo.speciesToUse for jSpecies in speciesInfo.speciesToUse if jSpecies >= iSpecies]
    have_triangular = [(os.path.exists(fn) or os.path.exists(fn + ".gz")) for fn in blast_fns_triangular]
    for qHave, fn in zip(have_triangular, blast_fns_triangular):
        if not qHave: print("BLAST results file is missing: %s" % fn)
    
    if qDoubleBlast:
        blast_fns_remainder = [files.FileHandler.GetBlastResultsFN(iSpecies, jSpecies) for iSpecies in speciesInfo.speciesToUse for jSpecies in speciesInfo.speciesToUse if jSpecies < iSpecies]
        have_remainder = [(os.path.exists(fn) or os.path.exists(fn + ".gz")) for fn in blast_fns_remainder]
        if not (all(have_triangular) and all(have_remainder)):
            for qHave, fn in zip(have_remainder, blast_fns_remainder):
                if not qHave: print("BLAST results file is missing: %s" % fn)
            if not all(have_triangular):
                files.FileHandler.LogFailAndExit()
            else:
                # would be able to do it using just one-way blast
                files.FileHandler.LogFailAndExit("ERROR: Required BLAST results files are present for using the one-way sequence search option (default) but not the double BLAST search ('-d' option)")
    else:
        if not all(have_triangular):
            files.FileHandler.LogFailAndExit()
                            
    # check SequenceIDs.txt and SpeciesIDs.txt files are present
    if not os.path.exists(files.FileHandler.GetSequenceIDsFN()):
        files.FileHandler.LogFailAndExit("ERROR: %s file must be provided if using previous calculated BLAST results" % files.FileHandler.GetSequenceIDsFN())
    
    file_ok, err_line = IDsFileOK(files.FileHandler.GetSequenceIDsFN())
    if not file_ok: 
        files.FileHandler.LogFailAndExit("ERROR: %s file contains a blank accession. Line:\n %s" % (files.FileHandler.GetSequenceIDsFN(), err_line))
    return speciesInfo, speciesToUse_names

# 6
def CreateSearchDatabases(seqsInfoObj, options, prog_caller):
    nDB = max(seqsInfoObj.speciesToUse) + 1
    for iSp in range(nDB):
        if options.search_program == "blast":
            command = " ".join(["makeblastdb", "-dbtype", "prot", "-in", files.FileHandler.GetSpeciesFastaFN(iSp), "-out", files.FileHandler.GetSpeciesDatabaseN(iSp)])
            util.PrintTime("Creating Blast database %d of %d" % (iSp + 1, nDB))
            RunBlastDBCommand(command) 
        else:
            command = prog_caller.GetSearchMethodCommand_DB(options.search_program, files.FileHandler.GetSpeciesFastaFN(iSp), files.FileHandler.GetSpeciesDatabaseN(iSp, options.search_program))
            util.PrintTime("Creating %s database %d of %d" % (options.search_program, iSp + 1, nDB))
            ret_code = parallel_task_manager.RunCommand(command, qPrintOnError=True, qPrintStderr=False)
            if ret_code != 0:
                files.FileHandler.LogFailAndExit("ERROR: diamond makedb failed")

# 7
def RunSearch(options, speciessInfoObj, seqsInfo, prog_caller):
    name_to_print = "BLAST" if options.search_program == "blast" else options.search_program
    if options.qStopAfterPrepare:
        util.PrintUnderline("%s commands that must be run" % name_to_print)
    else:        
        util.PrintUnderline("Running %s all-versus-all" % name_to_print)
    commands = GetOrderedSearchCommands(seqsInfo, speciessInfoObj, options.qDoubleBlast, options.search_program, prog_caller)
    if options.qStopAfterPrepare:
        for command in commands:
            print(command)
        util.Success()
    print("Using %d thread(s)" % options.nBlast)
    util.PrintTime("This may take some time....")  
    cmd_queue = mp.Queue()
    for iCmd, cmd in enumerate(commands):
        cmd_queue.put((iCmd+1, cmd))           
    runningProcesses = [mp.Process(target=parallel_task_manager.Worker_RunCommand, args=(cmd_queue, options.nBlast, len(commands), True)) for i_ in range(options.nBlast)]
    for proc in runningProcesses:
        proc.start()#
    for proc in runningProcesses:
        while proc.is_alive():
            proc.join()
    # remove BLAST databases
    util.PrintTime("Done all-versus-all sequence search")
    if options.search_program == "blast":
        for f in glob.glob(files.FileHandler.GetWorkingDirectory1_Read()[0] + "BlastDBSpecies*"):
            os.remove(f)
    if options.search_program == "mmseqs":
        for i in range(speciessInfoObj.nSpAll):
            for j in range(speciessInfoObj.nSpAll):
                tmp_dir = "/tmp/tmpBlast%d_%d.txt" % (i,j)
                if os.path.exists(tmp_dir):
                    try:
                        shutil.rmtree(tmp_dir)
                    except OSError:
                        time.sleep(1)
                        shutil.rmtree(tmp_dir, True)  # shutil / NFS bug - ignore errors, it's less crucial that the files are deleted

# 9
def GetOrthologues(speciesInfoObj, options, prog_caller):
    util.PrintUnderline("Analysing Orthogroups", True)

    orthologues.OrthologuesWorkflow(speciesInfoObj.speciesToUse, 
                                    speciesInfoObj.nSpAll, 
                                    prog_caller,
                                    options.msa_program,
                                    options.tree_program,
                                    options.recon_method,
                                    options.nBlast,
                                    options.nProcessAlg,
                                    options.qDoubleBlast,
                                    options.qAddSpeciesToIDs,
                                    options.qTrim,
                                    options.speciesTreeFN, 
                                    options.qStopAfterSeqs,
                                    options.qStopAfterAlignments,
                                    options.qStopAfterTrees,
                                    options.qMSATrees,
                                    options.qPhyldog,
                                    options.name,
                                    options.qSplitParaClades)
    util.PrintTime("Done orthologues")

def GetOrthologues_FromTrees(options):
    orthologues.OrthologuesFromTrees(options.recon_method, options.nBlast, options.nProcessAlg, options.speciesTreeFN, options.qAddSpeciesToIDs, options.qSplitParaClades)
 
def ProcessesNewFasta(fastaDir, q_dna, speciesInfoObj_prev = None, speciesToUse_prev_names=[]):
    """
    Process fasta files and return a Directory object with all paths completed.
    """
    # Check files present
    qOk = True
    if not os.path.exists(fastaDir):
        print("\nDirectory does not exist: %s" % fastaDir)
        util.Fail()
    files_in_directory = sorted([f for f in os.listdir(fastaDir) if os.path.isfile(os.path.join(fastaDir,f))])
    originalFastaFilenames = []
    excludedFiles = []
    for f in files_in_directory:
        if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in fastaExtensions and not f.startswith("._"):
            originalFastaFilenames.append(f)
        else:
            excludedFiles.append(f)
    if len(excludedFiles) != 0:
        print("\nWARNING: Files have been ignored as they don't appear to be FASTA files:")
        for f in excludedFiles:
            print(f)
        print("OrthoFinder expects FASTA files to have one of the following extensions: %s" % (", ".join(fastaExtensions)))
    speciesToUse_prev_names = set(speciesToUse_prev_names)
    if len(originalFastaFilenames) + len(speciesToUse_prev_names) < 2:
        print("ERROR: At least two species are required")
        util.Fail()
    if any([fn in speciesToUse_prev_names for fn in originalFastaFilenames]):
        print("ERROR: Attempted to add a second copy of a previously included species:")
        for fn in originalFastaFilenames:
            if fn in speciesToUse_prev_names: print(fn)
        print("")
        util.Fail()
    if len(originalFastaFilenames) == 0:
        print("\nNo fasta files found in supplied directory: %s" % fastaDir)
        util.Fail()
    if speciesInfoObj_prev == None:
        # Then this is a new, clean analysis 
        speciesInfoObj = files.SpeciesInfo()
    else:
        speciesInfoObj = speciesInfoObj_prev
    iSeq = 0
    iSpecies = 0
    # If it's a previous analysis:
    if len(speciesToUse_prev_names) != 0:
        with open(files.FileHandler.GetSpeciesIDsFN(), 'r') as infile:
            for line in infile: pass
        if line.startswith("#"): line = line[1:]
        iSpecies = int(line.split(":")[0]) + 1
    speciesInfoObj.iFirstNewSpecies = iSpecies
    newSpeciesIDs = []
    with open(files.FileHandler.GetSequenceIDsFN(), 'a') as idsFile, open(files.FileHandler.GetSpeciesIDsFN(), 'a') as speciesFile:
        for fastaFilename in originalFastaFilenames:
            newSpeciesIDs.append(iSpecies)
            outputFasta = open(files.FileHandler.GetSpeciesFastaFN(iSpecies, qForCreation=True), 'w')
            fastaFilename = fastaFilename.rstrip()
            speciesFile.write("%d: %s\n" % (iSpecies, fastaFilename))
            baseFilename, extension = os.path.splitext(fastaFilename)
            mLinesToCheck = 100
            qHasAA = False
            with open(fastaDir + os.sep + fastaFilename, 'r') as fastaFile:
                for iLine, line in enumerate(fastaFile):
                    if line.isspace(): continue
                    if len(line) > 0 and line[0] == ">":
                        newID = "%d_%d" % (iSpecies, iSeq)
                        acc = line[1:].rstrip()
                        if len(acc) == 0:
                            print("ERROR: %s contains a blank accession line on line %d" % (fastaDir + os.sep + fastaFilename, iLine+1))
                            util.Fail()
                        idsFile.write("%s: %s\n" % (newID, acc))
                        outputFasta.write(">%s\n" % newID)    
                        iSeq += 1
                    else:
                        line = line.upper()    # allow lowercase letters in sequences
                        if not qHasAA and (iLine < mLinesToCheck):
#                            qHasAA = qHasAA or any([c in line for c in ['D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y']])
                            qHasAA = qHasAA or any([c in line for c in ['E','F','I','L','P','Q']]) # AAs minus nucleotide ambiguity codes
                        outputFasta.write(line)
                outputFasta.write("\n")
            if (not qHasAA) and (not q_dna):
                qOk = False
                print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences. Use '-d' option" % fastaFilename)
            iSpecies += 1
            iSeq = 0
            outputFasta.close()
        if not qOk:
            util.Fail()
    if len(originalFastaFilenames) > 0: outputFasta.close()
    speciesInfoObj.speciesToUse = speciesInfoObj.speciesToUse + newSpeciesIDs
    speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1      # will be one of the new species
    return speciesInfoObj

def DeleteDirectoryTree(d):
    if os.path.exists(d): 
        try:
            shutil.rmtree(d)
        except OSError:
            time.sleep(1)
            shutil.rmtree(d, True)   

def CheckOptions(options, speciesToUse):
    """Check any optional arguments are valid once we know what species are in the analysis
    - user supplied species tree
    """
    if options.speciesTreeFN:
        expSpecies = list(SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN()).values())
        orthologues.CheckUserSpeciesTree(options.speciesTreeFN, expSpecies)
        
    if options.qStopAfterSeqs and (not options.qMSATrees):
        print("ERROR: Must use '-M msa' option to generate sequence files for orthogroups")
        util.Fail()
    if options.qStopAfterAlignments and (not options.qMSATrees):
        print("ERROR: Must use '-M msa' option to generate sequence files and infer multiple sequence alignments for orthogroups")
        util.Fail()

    # check can open enough files
    n_extra = 50
    q_do_orthologs = not any((options.qStopAfterPrepare, options.qStopAfterGroups, options.qStopAfterSeqs, options.qStopAfterAlignments, options.qStopAfterTrees))
    if q_do_orthologs and not options.qStartFromTrees:
        n_sp = len(speciesToUse)
        wd = files.FileHandler.GetWorkingDirectory_Write()
        wd_files_test = wd + "Files_test/"
        fh = []
        try:
            if not os.path.exists(wd_files_test):
                os.mkdir(wd_files_test)
            for i_sp in range(n_sp):
                di = wd_files_test + "Sp%d/" % i_sp
                if not os.path.exists(di):
                    os.mkdir(di)
                for j_sp in range(n_sp):
                    fnij = di + "Sp%d.txt" % j_sp
                    fh.append(open(fnij, 'w'))
            # create a few extra files to be safe
            for i_extra in range(n_extra):
                fh.append(open(wd_files_test + "Extra%d.txt" % i_extra, 'w'))
            # close the files again and delete
            for fhh in fh:
                fhh.close()
            DeleteDirectoryTree(wd_files_test)
        except IOError as e:
            if str(e).startswith("[Errno 24] Too many open files"):
                util.number_open_files_exception_advice(len(speciesToUse), False)
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                util.Fail()
            else:
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                print("ERROR: Attempted to open required files for OrthoFinder run but an unexpected error occurred. \n\nStacktrace:")
                raise
    return options

def main(args=None):    
    try:
        if args is None:
            args = sys.argv[1:]
        # Create PTM right at start
        ptm_initialised = parallel_task_manager.ParallelTaskManager_singleton()
        print("")
        print(("OrthoFinder version %s Copyright (C) 2014 David Emms\n" % util.version))
        prog_caller = GetProgramCaller()
        
        options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault = ProcessArgs(prog_caller, args)  
        
        files.InitialiseFileHandler(options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault)     
                    
        CheckDependencies(options, prog_caller, files.FileHandler.GetWorkingDirectory1_Read()[0]) 
            
        # if using previous Trees etc., check these are all present - Job for orthologues
        if options.qStartFromBlast and options.qStartFromFasta:
            # 0. Check Files
            speciesInfoObj, speciesToUse_names = ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast)
            print("\nAdding new species in %s to existing analysis in %s" % (fastaDir, continuationDir))
            # 3. 
            speciesInfoObj = ProcessesNewFasta(fastaDir, options.dna, speciesInfoObj, speciesToUse_names)
            files.FileHandler.LogSpecies()
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4.
            seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
            # 5.
            if options.speciesXMLInfoFN:   
                speciesXML = GetXMLSpeciesInfo(speciesInfoObj, options)
            # 6.    
            util.PrintUnderline("Dividing up work for BLAST for parallel processing")
            CreateSearchDatabases(speciesInfoObj, options, prog_caller)
            # 7.  
            RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
            # 8.
            DoOrthogroups(options, speciesInfoObj, seqsInfo)
            # 9.
            if not options.qStopAfterGroups:
                GetOrthologues(speciesInfoObj, options, prog_caller)   
        elif options.qStartFromFasta:
            # 3. 
            speciesInfoObj = None
            speciesInfoObj = ProcessesNewFasta(fastaDir, options.dna)
            files.FileHandler.LogSpecies()
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4
            seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
            # 5.
            if options.speciesXMLInfoFN:   
                speciesXML = GetXMLSpeciesInfo(speciesInfoObj, options)
            # 6.    
            util.PrintUnderline("Dividing up work for BLAST for parallel processing")
            CreateSearchDatabases(speciesInfoObj, options, prog_caller)
            # 7. 
            RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
            # 8.  
            DoOrthogroups(options, speciesInfoObj, seqsInfo)    
            # 9. 
            if not options.qStopAfterGroups:
                GetOrthologues(speciesInfoObj, options, prog_caller)
        elif options.qStartFromBlast:
            # 0.
            speciesInfoObj, _ = ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast)
            files.FileHandler.LogSpecies()
            print("Using previously calculated BLAST results in %s" % (files.FileHandler.GetWorkingDirectory1_Read()[0]))
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4.
            seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
            # 5.
            if options.speciesXMLInfoFN:   
                speciesXML = GetXMLSpeciesInfo(speciesInfoObj, options)
            # 8        
            DoOrthogroups(options, speciesInfoObj, seqsInfo)    
            # 9
            if not options.qStopAfterGroups:
                GetOrthologues(speciesInfoObj, options, prog_caller)
        elif options.qStartFromGroups:
            # 0.  
            speciesInfoObj, _ = ProcessPreviousFiles(continuationDir, options.qDoubleBlast)
            files.FileHandler.LogSpecies()
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            # 9
            GetOrthologues(speciesInfoObj, options, prog_caller)
        elif options.qStartFromTrees:
            speciesInfoObj, _ = ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast)
            files.FileHandler.LogSpecies()
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            GetOrthologues_FromTrees(options)
        else:
            raise NotImplementedError
            ptm = parallel_task_manager.ParallelTaskManager_singleton()
            ptm.Stop()
        d_results = os.path.normpath(files.FileHandler.GetResultsDirectory1()) + os.path.sep
        print("\nResults:\n    %s" % d_results)
        util.PrintCitation(d_results)
        files.FileHandler.WriteToLog("OrthoFinder run completed\n", True)
    except Exception as e:
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()
        raise
    ptm = parallel_task_manager.ParallelTaskManager_singleton()
    ptm.Stop()
        
        
     

    

    

