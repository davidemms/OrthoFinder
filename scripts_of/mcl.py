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
from __future__ import absolute_import
import sys
import csv
from typing import List, Set

from . import parallel_task_manager

from collections import defaultdict, Counter
import xml.etree.ElementTree as ET              # Y
from xml.etree.ElementTree import SubElement    # Y
from xml.dom import minidom

from . import util, files


def GetPredictedOGs(clustersFilename):
    """
    Returns:
        ogs: List[Set[str]] - List of sets of genes ids (as strings)
    """
    predictedOGs = []
    nOGsString = ""
    qContainsProfiles = False
    with open(clustersFilename, 'r') as clusterFile:
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
    return predictedOGs
    
def GetSingleID(speciesStartingIndices, seq, speciesToUse): 
    a, b = seq.split("_")
    iSpecies = int(a)
    iSeq = int(b)
    offset = speciesStartingIndices[speciesToUse.index(iSpecies)]
    return iSeq + offset  

def GetIDPair(speciesStartingIndices, singleID, speciesToUse):   
    for i, startingIndex in enumerate(speciesStartingIndices):
        if startingIndex > singleID:
            return "%d_%d" % (speciesToUse[i-1], singleID - speciesStartingIndices[i-1])
    return "%d_%d" % (speciesToUse[-1], singleID - speciesStartingIndices[len(speciesStartingIndices)-1]) 

def ConvertSingleIDsToIDPair(seqsInfo, clustersFilename, newFilename, q_unassigned=False):
    with open(clustersFilename, 'r') as clusterFile, open(newFilename, "w") as output:
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
                if q_unassigned and not line.startswith(" ") and len(line.split()) == 3:
                    # This is a single gene, not necessarily from this MCL clustering but here implicitly
                    continue
                if line[-2] == "$":
                    line = line[:-3]
                    appendDollar = True
                if line[0] != " ":
                    initialText, line = line.split(None, 1)
                # continuation of group
                ids = line.split()
                for id in ids:
                    idsString += GetIDPair(seqsInfo.seqStartingIndices, int(id), seqsInfo.speciesToUse) + " "
                output.write(initialText + "      " + idsString)
                if appendDollar:
                    output.write("$\n")
                else:
                    output.write("\n")


def write_updated_clusters_file(ogs: List[Set[str]], clustersFilename_pairs):
    with open(clustersFilename_pairs, 'w') as outfile:
        outfile.write("(mclmatrix")
        outfile.write("begin\n")
        for iog, og in enumerate(ogs):
            outfile.write("%d      %s $\n" % (iog, " ".join(sorted(og))))
        outfile.write(")\n")


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
        with open(outputFilename, util.csv_write_mode) as outputFile, open(singleGeneFilename, util.csv_write_mode) as singleGeneFile, open(outputFilename_counts, util.csv_write_mode) as outFile_counts:
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
