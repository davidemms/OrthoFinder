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

def GetPredictedOGs(clustersFilename):
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
    if not qContainsProfiles:
        assert(len(predictedOGs) == int(nOGsString) + 1)
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

def ConvertSingleIDsToIDPair(seqsInfo, clustersFilename, newFilename):
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
                    
