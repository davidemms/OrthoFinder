#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import os
import glob

#sys.path.append(os.path.split(__file__)[0] + os.sep + "..")

if __name__ == "__main__" and __package__ is None:   
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#    __package__ = "orthofinder"

#import orthofinder.scripts.tree
from scripts import tree

def ReplaceStringWithNewIDs(idMap, treeString):
    # Load species names
    sequences = idMap
    start = 0
    regex = re.compile('\d+_\d+')
    m = regex.search(treeString)
    while m:
        id = treeString[m.start():m.end()]
        newText = treeString[:m.start()] + sequences[id]
        start = len(newText)
        treeString = newText + treeString[m.end():]
        m = regex.search(treeString, start)
    return treeString
    
def ReplaceString(sequenceIDsFilename, treeString):
    # Load species names
    sequences = dict()
    with open(sequenceIDsFilename, 'rb') as speciesFile:
        for line in speciesFile:
            id, name = line.rstrip().split(": ", 1)
            name = name.split("|", 1)[0]
            sequences[id] = name
    return ReplaceStringWithNewIDs(sequences, treeString)  
    
def ReplaceFileWithNewIDs_ete(idsMap, treeFilename, newTreeFilename):     
    #with open(treeFilename, "rb") as inputTree: treeString = inputTree.next()
    t = tree.Tree(treeFilename, format=1)
    for node in t.get_leaves():
        node.name = idsMap[node.name]
    outputTreeString = t.write(format=1)
    with open(newTreeFilename, 'wb') as outFile: outFile.write(outputTreeString)          
        
def ReplaceAlignmnetFileWithNewIDs(idsMap, inFilename, outFilename):     
    raise Exception("Not implemented yet")
    
def ReplaceInFasta(inputFilename, idsDict):
    pathfilename, ext = os.path.splitext(inputFilename)
    newFilename = pathfilename + "_names" + ext
    with open(inputFilename, 'rb') as infile, open(newFilename, 'wb') as outfile:
        for line in infile:
            if len(line) < 1 or line[0] != ">":
                outfile.write(line)
            else:
                outfile.write(">%s\n" % idsDict[line[1:].rstrip()])
    
def ReplaceFileWithNewIDs(idsMap, treeFilename, newTreeFilename):     
    treeString = ""
    with open(treeFilename, 'rb') as treeFile, open(newTreeFilename, 'wb') as newTreeFile:
        for treeString in treeFile:
            newTreeFile.write(ReplaceStringWithNewIDs(idsMap, treeString))
        
def ReplaceIDsInFile(idsDict, treeFilename, newTreeFilename):
    treeString = ""
    with open(treeFilename, 'rb') as treeFile, open(newTreeFilename, 'wb') as newTreeFile:
        for treeString in treeFile:
            newTreeFile.write(ReplaceStringWithNewIDs(idsDict, treeString))
            
def ReplaceFilePipeline(sequenceIDsFilename, treeFilename):
    pathfilename, ext = os.path.splitext(treeFilename)
    path, filename = os.path.split(pathfilename)
#    newTreeFilename= filename + "_acc" + ext
    geneIdParts = filename.split("_", 2)
    newTreeFilename = geneIdParts[0] + "_" + geneIdParts[1]
    newTreeFilename = ReplaceString(sequenceIDsFilename, newTreeFilename) + "_" + geneIdParts[2] + ext
    if path != "":
        newTreeFilename = path + "/" + newTreeFilename
#    print(newTreeFilename)
    idsDict = dict()
    with open(sequenceIDsFilename, 'rb') as speciesFile:
        for line in speciesFile:
            id, name = line.rstrip().split(": ", 1)
            name = name.split("|", 1)[0]
            idsDict[id] = name
    ReplaceIDsInFile(idsDict, treeFilename, newTreeFilename)
  
def GetSpeciesSequenceIDsDict(sequenceIDsFilename, speciesIDsFN = None):
    idsDict = dict()
    with open(sequenceIDsFilename, 'rb') as speciesFile:
        for line in speciesFile:
            id, name = line.rstrip().split(": ", 1)
#            name = name.split("|", 1)[0]
            idsDict[id] = name
    if speciesIDsFN != None:
        speciesDict = dict()
        with open(speciesIDsFN, 'rb') as speciesFile:
            for line in speciesFile:
                id, name = line.rstrip().split(": ", 1)
                if id.startswith('#'): id = id[1:] 
                speciesDict[id] = name
#        idsDict = {id: (speciesDict[id.split("_")[0]].replace(". ", "_").replace(" ", "_") + "_" + seqName.split()[0].split("|")[0].replace(":", "_")) for id, seqName in idsDict.items()}
        idsDict = {id: (speciesDict[id.split("_")[0]].replace(". ", "_").replace(" ", "_") + "_" + seqName.split()[0].replace(":", "_")) for id, seqName in idsDict.items()}
    return idsDict
          
if __name__ == "__main__":
    print("Usage:")
    print("  ReplaceSequenceNamesInTrees.py sequenceIDsFilename treeFilename [speciesIDsFilename]")
    qDoAll = False
    if len(sys.argv) < 3:
        sys.exit()
 
    treeFilename = sys.argv[2]         
    sequenceIDsFilename = sys.argv[1]
    idsDict = GetSpeciesSequenceIDsDict(sequenceIDsFilename, sys.argv[3] if len(sys.argv) == 4 else None)
            
    if qDoAll:
        filesToDo = [f for f in glob.glob(treeFilename + "/*")]
    else:
        filesToDo = [treeFilename]
    for treeFilename in filesToDo:
        pathfilename, ext = os.path.splitext(treeFilename)
        newFilename = pathfilename + "_names" + ext
        ReplaceFileWithNewIDs_ete(idsDict, treeFilename, newFilename)

            
    
