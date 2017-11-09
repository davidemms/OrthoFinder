#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import glob
import argparse

if __name__ == "__main__" and __package__ is None:   
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import tree
        
def ReplaceFileWithNewIDs(idsMap, treeFilename, newTreeFilename):     
    t = tree.Tree(treeFilename, format=1)
    for node in t.get_leaves():
        node.name = idsMap[node.name]
    outputTreeString = t.write(format=1)
    with open(newTreeFilename, 'wb') as outFile: outFile.write(outputTreeString)          
                             
def GetSpeciesSequenceIDsDict(sequenceIDsFilename, speciesIDsFN = None):
    idsDict = dict()
    with open(sequenceIDsFilename, 'rb') as speciesFile:
        for line in speciesFile:
            id, name = line.rstrip().split(": ", 1)
            idsDict[id] = name
    if speciesIDsFN != None:
        speciesDict = dict()
        with open(speciesIDsFN, 'rb') as speciesFile:
            for line in speciesFile:
                id, name = line.rstrip().split(": ", 1)
                if id.startswith('#'): id = id[1:] 
                speciesDict[id] = name
        idsDict = {id: (speciesDict[id.split("_")[0]].replace(". ", "_").replace(" ", "_") + "_" + seqName.split()[0].replace(":", "_")) for id, seqName in idsDict.items()}
    return idsDict
          
if __name__ == "__main__":
    parser = argparse.ArgumentParser("tree_ids_to_accessions takes a tree with OrthoFinder IDs and outputs a tree with gene accessions")
    parser.add_argument("SequenceIDsFilename")
    parser.add_argument("Tree")
    parser.add_argument('SpeciesIDsFilename', nargs='?', default=None)
    parser.add_argument('-d', "--directory", action="store_true", help="Process all trees in 'Tree' directory")

    args = parser.parse_args()
    idsDict = GetSpeciesSequenceIDsDict(args.SequenceIDsFilename, args.SpeciesIDsFilename)
        
    if args.directory:
        filesToDo = [f for f in glob.glob(args.Tree + "/*")]
    else:
        filesToDo = [args.Tree]
        
    for treeFilename in filesToDo:
        try:
            sys.stdout.write(treeFilename)
            pathfilename, ext = os.path.splitext(treeFilename)
            newFilename = pathfilename + "_accessions" + ext
            ReplaceFileWithNewIDs(idsDict, treeFilename, newFilename)
        except:
            sys.stdout.write(" - skipped")
        print("")