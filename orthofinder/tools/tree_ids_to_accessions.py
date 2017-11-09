#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse

if __name__ == "__main__" and __package__ is None:   
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import tree, util
        
def ReplaceFileWithNewIDs(idsMap, treeFilename, newTreeFilename):     
    t = tree.Tree(treeFilename, format=1)
    for node in t.get_leaves():
        node.name = idsMap[node.name]
    outputTreeString = t.write(format=1)
    with open(newTreeFilename, 'wb') as outFile: outFile.write(outputTreeString)          
                             
def GetSpeciesSequenceIDsDict(sequenceIDsFilename, speciesIDsFN = None):
    try:
        extract = util.FirstWordExtractor(sequenceIDsFilename)
    except RuntimeError as error:
        print(error.message)
        if error.message.startswith("ERROR"): 
            util.Fail()
        else:
            print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup\nmore concisely but these were not unique. The full accession line will be used instead.\n")     
            extract = util.FullAccession(sequenceIDsFilename)
    idsDict = extract.GetIDToNameDict()    
    if speciesIDsFN != None:
        speciesDict = util.FullAccession(speciesIDsFN).GetIDToNameDict()
        speciesDict = {k:v.rsplit(".",1)[0].replace(".", "_").replace(" ", "_") for k,v in speciesDict.items()}
        idsDict = {seqID:speciesDict[seqID.split("_")[0]] + "_" + name for seqID, name in idsDict.items()}
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