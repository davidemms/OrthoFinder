#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script writes out gene sequence FASTA files for the phylogenetically determined
orthogroups at a user-specified node in the species tree. These sets of orthogroups
are in the files "Phylogenetic_Hierarchical_Orthogroups/N*tsv".

Background:
OrthoFinder infers the set of Hierarchical Orthogroups (HOGs) for each node in the 
species tree, e.g. N3. A single HOG at the N3 level is the set of genes descended 
from a single gene in the ancestral species at that node. 

OrthoFinder writes out gene sequence FASTA files for each of the orthogroups in 
"Orthogroups/Orthogroups.tsv" but not for each of the sets of the orthogroups in 
the "Phylogenetic_Hierarchical_Orthogroups/N*tsv" files. This is due to the disk
space and the number of files that would be required, since most users will only 
be interested in the orthogroups at a particular level or may not even require the
sequence files.

Usage:

"""

import os
import sys
import csv
import glob
import argparse
from collections import Counter

if __name__ == "__main__" and __package__ is None:   
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts_of import util, trees_msa, orthologues, files


def process_log(logFN):
    """
    Process to Log.txt file to get wd_base and the species list
    Args:
        logFN - path for Log.txt
    Returns:
        wd_base_list - list of paths to wd_base directories
        species_ids - list of ints
    Notes:
    Should work with relevant paths to allow directory to move
    """
    with open(logFN, 'r') as infile:
        for line in infile:
            if line.startswith("Species used:"):
                species_ids_lines = ""
                line = next(infile)
                while line.rstrip() != "":
                    species_ids_lines += line
                    line = next(infile)
            wd_base_str = "WorkingDirectory_Base: "
            if line.startswith(wd_base_str): 
                wd_base_anchor = line.rstrip()[len(wd_base_str):]
                if not os.path.exists(wd_base_anchor):
                    # try to see if it's a relative directory to current one
                    path, d_wd = os.path.split(wd_base_anchor[:-1])
                    path, d_res = os.path.split(path)
                    wd_base_anchor = os.path.split(logFN)[0] + ("/../%s/%s/" % (d_res, d_wd))
                    if not os.path.exists(wd_base_anchor):
                        print("ERROR: Missing directory: %s" % wd_base_anchor)
                        util.Fail()
                wd_base_prev = files.PreviousFilesLocator_new.GetWDBaseChain(wd_base_anchor)
    species_ids_lines = species_ids_lines.split("\n")
    i_species = [int(l.split(":")[0]) for l in species_ids_lines if l != ""]
    return wd_base_prev, i_species


def read_hierarchical_orthogroup(fn, i_skip = 3):
    """
    Read genes per species, this gives a more robust way of identifying the correct
    sequence in the case where gene names are repeated across multiple species.
    Args:
        fn - HOGs fn
    Returns:
        ogs - List of orthogroups, orthogroup is list of species, species is list 
              of genes
    """
    with open(fn, csv_read_mode) as infile:
        reader = csv.reader(infile, delimiter="\t")
        ogs = []
        next(reader)   # header
        for line in reader:
            genes = [[g for g in s.split(", ") if g != ""] for s in line[i_skip:]]
            ogs.append(genes)
    return ogs


def create_files_for_node(dres, node_name, dout):
    # Check required files are present
    if not os.path.exists(dres):
        print("ERROR: %s does not exist" % dres)
        util.Fail()
    # HOG files
    d_hogs = dres + "Phylogenetic_Hierarchical_Orthogroups/"
    fn_hogs = d_hogs + node_name + ".tsv"
    if not os.path.exists(fn_hogs):
        print("ERROR: HOGs file '%s' does not exist" % fn_hogs)
        util.Fail()
    # Check can write to output directory
    d_out_fasta = create_output_directories(dout, node_name)
    ogs = read_hierarchical_orthogroup(fn_hogs)
    # Log file - to get Working directory base
    fn_log = dres + "Log.txt"
    if not os.path.exists(fn_log):
        print("ERROR: %s does not exist" % fn_log)
        util.Fail()
    wd_base_prev, i_species = process_log(fn_log)
    # Sequence files
    for isp in i_species:
        if not any(os.path.exists(d + "/Species%d.fa" % isp) for d in wd_base_prev):
            print("ERROR: Can't find Species%d.fa in any of the Working Directories:" % isp)
            print("\n".join(wd_base_prev))
            util.Fail()
    fw = trees_msa.FastaWriter(wd_base_prev, i_species)
    # convert IDs - do on a per species basis to increase likelihood of success
    ids, ids_rev = GetIDs(wd_base_prev, i_species)
    # Convert one at a time in case of error, can continue with rest
    q_warning = False
    for iog, og in enumerate(ogs):
        try:
            og_ids = [ids_rev[g] for sp in og for g in sp]    # don't care what species they're from (may need to later if this helps identify sequences)
            og_ids = [orthologues.Seq(g) for g in og_ids if g != "inf_inf"]
            fn = d_out_fasta + node_name + (".HOG%07d.fa" % iog)
            fw.WriteSeqsToFasta_withNewAccessions(og_ids, fn, ids)
        except KeyError as e:
            q_warning = True
            print("WARNING: %s.HOG%07d.fa. Problem with sequences:" % (node_name, iog))
            errors = []
            for g in og:
                if g not in ids_rev:
                    errors.append(g)
                    continue
                seq_id = ids_rev[g]
                if seq_id not in ids:
                    errors.append(g)
            print("; ".join(errors))    
    if q_warning:
        print("\nWARNING some sequence IDs were not unique and so the correct gene sequence couldn't be identified. These have been skipped")


def create_output_directories(dout, node_name):
    """
    Attempts to create the output directories, otherwise prints error and calls
    util.Fail()
    Args:
        dout - user-supplied base output directory
        node_name - HOG node name
    Returns:
        d_out_fasta - output directory for the FASTA sequence files
    """
    if not dout.endswith("/"):
        dout += "/"
    if not os.path.exists(dout):
        try:
            os.mkdir(dout)
        except OSError as e:
            print("Could not create output directory: %s" % dout)
            print(str(e))
            util.Fail()
    d_out_node = dout + node_name + "/"
    if not os.path.exists(d_out_node):
        try:
            os.mkdir(d_out_node)
        except OSError as e:
            print("Could not create output directory: %s" % d_out_node)
            print(str(e))
            util.Fail()
    d_out_fasta = d_out_node + "HOG_Sequences/"
    if not os.path.exists(d_out_fasta):
        try:
            os.mkdir(d_out_fasta)
        except OSError as e:
            print("Could not create output directory: %s" % d_out_fasta)
            print(str(e))
            util.Fail()
    return d_out_fasta


def GetIDs(wd_base_list, i_species):
    """
    Get the sequence IDs including ambiguous accessions
    Args:
        wd_base_list - list of working directories
        i_species - list of ints for species
    Returns:
        ids - OrthoFinder id to user ID
        ids_rev - User ID to OrthoFinder ID
    """
    nSpAll = max(i_species) + 1
    try:
        ids = util.FirstWordExtractor(wd_base_list[0] + "SequenceIDs.txt").GetIDToNameDict()
    except RuntimeError as error:
        ids = util.FullAccession(wd_base_list[0] + "SequenceIDs.txt").GetIDToNameDict()
    species_start = ["%d_" % isp for isp in i_species]
    ids = {k:v for k, v in ids.items() if any(k.startswith(sp_start) for sp_start in species_start)}
    # replace any ambiguous sequences
    ids_rev = {v:k for k,v in ids.items()}
    accessions_count = Counter(ids.values())
    for acc, n in accessions_count.most_common():
        if n == 1: 
            break
        ids_rev[acc] = "inf_inf"
    ids["inf_inf"] = "AMBIGUOUS_ID"
    return ids, ids_rev


def create_files(orthofinder_results_dir, node, output_directory):
    if node == "all":
        d_hogs = orthofinder_results_dir + "Phylogenetic_Hierarchical_Orthogroups/"
        nodes = glob.glob(d_hogs + "*.tsv")
        nodes = [os.path.splitext(os.path.basename(fn))[0] for fn in nodes]
        print("Writing files for orthogroups at these levels:\n%s\n" % ", ".join(nodes))
    else:
        nodes = [node, ]
    print("Processing node:")
    for node in nodes:
        print(node)
        create_files_for_node(orthofinder_results_dir, node, output_directory)
    print("Done")


def main_create_files():
    with util.Finalise():
        parser = argparse.ArgumentParser(description="Create files for the orthogroups at a particular phylogenetic level")
        parser.add_argument("orthofinder_results", help="Input directory containing OrthoFinder results")
        parser.add_argument("output_directory", help="Output directory where new files will be written")
        parser.add_argument("node_name", help="Node name, e.g. 'N0', or 'all' for HOGs at all nodes")
        args = parser.parse_args()
        orthofinder_results_dir = args.orthofinder_results
        if not orthofinder_results_dir.endswith("/"):
            orthofinder_results_dir += "/"
        create_files(orthofinder_results_dir, args.node_name, args.output_directory)

if __name__ == "__main__":
    main_create_files()
