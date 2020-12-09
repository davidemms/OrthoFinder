#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import numpy as np


class MSA(object):
    def __init__(self, msa_dict):
        self.seqs = msa_dict
        self.length = len(list(msa_dict.values())[0]) if len(msa_dict) > 0 else 0
        self.n = len(msa_dict)


def main(infn, outfn, f, n_min, c):
    """
    Trim the alignment file. Auto reduce f if required
    Args:
        infn - input FASTA filename
        outfn - output FASTA filename
        f - Lower limit for fraction of meaningful characters
        n_min - The min trimmed length
        c - The min fraction of non-gap characters to conserve
    Notes:
        If f is too high then it will reduce f s.t. len >= n_min or the full 
        alignment is kept (if less than n_min).
        In some way fraction of gaps is unsatisfactory since as the number of sequences
        in the MSA gets larger a column could contain a high fraction of gaps and 
        yet still contain useful phylogenetic signal for (many) sequences that are
        closely related and don't have a gap in that column. For this reason the 
        c parameter is included which specifies the minimum total fraction of the 
        original non-gap characters that must be preserved.
    """
    msa = ReadAlignment(infn)   
    # print(infn)
    if msa.length <= n_min:
        copy_input_to_output(infn, outfn)
        # print("Already below min length: %d" % msa.length)
        return
    # vectorised is far quicker
    names = list(msa.seqs.keys())
    M = np.array([list(seq) for seq in msa.seqs.values()])
    maxGap = (1.-f)*msa.n
    gap_counts = (M == "-").sum(axis=0)
    aa_counts = msa.n - gap_counts
    aa_before = sum(aa_counts)
    i_keep = np.where(gap_counts <= maxGap)
    n_keep = i_keep[0].size     # it's an I, J tuple
    aa_after = sum(aa_counts[i_keep])
    # Check we have kept enough columns and characters
    if n_keep < n_min or aa_after < c*aa_before:
        # print("%0.3f: %d columns, %0.1f%% of characters retained" % (f, n_keep, 100.*aa_after/aa_before))
        f, i_keep = get_satifactory_f(gap_counts, aa_counts, msa.n, f, n_min, c)
    n_keep = i_keep[0].size
    if n_keep == msa.length:
        copy_input_to_output(infn, outfn)
    M = M[:, i_keep[0]]
    s,t = M.shape
    aa_after = s * t - (M == '-').sum()
    write_msa(M, names, outfn)
    # print("%0.3f: %d->%d, %0.1f%% characters retained. Trimmed %s" % (f, msa.length, i_keep[0].size, 100.*aa_after/aa_before, infn))


def get_satifactory_f(gap_counts, aa_counts, N, f_orig, n_min, c, tol = 0.001):
    """
    The f used was too large, get an f that gives n_min columns
    Args:
        gap_counts = np array of gaps per column
        aa_counts = np array of amino acids per column
        N - number of sequences
        f_orig - starting f for the search
        n_min - minimum alignment length
        c - minimum fraction of amino acids to keep
        tol - convergence tolerance
    Implementation:
        down - lower bound on non-gap percentage, giving larger alignment
        up - upper bound on non-gap character percentage
    """
    # binary search 
    aa_before = sum(aa_counts)
    up = f_orig
    down = 0.
    x_prev = f_orig
    x = 0.5*(down+up)
    while abs(x-x_prev) > tol:
        max_gap = (1.-x)*N
        i_keep = np.where(gap_counts <= max_gap)
        n_keep = i_keep[0].size
        aa_new = sum(aa_counts[i_keep])
        # print("%0.3f: %d columns, %0.1f%% of characters retained" % (x, n_keep, 100.*aa_new/aa_before))
        if n_keep < n_min or aa_new < c*aa_before:
            up = x
        elif (n_keep == n_min and aa_new >= c*aa_before) or (n_keep >= n_min and aa_new == c*aa_before):
            # one on the limit, the other satisfied - stop exactly here
            break
        else:
            down = x
        x_prev = x
        x = 0.5*(down+up)
    # Use the value which gave just above the min length
    x = down
    max_gap = (1.-x)*N
    i_keep = np.where(gap_counts <= max_gap)
    aa_new = sum(aa_counts[i_keep])
    # print("%0.3f: %d columns, %0.1f%% characters retained" % (x, n_keep, 100.*aa_new/aa_before))
    return x, i_keep


def copy_input_to_output(infn, outfn):        
    if (not os.path.exists(outfn)) or (not os.path.samefile(infn, outfn)):
        shutil.copy(infn, outfn)


def write_msa(M, names, outfn, nChar = 80):
    with open(outfn, 'w') as outfile:
        for iSeq, name in enumerate(names):
            outfile.write(">%s\n" % name)
            seq = M[iSeq,:].tolist()
            for i in range(0, len(seq), nChar):
                outfile.write("".join(seq[i:i+nChar]) + "\n")


def ReadAlignment(fn):
    msa = dict()
    accession = None
    length = None
    seq = ""
    with open(fn, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith(">"):
                if accession != None:
                    if length != None and len(seq) != length:
                        print("ERROR: Sequence length mismatch in MSA: %s & %d" % (length, len(seq)))
                        sys.exit()
                    msa[accession] = seq
                accession = line[1:]
                seq = ""
            else:
                seq += line.replace('*', '-')
        if accession != None:
            if length != None and len(seq) != length:
                print("Error: Sequence length mismatch in MSA: %s & %d" % (length, len(seq)))
                sys.exit()
            msa[accession] = seq
    return MSA(msa)    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infn", help="Input alignment fasta file")
    parser.add_argument("outfn", help="Output alignment fasta file")
    parser.add_argument("gt", 
                        type=float,
                        help="Lower limit for fraction of non-gap characters")
    parser.add_argument("min", type=int, help="Minimum number of columns to keep") 
    parser.add_argument("-c", "--conserve", type=float, 
                        default=0.0,
                        help="Conserve at least this fraction of non-gap characters") 
    args = parser.parse_args()
    main(args.infn, args.outfn, args.gt, args.min, args.conserve)