#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import itertools
import numpy as np
import scipy.sparse


class MSA(object):
    def __init__(self, fn):
        self.names = []         # list of names
        self.non_gap_pos = []   # list of list of non-gaps positions
        self.non_gaps = []      # list of list of non-gaps 
        self.n = 0              # number of sequences   
        self.length = None         # number of columns in MSA
        current_length = 0
        with open(fn, 'r') as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith(">"):
                    if len(self.names) > 1 and current_length != self.length:
                        print("Error: sequence %d has length %d, previous was %d" % (len(self.names), current_length, self.length))
                        sys.exit()
                    self.length = current_length
                    current_length = 0   
                    self.names.append(line[1:])
                    self.non_gap_pos.append([])
                    self.non_gaps.append([])
                else:
                    # process this bit of sequence
                    # the index arithmetic current_length + i works for base 0 indexing
                    self.non_gap_pos[-1].extend([current_length + i for i,c in enumerate(line) if (c != "*" and c != "-")])
                    self.non_gaps[-1].extend([c for c in line if (c != "*" and c != "-")])
                    current_length += len(line)  
        self.n = len(self.names)
        if self.n >= 2 and current_length != self.length:
            print("Error: Last sequence length is %d, previous was %d" % (current_length, self.length))
            sys.exit()
        self.length = current_length   # could just be one sequence in file
        self.length = 0 if self.length is None else self.length
        # create a sparse matrix so less of the subsequent analysis code has to change
        # Quite possibly this is the fastest way to do it too
        row_ind = [i_seq for i_seq, ngp in enumerate(self.non_gap_pos) for _ in range(len(ngp))]
        self.n_non_gaps = sum([len(ngp) for ngp in self.non_gap_pos])
        col_ind = list(itertools.chain.from_iterable(self.non_gap_pos))
        data = [1 for _ in range(self.n_non_gaps)]
        self.M = scipy.sparse.csr_matrix((data, (row_ind, col_ind)), shape=(self.n, self.length))
    
    def write_msa(self, i_cols, outfn, nChar = 80):
        with open(outfn, 'w') as outfile:
            for name, posn, chars in zip(self.names, self.non_gap_pos, self.non_gaps):
                outfile.write(">" + name + "\n")
                seq = ("-" * posn[0]) + chars[0]
                for ipos, ipos_m1, c in zip(posn[1:], posn[:-1], chars[1:]):
                    seq += "-" * (ipos - ipos_m1 - 1) + c
                seq += "-" * (self.length - posn[-1] -1)
                seq = "".join([seq[i] for i in i_cols])
                for i in range(0, len(seq), nChar):
                    outfile.write(seq[i:i+nChar] + "\n")

def main(infn, outfn, f=0.1, n_min=500, c=0.75):
    """
    Trim the alignment file. Auto reduce f parameter if required. OrthoFinder default
    parameters are used by default.
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
    if not os.path.exists(infn):
        print("ERROR, input file does not exist: %s" % infn)
        sys.exit()
    try:
        # BEWARE: Input and output files may be the same
        with open(outfn, 'a') as outfile:
            outfile.write(" ")
    except:
        print("ERROR, cannot write to output file: %s" % outfn)
        sys.exit()
    msa = MSA(infn)   
    # print(infn)
    if msa.length <= n_min:
        copy_input_to_output(infn, outfn)
        # print("Already below min length: %d" % msa.length)
        return
    # vectorised is far quicker
    # M = np.array([list(seq) for seq in msa.seqs.values()])
    # drop msa at this point
    n = msa.n
    length = msa.length
    maxGap = (1.-f)*n
    aa_counts = msa.M.sum(axis=0)[0]   # per column
    aa_counts = np.squeeze(np.asarray(aa_counts))
    gap_counts = n- aa_counts
    aa_before = msa.M.nnz
    i_keep = np.where(gap_counts <= maxGap)
    n_keep = i_keep[0].size     # it's an I, J tuple
    aa_after = sum(aa_counts[i_keep])
    # Check we have kept enough columns and characters
    if n_keep < n_min or aa_after < c*aa_before:
        # print("%0.3f: %d columns, %0.1f%% of characters retained" % (f, n_keep, 100.*aa_after/aa_before))
        f, i_keep = get_satifactory_f(gap_counts, aa_counts, n, f, n_min, c)
    n_keep = i_keep[0].size
    if n_keep == length:
        copy_input_to_output(infn, outfn)
    # M = M[:, i_keep[0]]
    # s,t = M.shape
    # aa_after = s * t - (M == '-').sum()
    msa.write_msa(i_keep[0], outfn)
    # write_msa(M, names, outfn)
    # print("%0.3f: %d->%d, %0.1f%% characters retained. Trimmed %s" % (f, length, i_keep[0].size, 100.*aa_after/aa_before, infn))


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
    # if the outfn doesn't exist doesn't the check for being the same is redundant  
    if (not os.path.exists(outfn)) or (not os.path.samefile(infn, outfn)):
        shutil.copy(infn, outfn)

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