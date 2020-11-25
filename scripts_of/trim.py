import os
import shutil
import argparse
import numpy as np


class MSA(object):
    def __init__(self, msa_dict):
        self.seqs = msa_dict
        self.length = len(list(msa_dict.values())[0])
        self.n = len(msa_dict)


def main(infn, outfn, f, n_min):
    """
    Trim the alignment file. Auto reduce f if required
    Args:
        infn - input FASTA filename
        outfn - output FASTA filename
        f - Lower limit for fraction of meaningful characters
    Notes:
        If f is too high then it will reduce f s.t. len >= n_min or the full 
        alignment is kept (if less than n_min).
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
    # M = np.array([list(msa.seqs[name]) for name in names])
    maxGap = (1.-f)*msa.n
    gap_counts = sum(M == "-")
    i_keep = np.where(gap_counts <= maxGap)
    n_keep = i_keep[0].size
    # Check we have kept enough columns
    if n_keep < n_min:
        # print("%0.3f: %d columns" % (f, n_keep))
        x, i_keep = get_satifactory_f(gap_counts, msa.n, f, n_min)
    n_keep = i_keep[0].size
    if n_keep == msa.length:
        copy_input_to_output(infn, outfn) 
    M = M[:, i_keep]
    write_msa(M, names, outfn)
    # print(msa.length)
    # print(i_keep[0].size)
    # print("%d->%d: Trimmed %s" % (msa.length, i_keep[0].size, infn))


def copy_input_to_output(infn, outfn):        
    if not os.path.samefile(infn, outfn):
        shutil.copy(infn, outfn)


def write_msa(M, names, outfn, nChar = 80):
    with open(outfn, 'w') as outfile:
        for iSeq, name in enumerate(names):
            outfile.write(">%s\n" % name)
            seq = M[iSeq,:].tolist()[0]
            for i in range(0, len(seq), nChar):
                outfile.write("".join(seq[i:i+nChar]) + "\n")


def get_satifactory_f(gap_counts, N, f_orig, n_min, tol = 0.001):
    """
    The f used was too large, get an f that gives n_min columns
    Args:
        gap_counts = np vector of gaps per column
        N - number of sequences
        f_orig - starting f for the search
        n_min - minimum alignment length
        tol - convergence tolerance
    Implementation:
        down - small meaningful character percentage, giving larger alignment
        up - large meaningful character percentage
    """
    # binary search 
    up = f_orig
    down = 0.
    x_prev = f_orig
    x = 0.5*(down+up)
    while abs(x-x_prev) > tol:
        max_gap = (1.-x)*N
        i_keep = np.where(gap_counts <= max_gap)
        n_keep = i_keep[0].size
        # print("%0.3f: %d columns" % (x, n_keep))
        if n_keep < n_min:
            up = x
        elif n_keep > n_min:
            down = x
        else:
            break
        x_prev = x
        x = 0.5*(down+up)
    # Use the value which gave just above the min length
    x = down
    max_gap = (1.-x)*N
    i_keep = np.where(gap_counts <= max_gap)
    # print("%0.3f: %d columns" % (x, n_keep))
    return x, i_keep

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
                        text  = "ERROR: Sequence length mismatch in MSA: %s & %d" % (length, len(seq))
                        files.FileHandler.LogFailAndExit(text)
                    msa[accession] = seq
                accession = line[1:]
                seq = ""
            else:
                seq += line.replace('*', '-')
        if accession != None:
            if length != None and len(seq) != length:
                text = "Error: Sequence length mismatch in MSA: %s & %d" % (length, len(seq))
                files.FileHandler.LogFailAndExit(text)
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
    args = parser.parse_args()
    main(args.infn, args.outfn, args.gt, args.min)