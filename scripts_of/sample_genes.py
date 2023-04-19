"""
Sample a set of genes to use in a representative tree that aims to capture all 
parts of interest.
What to try?
- k-mer coverage

"""
import os
import random
import argparse
import warnings
from collections import defaultdict

import sklearn
from sklearn import cluster
import numpy as np
from Bio import AlignIO

from . import fasta_writer, util

use_n_auto = util.version_parse_simple(sklearn.__version__) >= util.version_parse_simple('1.2.0')


def get_embedding_vectors():
    d = defaultdict(int)  # anything not found gets mapped to 0 (e.g. U, *, O, X)
    amino_acids = [b'A', b'R', b'N', b'D', b'B', b'C', b'E', b'Q', b'Z', b'G', b'H', b'I', b'L', b'K', b'M', b'F', b'P', b'S', b'T', b'W', b'Y', b'V']
    for i,a in enumerate(amino_acids):
        d[a] = i+1
    x = [0.0, 0.7616408720368554, -0.8992027785852091, -0.07873679384886653, -0.14037322274105996, -0.3299593767217029, 0.8310546727955828, -0.5399033158394786, -0.4815946221391955, -0.593800866335506, 0.797220717024523, -0.579087697236552, 0.6228711220468686, 0.35829799783893695, -0.7688731182964063, 0.20340513326201978, -0.2699133888307127, 0.487114705529672, 0.2494989007418178, 0.3971369000035997, -0.21447784861556976, -0.48582261696986334, 0.6735046248802461]
    y = [0.0, 0.12031822087625436, -0.006548892920800078, 0.6148427092374856, 0.8999931109979885, 0.9159034615116097, -0.2248776341570619, 0.5640041523838193, 0.24951671883889143, 0.527451024374676, 0.6520120008417274, -0.229131082802156, -0.878416213428115, -1.002773729034362, 0.20934299312386676, -0.6733718616949979, -0.9829046315103517, 0.7726524237395054, 0.39518870263656747, 0.027351205616309555, -0.5587327809552116, -0.7044914945390158, -0.6873284031366291]
    return x,y,d

def embed(c, ndim = 2):
    """
    Takes an MSA and returns a feature vector. It is recommended to trim first.
    Args:
        c - numpy char array (n_seqs, n_cols)
        ndim - the dimensionality of the AA embedding
    Returns:
        M - Matrix (n_seqs, ndim*n_cols)
    """
    if ndim != 2:
        raise Exception()
    nseqs, ncols = c.shape
    x0, x1, d = get_embedding_vectors()
    M = np.empty((nseqs, ncols*ndim))
    for i in range(nseqs):
        for j in range(ncols):
            k = d[c[i,j]]
            M[i, 2*j] = x0[k]
            M[i, 2*j+1] = x1[k]
    return M

def get_kmers(s, k):
    """
    Return a list of all kmers of length k in the sequences, s.
    Args:
        s - seqeunce (str)
        k - kmer length
    """
    kmers = [s[i:i+k] for i in range(len(s)-k+1)]
    return kmers


def kmers_in_seqs(accs, fw, k):
    """
    Args:
        accs - list of accessions
        fw - the FastaWriter object
        k - kmer length
        q_aligned - is this a MSA
    Returns 
        kmer_per_seq - list of sets of kmers in the seqs
        d - dict kmer -> list of accessions the kmer is present in

    """
    kmer_per_seq = []
    d = defaultdict(list)
    for i, acc in enumerate(accs):
        if i % 10000 == 0:
            print(i)
        s = fw.SeqLists[acc]
        kmers = get_kmers(s, k)
        kmer_per_seq.append(set(kmers))
        for kmer in kmers:
            d[kmer].append(acc)
    return kmer_per_seq, d


def select_from_unaligned(infn, k, n, nmax=50000):
    """
    Identify the representative sequences
    Args:
        infn - input filename
        k - k-mer length
        n - number to select
    Returns:
        seq_dict - dictionary: acc -> sequence of the selected sequences
        FastaWriter - The FastaWriter for all sequences in input file
    """
    fw = fasta_writer.FastaWriter(infn)
    print("Read sequences")
    accs = list(fw.SeqLists.keys())    # ordered list for the sequences
    if len(accs) > nmax:
        accs = np.random.choice(accs, nmax, replace=False)
    # get k-mers in each sequence
    kmer_per_seq, d = kmers_in_seqs(accs, fw, k)
    # create a vector space of sequences based on their kmer content
    K = list(d.keys()) # the basis
    lookup = {kmer:i for i, kmer in enumerate(K)}
    print("%d unique %d-mers found" % (len(K), k))
    L = list(map(len, kmer_per_seq))
    print("Mean %0.1f unique kmers per sequence" % np.mean(L))
    # S = np.array([[kmer in kmers for kmer in K] for kmers in kmer_per_seq])
    S = [np.zeros(len(K)) for _ in accs]
    for i, seq_kmers in enumerate(kmer_per_seq):
        for kmer in seq_kmers:
            # print((i, lookup[kmer]))
            S[i][lookup[kmer]] = 1
    print("Constructed features matrix")
    # Now just need to pick the most representative ones, use a heuristic k-means clustering
    if use_n_auto:
        kmeans = cluster.KMeans(n_clusters=n, random_state=0, n_init='auto').fit(S)
    else:
        kmeans = cluster.KMeans(n_clusters=n, random_state=0).fit(S)
    labels = kmeans.predict(S)
    # clusters = cluster.Birch(n_clusters=n).fit_predict(S)
    print(labels)
    print(kmeans.cluster_centers_)
    print((max(labels), len(kmeans.cluster_centers_)))
    cluster_representative = []
    for i_clust, centre in enumerate(kmeans.cluster_centers_):
        c = centre > 0.5
        reps = []
        similarity = []
        for i_sample, l in enumerate(labels):
            if l == i_clust:
                reps.append(i_sample)
                similarity.append(S[i_sample].dot(c))
        print("%d sequences for cluster" % len(similarity))
        if len(similarity) == 0:
            continue
        j = np.argmax(similarity)
        # print(reps[j])
        # print(np.where(c))
        cluster_representative.append(reps[j])
    cluster_representative = set(cluster_representative)
    print("Requested %d, %d were unique" % (n, len(cluster_representative)))
    rep_seqs = [accs[j] for j in cluster_representative]
    return rep_seqs, fw

def write_accession_list(infn, seqs):
    outfn = infn + ".selected.txt"
    with open(outfn, 'w') as outfile:
        outfile.write("\n".join(seqs))


def msa_biopython_matrix(fn):
    """
    Create a np.chararray matrix of the MSA
    Args:
        fn - input fasta fn for the MSA
    Returns:
        z - (n_seqs x msa_length) np.chararray
        accs - orders list of gene names
    """
    msa = AlignIO.read(open(fn), 'fasta')
    n = msa.get_alignment_length()
    m = len(msa)
    accs = [msa[i].name for i in range(m)]
    # keep all columns
    I = list(range(n))
    # print((m, len(I)))
    z = np.chararray((m, len(I)))
    # print("ok")
    for ii, i in enumerate(I):
        z[:, ii] = list(msa[:,i])
    return z, accs

def run_from_aligned(infn, n_sample):
    selected = select_from_aligned(infn, n_sample)
    outfn = infn + ".selected.txt"
    with open(outfn, 'w') as outfile:
        outfile.write("\n".join(selected))

def select_from_aligned(infn, n_sample, q_trim=True):
    """
    Args:
        infn - input FASTA MSA filename
        n_sample - the number of sequences to sample
    Post-condition:
        File is created: infn + ".selected.txt" with the names of the selected
        taxa, one per line.
    """
    # aln = AlignIO.read(open(infn), 'fasta')
    # we could get a distance matrix and cluster (e.g. spectral clustering), but
    # k-means is nice in that it allows the identification of a gene closest to the
    # centre of each cluster

    # # distance based
    # calculator = DistanceCalculator('blosum62')
    # dm = calculator.get_distance(aln)
    # # create np. matrix m from the triangular dm.matrix
    # labels = sklearn.cluster.spectral_clustering(m, n_sample=n, eigen_solver='arpack')

    # k-means
    # note, the requirement should really  trail off as the number of sequences increases, 
    # or better yet, the penalty for gaps should be increased for the MSA inference
    fn_trim_temp_out = "/tmp/" + os.path.basename(infn)
    fn_align_to_use = infn
    if q_trim:
        import trim
        trim.main(infn, fn_trim_temp_out)
        fn_align_to_use = fn_trim_temp_out
    m, accs = msa_biopython_matrix(fn_align_to_use)
    n_seqs, n_cols = m.shape
    if n_seqs <= n_sample:
        write_accession_list(infn, accs)
        return accs
    # remove any which are substantially gaps, they will form a cluster!!! We don't want those sequences!
    n_gaps = (m == b'-').sum(axis=1)
    med = np.percentile(n_gaps, 75)
    ikeep = np.where(n_gaps <= med)[0]
    if (len(ikeep) * 2 < len(accs) and len(accs) > 20):
        print("WARNING only %d are non-gappy out of %d: %s" % (len(ikeep), len(accs), infn))
    if len(ikeep) <= n_sample:
        write_accession_list(infn, [accs[i] for i in ikeep])
        return [accs[i] for i in ikeep]
    d_new_old = {new:old for new, old in enumerate(ikeep)}
    m = m[ikeep, ]
    M = embed(m)
    # print("embedding successful")
    # Cluster
    # Kmeans will warn if it finds fewer clusters than requested, ignore these warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if use_n_auto:
            kmeans = cluster.KMeans(n_clusters=n_sample, random_state=0, n_init='auto').fit(M)
        else:
            kmeans = cluster.KMeans(n_clusters=n_sample, random_state=0).fit(M)
    # print("clustering successful")
    labels = kmeans.predict(M)
    n_keep = len(ikeep)
    cluster_representative = []
    for i_clust, centre in enumerate(kmeans.cluster_centers_):
        c = centre > 0.5
        reps = []
        similarity = []
        for i_sample, l in enumerate(labels):
            if l == i_clust:
                reps.append(i_sample)
                similarity.append(M[i_sample].dot(c))
        if len(similarity) == 0:
            continue
        j = np.argmax(similarity)
        cluster_representative.append(reps[j])
    cluster_representative = list(set(cluster_representative))
    n_found = len(cluster_representative)
    n_extra = n_sample - n_found
    if n_extra > 0:
        # select some more to make it up to the total 
        not_used = set(range(n_keep)).difference(cluster_representative)
        cluster_representative.extend(random.sample(not_used, n_extra))
        # print("Found extra sequences. Have %d" % len(cluster_representative))
    # print(cluster_representative)
    selected = [accs[d_new_old[i]] for i in cluster_representative]
    # print("|".join(selected))
    # fw = fasta_writer.FastaWriter(infn)
    return selected


def run(infn, q_aligned, k, n):
    """
    Main function. Read file, select representatives, write.
    Args:
        infn - input filename
        k - k-mer length
        n - number to select
    """
    if q_aligned:
        run_from_aligned(infn, n)
    else:
        rep_seqs, fw = select_from_unaligned(infn, k, n)
        outfn = infn + ".selected.fa"
        fw.WriteSeqsToFasta(rep_seqs, outfn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Fasta file of unaligned sequences")
    parser.add_argument("-k", "--kmer", type=int, help="Length of k-mers to use", default=3)
    parser.add_argument("-n", "--number", type=int, help="Number of sequences to select", default=100)
    parser.add_argument("-a", "--aligned", action='store_true', help="Input is a MSA")
    args = parser.parse_args()
    run(args.fasta, args.aligned, args.kmer, args.number)
