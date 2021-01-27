from __future__ import absolute_import
from . import parallel_task_manager

import os
import numpy as np
import subprocess
from scipy import sparse
import time
import warnings
from collections import Counter, deque
import numpy.core.numeric as numeric     
from scipy.optimize import curve_fit         
import multiprocessing as mp
try: 
    import queue
except ImportError:
    import Queue as queue  


from . import util, files, blast_file_processor, matrices, mcl, orthologues, trees_msa, stats, wpgma
 

def WriteGraph_perSpecies(args):
    seqsInfo, graphFN, iSpec, d_pickle = args            
    # calculate the 2-way connections for one query species
    with open(graphFN + "_%d" % iSpec, 'w') as graphFile:
        connect2 = []
        for jSpec in range(seqsInfo.nSpecies):
            m1 = matrices.LoadMatrix("connect", iSpec, jSpec, d_pickle)
            m2tr = numeric.transpose(matrices.LoadMatrix("connect", jSpec, iSpec, d_pickle))
            connect2.append(m1 + m2tr)
            del m1, m2tr
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpec, d_pickle)
        B_connect = matrices.MatricesAnd_s(connect2, B)
        del B, connect2
        
        W = [b.sorted_indices().tolil() for b in B_connect]
        del B_connect
        for query in range(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]):
            offset = seqsInfo.seqStartingIndices[iSpec]
            graphFile.write("%d    " % (offset + query))
            for jSpec in range(seqsInfo.nSpecies):
                row = W[jSpec].getrowview(query)
                jOffset = seqsInfo.seqStartingIndices[jSpec]
                for j, value in zip(row.rows[0], row.data[0]):
                    graphFile.write("%d:%.3f " % (j + jOffset, value))
            graphFile.write("$\n")
        if iSpec == (seqsInfo.nSpecies - 1): graphFile.write(")\n")
        util.PrintTime("Written final scores for species %d to graph file" % iSpec)

"""
scnorm
-------------------------------------------------------------------------------
"""
class scnorm:
    @staticmethod
    def loglinear(x, a, b):
        return a*np.log10(x)+b     
    
    @staticmethod
    def GetLengthArraysForMatrix(m, len_i, len_j):
        I, J = m.nonzero()
        scores = [v for row in m.data for v in row]     # use fact that it's lil
        Li = np.array(len_i[I])
        Lj = np.array(len_j[J])
        return Li, Lj, scores
        
    @staticmethod
    def GetTopPercentileOfScores(L, S, percentileToKeep):
        # Get the top x% of hits at each length
        nScores = len(S)
        t_sort = sorted(zip(L, range(nScores)))
        indices = [j for i, j in t_sort]
        s_sorted = [S[i] for i in indices]
        l_sorted = [L[i] for i in indices]
        if nScores < 100:
            # then we can't split them into bins, return all for fitting
            return l_sorted, s_sorted
        nInBins = 1000 if nScores > 5000 else (200 if nScores > 1000 else 20)
        nBins, remainder = divmod(nScores, nInBins)
        topScores = []
        topLengths = []
        for i in range(nBins):
            first = i*nInBins
            last = min((i+1)*nInBins-1, nScores - 1)
            theseLengths = l_sorted[first:last+1]
            theseScores = s_sorted[first:last+1]
            cutOff = np.percentile(theseScores, percentileToKeep)
            lengthsToKeep = [thisL for thisL, thisScore in zip(theseLengths, theseScores) if thisScore >= cutOff]
            topLengths.extend(lengthsToKeep)
            topScores.extend([thisScore for thisL, thisScore in zip(theseLengths, theseScores) if thisScore >= cutOff])
        return topLengths, topScores
        
    @staticmethod
    def CalculateFittingParameters(Lf, S):
        pars,covar =  curve_fit(scnorm.loglinear, Lf, np.log10(S))
        return pars
           
    @staticmethod   
    def NormaliseScoresByLogLengthProduct(b, Lq, Lh, params): 
        rangeq = list(range(len(Lq)))
        rangeh = list(range(len(Lh)))
        li_vals = Lq**(-params[0])
        lj_vals = Lh**(-params[0])
        li_matrix = sparse.csr_matrix((li_vals, (rangeq, rangeq)))
        lj_matrix = sparse.csr_matrix((lj_vals, (rangeh, rangeh)))
        return sparse.lil_matrix(10**(-params[1]) * li_matrix * b * lj_matrix)

"""
WaterfallMethod
-------------------------------------------------------------------------------
"""       
class WaterfallMethod:  
    @staticmethod
    def GetSequenceLengths(seqsInfo):                
        sequenceLengths = []
        for iSpecies, iFasta in enumerate(seqsInfo.speciesToUse):
            sequenceLengths.append(np.zeros(seqsInfo.nSeqsPerSpecies[iFasta]))
            fastaFilename = files.FileHandler.GetSpeciesFastaFN(iFasta)
            currentSequenceLength = 0
            iCurrentSequence = -1
            qFirstLine = True
            with open(fastaFilename) as infile:
                for row in infile:
                    if len(row) > 1 and row[0] == ">":    
                        if qFirstLine:
                            qFirstLine = False
                        else:
                            sequenceLengths[iSpecies][iCurrentSequence] = currentSequenceLength
                            currentSequenceLength = 0
                        _, iCurrentSequence = util.GetIDPairFromString(row[1:])
                    else:
                        currentSequenceLength += len(row.rstrip())
            sequenceLengths[iSpecies][iCurrentSequence] = currentSequenceLength
        return sequenceLengths

    @staticmethod
    def NormaliseScores(B, Lengths, iSpecies, jSpecies):    
        Li, Lj, scores = scnorm.GetLengthArraysForMatrix(B, Lengths[iSpecies], Lengths[jSpecies])
        Lf = Li * Lj     
        topLf, topScores = scnorm.GetTopPercentileOfScores(Lf, scores, 95)   
        if len(topScores) > 1:
            fittingParameters = scnorm.CalculateFittingParameters(topLf, topScores)  
            return scnorm.NormaliseScoresByLogLengthProduct(B, Lengths[iSpecies], Lengths[jSpecies], fittingParameters)
        else:
            print("WARNING: Too few hits between species %d and species %d to normalise the scores, these hits will be ignored" % (iSpecies, jSpecies))
            return sparse.lil_matrix(B.get_shape())
    
            
    @staticmethod
    def ProcessBlastHits(seqsInfo, blastDir_list, Lengths, iSpecies, d_pickle, qDoubleBlast):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            # process up to the best hits for each species
            Bi = []
            for jSpecies in range(seqsInfo.nSpecies):
                Bij = blast_file_processor.GetBLAST6Scores(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[jSpecies], qDoubleBlast=qDoubleBlast)  
                Bij = WaterfallMethod.NormaliseScores(Bij, Lengths, iSpecies, jSpecies)
                Bi.append(Bij)
            matrices.DumpMatrixArray("B", Bi, iSpecies, d_pickle)
            BH = WaterfallMethod.GetBH_s(Bi, seqsInfo, iSpecies)
            matrices.DumpMatrixArray("BH", BH, iSpecies, d_pickle)
            util.PrintTime("Initial processing of species %d complete" % iSpecies)
        
    @staticmethod 
    def Worker_ProcessBlastHits(cmd_queue, d_pickle, qDoubleBlast):
        while True:
            try:
                args = cmd_queue.get(True, 1)
                WaterfallMethod.ProcessBlastHits(*args, d_pickle=d_pickle, qDoubleBlast=qDoubleBlast)
            except queue.Empty:
                return 
            except Exception:
                seqsInfo, _, _, iSpecies = args
                i = seqsInfo.speciesToUse[iSpecies]
                print("ERROR: Error processing files Blast%d_*" % i)
                raise

    @staticmethod
    def GetBH_s(pairwiseScoresMatrices, seqsInfo, iSpecies, tol=1e-3):
        nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpecies]]
        bestHitForSequence = -1*np.ones(nSeqs_i)
        H = [None for i_ in range(seqsInfo.nSpecies)] # create array of Nones to be replace by matrices
        for j in range(seqsInfo.nSpecies):
            if iSpecies == j:
                # identify orthologs then come back to paralogs
                continue
            W = pairwiseScoresMatrices[j]
            I = []
            J = []
            for kRow in range(nSeqs_i):
                values=W.getrowview(kRow)
                if values.nnz == 0:
                    continue
                m = max(values.data[0])
                bestHitForSequence[kRow] = m if m > bestHitForSequence[kRow] else bestHitForSequence[kRow]
                # get all above this value with tolerance
                temp = [index for index, value in zip(values.rows[0], values.data[0]) if value > m - tol]
                J.extend(temp)
                I.extend(kRow * np.ones(len(temp), dtype=np.dtype(int)))
            H[j] = sparse.csr_matrix((np.ones(len(I)), (I, J)), shape=W.get_shape())
        # now look for paralogs
        I = []
        J = []
        W = pairwiseScoresMatrices[iSpecies]
        for kRow in range(nSeqs_i):
            values=W.getrowview(kRow)
            if values.nnz == 0:
                continue
            temp = [index for index, value in zip(values.rows[0], values.data[0]) if value > bestHitForSequence[kRow] - tol]
            J.extend(temp)
            I.extend(kRow * np.ones(len(temp), dtype=np.dtype(int)))
        H[iSpecies] = sparse.csr_matrix((np.ones(len(I)), (I, J)), shape=W.get_shape())
        return H

    @staticmethod
    def ConnectCognates(seqsInfo, iSpecies, d_pickle, gathering_version): 
        # calculate RBH for species i
        BHix = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle)
        BHxi = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle, row=False)
        RBHi = matrices.MatricesAndTr_s(BHix, BHxi)   # twice as much work as before (only did upper triangular before)
        del BHix, BHxi
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpecies, d_pickle)
        connect = WaterfallMethod.ConnectAllBetterThanAnOrtholog_s(RBHi, B, seqsInfo, iSpecies, gathering_version) 
        matrices.DumpMatrixArray("connect", connect, iSpecies, d_pickle)
            
    @staticmethod 
    def Worker_ConnectCognates(cmd_queue, d_pickle, gathering_version):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            while True:
                try:
                    args = cmd_queue.get(True, 1)
                    WaterfallMethod.ConnectCognates(*args, d_pickle=d_pickle, gathering_version=gathering_version)
                except queue.Empty:
                    return  
                                   
    @staticmethod
    def WriteGraphParallel(seqsInfo, nProcess):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            with open(files.FileHandler.GetGraphFilename(), 'w') as graphFile:
                graphFile.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n" % (seqsInfo.nSeqs, seqsInfo.nSeqs)) 
                graphFile.write("\n(mclmatrix\nbegin\n\n") 
            pool = mp.Pool(nProcess)
            graphFN = files.FileHandler.GetGraphFilename()
            pool.map(WriteGraph_perSpecies, [(seqsInfo, graphFN, iSpec, files.FileHandler.GetPickleDir()) for iSpec in range(seqsInfo.nSpecies)])
            for iSp in range(seqsInfo.nSpecies):
                subprocess.call("cat " + graphFN + "_%d" % iSp + " >> " + graphFN, shell=True)
                os.remove(graphFN + "_%d" % iSp)
            # Cleanup
            pool.close()
            # matrices.DeleteMatrices("B", files.FileHandler.GetPickleDir()) 
            matrices.DeleteMatrices("connect", files.FileHandler.GetPickleDir()) 
    
    @staticmethod
    def GetMostDistant_s(RBH, B, seqsInfo, iSpec):
        # most distant RBB - as cut-off for connecting to other genes
        mostDistant = numeric.transpose(np.ones(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]])*1e9)
        # Best hit in another species: species-specific paralogues will now be connected - closer than any gene in any other species
        bestHit = numeric.transpose(np.zeros(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]))
        for kSpec in range(seqsInfo.nSpecies):
            B[kSpec] = B[kSpec].tocsr()
            if iSpec == kSpec:
                continue
            bestHit = np.maximum(bestHit, matrices.sparse_max_row(B[kSpec]))
            I, J = RBH[kSpec].nonzero()
            if len(I) > 0:
                mostDistant[I] = np.minimum(B[kSpec][I, J], mostDistant[I])
        # anything that doesn't have an RBB, set to distance to the closest gene in another species. I.e. it will hit just that gene and all genes closer to it in the its own species
        I = mostDistant > 1e8
        mostDistant[I] = bestHit[I] + 1e-6   # to connect to one in it's own species it must be closer than other species. We can deal with hits outside the species later 
        return mostDistant
        
    
    @staticmethod
    def GetMostDistant_s_estimate_from_relative_rbhs(RBH, B, seqsInfo, iSpec, p=50):
        """
        Get the cut-off score for each sequence
        Args:
            RBH - array of 0-1 RBH lil matrices
            B - array of corresponding lil score matrices (Similarity ~ 1/Distance)
            p - the percentile of times the estimate will include the most distant RBH 
        Returns:
            mostDistant - the cut-off for each sequence in this species
        Implementation:
			Stage 1: conversion factors for rbh from iSpec to furthest rbh
            - Get the relative similarity scores for RBHs from iSpec to species j
              versus species k where the two RBHs exist
            - Find the p^th percentile of what the rbh would be in species j given 
              what the rbh is in species k. This is the conversion factor for 
              similarity score of an rbh in j vs k
            - Take the column min: I.e. the expected score for the furthest RBH
            Stage 2: For each rbh observed, estimate what this would imply would be the cut-off 
             for the furthest rbh
             - Take the furthest of these estimates (could potentially use a percentile, but I
               think the point is that we should take diverged genes as a clue we should look
               further and then allow the trees to sort it out.
        """
        # 1. Reduce to column vectors
        q_debug_gather = False
        RBH = [rbh.tocsr() for rbh in RBH]
        B = [b.tocsr() for b in B]
        RBH_B = [(rbh.multiply(b)).tolil() for i, (rbh, b) in enumerate(zip(RBH, B)) if i!= iSpec]
        # create a vector of these scores
        nseqi = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]
        # nsp-1 x nseqi
        Z = np.matrix([[rhb_b.data[i][0] if rhb_b.getrowview(i).nnz == 1 else 0. for i in range(nseqi)] for rhb_b in RBH_B])   # RBH if it exists else zero
        nsp_m1 = Z.shape[0]
        Zr = 1./Z
        rs = []
        n_pairs = []
        for isp in range(nsp_m1):
            rs.append([])
            n_pairs.append([])
            for jsp in range(nsp_m1):
                if isp == jsp:
                    rs[-1].append(1.0)
                    n_pairs[-1].append(np.count_nonzero(Z[isp,:]))
                    continue
                ratios = np.multiply(Z[isp,:], Zr[jsp,:])
                i_nonzeros = np.where(np.logical_and(np.isfinite(ratios), ratios > 0))  # only those which a score is availabe for both
                ratios = ratios[i_nonzeros].reshape(-1,1) 
                r = np.percentile(ratios, 100-p)
                rs[-1].append(r)
                n_pairs[-1].append(len(ratios))
                # calculate the n vectors and take the min
        if q_debug_gather: print("\nSpecies %d" % iSpec)
        if q_debug_gather: print("RBH pair:")
        if q_debug_gather: print([x for x in range(len(RBH)) if x != iSpec])
        C = np.matrix(rs)   # Conversion matrix: 
        # To convert from a hit in species j to one in species i multiply by Cij
        if q_debug_gather: print(C)
        # To convert from a hit in species j to the furthest hit, need to take 
        # the column-wise minimum
        C = np.amin(C, axis=0)    # conversion from a RBH to the estmate of what the most distant RBH should be
        if q_debug_gather: print(C)
        if q_debug_gather: print("Number of data points:")
        if q_debug_gather: print(np.matrix(n_pairs))
        mostDistant = numeric.transpose(np.ones(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]])*1e9)

        bestHit = numeric.transpose(np.zeros(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]))
        j_conversion = 0
        for kSpec in range(seqsInfo.nSpecies):
            B[kSpec] = B[kSpec].tocsr()
            if iSpec == kSpec:
                continue
            species_factor_to_most_distant = C[0,j_conversion]
            if q_debug_gather: print(species_factor_to_most_distant)
            j_conversion += 1
            bestHit = np.maximum(bestHit, matrices.sparse_max_row(B[kSpec]))
            I, J = RBH[kSpec].nonzero()
            if len(I) > 0:
                if q_debug_gather: print("%d updated" % np.count_nonzero(species_factor_to_most_distant * B[kSpec][I, J] < mostDistant[I]))
                mostDistant[I] = np.minimum(species_factor_to_most_distant * B[kSpec][I, J], mostDistant[I])
        # anything that doesn't have an RBB, set to distance to the closest gene in another species. I.e. it will hit just that gene and all genes closer to it in the its own species
        I = mostDistant > 1e8
        mostDistant[I] = bestHit[I] + 1e-6   # to connect to one in its own species it must be closer than other species. We can deal with hits outside the species later 
        return mostDistant

    @staticmethod
    def ConnectAllBetterThanCutoff_s(B, mostDistant, seqsInfo, iSpec):
        connect = []
        nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]
        for jSpec in range(seqsInfo.nSpecies):
            M=B[jSpec].tolil()
            if iSpec != jSpec:
                IIJJ = [(i,j) for i, (valueRow, indexRow) in enumerate(zip(M.data, M.rows)) for j, v in zip(indexRow, valueRow) if v >= mostDistant[i]]
            else:
                IIJJ = [(i,j) for i, (valueRow, indexRow) in enumerate(zip(M.data, M.rows)) for j, v in zip(indexRow, valueRow) if (i != j) and v >= mostDistant[i]]
            II = [i for (i, j) in IIJJ]
            JJ = [j for (i, j) in IIJJ]
            onesArray = np.ones(len(IIJJ))
            mat = sparse.csr_matrix( (onesArray,  (II, JJ)), shape=(nSeqs_i,  seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[jSpec]]))
            connect.append(mat)
        return connect
    
    @staticmethod
    def ConnectAllBetterThanAnOrtholog_s(RBH, B, seqsInfo, iSpec, gathering_version):     
        if gathering_version < (3,0):   
            mostDistant = WaterfallMethod.GetMostDistant_s(RBH, B, seqsInfo, iSpec) 
        else:
            mostDistant = WaterfallMethod.GetMostDistant_s_estimate_from_relative_rbhs(RBH, B, seqsInfo, iSpec) 
        connect = WaterfallMethod.ConnectAllBetterThanCutoff_s(B, mostDistant, seqsInfo, iSpec)
        return connect


"""
Phase 2 - ConnectClusters
-------------------------------------------------------------------------------
"""

def ConnectClusters(clusters, nSp):
    """
    Work out all the relationships between clusters & connect into OGs any that 
    have been splintered
    Args:
        clusters - list of sets of strings

    Implementation:
    - Matrices available:
        - B - transformed bit scores
        - mostDistant - sequence-by-sequence cut-off used

    Tasks should be:
    1. Identify all relationships between orthogroups, ideally as a tree
    2. Those that are phylogenetically limited and closely related to another cluster
    should be joined to it and the later gene tree should be used to determine the
    relationships.

    Similarity measure between clusters: (Between group similarity)/(Within group similarity)
    """
    S = GetClusterSimilarities(clusters, nSp)
    # clusters = ConnectClusters_v1(S, clusters)   # connect those better than threshold
    clusters = ConnectClusters_v2(S, clusters)    # calculate connected components  
    return clusters
        
def get_latest_location(clusts, i):
    i_orig = i
    while isinstance(clusts[i], int):
        i = clusts[i]
    if i != i_orig:
        clusts[i_orig] = i # short-cut a potential series of links
    return i

def GetClusterSimilarities(clusters, nSp):
    """
    Create the similarity matrix between clusters - average score for pairs of 
    sequences across the two clusters.
    Args:
        clusters - list of lists of strings describing the genes in each cluster 
        nSp - number of species
    Returns
        S - lil sparse matrix of average similarity scores between the genes in 
            each pair of clusters
    """
    N = len(clusters)
    S = sparse.csr_matrix((N, N))
    d_pickle = files.FileHandler.GetPickleDir()
    g_to_clust = order_cluster_data(clusters, nSp)
    for iSpec in range(nSp):
        for jSpec in range(nSp):
            S += sum_cluster_hits(iSpec, jSpec, d_pickle, g_to_clust, N, "B") # CSR

    # Now normalise the entries by the number of gene pairs between each orthogroup
    n = np.array([len(c) for c in clusters])
    S = S.tolil()   # can optimise the choice for production
    for i, (Js, Vs) in enumerate(zip(S.rows, S.data)):
        ni = n[i]
        # only n*n-n pairs for self-self
        S.data[i] = [v/(ni*n[j]-ni*(i==j)) for j, v in zip(Js, Vs)]
    return S

def ConnectClusters_v2(S, clusters, threshold=0.5):
    """
    V2 - Get teh number of connected components
    Args:
        S - lil matrix of similarity scores
        clusters - list of sets of genes (strings)
    Returns:
        clusters - list of sets of genes z
    """
    n = list(map(len, clusters))
    i_single = n.index(1)
    S = normalise_cluster_scores(S, i_single)
    S = 0.5*(S + S.transpose())   # symmetrise  (or should it be one or the other is greater than the threshold?)
    print(S.shape)
    print(S.nnz)
    # set best hit for each singleton to be above the threshold
    print("%d genes" % sum(n))
    t2 = threshold*1.1
    S.setdiag(0)
    j_largest = S.argmax(0)
    n_clust = S.get_shape()[0]
    print("First singleton: %d" % i_single)
    for i in range(i_single, n_clust):
        # correct order of i & j if we used argmax(axis=0)
        if S[j_largest[0,i], i] > 0:
            S[j_largest[0,i], i] = t2     # if all entries were zero, would still return an argmax so check it is non-zero    

    # threshold the edges in the graph/matrix
    S_thresh = S.multiply(S>threshold)
    S_thresh = S_thresh.tolil()
    # print(S)
    C, n_comps = ConnectedComponents(S_thresh)
    print("%d components" % n_comps)

    # # what is the size of these?
    # n_genes = [0 for _ in range(1+max(C))]
    # for i_cluster, assignment in enumerate(C):
    #     n_genes[assignment] += n[i_cluster]
    # print("%d singletons" % n_genes.count(1))

    # import matplotlib.pyplot as plt
    # fig, [ax0, ax1] = plt.subplots(2,1)
    # ax0.set_yscale('log')
    # ax1.set_yscale('log')
    # ax0.hist(n, 50)
    # ax1.hist(n_genes, 50)
    # ax0.set_xlim([0,2500])
    # ax1.set_xlim([0,2500])
    # plt.show()
    # util.Fail()

    """    
            ************* Next Steps *************    
    - Get the phylogenetic extent of each cluster
    - Load the species tree
    - For each cluster do a leaf to branch traverse, combining as appropriate:
        - Use phylogenetic extent
        - Consider proportion of sequences that hit between the two clusters (note
          that this could be distorted by the limit on hits. Perhaps we need to only
          average over those ones we got hits for? Perhaps we should record which 
          sequences maxed out on the diamond search)
        - Do we need to look again at the scores?
    """

    # For each connected component create the WPGMA tree
    comps = [[] for _ in range(n_comps)]
    for i_clust, i_assign in enumerate(C):
        comps[i_assign].append(i_clust)
    Scoo = S.tocoo()
    for connect_comp in comps:
        if len(connect_comp) < 3:
            continue
        T = DoTreeOfClusters(Scoo, connect_comp)
        if len(str(T)) > 6:
            print(T)
        # create OGs

    return clusters

def DoTreeOfClusters(S, connect_comp):
    if len(connect_comp) < 2:
        return connect_comp
    # create a dense matrix
    connect_comp = np.array(connect_comp)
    M = matrices.coo_submatrix_pull(S, connect_comp, connect_comp)
    M = M.todense()
    # do tree
    try:
        return wpgma.wpgma(M, list(connect_comp))
    except:
        print(connect_comp)
        print(M)
        raise

def ConnectClusters_v1(S, clusters):
    """
    V1 - Connect clusters if they have a between score >= 0.5 of the within cluster
    score or they are a singleton
    Args:
        S - lil matrix of similarity scores
        clusters - list of sets of genes (strings)
    Returns:
        clusters - list of sets of genes z
    """
    sizes = list(map(len, clusters))
    i_single = sizes.index(1)
    S = normalise_cluster_scores(S, i_single)
    # What are the largest non-diagonal entries
    # Hits (between clusters)
    H = []
    IJs = []
    # can be more efficient
    for i, (Js, Vs) in enumerate(zip(S.rows, S.data)):
        H.extend([v for j, v in zip(Js, Vs) if i != j])
        IJs.extend([(i,j) for j in Js if i != j])

    # indices = list(range(len(H)))

    print("i_single: %d" % i_single)
    top_hits = sorted(zip(H, IJs), reverse=True)
    print([xx[0] for xx in top_hits])
    n_print = 100
    for x, ij in top_hits:
        if ij[0] >= i_single or ij[1] >= i_single:
            continue
        print((x, ij))
        n_print -= 1
        if n_print == 0:
            break

    print("# Cut largest clusters down first so that we know what the computational cost currently is")

    # What criteria could be used to determine the clusters should be joined?
    # Joined in the OrthoFinder graph
    # Singleton gene (especially since they can be ejected due to number of hits exceeding max requested) 
    # It joins together the complete set of sequences that are hit by any of them
    # It makes minimal difference to the cost of the trees (i.e. small are fine, after we've got the largest OGs down)
    # If they both contain all species we don't want to join
    # ... but where do we cut?
    # If there is significant overlap between the best between hist and the worst within hits

    # join up, sort out
    # Join all single genes
    # Join all ...
    
    # Join all singletons

    # Join all others if between score >= 0.5 * within score and not all species present in both
    trees_before = np.array([len(c) for c in clusters if len(c) >= 4])
    # join them in the smallest index, replace larger index with pointer to joined clusters
    print("Put a better check in than looking to see if it's OG 0")
    for x, (i,j) in top_hits:
        if i < 5 or j < 5:
            continue   # don't do this, it is dangerous if the OG 0 genes aren't in this dataset
        if x < 0.5 and i < i_single and j < i_single:
            continue
        i_new = get_latest_location(clusters, i)
        j_new = get_latest_location(clusters, j)
        if i_new == j_new:
            # these have already been joined
            continue
        i_low = min(i_new, j_new)
        j_high = max(i_new, j_new)
        print((i_low, j_high))
        # print(clusters[i_low])
        clusters[i_low].update(clusters[j_high])
        # print(clusters[i_low])
        clusters[j_high] = i_low

    # remove the pointers to the combined trees
    clusters = [c for c in clusters if not isinstance(c, int)]
    print("%d clusters before, cost %d" % (len(trees_before), sum(trees_before*trees_before)))
    trees_after = np.array([len(c) for c in clusters if len(c) >= 4])
    print("%d clusters after, cost %d" % (len(trees_after), sum(trees_after*trees_after)))
    return clusters

def normalise_cluster_scores(S, i_single):
    """
    Divide non-zero entries of a square lil matrix by the diagonal entries
    Args:
        S - nxn lil sparse matrix
    Returns:
        S - nxn lil sparse matrix
    """
    # What is the biggest ratio of scores between vs within
    # (note, I will look at this compared to the score for each of the individual 
    # clusters, either of these comparisons might indicate a connection should be made)
    # Since S is symmetrical, ony need to do this once and then look at Sij and Sji
    # with np.errstate(divide='ignore'):
    #     within_m1 = 1./S.diagonal()   # 1/(within cluster similarity)
    S0=S.copy()
    d = S.diagonal()
    d[i_single:] = 1.0
    within_m1 = 1./d
    for i, (Js, Vs) in enumerate(zip(S.rows, S.data)):
        # only n*n-n pairs for self-self
        S.data[i] = [v*within_m1[i] for j, v in zip(Js, Vs)]
    return S

def order_cluster_data(clusters, nSp):
    """
    Return the clusters sorted by species and the cluster lookup dict
    Args:
        clusters - list of lists of strings describing the genes in each cluster 
    Returns:
        C - clusters as list (clusters) of lists (species within clusters) of lists (genes from that species) of ints
        Singletons - list of genes as (iSp, iSeq)
        g_to_clust - dist from (iSp, iSeq) to iCluster
    """
    # C = []
    # Singletons = []
    g_to_clust = dict()
    for iClus, clust in enumerate(clusters):
        if len(clust) == 1:
            break
        # c = [[] for _ in nSp]
        for g in clust:
            iSp, iSeq = map(int, g.split("_"))
            # c[iSp].append(iSeq)
            g_to_clust[(iSp, iSeq)] = iClus
        # C.append(c)
    # Now process the singletons
    # print(iClus, clust)
    nMulti = iClus   # since we are one past the end
    for iSingleton, clust in enumerate(clusters[nMulti:]):
        # print(nMulti + iSingleton, clust)
        iSp, iSeq = map(int,next(g for g in clust).split("_"))
        # Singletons.append(map(int, g))
        g_to_clust[(iSp, iSeq)] = nMulti + iSingleton
    return g_to_clust

def ConnectedComponents(S):
    """
    Return a list of the connected components of sparse lil matrix S
    Args:
        S - lil sparse matrix (0,1)
    Returns:
        C - list of components
    """
    # A = []   # adjacency list
    # for i, (Js, Vs) in enumerate(zip(S.rows, S.data)):
    #     A.extend([(i,j) for j, v in zip(Js, Vs) if i<j and v>0])
    n = S.get_shape()[0]
    print("%d vertices, %d edges" % (n, S.nnz/2))
    explored = [False for _ in range(n)]
    i_comp = -1
    assignments = [None for _ in range(n)]
    q = deque()
    for u in range(n):
        if explored[u]:
            continue
        i_comp += 1
        q.append(u)
        explored[u] = True
        while q:
            v = q.popleft()
            assignments[v] = i_comp
            for w in S.rows[v]:
                if not explored[w]:
                    q.append(w)
                    explored[w] = True
    n_comps = i_comp + 1
    return assignments, n_comps


def sum_cluster_hits(iSpec, jSpec, d_pickle, g_to_clust, N, input_matrix_name):
    """
    Sum up all the hits between (or within) clusters in the matrix B_{iSpec, jSpec}
    Args:
        iSpec - iSpecies
        jSpec - jSpecies
        d_pickle - Pickle directory
        g_to_clust - dict from a gene to its cluster
        N - Number of clusters
    Returns:
        S - (nClust x nClust) sum of hits as a CSR sparse matrix 
    """
    # Iterate through the matrix and add the scores 
    B = matrices.LoadMatrix(input_matrix_name, iSpec, jSpec, d_pickle)    # lil_matrix
    iClusts = []
    jClusts = []
    Vals = []
    q_same_species = iSpec == jSpec
    for i, (Js, Vs) in enumerate(zip(B.rows, B.data)):
        for j, v in zip(Js, Vs):
            # don't count self hits
            if q_same_species and i == j:
                continue
            iClusts.append(g_to_clust[(iSpec, i)])
            jClusts.append(g_to_clust[(jSpec, j)])
            Vals.append(v)
    S = sparse.coo_matrix((Vals, (iClusts, jClusts)), shape=(N,N))
    return S.tocsr()

    

"""
Orthogroups
-------------------------------------------------------------------------------
"""

def DoOrthogroups(options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML=None):
    # Run Algorithm, cluster and output cluster files with original accessions
    util.PrintUnderline("Running OrthoFinder algorithm")
    # it's important to free up the memory from python used for processing the genomes
    # before launching MCL becuase both use sizeable ammounts of memory. The only
    # way I can find to do this is to launch the memory intensive python code 
    # as separate process that exits before MCL is launched.
    Lengths = WaterfallMethod.GetSequenceLengths(seqsInfo)
    
    # Process BLAST hits
    util.PrintTime("Initial processing of each species")
    cmd_queue = mp.Queue()
    blastDir_list = files.FileHandler.GetBlastResultsDir()
    for iSpecies in range(seqsInfo.nSpecies):
        cmd_queue.put((seqsInfo, blastDir_list, Lengths, iSpecies))
    files.FileHandler.GetPickleDir()     # create the pickle directory before the parallel processing to prevent a race condition
    runningProcesses = [mp.Process(target=WaterfallMethod.Worker_ProcessBlastHits, args=(cmd_queue, files.FileHandler.GetPickleDir(), options.qDoubleBlast)) for i_ in range(options.nProcessAlg)]
    for proc in runningProcesses:
        proc.start()
    parallel_task_manager.ManageQueue(runningProcesses, cmd_queue)
    
    # if options
    cmd_queue = mp.Queue()
    for iSpecies in range(seqsInfo.nSpecies):
        cmd_queue.put((seqsInfo, iSpecies))
    runningProcesses = [mp.Process(target=WaterfallMethod.Worker_ConnectCognates, args=(cmd_queue, files.FileHandler.GetPickleDir(), options.gathering_version)) for i_ in range(options.nProcessAlg)]
    for proc in runningProcesses:
        proc.start()
    parallel_task_manager.ManageQueue(runningProcesses, cmd_queue)
    
    util.PrintTime("Connected putative homologues") 
    WaterfallMethod.WriteGraphParallel(seqsInfo, options.nProcessAlg)
    
    # 5b. MCL     
    clustersFilename, clustersFilename_pairs = files.FileHandler.CreateUnusedClustersFN(options.mclInflation) 
    graphFilename = files.FileHandler.GetGraphFilename() 
    mcl.MCL.RunMCL(graphFilename, clustersFilename, options.nProcessAlg, options.mclInflation)
    mcl.ConvertSingleIDsToIDPair(seqsInfo, clustersFilename, clustersFilename_pairs)   
    
    util.PrintUnderline("Writing orthogroups to file")
    ogs = mcl.GetPredictedOGs(clustersFilename_pairs)
    
    # 5c. Connect clusters of limited phylogenetic extent
    ogs = ConnectClusters(ogs, seqsInfo.nSpecies)
    
    mcl.RewriteMCLpairsFile(ogs, clustersFilename_pairs)

    resultsBaseFilename = files.FileHandler.GetOrthogroupResultsFNBase()
    idsDict = mcl.MCL.WriteOrthogroupFiles(ogs, [files.FileHandler.GetSequenceIDsFN()], resultsBaseFilename, clustersFilename_pairs)
    mcl.MCL.CreateOrthogroupTable(ogs, idsDict, speciesNamesDict, speciesInfoObj.speciesToUse, resultsBaseFilename)
    
    # Write Orthogroup FASTA files    
    ogSet = orthologues.OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll, options.qAddSpeciesToIDs, idExtractor = util.FirstWordExtractor)
    treeGen = trees_msa.TreesForOrthogroups(None, None, None)
    fastaWriter = trees_msa.FastaWriter(files.FileHandler.GetSpeciesSeqsDir(), speciesInfoObj.speciesToUse)
    d_seqs = files.FileHandler.GetResultsSeqsDir()
    if not os.path.exists(d_seqs): os.mkdir(d_seqs)
    treeGen.WriteFastaFiles(fastaWriter, ogSet.OGs(qInclAll=True), idsDict, False)
    
    stats.Stats(ogs, speciesNamesDict, speciesInfoObj.speciesToUse, files.FileHandler.iResultsVersion)
    if options.speciesXMLInfoFN:
        mcl.MCL.WriteOrthoXML(speciesXML, ogs, seqsInfo.nSeqsPerSpecies, idsDict, resultsBaseFilename + ".orthoxml", speciesInfoObj.speciesToUse)
    print("")
    util.PrintTime("Done orthogroups")
    files.FileHandler.LogOGs()