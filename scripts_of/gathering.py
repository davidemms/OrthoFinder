from __future__ import absolute_import
from . import parallel_task_manager

import os
import numpy as np
import subprocess
from scipy import sparse
import random
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

from . import util, files, blast_file_processor, matrices, mcl, orthologues, trees_msa, stats


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

def WriteGraph_perSpecies_homology(args):
    seqsInfo, graphFN, iSpec, d_pickle = args
    # calculate the 2-way connections for one query species
    # W = [matrices.LoadMatrix("B", iSpec, jSpec, d_pickle).tolil() for jSpec in range(seqsInfo.nSpecies)]
    W = []
    for jSpec in range(seqsInfo.nSpecies):
        w1 = matrices.LoadMatrix("B", iSpec, jSpec, d_pickle)
        matrices.DumpMatrix("H", (w1>0).tolil(), iSpec, jSpec, d_pickle)
        w2tr = numeric.transpose(matrices.LoadMatrix("B", jSpec, iSpec, d_pickle))
        W.append((w1 + w2tr > 0).tolil())  # symmetrise
    # matrices.DumpMatrixArray("H", W, iSpec, d_pickle)
    with open(graphFN + "_%d" % iSpec, 'w') as graphFile:
        for query in range(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]):
            offset = seqsInfo.seqStartingIndices[iSpec]
            graphFile.write("%d    " % (offset + query))
            for jSpec in range(seqsInfo.nSpecies):
                row = W[jSpec].getrowview(query)
                jOffset = seqsInfo.seqStartingIndices[jSpec]
                for j in row.rows[0]:
                    graphFile.write("%d:%.3f " % (j + jOffset, 1.0))
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
    def NormaliseScores(B, Lengths, iSpecies, jSpecies):
        Li, Lj, scores = scnorm.GetLengthArraysForMatrix(B, Lengths[iSpecies], Lengths[jSpecies])
        Lf = Li * Lj
        topLf, topScores = scnorm.GetTopPercentileOfScores(Lf, scores, 95)
        if len(topScores) > 1:
            fittingParameters = scnorm.CalculateFittingParameters(topLf, topScores)
            return scnorm.NormaliseScoresByLogLengthProduct(B, Lengths[iSpecies], Lengths[jSpecies], fittingParameters)
        elif iSpecies == jSpecies and B.nnz == 0:
            print(
                "WARNING: THIS IS UNCOMMON, there are zero hits when searching the genes in species %d against itself. Check the input proteome contains all the genes from that species and check the search program is working (default is diamond)." % iSpecies)
            return sparse.lil_matrix(B.get_shape())
        else:
            print(
                "WARNING: Too few hits between species %d and species %d to normalise the scores, these hits will be ignored" % (
                iSpecies, jSpecies))
            return sparse.lil_matrix(B.get_shape())
    @staticmethod
    def NormalisedBitScore(B, Lengths, iSpecies, jSpecies):
        """
        Args:
            B - LIL matrix
        Returns
            B' - LIL matrix
        """
        if B.nnz == 0:
            return B
        Lq = Lengths[iSpecies]
        Lh = Lengths[jSpecies]
        rangeq = list(range(len(Lq)))
        rangeh = list(range(len(Lh)))
        li_vals = Lq**(-0.5)
        lj_vals = Lh**(-0.5)
        li_matrix = sparse.csr_matrix((li_vals, (rangeq, rangeq)))
        lj_matrix = sparse.csr_matrix((lj_vals, (rangeh, rangeh)))
        return sparse.lil_matrix(li_matrix * B * lj_matrix)

    @staticmethod
    def ProcessBlastHits(seqsInfo, blastDir_list, Lengths, iSpecies, d_pickle, qDoubleBlast, v2_scores, q_allow_empty=False):
        """
        iLimitNewSpecies: int - Only process fasta files with OrthoFidner ID >= iLimitNewSpecies
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # process up to the best hits for each species
            Bi = []
            for jSpecies in range(seqsInfo.nSpecies):
                Bij = blast_file_processor.GetBLAST6Scores(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies],
                                                           seqsInfo.speciesToUse[jSpecies], qDoubleBlast=qDoubleBlast, q_allow_empty=q_allow_empty)
                if v2_scores:
                    Bij = WaterfallMethod.NormalisedBitScore(Bij, Lengths, iSpecies, jSpecies)
                else:
                    Bij = WaterfallMethod.NormaliseScores(Bij, Lengths, iSpecies, jSpecies)
                Bi.append(Bij)
            matrices.DumpMatrixArray("B", Bi, iSpecies, d_pickle)
            BH = WaterfallMethod.GetBH_s(Bi, seqsInfo, iSpecies)
            matrices.DumpMatrixArray("BH", BH, iSpecies, d_pickle)
            util.PrintTime("Initial processing of species %d complete" % iSpecies)

    @staticmethod
    def Worker_ProcessBlastHits(seqsInfo, blastDir_list, Lengths, cmd_queue, d_pickle, qDoubleBlast, v2_scores, q_allow_empty):
        while True:
            try:
                iSpecies = cmd_queue.get(True, 1)
                WaterfallMethod.ProcessBlastHits(seqsInfo, blastDir_list, Lengths, iSpecies, d_pickle=d_pickle,
                                                 qDoubleBlast=qDoubleBlast, v2_scores=v2_scores, q_allow_empty=q_allow_empty)
            except queue.Empty:
                return
            except Exception:
                i = seqsInfo.speciesToUse[iSpecies]
                print("ERROR: Error processing files Blast%d_*" % i)
                raise

    @staticmethod
    def GetBH_s(pairwiseScoresMatrices, seqsInfo, iSpecies, tol=1e-3):
        nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpecies]]
        bestHitForSequence = -1 * np.ones(nSeqs_i)
        H = [None for i_ in range(seqsInfo.nSpecies)]  # create array of Nones to be replace by matrices
        for j in range(seqsInfo.nSpecies):
            if iSpecies == j:
                # identify orthologs then come back to paralogs
                continue
            W = pairwiseScoresMatrices[j]
            I = []
            J = []
            for kRow in range(nSeqs_i):
                values = W.getrowview(kRow)
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
            values = W.getrowview(kRow)
            if values.nnz == 0:
                continue
            temp = [index for index, value in zip(values.rows[0], values.data[0]) if
                    value > bestHitForSequence[kRow] - tol]
            J.extend(temp)
            I.extend(kRow * np.ones(len(temp), dtype=np.dtype(int)))
        H[iSpecies] = sparse.csr_matrix((np.ones(len(I)), (I, J)), shape=W.get_shape())
        return H

    @staticmethod
    def ConnectCognates(seqsInfo, iSpecies, d_pickle, v2_scores=False):
        # calculate RBH for species i
        BHix = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle)
        BHxi = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle, row=False)
        RBHi = matrices.MatricesAndTr_s(BHix, BHxi)  # twice as much work as before (only did upper triangular before)
        del BHix, BHxi
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpecies, d_pickle)
        connect = WaterfallMethod.ConnectAllBetterThanAnOrtholog_s(RBHi, B, seqsInfo, iSpecies, v2_scores)
        matrices.DumpMatrixArray("connect", connect, iSpecies, d_pickle)

    @staticmethod
    def Worker_ConnectCognates(cmd_queue, d_pickle, v2_scores=True):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            while True:
                try:
                    args = cmd_queue.get(True, 1)
                    WaterfallMethod.ConnectCognates(*args, d_pickle=d_pickle, v2_scores=v2_scores)
                except queue.Empty:
                    return

    @staticmethod
    def WriteGraphParallel(seqsInfo, nProcess, i_unassigned=None, func=WriteGraph_perSpecies):
        graphFN = files.FileHandler.GetGraphFilename(i_unassigned)
        # Should use PTM?
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with open(graphFN, 'w') as graphFile:
                graphFile.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n" % (seqsInfo.nSeqs, seqsInfo.nSeqs))
                graphFile.write("\n(mclmatrix\nbegin\n\n")
            pool = mp.Pool(nProcess)
            pool.map(func, [(seqsInfo, graphFN, iSpec, files.FileHandler.GetPickleDir()) for iSpec in
                                             range(seqsInfo.nSpecies)])
            for iSp in range(seqsInfo.nSpecies):
                subprocess.call("cat " + graphFN + "_%d" % iSp + " >> " + graphFN, shell=True)
                os.remove(graphFN + "_%d" % iSp)
            # Cleanup
            pool.close()
            matrices.DeleteMatrices("B", files.FileHandler.GetPickleDir())
            matrices.DeleteMatrices("connect", files.FileHandler.GetPickleDir())
        return graphFN

    @staticmethod
    def GetMostDistant_s(RBH, B, seqsInfo, iSpec):
        # most distant RBB - as cut-off for connecting to other genes
        mostDistant = numeric.transpose(np.ones(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]) * 1e9)
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
        mostDistant[I] = bestHit[I] + 1e-6  # to connect to one in it's own species it must be closer than other species. We can deal with hits outside the species later
        return mostDistant


    @staticmethod
    def GetMostDistant_s_estimate_from_relative_rbhs(RBH, B, seqsInfo, iSpec, p=90):
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
        for isp in range(nsp_m1):
            rs.append([])
            for jsp in range(nsp_m1):
                if isp == jsp:
                    rs[-1].append(1.0)
                    continue
                ratios = np.multiply(Z[isp,:], Zr[jsp,:])
                i_nonzeros = np.where(np.logical_and(np.isfinite(ratios), ratios > 0))  # only those which a score is availabe for both
                if i_nonzeros[0].size == 0:
                    rs[-1].append(1.)  # no conversion between this pair
                else:
                    ratios = ratios[i_nonzeros]
                    r = np.percentile(np.asarray(ratios), 100-p)
                    rs[-1].append(r)
                # calculate the n vectors and take the min
        C = np.matrix(rs)   # Conversion matrix:
        # To convert from a hit in species j to one in species i multiply by Cij
        # To convert from a hit in species j to the furthest hit, need to take
        # the column-wise minimum
        C = np.amin(C, axis=0)    # conversion from a RBH to the estmate of what the most distant RBH should be
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
            M = B[jSpec].tolil()
            if iSpec != jSpec:
                IIJJ = [(i, j) for i, (valueRow, indexRow) in enumerate(zip(M.data, M.rows)) for j, v in
                        zip(indexRow, valueRow) if v >= mostDistant[i]]
            else:
                IIJJ = [(i, j) for i, (valueRow, indexRow) in enumerate(zip(M.data, M.rows)) for j, v in
                        zip(indexRow, valueRow) if (i != j) and v >= mostDistant[i]]
            II = [i for (i, j) in IIJJ]
            JJ = [j for (i, j) in IIJJ]
            onesArray = np.ones(len(IIJJ))
            mat = sparse.csr_matrix((onesArray, (II, JJ)),
                                    shape=(nSeqs_i, seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[jSpec]]))
            connect.append(mat)
        return connect

    @staticmethod
    def ConnectAllBetterThanAnOrtholog_s(RBH, B, seqsInfo, iSpec, v2_scores):
        if v2_scores:
            mostDistant = WaterfallMethod.GetMostDistant_s_estimate_from_relative_rbhs(RBH, B, seqsInfo, iSpec)
        else:
            mostDistant = WaterfallMethod.GetMostDistant_s(RBH, B, seqsInfo, iSpec)
        connect = WaterfallMethod.ConnectAllBetterThanCutoff_s(B, mostDistant, seqsInfo, iSpec)
        return connect


def BitScore(sequence):
    pass

def GetMaxBitscores(seqsInfo):
    bit_scores = []
    for iSpecies, iFasta in enumerate(seqsInfo.speciesToUse):
        bit_scores.append(np.zeros(seqsInfo.nSeqsPerSpecies[iFasta]))
        fastaFilename = files.FileHandler.GetSpeciesFastaFN(iFasta)
        current_sequence = ""
        qFirstLine = True
        with open(fastaFilename) as infile:
            for row in infile:
                if len(row) > 1 and row[0] == ">":
                    if qFirstLine:
                        qFirstLine = False
                    else:
                        bit_scores[iSpecies][iCurrentSequence] = BitScore(current_sequence)
                        current_sequence = ""
                    _, iCurrentSequence = util.GetIDPairFromString(row[1:])
                else:
                    current_sequence += row.rstrip()
        bit_scores[iSpecies][iCurrentSequence] = BitScore(current_sequence)
    return bit_scores

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


def DoOrthogroups(options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML=None, i_unassigned=None):
    # Run Algorithm, cluster and output cluster files with original accessions
    q_unassigned = i_unassigned is not None
    util.PrintUnderline("Running OrthoFinder algorithm" + (" for clade-specific genes" if q_unassigned else ""))
    # it's important to free up the memory from python used for processing the genomes
    # before launching MCL because both use sizeable amounts of memory. The only
    # way I can find to do this is to launch the memory intensive python code
    # as separate process that exits before MCL is launched.

    Lengths = GetSequenceLengths(seqsInfo)  # Alternatively, self-self bit scores, but it amounts to the same thing

    # Process BLAST hits
    util.PrintTime("Initial processing of each species")
    cmd_queue = mp.Queue()
    blastDir_list = files.FileHandler.GetBlastResultsDir()
    if q_unassigned:
        blastDir_list = blastDir_list[:1]  # only use latet directory with unassigned gene searches
    for iSpeciesJob in range(seqsInfo.nSpecies):  # The i-th job, not the OrthoFinder species ID
        cmd_queue.put(iSpeciesJob)
    files.FileHandler.GetPickleDir()  # create the pickle directory before the parallel processing to prevent a race condition
    # Should use PTM?
    # args_list = [(seqsInfo, blastDir_list, Lengths, cmd_queue, files.FileHandler.GetPickleDir(), options.qDoubleBlast, options.v2_scores, q_unassigned)
    #              for i_ in range(options.nProcessAlg)]
    # parallel_task_manager.RunParallelMethods(WaterfallMethod.Worker_ProcessBlastHits, args_list, options.nProcessAlg)
    runningProcesses = [mp.Process(target=WaterfallMethod.Worker_ProcessBlastHits,
                                   args=(seqsInfo, blastDir_list, Lengths, cmd_queue, files.FileHandler.GetPickleDir(), options.qDoubleBlast, options.v2_scores, q_unassigned))
                        for i_ in range(options.nProcessAlg)]
    for proc in runningProcesses:
        proc.start()
    parallel_task_manager.ManageQueue(runningProcesses, cmd_queue)

    if options.gathering_version < (3, 0):
        cmd_queue = mp.Queue()
        for iSpecies in range(seqsInfo.nSpecies):
            cmd_queue.put((seqsInfo, iSpecies))
        # args_list = [(cmd_queue, files.FileHandler.GetPickleDir(), options.v2_scores) for i_ in range(options.nProcessAlg)]
        # parallel_task_manager.RunParallelMethods(WaterfallMethod.Worker_ConnectCognates, args_list, options.nProcessAlg)
        runningProcesses = [
            mp.Process(target=WaterfallMethod.Worker_ConnectCognates, args=(cmd_queue, files.FileHandler.GetPickleDir(), options.v2_scores))
            for i_ in range(options.nProcessAlg)]
        for proc in runningProcesses:
            proc.start()
        parallel_task_manager.ManageQueue(runningProcesses, cmd_queue)

        util.PrintTime("Connected putative homologues")
        graphFilename = WaterfallMethod.WriteGraphParallel(seqsInfo, options.nProcessAlg, i_unassigned)

        # 5b. MCL
        clustersFilename, clustersFilename_pairs = files.FileHandler.CreateUnusedClustersFN(
            "_I%0.1f" % options.mclInflation, i_unassigned)
        mcl.MCL.RunMCL(graphFilename, clustersFilename, options.nProcessAlg, options.mclInflation)
        # If processing unassigned, then ignore all 'unclustered' genes - they will include any genes not included in this search
        mcl.ConvertSingleIDsToIDPair(seqsInfo, clustersFilename, clustersFilename_pairs, q_unassigned)
    elif options.gathering_version == (3, 2):
        graphFilename = WaterfallMethod.WriteGraphParallel(seqsInfo, options.nProcessAlg, i_unassigned, WriteGraph_perSpecies_homology)
        clustersFilename, clustersFilename_pairs = files.FileHandler.CreateUnusedClustersFN("_I%0.1f" % options.mclInflation, i_unassigned)
        mcl.MCL.RunMCL(graphFilename, clustersFilename, options.nProcessAlg, options.mclInflation)
        mcl.ConvertSingleIDsToIDPair(seqsInfo, clustersFilename, clustersFilename_pairs, q_unassigned)
    if not q_unassigned:
        post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo, speciesNamesDict, options, speciesXML)
    return clustersFilename_pairs


def post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo, speciesNamesDict, options, speciesXML, q_incremental=False):
    """
    Write OGs & statistics to results files, write Fasta files.
    Args:
        q_incremental - These are not the final orthogroups, don't write results
    """
    ogs = mcl.GetPredictedOGs(clustersFilename_pairs)
    resultsBaseFilename = files.FileHandler.GetOrthogroupResultsFNBase()
    util.PrintUnderline("Writing orthogroups to file")
    idsDict = mcl.MCL.WriteOrthogroupFiles(ogs, [files.FileHandler.GetSequenceIDsFN()], resultsBaseFilename,
                                           clustersFilename_pairs)
    if not q_incremental:
        mcl.MCL.CreateOrthogroupTable(ogs, idsDict, speciesNamesDict, speciesInfoObj.speciesToUse, resultsBaseFilename)

    # Write Orthogroup FASTA files
    ogSet = orthologues.OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse,
                                       speciesInfoObj.nSpAll, options.qAddSpeciesToIDs,
                                       idExtractor=util.FirstWordExtractor)
    treeGen = trees_msa.TreesForOrthogroups(None, None, None)
    fastaWriter = trees_msa.FastaWriter(files.FileHandler.GetSpeciesSeqsDir(), speciesInfoObj.speciesToUse)
    d_seqs = files.FileHandler.GetResultsSeqsDir()
    if not os.path.exists(d_seqs): os.mkdir(d_seqs)
    treeGen.WriteFastaFiles(fastaWriter, ogSet.OGsAll(), idsDict, False)

    if not q_incremental:
        # stats.Stats(ogs, speciesNamesDict, speciesInfoObj.speciesToUse, files.FileHandler.iResultsVersion)
        if options.speciesXMLInfoFN:
            mcl.MCL.WriteOrthoXML(speciesXML, ogs, seqsInfo.nSeqsPerSpecies, idsDict, resultsBaseFilename + ".orthoxml",
                                  speciesInfoObj.speciesToUse)
        print("")
        util.PrintTime("Done orthogroups")
        files.FileHandler.LogOGs()