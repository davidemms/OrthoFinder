from __future__ import absolute_import
from . import parallel_task_manager

import os
import numpy as np
import subprocess
from scipy import sparse
import time
import warnings
from collections import Counter
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
    species_similarities = matrices.LoadMatrix("species_similarities", iSpec, 0, d_pickle)
    with open(graphFN + "_%d" % iSpec, 'w') as graphFile:
        connect2 = []
        for jSpec in range(seqsInfo.nSpecies):
            m1 = matrices.LoadMatrix("connect", iSpec, jSpec, d_pickle)
            m2tr = numeric.transpose(matrices.LoadMatrix("connect", jSpec, iSpec, d_pickle))
            connect2.append(m1 + m2tr)
            del m1, m2tr
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpec, d_pickle)
        B_connect = matrices.MatricesAnd_s(connect2, B)
        # normalise the scores with the mean species-pair similarities
        n = [1.0/s for s in species_similarities]
        B_connect = [BB.multiply(nn) for BB,nn in zip(B_connect, n)]
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
    def GetD0(Bii, li):
        """
        Gets the bitscore/length for an identical sequence i.e. distance zero (~zero substitutions per site)
        Args:
            Bii - lil matrix of the bit-scores for hits self-self hits 
            li - the numpy array of corresponding sequence lengths
        """
        S = Bii.diagonal()
        d0 = S/li
        # If there is no self-self hit replace with usual value for a self self
        d0[d0 < 0.5] = 2.  
        return d0

    @staticmethod
    def EstimateReciprocalDistances(Bij, d0, lj):
        """
        Estimate an approximation proportional to 1/(substitutions per site), which is itself a proxy for evolutionary distance
        Args:
            Bij - lil matrix of the bit-scores for hits between iSpecies & jSpecies 
            d0 -
            lj - 
        Returns:
            Dij - 1/(the estimated distance measure) - since with sparse matrixes it is easier to use similarity scores rather than distance matrixes
        Implementation:
            From a branch-length (substitutions per site) between approximately 0 and 1.5 there is an approximately (negative) linear
            relationship between bitscore/genelength corresponding to a bitscore/length of about 0.6. Above this it approaches an
            asymptote and this behaviour will need to be more closely modelled. In the initial implementation it is assumed that the
            genes are within this linear regime.
            Take 0 branch length to be self-self bitscore / length (or 2.0 if no self-self hit) so distance is B0/L0 - B1/L1 
        """
        # Diagonal entries for Bii
        q, h = Bij.get_shape()
        rangeh = list(range(h))
        Lj = sparse.csr_matrix((1./lj, (rangeh, rangeh)))
        Dij = Bij*Lj
        drep = np.repeat(d0, np.diff(Dij.indptr))
        Dij.data = 1./(drep - Dij.data)
        return Dij.tolil()  

    @staticmethod
    def EstimateSimilarity1(Bij, lj):
        # Diagonal entries for Bii
        q, h = Bij.get_shape()
        rangeh = list(range(h))
        Lj_recip = sparse.csr_matrix((1./lj, (rangeh, rangeh)))
        Dij = Bij*Lj_recip
        return Dij.tolil() 
            
    @staticmethod
    def ProcessBlastHits(seqsInfo, blastDir_list, Lengths, iSpecies, d_pickle, qDoubleBlast, qNewMeasure=True):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            # process up to the best hits for each species
            Bi = []
            if qNewMeasure:
                version = "2a"
                if version == "0":
                    # The cover-trees version (which has a max distance of approx 2.0)
                    # Do self-self first so can use to calculate the distance measure
                    Bii = blast_file_processor.GetBLAST6Scores(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[iSpecies], qDoubleBlast=qDoubleBlast) 
                    D0 = WaterfallMethod.GetD0(Bii, Lengths[iSpecies])
                    for jSpecies in range(seqsInfo.nSpecies):
                        if jSpecies == iSpecies:
                            Bij = Bii.copy()
                        else:
                            Bij = blast_file_processor.GetBLAST6Scores(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[jSpecies], qDoubleBlast=qDoubleBlast)
                        Dij = WaterfallMethod.EstimateReciprocalDistances(Bij, D0, Lengths[jSpecies])
                        Bi.append(Dij)
                elif version == "1":
                    # B/L, where L is the length of the hit sequence
                    Bii = blast_file_processor.GetBLAST6Scores(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[iSpecies], qDoubleBlast=qDoubleBlast) 
                    for jSpecies in range(seqsInfo.nSpecies):
                        if jSpecies == iSpecies:
                            Bij = Bii.copy()
                        else:
                            Bij = blast_file_processor.GetBLAST6Scores(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[jSpecies], qDoubleBlast=qDoubleBlast)
                        Dij = WaterfallMethod.EstimateSimilarity1(Bij, Lengths[jSpecies])
                        Bi.append(Dij)
                elif version == "2a":
                    # B/l, where l is the length of the alignment
                    Bii = blast_file_processor.GetBLAST6Scores_by_length(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[iSpecies], qDoubleBlast=qDoubleBlast) 
                    for jSpecies in range(seqsInfo.nSpecies):
                        if jSpecies == iSpecies:
                            Bij = Bii.copy()
                        else:
                            Bij = blast_file_processor.GetBLAST6Scores_by_length(seqsInfo, blastDir_list, seqsInfo.speciesToUse[iSpecies], seqsInfo.speciesToUse[jSpecies], qDoubleBlast=qDoubleBlast)
                        Bi.append(Bij)
                else:
                    raise NotImplementedError("Distance measure selected has not been implemented")
            else:
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
    def ConnectCognates(seqsInfo, iSpecies, d_pickle, gather_version): 
        # calculate RBH for species i
        BHix = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle)
        BHxi = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle, row=False)
        RBHi = matrices.MatricesAndTr_s(BHix, BHxi)   # twice as much work as before (only did upper triangular before)
        del BHix, BHxi
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpecies, d_pickle)
        connect, species_similarities = WaterfallMethod.ConnectAllBetterThanAnOrtholog_s(RBHi, B, seqsInfo, iSpecies, gather_version) 
        matrices.DumpMatrixArray("connect", connect, iSpecies, d_pickle)
        matrices.DumpMatrix("species_similarities", species_similarities, iSpecies, 0, d_pickle)
            
    @staticmethod 
    def Worker_ConnectCognates(cmd_queue, d_pickle, gather_version):
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            while True:
                try:
                    args = cmd_queue.get(True, 1)
                    WaterfallMethod.ConnectCognates(*args, d_pickle=d_pickle, gather_version=gather_version)
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
            matrices.DeleteMatrices("B", files.FileHandler.GetPickleDir()) 
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
            species_similarities - the mean species similarity scores (with the max 
                                   value for its own hits - note, could be higher)
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
        # Could just get all pairs (in fact, the ratio)
        # For each gene, what is the most distant RBH
        # For each RBH what factor would we apply to get the RBH
        # 1. Reduce to column vectors
        q_debug_gather = False
        RBH = [rbh.tocsr() for rbh in RBH]
        B = [b.tocsr() for b in B]
        RBH_B = [(rbh.multiply(b)).tolil() for i, (rbh, b) in enumerate(zip(RBH, B)) if i!= iSpec]
        # create a vector of these scores
        nseqi = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]
        # nsp-1 x nseqi
        Z = np.matrix([[rhb_b.data[i][0] if rhb_b.getrowview(i).nnz == 1 else 0. for i in range(nseqi)] for rhb_b in RBH_B])   # RBH if it exists else zero
        # Get the mean similarity scores against each other species
        species_similarities = Z.sum(1)/(Z!=0).sum(1).astype(float)
        # this is a matrix, convert it to a list
        species_similarities = species_similarities.flatten().tolist()[0]
        # add an entry for this species self-self hits
        species_similarities = species_similarities[:iSpec] + [max(species_similarities)] + species_similarities[iSpec:]
        # for each pair of these species we want the entires that are both non zero
        # print(species_similarities)
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
                ratios = ratios[i_nonzeros]   # Dj = r Di  so (1/Di) = r (1/Dj)
                # we will look for genes with 1/D > cut-off
                # instead of 1/Di will use r/Dj
                # so want the value of r such that 95% are bigger than it
                ratios = ratios.reshape(-1,1)
                # import matplotlib.pyplot as plt
                # plt.hist(ratios)
                # plt.show()
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
        # So, if X has an RBH to a gene in species A with score a
        # And if r = 0.63 
        # Then 95% of the time the RBH in species B would be 0.63a
        # We should use the lowest of these estimates as the cut-off for inclusion
        # Now we can estimate the mostDistant hit to accept for each sequence based 
        # on the RBHs seen in any species

        # most distant RBB - as cut-off for connecting to other genes
        mostDistant = numeric.transpose(np.ones(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]])*1e9)
        # now go through the RBHs for each sequence and whenever there are multiple, 
        # calculate the expected distance to the furthest in the orthogroup from each one 
        # and store the largest number

        # do this by going through each of the other species and if one of those has an
        # RBH then estimate from it the cut-off and if it's less than the previous then update it

        # Best hit in another species: species-specific paralogues will now be connected - closer than any gene in any other species
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
        return mostDistant, species_similarities

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
    def ConnectAllBetterThanAnOrtholog_s(RBH, B, seqsInfo, iSpec, gather_version):     
        if gather_version < (3,0):   
            mostDistant = WaterfallMethod.GetMostDistant_s(RBH, B, seqsInfo, iSpec) 
        else:
            mostDistant, species_similarities = WaterfallMethod.GetMostDistant_s_estimate_from_relative_rbhs(RBH, B, seqsInfo, iSpec) 
        connect = WaterfallMethod.ConnectAllBetterThanCutoff_s(B, mostDistant, seqsInfo, iSpec)
        return connect, species_similarities


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