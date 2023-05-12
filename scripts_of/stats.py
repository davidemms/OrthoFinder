from __future__ import absolute_import

import os
import csv
import shutil
import datetime
import numpy as np
from collections import Counter

from . import util, files


def OrthogroupsMatrix(iSpecies, properOGs):
    speciesIndexDict = {iSp: iCol for iCol, iSp in enumerate(iSpecies)}
    nSpecies = len(iSpecies)
    nGroups = len(properOGs)
    # (i, j)-th entry of ogMatrix gives the number of genes from i in orthogroup j
    ogMatrix = np.zeros((nGroups, nSpecies))
    for i_og, og in enumerate(properOGs):
        for species, _ in og:
            ogMatrix[i_og, speciesIndexDict[species]] += 1
    return ogMatrix


def Stats_SpeciesOverlaps(fn, speciesNamesDict, iSpecies, speciesPresence):
    """ Number of orthogroups in which each species-pair is present. Called by Stats"""
    with open(fn, util.csv_write_mode) as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow([""] + [speciesNamesDict[index] for index in iSpecies])
        for iSp in iSpecies:
            overlap = [len([1 for og in speciesPresence if (iSp in og and jSp in og)]) for jSp in iSpecies]
            writer.writerow([speciesNamesDict[iSp]] + overlap)


def Stats_SizeTable(writer_sum, writer_sp, properOGs, allGenesCounter, iSpecies, speciesPresence):
    """ Overall and per-species histogram tables of orthogroup sizes. Called by Stats"""
    bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 21, 51, 101, 151, 201, 501, 1001, 9e99]
    writer_sum.writerow([])
    writer_sp.writerow([])
    nGenesPerOG = [len(og) for og in properOGs]
    nGenesPerOGPerSpecies = [[len([1 for g in og if g[0] == iSp]) for og in properOGs] for iSp in iSpecies]
    counters_GenesPerOGPerSpecies = [Counter(spCount) for spCount in nGenesPerOGPerSpecies]
    counters_GenesPerOG = Counter(nGenesPerOG)
    nSp = len(iSpecies)
    nOGs = len(properOGs)
    nGenesTotal = sum([len(og) for og in properOGs])
    percentFormat = "%0.1f"
    writer_sum.writerow(
        ["Average number of genes per-species in orthogroup", "Number of orthogroups", "Percentage of orthogroups",
         "Number of genes", "Percentage of genes"])
    # Per-species tables
    table_NO = [["Number of genes per-species in orthogroup"] + ["Number of orthogroups" for _ in
                                                                 iSpecies]]  # number of orthogroups
    table_PO = [["Number of genes per-species in orthogroup"] + ["Percentage of orthogroups" for _ in
                                                                 iSpecies]]  # percentage of orthogroups
    table_NG = [
        ["Number of genes per-species in orthogroup"] + ["Number of genes" for _ in iSpecies]]  # Number of genes
    table_PG = [["Number of genes per-species in orthogroup"] + ["Percentage of genes" for _ in
                                                                 iSpecies]]  # percentage of genes
    for start, end in zip(bins, bins[1:]):
        binName = "<1" if start == 0 else (
                    "'%d" % start) if start + 1 == end else "'%d+" % start if end == 9e99 else "%d-%d" % (
        start, end - 1)
        nOrthogroups = sum([count for size, count in counters_GenesPerOG.items() if start * nSp <= size < end * nSp])
        nGenes = sum([size * count for size, count in counters_GenesPerOG.items() if start * nSp <= size < end * nSp])
        row_sum = [binName, nOrthogroups, percentFormat % (100. * nOrthogroups / nOGs), nGenes,
                   percentFormat % (100. * nGenes / nGenesTotal)]
        writer_sum.writerow(row_sum)
        # Per-species stats
        if binName == "<1": binName = "'0"
        nOrthogroups_ps = [sum([number for size, number in c.items() if start <= size < end]) for c in
                           counters_GenesPerOGPerSpecies]
        nGenes_ps = [sum([size * number for size, number in c.items() if start <= size < end]) for c in
                     counters_GenesPerOGPerSpecies]
        table_NO.append([binName] + nOrthogroups_ps)
        table_PO.append([binName] + [percentFormat % (100. * n / nOGs) for n in nOrthogroups_ps])
        table_NG.append([binName] + nGenes_ps)
        table_PG.append(
            [binName] + [percentFormat % (100. * n / allGenesCounter[iSp]) for iSp, n in zip(iSpecies, nGenes_ps)])
    writer_sp.writerow([])
    for r in table_NO: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_PO: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_NG: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_PG: writer_sp.writerow(r)

    # Species presence
    n = list(map(len, speciesPresence))
    writer_sum.writerow([])
    writer_sum.writerow(["Number of species in orthogroup", "Number of orthogroups"])
    for i in range(1, nSp + 1):
        writer_sum.writerow([i, n.count(i)])


def Stats(ogs, speciesNamesDict, iSpecies, iResultsVersion, fastaWriter, ids_dict, q_fast_add=False):
    """ Top-level method for calculation of stats for the orthogroups"""
    allOgs = [[list(map(int, g.split("_"))) for g in og] for og in ogs]
    properOGs = [og for og in allOgs if len(og) > 1]
    iogs_properOGs = [iog for iog, og in enumerate(allOgs) if len(og) > 1]
    allGenes = [g for og in allOgs for g in og]
    ogStatsResultsDir = files.FileHandler.GetOGsStatsResultsDirectory()
    filename_sp = ogStatsResultsDir + "Statistics_PerSpecies" + (
        "" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".tsv"
    filename_sum = ogStatsResultsDir + "Statistics_Overall" + (
        "" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".tsv"
    filename_overlap = ogStatsResultsDir + "Orthogroups_SpeciesOverlaps" + (
        "" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".tsv"
    filename_single_copy = files.FileHandler.GetOrthogroupResultsFNBase() + "_SingleCopyOrthologues.txt"
    percentFormat = "%0.1f"
    with open(filename_sp, util.csv_write_mode) as outfile_species, open(filename_sum, util.csv_write_mode) as outfile_sum:
        writer_sp = csv.writer(outfile_species, delimiter="\t")
        writer_sum = csv.writer(outfile_sum, delimiter="\t")
        # header
        writer_sp.writerow([""] + [speciesNamesDict[index] for index in iSpecies])

        # Number of genes
        allGenesCounter = Counter([g[0] for g in allGenes])
        nGenes = sum(allGenesCounter.values())
        writer_sum.writerow(["Number of species", len(iSpecies)])
        writer_sp.writerow(["Number of genes"] + [allGenesCounter[iSp] for iSp in iSpecies])
        writer_sum.writerow(["Number of genes", nGenes])

        # Number of assigned/unassigned genes
        assignedGenesCounter = Counter([g[0] for og in properOGs for g in og])
        nAssigned = sum(assignedGenesCounter.values())
        writer_sp.writerow(["Number of genes in orthogroups"] + [assignedGenesCounter[iSp] for iSp in iSpecies])
        writer_sum.writerow(["Number of genes in orthogroups"] + [nAssigned])
        writer_sp.writerow(
            ["Number of unassigned genes"] + [allGenesCounter[iSp] - assignedGenesCounter[iSp] for iSp in iSpecies])
        writer_sum.writerow(["Number of unassigned genes"] + [nGenes - nAssigned])
        # Percentages
        pAssigned = 100. * nAssigned / nGenes
        writer_sp.writerow(["Percentage of genes in orthogroups"] + [
            percentFormat % (100. * assignedGenesCounter[iSp] / allGenesCounter[iSp]) for iSp in iSpecies])
        writer_sum.writerow(["Percentage of genes in orthogroups", percentFormat % pAssigned])
        writer_sp.writerow(["Percentage of unassigned genes"] + [
            percentFormat % (100 * (1. - (float(assignedGenesCounter[iSp]) / allGenesCounter[iSp]))) for iSp in
            iSpecies])
        writer_sum.writerow(
            ["Percentage of unassigned genes", percentFormat % (100 * (1. - (float(nAssigned) / nGenes)))])

        # Number of Orthogroups
        speciesPresence = [set([g[0] for g in og]) for og in properOGs]
        nOgs = len(properOGs)
        writer_sum.writerow(["Number of orthogroups", nOgs])
        writer_sp.writerow(
            ["Number of orthogroups containing species"] + [sum([iSp in og_sp for og_sp in speciesPresence]) for iSp in
                                                            iSpecies])
        writer_sp.writerow(["Percentage of orthogroups containing species"] + [percentFormat % (
            (100. * sum([iSp in og_sp for og_sp in speciesPresence]) / len(properOGs)) if len(properOGs) > 0 else 0.)
                                                                               for iSp in iSpecies])

        # Species specific orthogroups - orthogroups-based
        speciesSpecificOGsCounter = Counter([next(iter(og_sp)) for og_sp in speciesPresence if len(og_sp) == 1])
        writer_sp.writerow(
            ["Number of species-specific orthogroups"] + [speciesSpecificOGsCounter[iSp] for iSp in iSpecies])
        writer_sum.writerow(["Number of species-specific orthogroups", sum(speciesSpecificOGsCounter.values())])

        # Species specific orthogroups - gene-based
        iSpeciesSpecificOGs = [i for i, og_sp in enumerate(speciesPresence) if len(og_sp) == 1]
        iSpSpecificOGsGeneCounts = [
            sum([len(properOGs[iog]) for iog in iSpeciesSpecificOGs if properOGs[iog][0][0] == iSp]) for iSp in
            iSpecies]
        writer_sp.writerow(["Number of genes in species-specific orthogroups"] + iSpSpecificOGsGeneCounts)
        writer_sum.writerow(["Number of genes in species-specific orthogroups", sum(iSpSpecificOGsGeneCounts)])
        writer_sp.writerow(["Percentage of genes in species-specific orthogroups"] + [
            percentFormat % (100. * n_ss / allGenesCounter[iSp]) for n_ss, iSp in
            zip(iSpSpecificOGsGeneCounts, iSpecies)])
        writer_sum.writerow(["Percentage of genes in species-specific orthogroups",
                             percentFormat % (100. * sum(iSpSpecificOGsGeneCounts) / nGenes)])

        # 'averages'
        l = list(sorted(list(map(len, properOGs))))
        writer_sum.writerow(["Mean orthogroup size", "%0.1f" % np.mean(l)])
        writer_sum.writerow(["Median orthogroup size", np.median(l)])
        L = np.cumsum(l)
        j, _ = next((i, x) for i, x in enumerate(L) if x > nAssigned / 2)
        writer_sum.writerow(["G50 (assigned genes)", l[j]])
        l2 = list(reversed(list(map(len, ogs))))
        L2 = np.cumsum(l2)
        j2, _ = next((i, x) for i, x in enumerate(L2) if x > nGenes / 2)
        G50 = l2[j2]
        writer_sum.writerow(["G50 (all genes)", G50])
        writer_sum.writerow(["O50 (assigned genes)", len(l) - j])
        O50 = len(l2) - j2
        writer_sum.writerow(["O50 (all genes)", O50])

        # Single-copy orthogroups
        ogMatrix = OrthogroupsMatrix(iSpecies, properOGs)  # use iogs_properOGs for indexing
        nSpecies = len(iSpecies)
        nPresent = (ogMatrix > np.zeros((1, nSpecies))).sum(1)
        nCompleteOGs = list(nPresent).count(nSpecies)
        singleCopyOGs = (ogMatrix == np.ones((1, nSpecies))).all(1).nonzero()[0]
        nSingleCopy = len(singleCopyOGs)
        writer_sum.writerow(["Number of orthogroups with all species present", nCompleteOGs])
        writer_sum.writerow(["Number of single-copy orthogroups", nSingleCopy])
        with open(filename_single_copy, 'w') as outfile_singlecopy:
            outfile_singlecopy.write("\n".join(["N0.HOG%07d" % iogs_properOGs[i_] for i_ in singleCopyOGs]))
        # Link single-copy orthologues
        g_fmt = files.FileHandler.GetResultsSeqsDir_SingleCopy() +  "N0.HOG%07d.fa"
        for i in singleCopyOGs:
            out_fn = g_fmt % iogs_properOGs[i]
            og = ["%d_%d" % tuple(g) for g in properOGs[i]]
            fastaWriter.WriteSeqsToFasta_withNewAccessions(og, out_fn, ids_dict)

        # Results filenames
        writer_sum.writerow(["Date", str(datetime.datetime.now()).split()[0]])
        writer_sum.writerow(
            ["Orthogroups file", "Orthogroups" + ("" if iResultsVersion == 0 else "_%d" % iResultsVersion) + ".tsv"])
        writer_sum.writerow(["Unassigned genes file", "Orthogroups" + (
            "" if iResultsVersion == 0 else "_%d" % iResultsVersion) + "_UnassignedGenes.tsv"])
        writer_sum.writerow(["Per-species statistics", os.path.split(filename_sp)[1]])
        writer_sum.writerow(["Overall statistics", os.path.split(filename_sum)[1]])
        writer_sum.writerow(["Orthogroups shared between species", os.path.split(filename_overlap)[1]])

        # Sizes
        Stats_SizeTable(writer_sum, writer_sp, properOGs, allGenesCounter, iSpecies, speciesPresence)
        Stats_SpeciesOverlaps(filename_overlap, speciesNamesDict, iSpecies, speciesPresence)

    summaryText = """OrthoFinder assigned %d genes (%0.1f%% of total) to %d orthogroups. Fifty percent of genes were in orthogroups with %d or more genes (G50 was %d) and were contained in the largest %d orthogroups (O50 was %d). There were %d orthogroups with all species present and %d of these consisted entirely of single-copy genes.""" % (
    nAssigned, pAssigned, nOgs, G50, G50, O50, O50, nCompleteOGs, nSingleCopy)
    if q_fast_add:
        print("The majority of genes have been assigned to existing orthogroups, however, the remaining clade-specific genes not seen in the core species were also analysed with the following results:")
    print(summaryText)


def add_unassigned_genes(ogs, all_seq_ids):
    """Extend OGs with unassigned genes as singletons"""
    all_assigned = set([g for og in ogs for g in og])
    unassigned = set(all_seq_ids).difference(all_assigned)
    ogs.extend([{g,} for g in unassigned])
    return ogs

