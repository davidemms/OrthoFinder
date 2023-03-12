import csv
import glob
import gzip
import os
import string
import random
import subprocess
from collections import defaultdict
from typing import Optional
import multiprocessing as mp

import numpy as np

from . import util, files, parallel_task_manager, mcl, orthologues, sample_genes, fasta_writer


class XcelerateConfig(object):
    def __init__(self):
        # use kmeans clustering to select 1/kmeans_divide genes per orthogroup (to a minimum of 20). Otherwise, use all.
        self.kmeans_divide: Optional[int] = 10  # 10 is a good value if using this option


xcelerate_config = XcelerateConfig()


def check_for_orthoxcelerate(input_dir):
    # Fill in details as developed
    return True


def prepare_accelerate_database(input_dir, wd_list, nSpAll):
    if xcelerate_config.kmeans_divide is None:
        # create_shoot_db.create_full_database(input_dir, q_ids=True, subtrees_dir="")
        fn_diamond_db = create_profiles_database(input_dir, wd_list, nSpAll, selection="all", q_ids=True, subtrees_dir="")
    else:
        fn_diamond_db = create_profiles_database(input_dir, wd_list, nSpAll, selection="kmeans", q_ids=True, divide=xcelerate_config.kmeans_divide, subtrees_dir="")
    return fn_diamond_db


def RunSearch(options, speciessInfoObj, fn_diamond_db, prog_caller):
    name_to_print = "BLAST" if options.search_program == "blast" else options.search_program
    util.PrintUnderline("Running %s profiles search" % name_to_print)
    commands, results_files = GetOrderedSearchCommands(speciessInfoObj, fn_diamond_db, prog_caller)
    if options.qStopAfterPrepare:
        for command in commands:
            print(command)
        util.Success()
    print("Using %d thread(s)" % options.nBlast)
    cmd_queue = mp.Queue()
    for iCmd, cmd in enumerate(commands):
        cmd_queue.put((iCmd+1, cmd))
    runningProcesses = [mp.Process(target=parallel_task_manager.Worker_RunCommand, args=(cmd_queue, options.nBlast, len(commands), True)) for _ in range(options.nBlast)]
    for proc in runningProcesses:
        proc.start()
    for proc in runningProcesses:
        while proc.is_alive():
            proc.join()
    util.PrintTime("Done profiles search")
    return results_files


def GetOrderedSearchCommands(speciesInfoObj, diamond_db, prog_caller, search_program="diamond"):
    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing
    the BLAST commands.
    """
    iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    commands = [prog_caller.GetSearchMethodCommand_Search(
        search_program,
        files.FileHandler.GetSpeciesFastaFN(iFasta),
        diamond_db,
        files.FileHandler.GetBlastResultsFN(iFasta, 0, qForCreation=True))
        for iFasta in iSpeciesNew
    ]
    results_files = [files.FileHandler.GetBlastResultsFN(iFasta, 0, qForCreation=True) + ".gz" for iFasta in iSpeciesNew]
    return commands, results_files


def ogs_from_diamond_results(fn_og_results_out, q_ignore_sub=False):
    """
    Get the OG based on the DIAMOND results
    Args:
        fn_og_results_out - fn of compressed DIAMOND results
        q_ignore_sub - ignore any subtrees and just look at overall OGs
    Returns:
        iog        : List[str] "int.int" or "int" of ordered, ambiguous assignments
    Info:
        Hits to other groups are returned if np.log10 difference is less than 10
        and -np.log10 score > difference.
    """
    ogs_all_genes = defaultdict(list)
    scores_all_genes = defaultdict(list)
    with gzip.open(fn_og_results_out, 'rt') as infile:
        reader = csv.reader(infile, delimiter="\t")
        for line in reader:
            gene = line[0]
            og = line[1].split("_")[0]
            ogs_all_genes[gene].append(og.split(".", 1)[0] if q_ignore_sub else og)
            scores_all_genes[gene].append(float(line[-2]))
    all_genes = list(ogs_all_genes.keys())
    for gene in all_genes:
        ogs = ogs_all_genes[gene]
        scores = scores_all_genes[gene]
        sortedTuples = sorted(zip(scores, ogs))
        scores = [i for i, j in sortedTuples]
        ogs = [j for i, j in sortedTuples]
        unique = [i for i in range(len(ogs)) if ogs[i] not in ogs[:i]]
        ogs = [ogs[i] for i in unique]
        scores = [scores[i] for i in unique]
        if len(ogs) > 0:
            # filter out all ogs other than those which are potentially worth considering
            sliver = np.nextafter(0, 1)
            scores_ml10 = [-np.log10(s+sliver) for s in scores]
            s0 = scores_ml10[0]
            # in a test of 15k sequences only 12 passed the first test but failed s>(s0-s)
            # it is not worth arguing over whether it's a good second criteria
            # scores = [s for s in scores if s0-s<10 and s>(s0-s)]
            scores_ml10 = [s for s in scores_ml10 if s0-s<10]
            # ogs_all_genes[gene] = ogs[:len(scores_ml10)]
            # scores_all_genes[gene] = scores[:len(scores_ml10)]
            ogs_all_genes[gene] = ogs[0]
            # scores_all_genes[gene] = scores[0]
    return ogs_all_genes


def assign_genes(results_files):
    wd_input_clusters = files.FileHandler.GetWorkingDirectory1_Read()[1]  # first is the current working directory, second is the most recent of the input directories
    fn_clusters = glob.glob(wd_input_clusters + "clusters*_id_pairs.txt")
    if len(fn_clusters) == 0:
        print("ERROR: Couldn't find previous orthogroups in %s" % wd_input_clusters)
        util.Fail()
    elif len(fn_clusters) > 1:
        print("WARNING: Found multiple orthogroup files: %s" % fn_clusters)
    fn_clusters = fn_clusters[0]
    ogs = mcl.GetPredictedOGs(fn_clusters)
    for fn in results_files:
        ogs_all_genes = ogs_from_diamond_results(fn)
        for gene, og in ogs_all_genes.items():
            ogs[int(og)].add(gene)
    clustersFilename, clustersFilename_pairs = files.FileHandler.CreateUnusedClustersFN()
    mcl.write_updated_clusters_file(ogs, clustersFilename_pairs)
    return clustersFilename_pairs


def create_profiles_database(din, wd_list, nSpAll, selection="kmeans", min_for_profile=20, q_ids=True, divide=10, subtrees_dir=""):
    """
    Create a fasta file with profile genes from each orthogroup
    Args:
        din - Input OrthoFinder results directory
        selection - "kmeans|random|all"  - method to use for sampling genes from orthogroups
        min_for_profile - The min number of genes to use from each orthogroup, when available
        divide - Select `|genes| / divide` genes per orthogroup, as long as the number is greater than min_for_profile
        q_ids - Convert subtrees (with user gene accessions) back to IDs for profiles database
    Notes:
    If the trees have been split into subtrees then profiles will be created for
    the subtrees instead.
    """
    if selection not in ("kmeans", "random", "all"):
        raise RuntimeError("selection method '%s' is not defined" % selection)
    wd = din + "WorkingDirectory/"
    if subtrees_dir:
        subtrees_label = "." + os.path.split(subtrees_dir)[1]
        pat_super = din + subtrees_dir + "/super/OG%07d.super.tre"
        pat_sub_msa_glob = din + subtrees_dir + "/msa_sub/OG%07d.*.fa"
    else:
        subtrees_label = ""
    ogs = mcl.GetPredictedOGs(wd + "clusters_OrthoFinder_I1.5.txt_id_pairs.txt")
    fw = fasta_writer.FastaWriter(wd + "Species*fa", qGlob=True)
    seq_write = []
    seq_convert = dict()
    # print("WARNING: Check all gene names, can't start with '__'")
    # If there are subtrees then we need to convert their IDs in the profile file
    # back to internal IDs
    if q_ids:
        og_set = orthologues.OrthoGroupsSet(wd_list, list(range(nSpAll)), nSpAll, True)
        ids = og_set.Spec_SeqDict()
        ids_rev = {v: k for k, v in ids.items()}
    nToDo = len(ogs)
    for iog, og in enumerate(ogs):
        if iog >= 0 and divmod(iog, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
            print("Done %d of %d" % (iog, nToDo))
        og_id = "%07d" % iog
        q_subtrees = subtrees_dir and os.path.exists(pat_super % iog)
        fn_msa = wd + "Alignments_ids/OG%07d.fa" % iog
        if q_subtrees:
            print("Subtrees: %d" % iog)
            fns_msa = list(glob.glob(pat_sub_msa_glob % iog))
        elif os.path.exists(fn_msa):
            fns_msa = [fn_msa, ]
        else:
            fns_msa = [wd + "Sequences_ids/OG%07d.fa" % iog, ]
        for fn in fns_msa:
            if not os.path.exists(fn):
                print("File does not exist, skipping: %s" % fn)
                continue
            i_part = os.path.basename(fn).rsplit(".", 2)[1] if q_subtrees else None
            fw_temp = fasta_writer.FastaWriter(fn)
            n_in_og = len(fw_temp.SeqLists)
            n_for_profile = n_in_og // divide
            if n_for_profile < min_for_profile:
                n_for_profile = min_for_profile
            if selection == "kmeans" and len(fw_temp.SeqLists) > n_for_profile:
                if q_subtrees:
                    # MSA needs to be modified
                    letters = string.ascii_lowercase
                    fn_temp = "/tmp/shoot_db_create" + "".join(random.choice(letters) for i in range(6)) + os.path.basename(fn)
                    fw_temp.WriteSeqsToFasta([g for g in fw_temp.SeqLists if not g.startswith("SHOOTOUTGROUP_")],
                                            fn_temp)
                    fn = fn_temp
                # Don't trim as OrthoFinder has already trimmed by default
                s = sample_genes.select_from_aligned(fn, n_for_profile, q_trim=False)
                if q_subtrees:
                    os.remove(fn_temp)
            elif selection == "kmeans" or selection == "random":
                if q_subtrees:
                    fw_temp = fasta_writer.FastaWriter(fn)
                    og = [g for g in fw_temp.SeqLists if not g.startswith("SHOOTOUTGROUP_")]
                s = sample_random(og, n_for_profile)
            else:
                s = [g for g in fasta_writer.FastaWriter(fn).SeqLists.keys() if not g.startswith("SHOOTOUTGROUP_")]
            if q_ids and q_subtrees:
                s = [ids_rev[ss] for ss in s]
            if q_subtrees:
                og_id_full = og_id + "." + i_part
            else:
                og_id_full = og_id
            seq_write.extend(s)
            for ss in s:
                seq_convert[ss] = og_id_full + "_" + ss
    if selection == "all":
        fn_fasta = din + "profile_sequences%s.all.fa" % subtrees_label
    else:
        fn_fasta = din + "profile_sequences%s.d%d_min%d_%s.fa" % (subtrees_label, divide, min_for_profile, selection)
    fw.WriteSeqsToFasta_withNewAccessions(seq_write, fn_fasta, seq_convert)
    fn_diamond_db = fn_fasta + ".db"
    subprocess.call(["diamond", "makedb", "--in", fn_fasta, "-d", fn_diamond_db])
    return fn_diamond_db


def sample_random(og, n_max):
    """
    Sample min(n_max, |og|) genes randomly from clade
    Args:
        og - set of strings
        n_max - max number of genes to sample
    Returns:
        genes - list of genes
    """
    return random.sample(og, min(n_max, len(og)))
