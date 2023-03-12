import csv
import glob
import gzip
from collections import defaultdict
from typing import Optional
import multiprocessing as mp

import numpy as np

from shoot import create_shoot_db

from . import util, files, parallel_task_manager, mcl


class XcelerateConfig(object):
    def __init__(self):
        # use kmeans clustering to select 1/kmeans_divide genes per orthogroup (to a minimum of 20). Otherwise, use all.
        self.kmeans_divide: Optional[int] = 10  # 10 is a good value if using this option


xcelerate_config = XcelerateConfig()


def check_for_orthoxcelerate(input_dir):
    # Fill in details as developed
    return True


def prepare_accelerate_database(input_dir):
    if xcelerate_config.kmeans_divide is None:
        # create_shoot_db.create_full_database(input_dir, q_ids=True, subtrees_dir="")
        fn_diamond_db = create_shoot_db.create_profiles_database(input_dir, selection="all", q_ids=True, subtrees_dir="")
    else:
        fn_diamond_db = create_shoot_db.create_profiles_database(input_dir, selection="kmeans", q_ids=True, divide=xcelerate_config.kmeans_divide, subtrees_dir="")
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

