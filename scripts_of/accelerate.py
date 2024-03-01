import csv
import glob
import gzip
import os
import string
import random
from collections import defaultdict
from typing import Optional
import multiprocessing as mp

import numpy as np

from . import util, files, parallel_task_manager, mcl, orthologues, sample_genes, fasta_writer, tree, program_caller


class XcelerateConfig(object):
    def __init__(self):
        self.n_for_profiles: Optional[int] = 10  # 10 is a good value if using this option


xcelerate_config = XcelerateConfig()


def check_for_orthoxcelerate(input_dir, speciesInfoObj):
    # Add any specific checks required here
    if speciesInfoObj.speciesToUse != list(range(speciesInfoObj.nSpAll)):
        print("ERROR: Removing species from 'core' results directory is not supported for an --assign analysis.")
        return False
    return True


def prepare_accelerate_database(input_dir, wd_list, nSpAll):
    if xcelerate_config.n_for_profiles is None:
        # create_shoot_db.create_full_database(input_dir, q_ids=True, subtrees_dir="")
        fn_diamond_db, q_hogs = create_profiles_database(input_dir, wd_list, nSpAll, selection="all", q_ids=True, subtrees_dir="")
    else:
        fn_diamond_db, q_hogs = create_profiles_database(input_dir, wd_list, nSpAll, selection="kmeans", q_ids=True, n_for_profile=xcelerate_config.n_for_profiles, subtrees_dir="")
    return fn_diamond_db, q_hogs


def RunSearch(options, speciessInfoObj, fn_diamond_db, prog_caller, q_one_query=False):
    name_to_print = "BLAST" if options.search_program == "blast" else options.search_program
    util.PrintUnderline("Running %s profiles search" % name_to_print)
    commands, results_files = GetOrderedSearchCommands(speciessInfoObj, fn_diamond_db, prog_caller, q_one_query=q_one_query, threads=options.nBlast)
    if q_one_query:
        return_code = parallel_task_manager.RunCommand(commands[0], qPrintOnError=True)
        if return_code != 0:
            print("ERROR: DIAMOND search failed, see messages above")
            util.Fail()
        util.PrintTime("Done profiles search\n")
        return results_files
    program_caller.RunParallelCommands(options.nBlast, commands, qListOfList=False, q_print_on_error=True)
    util.PrintTime("Done profiles search")
    return results_files


def GetOrderedSearchCommands(speciesInfoObj, diamond_db, prog_caller, q_one_query, threads=1, search_program="diamond",):
    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing
    the BLAST commands.
    """
    iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    if q_one_query:
        wd = files.FileHandler.GetSpeciesSeqsDir()[0]
        fn_single_fasta = "%sall_sequences.fa" % wd
        with open(fn_single_fasta, 'w') as outfile:
            for iFasta in iSpeciesNew:
                with open(files.FileHandler.GetSpeciesFastaFN(iFasta), 'r') as infile:
                    for line in infile:
                        outfile.write(line)
                    outfile.write("\n")
        results = wd + "Blast_all_sequences.txt"
        # commands = [prog_caller.GetSearchMethodCommand_Search(search_program, fn_single_fasta, diamond_db, results)]
        commands = ["diamond blastp --ignore-warnings -d %s -q %s -o %s --more-sensitive -p %d --quiet -e 0.001 --compress 1" % (diamond_db, fn_single_fasta, results, threads)]
        results_files = [results + ".gz"]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(
            search_program,
            files.FileHandler.GetSpeciesFastaFN(iFasta),
            diamond_db,
            files.FileHandler.GetBlastResultsFN(iFasta, -1, qForCreation=True))
            for iFasta in iSpeciesNew
        ]
        results_files = [files.FileHandler.GetBlastResultsFN(iFasta, -1, qForCreation=True) + ".gz" for iFasta in iSpeciesNew]
    return commands, results_files


def ogs_from_diamond_results(fn_og_results_out, q_ignore_sub=False):
    """
    Get the OG based on the DIAMOND results
    Args:
        fn_og_results_out - fn of compressed DIAMOND results
        q_ignore_sub - ignore any subtrees and just look at overall OGs
    Returns:
        iog        : List[str] "int.int" or "int" of ordered, ambiguous assignments
        species_closest_hits : Dict[str, Dict[str, int]] search_species -> Dict[hit_species, n_closest_hits]
    Info:
        Hits to other groups are returned if np.log10 difference is less than 10
        and -np.log10 score > difference.
    """
    ogs_sp_hits = defaultdict(list)
    scores_all_genes = defaultdict(list)
    with gzip.open(fn_og_results_out, 'rt') as infile:
        reader = csv.reader(infile, delimiter="\t")
        for line in reader:
            gene = line[0]
            og, sp, _ = line[1].split("_")
            ogs_sp_hits[gene].append((og.split(".", 1)[0] if q_ignore_sub else og, sp))
            scores_all_genes[gene].append(float(line[-2]))
    all_genes = list(ogs_sp_hits.keys())
    og_assignments = defaultdict()
    species_closest_hits = defaultdict(lambda: defaultdict(int))
    for gene in all_genes:
        ogs_sp = ogs_sp_hits[gene]
        scores = scores_all_genes[gene]
        sortedTuples = sorted(zip(scores, ogs_sp))
        scores = [i for i, j in sortedTuples]
        ogs = [j[0] for i, j in sortedTuples]
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
            scores_ml10ml10 = [s for s in scores_ml10 if s0-s<10]
            # ogs_all_genes[gene] = ogs[:len(scores_ml10)]
            # scores_all_genes[gene] = scores[:len(scores_ml10)]
            og_assignments[gene] = ogs[0]
            # print(sortedTuplesples)
            sp_hit = sortedTuples[0][1][1]
            # sys.exit()
            search_sp = gene.split("_")[0]
            species_closest_hits[search_sp][sp_hit] += 1
            # scores_all_genes[gene] = scores[0]
    return og_assignments, species_closest_hits


def get_original_orthogroups():
    wd_input_clusters = files.FileHandler.GetWorkingDirectory1_Read()[1]  # first is the current working directory, second is the most recent of the input directories
    fn_clusters = glob.glob(wd_input_clusters + "clusters*_id_pairs.txt")
    if len(fn_clusters) == 0:
        print("ERROR: Couldn't find previous orthogroups in %s" % wd_input_clusters)
        util.Fail()
    elif len(fn_clusters) > 1:
        print("WARNING: Found multiple orthogroup files: %s" % fn_clusters)
    fn_clusters = fn_clusters[0]
    ogs = mcl.GetPredictedOGs(fn_clusters)
    return ogs


def assign_genes(results_files):
    """
    Returns OGs with the new species added + species group: Dict[query_species, closest_species]
    """
    ogs = defaultdict(set)
    species_closest_hits_totals = defaultdict(lambda: defaultdict(int))
    for fn in results_files:
        ogs_all_genes, species_closest_hits = ogs_from_diamond_results(fn)
        for query_species, hits in species_closest_hits.items():
            for hit_species, count in hits.items():
                species_closest_hits_totals[query_species][hit_species] += count
        for gene, og in ogs_all_genes.items():
            ogs[int(og)].add(gene)
    species_group = dict()
    for query_species, hits in species_closest_hits_totals.items():
        species_group[query_species] = max(hits, key=hits.get)
        # print((query_species, max(hits, key=hits.get), hits))
    return ogs, species_group


# def write_all_orthogroups(ogs: List[Set[str]], ogs_new_species: Dict[int, Set[str]], ogs_clade_specific: List[List[Set[str]]]):
def write_all_orthogroups(ogs, ogs_new_species, ogs_clade_specific_lists):
    for iog, genes in ogs_new_species.items():
        ogs[iog].update(genes)
    for ogs_clade_specific in ogs_clade_specific_lists:
        for og in ogs_clade_specific:
            ogs.append(og)
    clustersFilename, clustersFilename_pairs = files.FileHandler.CreateUnusedClustersFN()
    mcl.write_updated_clusters_file(ogs, clustersFilename_pairs)
    return clustersFilename_pairs


def create_profiles_database(din, wd_list, nSpAll, selection="kmeans", n_for_profile=20, q_ids=True,
                             subtrees_dir="", q_hogs=False):
    """
    Create a fasta file with profile genes from each orthogroup
    Args:
        din - Input OrthoFinder results directory
        selection - "kmeans|random|all"  - method to use for sampling genes from orthogroups
        n_for_profile - The number of genes to use from each orthogroup, when available
        q_ids - Convert subtrees (with user gene accessions) back to IDs for profiles database
        q_hogs - The use of HOGs is currently not recommended due to the divergence of the resulting
                 orthogroups from an analysis with the original method
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
    fn_base = "profile_sequences.hogs" if q_hogs else "profile_sequences"
    if selection == "all":
        fn_fasta = wd + fn_base + ".%s.all.fa" % subtrees_label
    else:
        fn_fasta = wd + fn_base + ".%s.%d_%s.fa" % (subtrees_label, n_for_profile, selection)
    fn_diamond_db = fn_fasta + ".dmnd"
    if os.path.exists(fn_diamond_db):
        print("Profiles database already exists and will be reused: %s" % fn_diamond_db)
        return fn_diamond_db, q_hogs
    og_set = orthologues.OrthoGroupsSet(wd_list, list(range(nSpAll)), nSpAll, True)
    ids = og_set.Spec_SeqDict()
    ids_rev = {v: k for k, v in ids.items()}
    if q_hogs:
        try:
            ids_simple = og_set.SequenceDict()
            ids_simple_rev = {v: k for k, v in ids_simple.items()}
            ogs = read_hogs(din, "N0", ids_simple_rev)
        except RuntimeError:
            print("ERROR: Cannot read HOGs file, please report this error: https://github.com/davidemms/OrthoFinder/issues")
            q_hogs = False
            print("WARNING: Using MCL-based orthogroups as a fall-back")
            # util.Fail()
    if not q_hogs:
        clusters_filename = list(glob.glob(wd + "clusters_OrthoFinder*id_pairs.txt"))
        if len(clusters_filename) == 0:
            print("ERROR: Can't find %s" % wd + "clusters_OrthoFinder*id_pairs.txt")
        ogs = mcl.GetPredictedOGs(clusters_filename[0])
        fn_fasta = fn_fasta[:-7] + ".fa"
    fw = fasta_writer.FastaWriter(wd + "Species*fa", qGlob=True)
    seq_write = []
    seq_convert = dict()
    # print("WARNING: Check all gene names, can't start with '__'")
    # If there are subtrees then we need to convert their IDs in the profile file
    # back to internal IDs
    nToDo = len(ogs)
    for iog, og in enumerate(ogs):
        if iog >= 0 and divmod(iog, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
            util.PrintTime("Done %d of %d" % (iog, nToDo))
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
            n_for_profile = min(n_in_og, n_for_profile)
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
    print("")
    fw.WriteSeqsToFasta_withNewAccessions(seq_write, fn_fasta, seq_convert)
    parallel_task_manager.RunCommand(" ".join(["diamond", "makedb", "--in", fn_fasta, "-d", fn_diamond_db]), qPrintOnError=True, qPrintStderr=False)
    return fn_diamond_db, q_hogs


class DummyIDs:
    def __getitem__(self, arg):
        return arg


def read_hogs(din, hog_name, ids_rev=None):
    fn_ids = din + "WorkingDirectory/%s.ids.tsv" % hog_name
    if os.path.exists(fn_ids):
        fn = fn_ids
        ids_rev = DummyIDs()
    elif ids_rev is None:
        return []
    else:
        fn = din + "Phylogenetic_Hierarchical_Orthogroups/%s.tsv" % hog_name
    if not os.path.exists(fn):
        print("ERROR: %s does not exist" % fn)
        raise RuntimeError
    ogs = []
    with open(fn, util.csv_read_mode) as infile:
        reader = csv.reader(infile, delimiter="\t")
        next(reader)  # header
        for line in reader:
            ogs.append([])
            for species in line[3:]:
                genes = species.split(", ")
                genes = [ids_rev[g] for g in genes if g != '']
                ogs[-1].extend(genes)
            ogs[-1] = set(ogs[-1])
    return ogs


def write_unassigned_fasta(ogs_orig_list, ogs_new_genes, speciesInfoObj):
    if ogs_new_genes is not None:
        assigned_genes = set.union(*ogs_new_genes.values())
        assigned_genes.update(set.union(*[og for og in ogs_orig_list if len(og) > 1]))
    else:
        assigned_genes = set.union(*[og for og in ogs_orig_list if len(og) > 1])
    # Write out files for all unassigned genes
    # iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    # for iSp in iSpeciesNew:
    n_unassigned = []
    for iSp in range(speciesInfoObj.nSpAll):
        fw = fasta_writer.FastaWriter(files.FileHandler.GetSpeciesFastaFN(iSp))
        unassigned = set(fw.SeqLists.keys()).difference(assigned_genes)
        n_unassigned.append(len(unassigned))
        fw.WriteSeqsToFasta(unassigned, files.FileHandler.GetSpeciesUnassignedFastaFN(iSp, qForCreation=True))
    return n_unassigned


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


def get_new_species_clades(rooted_species_tree_fn, core_species_ids, n_core_species=2):
    """
    Args:
        rooted_species_tree_fn - ids format
        core_species_ids: Set[str]
        n_core_species - maximum number of core species in a 'new clade of species'. 1 is the minimum, but 2 makes more
                         allowances for gene loss in one core species & still recovering associated orthogroup
    Returns:
        a list of sorted lists of species IDs
    """
    core_species_ids = set(map(str, core_species_ids))
    t = tree.Tree(rooted_species_tree_fn, format=1)
    species_clades = []

    def is_new_clade(node):
        return len(core_species_ids.intersection(node.get_leaf_names())) <= n_core_species

    for n in t.get_leaves(is_leaf_fn=is_new_clade):
        species_group = list(sorted(map(int, n.get_leaf_names())))
        if set(species_group).difference(core_species_ids):
            species_clades.append(species_group)
    return species_clades
