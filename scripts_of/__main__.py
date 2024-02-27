#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2014 David Emms
#
# This program (OrthoFinder) is distributed under the terms of the GNU General Public License v3
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#  When publishing work that uses OrthoFinder please cite:
#      Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
#      improves orthogroup inference accuracy, Genome Biology 16:157
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.com
from __future__ import absolute_import

# first import parallel task manager to minimise RAM overhead for small processes
import multiprocessing as mp                    # optional  (problems on OpenBSD)
import platform                                 # Y
if __name__ == "__main__":
    if platform.system() == "Darwin":
        # https://github.com/davidemms/OrthoFinder/issues/570
        # https://github.com/davidemms/OrthoFinder/issues/663
        mp.set_start_method('fork')
    # else:
    #     # Should be more RAM efficient than fork and the time penalty
    #     # should be very small as we never try to create many processes
    #     mp.set_start_method('spawn')

from . import parallel_task_manager
ptm_initialised = parallel_task_manager.ParallelTaskManager_singleton()

import os                                       # Y
os.environ["OPENBLAS_NUM_THREADS"] = "1"    # fix issue with numpy/openblas. Will mean that single threaded options aren't automatically parallelised 

import sys                                      # Y
import copy                                     # Y
import subprocess                               # Y
import glob                                     # Y
import shutil                                   # Y
import time                                     # Y
import itertools                                # Y
import datetime                                 # Y
from collections import Counter                 # Y
from scipy.optimize import curve_fit            # install
import numpy as np                              # install
import csv                                      # Y
import scipy.sparse as sparse                   # install
import os.path                                  # Y
import numpy.core.numeric as numeric            # install
from collections import defaultdict             # Y
import xml.etree.ElementTree as ET              # Y
from xml.etree.ElementTree import SubElement    # Y
from xml.dom import minidom                     # Y
try: 
    import queue
except ImportError:
    import Queue as queue                       # Y
import warnings                                 # Y

from . import blast_file_processor, files, mcl, util, matrices, orthologues, program_caller, trees_msa, \
    split_ortholog_files, gathering, astral, trees2ologs_of, tree
from . import accelerate as acc

# Get directory containing script/bundle
if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
max_int = sys.maxsize
ok = False
while not ok:
    try:
        csv.field_size_limit(max_int)
        ok = True
    except OverflowError:
        max_int = int(max_int/10)
sys.setrecursionlimit(10**6)
    
fastaExtensions = {"fa", "faa", "fasta", "fas", "pep"}
# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f) 

my_env = os.environ.copy()
# use orthofinder supplied executables by preference
my_env['PATH'] = os.path.join(__location__, 'bin:') + my_env['PATH']
# Fix LD_LIBRARY_PATH when using pyinstaller 
if getattr(sys, 'frozen', False):
    if 'LD_LIBRARY_PATH_ORIG' in my_env:
        my_env['LD_LIBRARY_PATH'] = my_env['LD_LIBRARY_PATH_ORIG']  
    else:
        my_env['LD_LIBRARY_PATH'] = ''  
    if 'DYLD_LIBRARY_PATH_ORIG' in my_env:
        my_env['DYLD_LIBRARY_PATH'] = my_env['DYLD_LIBRARY_PATH_ORIG']  
    else:
        my_env['DYLD_LIBRARY_PATH'] = ''


def RunBlastDBCommand(command):
    capture = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env, shell=True)
    stdout, stderr = capture.communicate()
    try:
        stdout = stdout.decode()
        stderr = stderr.decode()
    except (UnicodeDecodeError, AttributeError):
        stdout = stdout.encode()
        stderr = stderr.encode()
    n_stdout_lines = stdout.count("\n")
    n_stderr_lines = stderr.count("\n")
    nLines_success= 12
    if n_stdout_lines > nLines_success or n_stderr_lines > 0 or capture.returncode != 0:
        print("\nWARNING: Likely problem with input FASTA files")
        if capture.returncode != 0:
            print("makeblastdb returned an error code: %d" % capture.returncode)
        else:
            print("makeblastdb produced unexpected output")
        print("Command: %s" % " ".join(command))
        print("stdout:\n-------")
        print(stdout)
        if len(stderr) > 0:
            print("stderr:\n-------")
            print(stderr)
            
def SpeciesNameDict(speciesIDsFN):
    # type: (str) -> Dict[int, str]
    speciesNamesDict = dict()
    with open(speciesIDsFN, 'r') as speciesNamesFile:
        for line in speciesNamesFile:
            if line.startswith("#"): continue
            line = line.rstrip()
            if not line: continue
            short, full = line.split(": ")
            speciesNamesDict[int(short)] = full.rsplit(".", 1)[0]
    return speciesNamesDict

"""
RunInfo
-------------------------------------------------------------------------------
"""
# Redundant?  
def GetNumberOfSequencesInFile(filename):
    count = 0
    with open(filename) as infile:
        for line in infile:
            if line.startswith(">"): count+=1
    return count


def GetOrderedSearchCommands(seqsInfo, speciesInfoObj, qDoubleBlast, search_program, prog_caller, n_genes_per_species=None, q_new_species_unassigned_genes=False):
    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands 
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing 
    the BLAST commands.
    n_genes_per_species: List[int] - number of genes per species
    """
    iSpeciesPrevious = list(range(speciesInfoObj.iFirstNewSpecies))
    iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    if q_new_species_unassigned_genes:
        if n_genes_per_species is not None:
            exclude = {isp for isp, n_genes in enumerate(n_genes_per_species) if n_genes == 0}
            iSpeciesNew = list(set(iSpeciesNew).difference(exclude))
        speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew)]
        taskSizes = [n_genes_per_species[i]*n_genes_per_species[j] for i,j in speciesPairs]
    else:
        speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew) if (qDoubleBlast or i <=j)] + \
                       [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesPrevious) if (qDoubleBlast or i <=j)] + \
                       [(i, j) for i, j in itertools.product(iSpeciesPrevious, iSpeciesNew) if (qDoubleBlast or i <=j)]
        taskSizes = [seqsInfo.nSeqsPerSpecies[i]*seqsInfo.nSeqsPerSpecies[j] for i,j in speciesPairs]
    taskSizes, speciesPairs = util.SortArrayPairByFirst(taskSizes, speciesPairs, True)
    if search_program == "blast":
        commands = [" ".join(["blastp", "-outfmt", "6", "-evalue", "0.001",
                              "-query", files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta) if q_new_species_unassigned_genes else files.FileHandler.GetSpeciesFastaFN(iFasta),
                              "-db", files.FileHandler.GetSpeciesDatabaseN(iDB),
                              "-out", files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)]) for iFasta, iDB in speciesPairs]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(search_program,
                        files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta) if q_new_species_unassigned_genes else files.FileHandler.GetSpeciesFastaFN(iFasta),
                        files.FileHandler.GetSpeciesDatabaseN(iDB, search_program),
                        files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)) for iFasta, iDB in speciesPairs]
    return commands


def GetOrderedSearchCommands_clades(seqsInfo, speciesInfoObj, qDoubleBlast, search_program, prog_caller, n_genes_per_species, species_clades):
    """
    Search all species
    """
    exclude = {isp for isp, n_genes in enumerate(n_genes_per_species) if n_genes == 0}
    speciesPairs = []
    for clade in species_clades:
        clade = list(set(clade).difference(exclude))
        speciesPairs.extend([(i, j) for i, j in itertools.product(clade, clade)])
    if search_program == "blast":
        commands = [" ".join(["blastp", "-outfmt", "6", "-evalue", "0.001",
                              "-query", files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta),
                              "-db", files.FileHandler.GetSpeciesDatabaseN(iDB),
                              "-out", files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)]) for iFasta, iDB in speciesPairs]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(search_program,
                        files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta),
                        files.FileHandler.GetSpeciesDatabaseN(iDB, search_program),
                        files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)) for iFasta, iDB in speciesPairs]
    return commands

"""
OrthoFinder
-------------------------------------------------------------------------------
"""   
g_mclInflation = 1.2

def CanRunBLAST():
    if parallel_task_manager.CanRunCommand("makeblastdb -help") and parallel_task_manager.CanRunCommand("blastp -help"):
        return True
    else:
        print("ERROR: Cannot run BLAST+")
        program_caller.ProgramCaller.PrintDependencyCheckFailure("makeblastdb -help\nblastp -help")
        print("Please check BLAST+ is installed and that the executables are in the system path\n")
        return False

def CanRunMCL():
    command = "mcl -h"
    if parallel_task_manager.CanRunCommand(command):
        return True
    else:
        print("ERROR: Cannot run MCL with the command \"%s\"" % command)
        program_caller.ProgramCaller.PrintDependencyCheckFailure(command)
        print("Please check MCL is installed and in the system path. See information above.\n")
        return False

def CanRunASTRAL():
    d_deps_test = files.FileHandler.GetDependenciesCheckDir()
    fn_test = d_deps_test + "astral.input.nwk"
    with open(fn_test, 'w') as outfile:
        for _ in range(20):
            outfile.write("(((A:1,B:1):2,(C:1,D:1):2,E:3):4);i\n")
    fn_output = d_deps_test + "astral.output.nwk"
    cmd = " ".join(["astral-pro", "-i", fn_test, "-o", fn_output, "-t", "2"])
    if parallel_task_manager.CanRunCommand(cmd, qAllowStderr=True, qRequireStdout=False, qCheckReturnCode=True):
        return True
    else:
        print("ERROR: Cannot run astral-pro")
        program_caller.ProgramCaller.PrintDependencyCheckFailure(cmd)
        print("Please check astral-pro is installed and that the executables are in the system path\n")
        return False
    
def GetProgramCaller():
    config_file = os.path.join(__location__, 'config.json') 
    pc = program_caller.ProgramCaller(config_file if os.path.exists(config_file) else None)
    config_file_user = os.path.expanduser("~/config_orthofinder_user.json")
    if os.path.exists(config_file_user):
        pc_user = program_caller.ProgramCaller(config_file_user)
        pc.Add(pc_user)
    return pc

def PrintHelp(prog_caller):  
    msa_ops = prog_caller.ListMSAMethods()
    tree_ops = prog_caller.ListTreeMethods()
    search_ops = prog_caller.ListSearchMethods()
    
    print("SIMPLE USAGE:") 
    print("Run full OrthoFinder analysis on FASTA format proteomes in <dir>")
    print("  orthofinder [options] -f <dir>")   
    # print("")
    # print("Add new species in <dir1> to previous run in <dir2> and run new analysis")
    # print("  orthofinder [options] -f <dir1> -b <dir2>")
    print("")
    print("To assign species from <dir1> to existing OrthoFinder orthogroups in <dir2>")
    print("  orthofinder [options] --assign <dir1> --core <dir2>")
    print("") 
      
    print("OPTIONS:")
    print(" -t <int>        Number of parallel sequence search threads [Default = %d]" % util.nThreadsDefault)
    print(" -a <int>        Number of parallel analysis threads")
    print(" -d              Input is DNA sequences")
    print(" -M <txt>        Method for gene tree inference. Options 'dendroblast' & 'msa'")
    print("                 [Default = msa]")
    print(" -S <txt>        Sequence search program [Default = diamond]")
    print("                 Options: " + ", ".join(['blast'] + search_ops))
    print(" -A <txt>        MSA program, requires '-M msa' [Default = mafft]")
    print("                 Options: " + ", ".join(msa_ops))
    print(" -T <txt>        Tree inference method, requires '-M msa' [Default = fasttree]")
    print("                 Options: " + ", ".join(tree_ops)) 
#    print(" -R <txt>        Tree reconciliation method [Default = of_recon]")
#    print("                 Options: of_recon, dlcpar, dlcpar_convergedsearch")
    print(" -s <file>       User-specified rooted species tree")
    print(" -I <int>        MCL inflation parameter [Default = %0.1f]" % g_mclInflation)
    print(" --save-space    Only create one compressed orthologs file per species")
    print(" -x <file>       Info for outputting results in OrthoXML format")
    print(" -p <dir>        Write the temporary pickle files to <dir>")
    print(" -1              Only perform one-way sequence search")
    print(" -X              Don't add species names to sequence IDs")
    print(" -y              Split paralogous clades below root of a HOG into separate HOGs")
    print(" -z              Don't trim MSAs (columns>=90% gap, min. alignment length 500)")
    print(" -n <txt>        Name to append to the results directory")  
    print(" -o <txt>        Non-default results directory")  
    print(" -h              Print this help text")

    print("")    
    print("WORKFLOW STOPPING OPTIONS:")   
    print(" -op             Stop after preparing input files for BLAST" )
    print(" -og             Stop after inferring orthogroups")
    print(" -os             Stop after writing sequence files for orthogroups")
    print("                 (requires '-M msa')")
    print(" -oa             Stop after inferring alignments for orthogroups")
    print("                 (requires '-M msa')")
    print(" -ot             Stop after inferring gene trees for orthogroups " )
   
    print("")   
    print("WORKFLOW RESTART COMMANDS:") 
    print(" -b  <dir>         Start OrthoFinder from pre-computed BLAST results in <dir>")   
    print(" -fg <dir>         Start OrthoFinder from pre-computed orthogroups in <dir>")
    print(" -ft <dir>         Start OrthoFinder from pre-computed gene trees in <dir>")
    
    print("")
    print("LICENSE:")
    print(" Distributed under the GNU General Public License (GPLv3). See License.md")
    util.PrintCitation() 
    
"""
Main
-------------------------------------------------------------------------------
"""   

def GetDirectoryArgument(arg, args):
    if len(args) == 0:
        print("Missing option for command line argument %s" % arg)
        util.Fail()
    directory = os.path.abspath(args.pop(0))
    if not os.path.isfile(directory) and directory[-1] != os.sep: 
        directory += os.sep
    if not os.path.exists(directory):
        print("Specified directory doesn't exist: %s" % directory)
        util.Fail()
    return directory

#def GetOrthogroupsDirectory(suppliedDir, options):
#    """
#    Possible directory structures
#    1. Default: 
#        FastaFiles/Results_<date>/                          <- Orthogroups spreadsheets                 
#        FastaFiles/Results_<date>/WorkingDirectory/         <- Sequence and BLAST files
#        FastaFiles/Results_<date>/Orthologues_<date>/       <- Orthologues
#        FastaFiles/Results_<date>/Orthologues_<date>/WorkingDirectory/, Trees/, Orthologues 
#    2. From BLAST: 
#        <MainDirectory>/                                    <- Orthogroups spreadsheets / Sequence and BLAST files
#        FastaFiles/Results_<date>/WorkingDirectory/
#        FastaFiles/Results_<date>/Orthologues_<date>/
#        FastaFiles/Results_<date>/Orthologues_<date>/WorkingDirectory/, Trees/, Orthologues 
#    """
   
# Control
class Options(object):#
    def __init__(self):
        self.nBlast = util.nThreadsDefault
        self.nProcessAlg = None
        self.qFastAdd = False  # Add species in near-linear time
        self.qStartFromBlast = False  # remove, just store BLAST to do
        self.qStartFromFasta = False  # local to argument checking
        self.qStartFromGroups = False
        self.qStartFromTrees = False
        self.qStopAfterPrepare = False
        self.qStopAfterGroups = False
        self.qStopAfterSeqs = False
        self.qStopAfterAlignments = False
        self.qStopAfterTrees = False
        self.qMSATrees = True  # Updated default
        self.qAddSpeciesToIDs = True
        self.qTrim = True
        self.gathering_version = (1, 0)    # < 3 is the original method
        self.search_program = "diamond"
        self.msa_program = "mafft"
        self.tree_program = "fasttree"
        self.recon_method = "of_recon"
        self.name = None   # name to identify this set of results
        self.qDoubleBlast = True
        self.qSplitParaClades = False
        self.qPhyldog = False
        self.speciesXMLInfoFN = None
        self.speciesTreeFN = None
        self.mclInflation = g_mclInflation
        self.dna = False
        self.fewer_open_files = True  # By default only open O(n) orthologs files at a time
        self.save_space = False  # On complete, have only one orthologs file per species
        self.v2_scores = False
        self.root_from_previous = False

    def what(self):
        for k, v in self.__dict__.items():
            if v == True:
                print(k)
                                 
def ProcessArgs(prog_caller, args):
    """ 
    Workflow
    | 1. Fasta Files | 2.  Prepare files    | 3.   Blast    | 4. Orthogroups    | 5.   Gene Trees     | 6.   Reconciliations/Orthologues   |

    Options
    Start from:
    -f: 1,2,..,6    (start from fasta files, --fasta)
    -b: 4,5,6       (start from blast results, --blast)
    -fg: 5,6         (start from orthogroups/do orthologue workflow, --from-groups)
    -ft: 6           (start from gene tree/do reconciliation, --from-trees)
    Stop at:
    -op: 2           (only prepare, --only-prepare)
    -og: 4           (orthogroups, --only-groups)
    """
    if len(args) == 0 or args[0] == "--help" or args[0] == "help" or args[0] == "-h":
        PrintHelp(prog_caller)
        util.Success() 

    options = Options()
    fastaDir = None
    continuationDir = None
    resultsDir_nonDefault = None
    pickleDir_nonDefault = None
    q_selected_msa_options = False
    q_selected_tree_options = False
    q_selected_search_option = False
    user_specified_M = False
    
    """
    -f: store fastaDir
    -b: store workingDir
    -fg: store orthologuesDir 
    -ft: store orthologuesDir 
    + xml: speciesXMLInfoFN
    """    
    
    while len(args) > 0:
        arg = args.pop(0)    
        if arg == "-f" or arg == "--fasta":
            if options.qStartFromFasta:
                print("Repeated argument: -f/--fasta\n")
                util.Fail()
            options.qStartFromFasta = True
            fastaDir = GetDirectoryArgument(arg, args)
        elif arg == "-b" or arg == "--blast":
            if options.qStartFromBlast:
                print("Repeated argument: -b/--blast\n")
                util.Fail()
            options.qStartFromBlast = True
            continuationDir = GetDirectoryArgument(arg, args)
        elif arg == "--assign":
            options.qFastAdd = True
            fastaDir = GetDirectoryArgument(arg, args)
        elif arg == "--core":
            options.qFastAdd = True
            continuationDir = GetDirectoryArgument(arg, args)
        elif arg == "-fg" or arg == "--from-groups":
            if options.qStartFromGroups:
                print("Repeated argument: -fg/--from-groups\n")
                util.Fail()
            options.qStartFromGroups = True
            continuationDir = GetDirectoryArgument(arg, args)
        elif arg == "-ft" or arg == "--from-trees":
            if options.qStartFromTrees:
                print("Repeated argument: -ft/--from-trees\n")
                util.Fail()
            options.qStartFromTrees = True
            continuationDir = GetDirectoryArgument(arg, args)
        elif arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nBlast = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s\n" % arg)
                util.Fail()    
        elif arg == "-a" or arg == "--algthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nProcessAlg = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s\n" % arg)
                util.Fail()   
        elif arg == "-1":
            options.qDoubleBlast = False
        elif arg == "-d" or arg == "--dna":
            options.dna = True
            if not q_selected_search_option:
                options.search_program = "blast_nucl"
        elif arg == "-X":
            options.qAddSpeciesToIDs = False
        elif arg == "-y":
            options.qSplitParaClades = True
        elif arg == "-z":
            options.qTrim = False
        elif arg == "-I" or arg == "--inflation":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.mclInflation = float(arg)
            except:
                print("Incorrect argument for MCL inflation parameter: %s\n" % arg)
                util.Fail()
        elif arg == "-c1":
            print("\nThe option 'c1' has been renamed '--c-homologs'")
            util.Fail()
        elif arg == "--c-homologs":
            options.gathering_version = (3, 2)
        elif arg == "--save-space":
            options.save_space = True
        elif arg == "-x" or arg == "--orthoxml":  
            if options.speciesXMLInfoFN:
                print("Repeated argument: -x/--orthoxml")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.speciesXMLInfoFN = args.pop(0)
        elif arg == "-n" or arg == "--name":  
            if options.name:
                print("Repeated argument: -n/--name")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.name = args.pop(0)
            while options.name.endswith("/"): options.name = options.name[:-1]
            if any([symbol in options.name for symbol in [" ", "/"]]): 
                print("Invalid symbol for command line argument %s\n" % arg)
                util.Fail()
        elif arg == "-o" or arg == "--output":  
            if resultsDir_nonDefault != None:
                print("Repeated argument: -o/--output")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            resultsDir_nonDefault = args.pop(0)
            while resultsDir_nonDefault.endswith("/"): resultsDir_nonDefault = resultsDir_nonDefault[:-1]
            resultsDir_nonDefault += "/"
            if os.path.exists(resultsDir_nonDefault):
                print("ERROR: non-default output directory already exists: %s\n" % resultsDir_nonDefault)
                util.Fail()
            if " " in resultsDir_nonDefault:
                print("ERROR: non-default output directory cannot include spaces: %s\n" % resultsDir_nonDefault)
                util.Fail()
            checkDirName = resultsDir_nonDefault
            while checkDirName.endswith("/"):
                checkDirName = checkDirName[:-1]
            path, newDir = os.path.split(checkDirName)
            if path != "" and not os.path.exists(path):
                print("ERROR: location '%s' for results directory '%s' does not exist.\n" % (path, newDir))
                util.Fail()
        elif arg == "-s" or arg == "--speciestree":  
            if options.speciesTreeFN:
                print("Repeated argument: -s/--speciestree")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.speciesTreeFN = args.pop(0)
        elif arg == "--scores-v2":
            options.v2_scores = True
        elif arg == "-S" or arg == "--search":
            choices = ['blast'] + prog_caller.ListSearchMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.search_program = arg
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()
        elif arg == "-M" or arg == "--method":
            arg_M_or_msa = arg
            user_specified_M = True
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg == "msa": 
                options.qMSATrees = True
            elif arg == "phyldog": 
                options.qPhyldog = True
                options.recon_method = "phyldog"
                options.qMSATrees = False
            elif arg == "dendroblast":
                options.qMSATrees = False
            else:
                print("Invalid argument for option %s: %s" % (arg_M_or_msa, arg))
                print("Valid options are 'dendroblast' and 'msa'\n")
                util.Fail()
        elif arg == "-A" or arg == "--msa_program":
            choices = ['mafft'] + prog_caller.ListMSAMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.msa_program = arg
                q_selected_msa_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()
        elif arg == "-T" or arg == "--tree_program":
            choices = ['fasttree'] + prog_caller.ListTreeMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.tree_program = arg
                q_selected_tree_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()
        elif arg == "-R" or arg == "--recon_method":
            choices = ['of_recon', 'dlcpar', 'dlcpar_convergedsearch', 'only_overlap']
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.recon_method = arg
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()
        elif arg == "-p":
            pickleDir_nonDefault = GetDirectoryArgument(arg, args)
        elif arg == "-op" or arg == "--only-prepare":
            options.qStopAfterPrepare = True
        elif arg == "-og" or arg == "--only-groups":
            options.qStopAfterGroups = True
        elif arg == "-os" or arg == "--only-seqs":
            options.qStopAfterSeqs = True
        elif arg == "-oa" or arg == "--only-alignments":
            options.qStopAfterAlignments = True
        elif arg == "-ot" or arg == "--only-trees":
            options.qStopAfterTrees = True
        elif arg == "-h" or arg == "--help":
            PrintHelp(prog_caller)
            util.Success()
        else:
            print("Unrecognised argument: %s\n" % arg)
            util.Fail()

    # set a default for number of algorithm threads
    if options.nProcessAlg is None:
        options.nProcessAlg = min(16, max(1, int(options.nBlast/8)))

    if options.nBlast < 1:
        print("ERROR: Number of '-t' threads cannot be fewer than 1, got %d" % options.nBlast)
        util.Fail()    

    if options.nProcessAlg < 1:
        print("ERROR: Number of '-a' threads cannot be fewer than 1, got %d" % options.nProcessAlg)
        util.Fail()  

    # check argument combinations       
    if not (options.qStartFromFasta or options.qStartFromBlast or options.qStartFromGroups or options.qStartFromTrees or options.qFastAdd):
        print("ERROR: Please specify the input directory for OrthoFinder using one of the options: '-f', '-b', '-fg' or '-ft', '--assign'.")
        util.Fail()

    if options.qFastAdd:
        if (options.qStartFromFasta or options.qStartFromBlast or options.qStartFromGroups or options.qStartFromTrees):
            print("ERROR: Incompatible options used with --assign, cannot accept: '-f', '-b', '-fg' or '-ft'")
            util.Fail()
        if fastaDir is None:
            print("ERROR: '--assign' option also requires '--core' directory to be specified")
            util.Fail()
        if continuationDir is None:
            print("ERROR: '--core' option also requires '--assign' directory to be specified")
            util.Fail()
        if not options.qMSATrees:
            print("ERROR: --assign requires MSA trees, option '-M dendroblast' is invalid")
            util.Fail()

    if options.qStartFromFasta and (options.qStartFromTrees or options.qStartFromGroups):
        print("ERROR: Incompatible arguments, -f (start from fasta files) and" + (" -fg (start from orthogroups)" if options.qStartFromGroups else " -ft (start from trees)"))
        util.Fail()
        
    if options.qStartFromBlast and (options.qStartFromTrees or options.qStartFromGroups):
        print("ERROR: Incompatible arguments, -b (start from pre-calcualted BLAST results) and" + (" -fg (start from orthogroups)" if options.qStartFromGroups else " -ft (start from trees)"))
        util.Fail()      

    if options.qStartFromTrees and options.qStartFromGroups:
        print("ERROR: Incompatible arguments, -fg (start from orthogroups) and -ft (start from trees)")
        util.Fail()    

    if options.qStopAfterSeqs and (not options.qMSATrees):
        print("ERROR: Argument '-os' (stop after sequences) also requires option '-M msa'")
        util.Fail()   

    if options.qStopAfterAlignments and (not options.qMSATrees):
        print("ERROR: Argument '-oa' (stop after alignments) also requires option '-M msa'")
        util.Fail()     

    if (q_selected_msa_options or q_selected_tree_options) and (not options.qMSATrees and not options.qPhyldog):
        print("ERROR: Argument '-A' or '-T' (multiple sequence alignment/tree inference program) also requires option '-M msa'")
        util.Fail()       
        
    if options.qPhyldog and (not options.speciesTreeFN):
        print("ERROR: Phyldog currently needs a species tree to be provided")
        util.Fail()          

    if resultsDir_nonDefault != None and ((not options.qStartFromFasta) or options.qStartFromBlast):
        print("ERROR: Incompatible arguments, -o (non-default output directory) can only be used with a new OrthoFinder run using option '-f'")
        util.Fail()       
        
    if options.search_program not in (prog_caller.ListSearchMethods() + ['blast']):
        print("ERROR: Search program (%s) not configured in config.json file" % options.search_program)
        util.Fail()

    util.PrintTime("Starting OrthoFinder %s" % util.version)    
    print("%d thread(s) for highly parallel tasks (BLAST searches etc.)" % options.nBlast)
    print("%d thread(s) for OrthoFinder algorithm\n" % options.nProcessAlg)

    if options.qFastAdd and not q_selected_msa_options:
        print("INFO: For --assign defaulting to 'mafft --memsave' to reduce RAM usage\n")
        options.msa_program = "mafft_memsave"

    if options.qFastAdd and not q_selected_tree_options:
        print("INFO: For --assign defaulting to 'FastTree -fastest' to reduce RAM usage\n")
        options.tree_program = "fasttree_fastest"

    return options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault, user_specified_M

def GetXMLSpeciesInfo(seqsInfoObj, options):
    # speciesInfo:  name, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
    util.PrintUnderline("Reading species information file")
    # do this now so that we can alert user to any errors prior to running the algorithm
    speciesXML = [[] for i_ in seqsInfoObj.speciesToUse]
    speciesNamesDict = SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
    speciesRevDict = {v:k for k,v in speciesNamesDict.items()}
    userFastaFilenames = [os.path.split(speciesNamesDict[i])[1] for i in seqsInfoObj.speciesToUse]
    with open(options.speciesXMLInfoFN, 'r') as speciesInfoFile:
        reader = csv.reader(speciesInfoFile, delimiter = "\t")
        for iLine, line in enumerate(reader):
            if len(line) != 5:
                # allow for an extra empty line at the end
                if len(line) == 0 and iLine == len(userFastaFilenames):
                    continue
                print("ERROR")
                print("Species information file %s line %d is incorrectly formatted." % (options.speciesXMLInfoFN, iLine + 1))
                print("File should be contain one line per species")
                print("Each line should contain 5 tab-delimited fields:")
                print("  fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseFastaFilename")
                print("See README file for more information.")
                util.Fail() 
            fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile = line
            try:
                iSpecies = speciesRevDict[os.path.splitext(fastaFilename)[0]]
            except KeyError:
                print("Skipping %s from line %d as it is not being used in this analysis" % (fastaFilename, iLine+1))
                continue
            speciesXML[seqsInfoObj.speciesToUse.index(iSpecies)] = line   
    # check information has been provided for all species
    speciesMissing = False        
    for iPos, iSpecies in enumerate(seqsInfoObj.speciesToUse):
        if speciesXML[iPos] == []:
            if not speciesMissing:
                print("ERROR")
                print("Species information file %s does not contain information for all species." % options.speciesXMLInfoFN)
                print("Information is missing for:") 
                speciesMissing = True
            print(speciesNamesDict[iSpecies])
    if speciesMissing:
        util.Fail()
    return speciesXML

def IDsFileOK(filename):
    """
    It is best to detect any issues with input files at start, perform all required checks here
    """
    with open(filename, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if len(line) == 0: continue
            tokens = line.split(": ", 1)
            if len(tokens) !=2 or len(tokens[1]) == 0:
                return False, line
    return True, None

def CheckDependencies(options, user_specified_m, prog_caller, dirForTempFiles):
    util.PrintUnderline("Checking required programs are installed")
    if not user_specified_m:
        print("Running with the recommended MSA tree inference by default. To revert to legacy method use '-M dendroblast'.\n")
    if options.qStartFromFasta or options.qFastAdd:
        if options.search_program == "blast":
            if not CanRunBLAST(): util.Fail()
        else:
            d_deps_check = files.FileHandler.GetDependenciesCheckDir()
            success, stdout, stderr, cmd = prog_caller.TestSearchMethod(options.search_program, d_deps_check)
            if not success:
                print("\nERROR: Cannot run %s" % options.search_program)
                prog_caller.PrintDependencyCheckFailure(cmd)
                print("Please check %s is installed and that the executables are in the system path\n" % options.search_program)
                util.Fail()
    if (options.qStartFromFasta or options.qStartFromBlast) and not CanRunMCL():
        util.Fail()
    if not (options.qStopAfterPrepare or options.qStopAfterSeqs or options.qStopAfterGroups or options.qStartFromTrees):
        if not orthologues.CanRunOrthologueDependencies(dirForTempFiles, 
                                                            options.qMSATrees, 
                                                            options.qPhyldog, 
                                                            options.qStopAfterTrees, 
                                                            options.msa_program, 
                                                            options.tree_program, 
                                                            options.recon_method,
                                                            prog_caller, 
                                                            options.qStopAfterAlignments):
            print("Dependencies have been met for inference of orthogroups but not for the subsequent orthologue inference.")
            print("Either install the required dependencies or use the option '-og' to stop the analysis after the inference of orthogroups.\n")
            util.Fail()
    if options.qFastAdd:
        if not CanRunASTRAL():
            util.Fail()


# 0
def ProcessPreviousFiles(workingDir_list, qDoubleBlast, check_blast=True):
    """Checks for:
    workingDir should be the WorkingDirectory containing Blast*.txt files
    
    SpeciesIDs.txt
    Species*.fa
    Blast*txt
    SequenceIDs.txt
    
    Checks which species should be included
    
    """
    # check BLAST results directory exists
    if not os.path.exists(workingDir_list[0]):
        err_text = "ERROR: Previous/Pre-calculated BLAST results directory does not exist: %s\n" % workingDir_list[0]
        files.FileHandler.LogFailAndExit(err_text)
        
    speciesInfo = util.SpeciesInfo()
    if not os.path.exists(files.FileHandler.GetSpeciesIDsFN()):
        err_text = "ERROR: %s file must be provided if using previously calculated BLAST results" % files.FileHandler.GetSpeciesIDsFN()
        files.FileHandler.LogFailAndExit(err_text)
    file_ok, err_line = IDsFileOK(files.FileHandler.GetSpeciesIDsFN())
    if not file_ok: 
        files.FileHandler.LogFailAndExit("ERROR: %s file contains a blank accession. Line:\n %s" % (files.FileHandler.GetSpeciesIDsFN(), err_line))
    speciesInfo.speciesToUse, speciesInfo.nSpAll, speciesToUse_names = util.GetSpeciesToUse(files.FileHandler.GetSpeciesIDsFN())
 
    # check fasta files are present 
    previousFastaFiles = files.FileHandler.GetSortedSpeciesFastaFiles()
    if len(previousFastaFiles) == 0:
        err_text = "ERROR: No processed fasta files in the supplied previous working directories:\n" + "\n".join(workingDir_list) + "\n"
        files.FileHandler.LogFailAndExit(err_text)
    tokens = previousFastaFiles[-1][:-3].split("Species")
    lastFastaNumberString = tokens[-1]
    iLastFasta = 0
    nFasta = len(previousFastaFiles)
    try:
        iLastFasta = int(lastFastaNumberString)
    except:
        files.FileHandler.LogFailAndExit("ERROR: Filenames for processed fasta files are incorrect: %s\n" % previousFastaFiles[-1])
    if nFasta != iLastFasta + 1:
        files.FileHandler.LogFailAndExit("ERROR: Not all expected fasta files are present. Index of last fasta file is %s but found %d fasta files.\n" % (lastFastaNumberString, len(previousFastaFiles)))
    
    # check BLAST files
    if check_blast:
        blast_fns_triangular = [files.FileHandler.GetBlastResultsFN(iSpecies, jSpecies, raise_exception=False) for iSpecies in speciesInfo.speciesToUse for jSpecies in speciesInfo.speciesToUse if jSpecies >= iSpecies]
        have_triangular = [(os.path.exists(fn) or os.path.exists(fn + ".gz")) for fn in blast_fns_triangular]
        for qHave, fn in zip(have_triangular, blast_fns_triangular):
            if not qHave: print("Required BLAST results file not found: %s" % fn)

        if qDoubleBlast:
            blast_fns_remainder = [files.FileHandler.GetBlastResultsFN(iSpecies, jSpecies, raise_exception=False) for iSpecies in speciesInfo.speciesToUse for jSpecies in speciesInfo.speciesToUse if jSpecies < iSpecies]
            have_remainder = [(os.path.exists(fn) or os.path.exists(fn + ".gz")) for fn in blast_fns_remainder]
            if not (all(have_triangular) and all(have_remainder)):
                for qHave, fn in zip(have_remainder, blast_fns_remainder):
                    if not qHave: print("Required BLAST results file not found: %s" % fn)
                if not all(have_triangular):
                    files.FileHandler.LogFailAndExit()
                else:
                    # would be able to do it using just one-way blast
                    files.FileHandler.LogFailAndExit("ERROR: Required BLAST results files are present for using the one-way sequence search option (default) but not the double BLAST search ('-d' option)")
        else:
            if not all(have_triangular):
                files.FileHandler.LogFailAndExit()
                            
    # check SequenceIDs.txt and SpeciesIDs.txt files are present
    if not os.path.exists(files.FileHandler.GetSequenceIDsFN()):
        files.FileHandler.LogFailAndExit("ERROR: %s file must be provided if using previous calculated BLAST results" % files.FileHandler.GetSequenceIDsFN())
    
    file_ok, err_line = IDsFileOK(files.FileHandler.GetSequenceIDsFN())
    if not file_ok: 
        files.FileHandler.LogFailAndExit("ERROR: %s file contains a blank accession. Line:\n %s" % (files.FileHandler.GetSequenceIDsFN(), err_line))
    return speciesInfo, speciesToUse_names

# 6
def CreateSearchDatabases(speciesInfoObj, options, prog_caller, q_unassigned_genes=False):
    iSpeciesToDo = range(max(speciesInfoObj.speciesToUse) + 1)
    for iSp in iSpeciesToDo:
        fn_fasta = files.FileHandler.GetSpeciesUnassignedFastaFN(iSp) if q_unassigned_genes else files.FileHandler.GetSpeciesFastaFN(iSp)
        if options.search_program == "blast":
            command = " ".join(["makeblastdb", "-dbtype", "prot", "-in", fn_fasta, "-out", files.FileHandler.GetSpeciesDatabaseN(iSp)])
            util.PrintTime("Creating Blast database %d of %d" % (iSp + 1, len(iSpeciesToDo)))
            RunBlastDBCommand(command) 
        else:
            command = prog_caller.GetSearchMethodCommand_DB(options.search_program, fn_fasta, files.FileHandler.GetSpeciesDatabaseN(iSp, options.search_program))
            util.PrintTime("Creating %s database %d of %d" % (options.search_program, iSp + 1, len(iSpeciesToDo)))
            ret_code = parallel_task_manager.RunCommand(command, qPrintOnError=True, qPrintStderr=False)
            if ret_code != 0:
                files.FileHandler.LogFailAndExit("ERROR: diamond makedb failed")

# 7
def RunSearch(options, speciessInfoObj, seqsInfo, prog_caller, q_new_species_unassigned_genes=False, n_genes_per_species=None, species_clades=None):
    """
    n_genes_per_species: List[int] - optional, for use with unassigned genes. If a species has zero unassigned genes, don't search with/agaisnt it.
    """
    name_to_print = "BLAST" if options.search_program == "blast" else options.search_program
    if options.qStopAfterPrepare:
        util.PrintUnderline("%s commands that must be run" % name_to_print)
    elif species_clades is not None:
        util.PrintUnderline("Running %s searches for new non-core species clades" % name_to_print)
    else:
        util.PrintUnderline("Running %s all-versus-all" % name_to_print)
    if species_clades is None:
        commands = GetOrderedSearchCommands(seqsInfo, speciessInfoObj, options.qDoubleBlast, options.search_program, prog_caller, n_genes_per_species, q_new_species_unassigned_genes=q_new_species_unassigned_genes)
    else:
        commands = GetOrderedSearchCommands_clades(seqsInfo, speciessInfoObj, options.qDoubleBlast, options.search_program, prog_caller, n_genes_per_species, species_clades)
    if options.qStopAfterPrepare:
        for command in commands:
            print(command)
        util.Success()
    print("Using %d thread(s)" % options.nBlast)
    util.PrintTime("This may take some time....")
    program_caller.RunParallelCommands(options.nBlast, commands, qListOfList=False,
                                       q_print_on_error=True, q_always_print_stderr=False)

    # remove BLAST databases
    util.PrintTime("Done all-versus-all sequence search")
    if options.search_program == "blast":
        for f in glob.glob(files.FileHandler.GetWorkingDirectory1_Read()[0] + "BlastDBSpecies*"):
            os.remove(f)
    if options.search_program == "mmseqs":
        for i in range(speciessInfoObj.nSpAll):
            for j in range(speciessInfoObj.nSpAll):
                tmp_dir = "/tmp/tmpBlast%d_%d.txt" % (i,j)
                if os.path.exists(tmp_dir):
                    try:
                        shutil.rmtree(tmp_dir)
                    except OSError:
                        time.sleep(1)
                        shutil.rmtree(tmp_dir, True)  # shutil / NFS bug - ignore errors, it's less crucial that the files are deleted

# 9
def GetOrthologues(speciesInfoObj, options, prog_caller, i_og_restart=0):
    util.PrintUnderline("Analysing Orthogroups", True)
    orthologues.OrthologuesWorkflow(speciesInfoObj.speciesToUse,
                                    speciesInfoObj.nSpAll,
                                    prog_caller,
                                    options.msa_program,
                                    options.tree_program,
                                    options.recon_method,
                                    options.nBlast,
                                    options.nProcessAlg,
                                    options.qDoubleBlast,
                                    options.qAddSpeciesToIDs,
                                    options.qTrim,
                                    options.fewer_open_files,
                                    options.speciesTreeFN,
                                    options.qStopAfterSeqs,
                                    options.qStopAfterAlignments,
                                    options.qStopAfterTrees,
                                    options.qMSATrees,
                                    options.qPhyldog,
                                    options.name,
                                    options.qSplitParaClades,
                                    options.save_space,
                                    root_from_previous = False,
                                    i_og_restart=i_og_restart)
    util.PrintTime("Done orthologues")


def BetweenCoreOrthogroupsWorkflow(continuationDir, speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files, q_hogs):
    """
    Infer clade-specific orthogroups for the new species clades
    n_unassigned: List[int] - number of unassigned genes per species
    """
    # Get current orthogroups - original orthogroups plus genes assigned to them
    if q_hogs:
        ogs = acc.read_hogs(continuationDir, "N0")
    else:
        ogs = acc.get_original_orthogroups()
    i_og_restart = 0
    ogs_new_species, _ = acc.assign_genes(results_files)
    clustersFilename_pairs = acc.write_all_orthogroups(ogs, ogs_new_species, [])  # this updates ogs

    ogSet = orthologues.OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse,
                                       speciesInfoObj.nSpAll, options.qAddSpeciesToIDs, idExtractor=util.FirstWordExtractor)

    if options.qStopAfterGroups and options.speciesTreeFN is None:
        # Can't infer clade-specific groups
        print("\nSpecies tree required for clade-speicfic orthogroups - skipping")
        return clustersFilename_pairs, i_og_restart

    n_unassigned = acc.write_unassigned_fasta(ogs, None, speciesInfoObj)

    # Get/Infer species tree
    if options.speciesTreeFN is None:
        # Infer gene trees
        # We write orthogroup & stats results files in the following code, which we should avoid & only do once all OGs are done.
        gathering.post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo, speciesNamesDict, options,
                                              speciesXML=None, q_incremental=True)

        orthologues.InferGeneAndSpeciesTrees(ogSet,
                                              prog_caller, options.msa_program, options.tree_program,
                                              options.nBlast, options.nProcessAlg, options.qDoubleBlast,
                                              options.qAddSpeciesToIDs,
                                              options.qTrim, userSpeciesTree=None, qStopAfterSeqs=False,
                                              qStopAfterAlign=False, qMSA=options.qMSATrees,
                                              qPhyldog=False, results_name=options.name,
                                              root_from_previous=True)

        # Infer species tree
        astral_fn = files.FileHandler.GetAstralFilename()
        astral.create_input_file(files.FileHandler.GetOGsTreeDir(), astral_fn)
        species_tree_unrooted_fn = files.FileHandler.GetSpeciesTreeUnrootedFN()
        parallel_task_manager.RunCommand(astral.get_astral_command(astral_fn, species_tree_unrooted_fn, options.nBlast))

        # Root it
        core_rooted_species_tree = tree.Tree(files.FileHandler.GetCoreSpeciesTreeIDsRootedFN(), format=1)
        species_to_speices_map = lambda x: x
        rooted_species_tree_ids, qHaveSupport = trees2ologs_of.CheckAndRootTree(species_tree_unrooted_fn, core_rooted_species_tree, species_to_speices_map)
        rooted_species_tree_fn = files.FileHandler.GetSpeciesTreeIDsRootedFN()
        rooted_species_tree_ids.write(outfile=rooted_species_tree_fn)

        spTreeUnrootedFN = files.FileHandler.GetSpeciesTreeResultsFN(None, True)
        util.RenameTreeTaxa(rooted_species_tree_ids, spTreeUnrootedFN, ogSet.SpeciesDict(), qSupport=qHaveSupport, qFixNegatives=True)

        labeled_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
        util.RenameTreeTaxa(rooted_species_tree_ids, labeled_tree_fn, ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N')
        i_og_restart = len(ogs)  # Need to process the clade-specific orthogroups only
    else:
        util.PrintUnderline("Using user-supplied species tree")
        spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
        orthologues.ConvertUserSpeciesTree(options.speciesTreeFN, ogSet.SpeciesDict(), spTreeFN_ids)
        rooted_species_tree_fn = spTreeFN_ids

    # Identify clades for clade-specific orthogroup inference
    iSpeciesCore = set(speciesInfoObj.get_original_species())
    species_clades = acc.get_new_species_clades(rooted_species_tree_fn, iSpeciesCore)
    util.PrintUnderline("Identifying clade-specific orthogroups for the following clades:", qHeavy=True)
    species_dict = ogSet.SpeciesDict()
    for i, clade in enumerate(species_clades):
        print(str(i) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
    print("")

    # Clade-specific orthogroup inference
    CreateSearchDatabases(speciesInfoObj, options, prog_caller, q_unassigned_genes=True)
    # provide list of clades, only run these searches (and only if the fasta files are non-empty)
    RunSearch(options, speciesInfoObj, seqsInfo, prog_caller, n_genes_per_species=n_unassigned, species_clades=species_clades)
    # process the results files - only if they are present and non-empty
    options.v2_scores = True
    n_clades = len(species_clades)
    clustersFilename_pairs_unassigned_all = []
    for i_clade, clade in enumerate(species_clades):
        util.PrintUnderline("OrthoFinder clutering on new species clade %d of %d" % (i_clade+1, n_clades))
        print(str(i) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
        speciesInfo_clade = copy.deepcopy(speciesInfoObj)
        speciesInfo_clade.speciesToUse = clade
        seqsInfo_clade = util.SeqsInfoRecompute(seqsInfo, clade)
        clustersFilename_pairs_unassigned = gathering.DoOrthogroups(options, speciesInfo_clade, seqsInfo_clade,
                                                         speciesNamesDict, speciesXML=None, i_unassigned=i_clade)
        clustersFilename_pairs_unassigned_all.append(clustersFilename_pairs_unassigned)
    ogs_clade_specific_list = [mcl.GetPredictedOGs(filename) for filename in clustersFilename_pairs_unassigned_all]

    # OGs have had assigned genes added to them already
    clustersFilename_pairs = acc.write_all_orthogroups(ogs, {}, ogs_clade_specific_list)
    return clustersFilename_pairs, i_og_restart


def clade_specific_orthogroups_v1(speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files):
    ogs = acc.get_original_orthogroups()
    i_og_restart = len(ogs)
    ogs_new_species, _ = acc.assign_genes(results_files)
    n_unassigned = acc.write_unassigned_fasta(ogs, ogs_new_species, speciesInfoObj)
    CreateSearchDatabases(speciesInfoObj, options, prog_caller, q_unassigned_genes=True)
    RunSearch(options, speciesInfoObj, seqsInfo, prog_caller, n_genes_per_species=n_unassigned, q_new_species_unassigned_genes=True)

    # Update info objects for clustering of only the new species
    speciesInfoObj_for_unassigned = copy.deepcopy(speciesInfoObj)
    speciesInfoObj_for_unassigned.speciesToUse = [iSp for iSp in speciesInfoObj.speciesToUse if iSp >= speciesInfoObj.iFirstNewSpecies]

    nSpecies_unassigned = len(speciesInfoObj_for_unassigned.speciesToUse)
    nSeqs_in_new_species = sum([seqsInfo.nSeqsPerSpecies[iSp] for iSp in speciesInfoObj_for_unassigned.speciesToUse])
    n_offset = seqsInfo.seqStartingIndices[-nSpecies_unassigned]
    starting_indicies = [n-n_offset for n in seqsInfo.seqStartingIndices[-nSpecies_unassigned:]]
    seqsInfo_for_unassigned = util.SequencesInfo(
        nSeqs = nSeqs_in_new_species,
        nSpecies = nSpecies_unassigned,
        speciesToUse = speciesInfoObj_for_unassigned.speciesToUse,
        seqStartingIndices = starting_indicies,
        nSeqsPerSpecies = seqsInfo.nSeqsPerSpecies
    )
    options.v2_scores = True
    clustersFilename_pairs_unassigned = gathering.DoOrthogroups(options, speciesInfoObj_for_unassigned, seqsInfo_for_unassigned,
                                                     speciesNamesDict, speciesXML=None, i_unassigned=0)
    ogs_clade_specific = mcl.GetPredictedOGs(clustersFilename_pairs_unassigned)
    clustersFilename_pairs = acc.write_all_orthogroups(ogs, ogs_new_species, [ogs_clade_specific])
    return clustersFilename_pairs, i_og_restart

def GetOrthologues_FromTrees(options):
    orthologues.OrthologuesFromTrees(options.recon_method, options.nBlast, options.nProcessAlg, options.speciesTreeFN,
                                     options.qAddSpeciesToIDs, options.qSplitParaClades, options.fewer_open_files)
 
def ProcessesNewFasta(fastaDir, q_dna, speciesInfoObj_prev = None, speciesToUse_prev_names=[]):
    """
    Process fasta files and return a Directory object with all paths completed.
    """
    # Check files present
    qOk = True
    if not os.path.exists(fastaDir):
        print("\nDirectory does not exist: %s" % fastaDir)
        util.Fail()
    files_in_directory = sorted([f for f in os.listdir(fastaDir) if os.path.isfile(os.path.join(fastaDir,f))])
    originalFastaFilenames = []
    excludedFiles = []
    for f in files_in_directory:
        if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in fastaExtensions and not f.startswith("._"):
            originalFastaFilenames.append(f)
        else:
            excludedFiles.append(f)
    if len(excludedFiles) != 0:
        print("\nWARNING: Files have been ignored as they don't appear to be FASTA files:")
        for f in excludedFiles:
            print(f)
        print("OrthoFinder expects FASTA files to have one of the following extensions: %s" % (", ".join(fastaExtensions)))
    speciesToUse_prev_names = set(speciesToUse_prev_names)
    if len(originalFastaFilenames) + len(speciesToUse_prev_names) < 2:
        print("ERROR: At least two species are required")
        util.Fail()
    if any([fn in speciesToUse_prev_names for fn in originalFastaFilenames]):
        print("ERROR: Attempted to add a second copy of a previously included species:")
        for fn in originalFastaFilenames:
            if fn in speciesToUse_prev_names: print(fn)
        print("")
        util.Fail()
    if len(originalFastaFilenames) == 0:
        print("\nNo fasta files found in supplied directory: %s" % fastaDir)
        util.Fail()
    if speciesInfoObj_prev == None:
        # Then this is a new, clean analysis 
        speciesInfoObj = util.SpeciesInfo()
    else:
        speciesInfoObj = speciesInfoObj_prev
    iSeq = 0
    iSpecies = 0
    # If it's a previous analysis:
    if len(speciesToUse_prev_names) != 0:
        with open(files.FileHandler.GetSpeciesIDsFN(), 'r') as infile:
            for line in infile: pass
        if line.startswith("#"): line = line[1:]
        iSpecies = int(line.split(":")[0]) + 1
    speciesInfoObj.iFirstNewSpecies = iSpecies
    newSpeciesIDs = []
    with open(files.FileHandler.GetSequenceIDsFN(), 'a') as idsFile, open(files.FileHandler.GetSpeciesIDsFN(), 'a') as speciesFile:
        for fastaFilename in originalFastaFilenames:
            newSpeciesIDs.append(iSpecies)
            outputFasta = open(files.FileHandler.GetSpeciesFastaFN(iSpecies, qForCreation=True), 'w')
            fastaFilename = fastaFilename.rstrip()
            speciesFile.write("%d: %s\n" % (iSpecies, fastaFilename))
            baseFilename, extension = os.path.splitext(fastaFilename)
            mLinesToCheck = 100
            qHasAA = False
            with open(fastaDir + os.sep + fastaFilename, 'r') as fastaFile:
                for iLine, line in enumerate(fastaFile):
                    if line.isspace(): continue
                    if len(line) > 0 and line[0] == ">":
                        newID = "%d_%d" % (iSpecies, iSeq)
                        acc = line[1:].rstrip()
                        if len(acc) == 0:
                            print("ERROR: %s contains a blank accession line on line %d" % (fastaDir + os.sep + fastaFilename, iLine+1))
                            util.Fail()
                        idsFile.write("%s: %s\n" % (newID, acc))
                        outputFasta.write(">%s\n" % newID)    
                        iSeq += 1
                    else:
                        line = line.upper()    # allow lowercase letters in sequences
                        if not qHasAA and (iLine < mLinesToCheck):
#                            qHasAA = qHasAA or any([c in line for c in ['D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y']])
                            qHasAA = qHasAA or any([c in line for c in ['E','F','I','L','P','Q']]) # AAs minus nucleotide ambiguity codes
                        outputFasta.write(line)
                outputFasta.write("\n")
            if (not qHasAA) and (not q_dna):
                qOk = False
                print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences. Use '-d' option" % fastaFilename)
            iSpecies += 1
            iSeq = 0
            outputFasta.close()
        if not qOk:
            util.Fail()
    if len(originalFastaFilenames) > 0: outputFasta.close()
    speciesInfoObj.speciesToUse = speciesInfoObj.speciesToUse + newSpeciesIDs
    speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1      # will be one of the new species
    return speciesInfoObj

def DeleteDirectoryTree(d):
    if os.path.exists(d): 
        try:
            shutil.rmtree(d)
        except OSError:
            time.sleep(1)
            shutil.rmtree(d, True)   

def CheckOptions(options, speciesToUse):
    """Check any optional arguments are valid once we know what species are in the analysis
    - user supplied species tree
    """
    if options.speciesTreeFN:
        expSpecies = list(SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN()).values())
        orthologues.CheckUserSpeciesTree(options.speciesTreeFN, expSpecies)

    # check can open enough files
    n_extra = 50
    q_do_orthologs = not any((options.qStopAfterPrepare, options.qStopAfterGroups, options.qStopAfterSeqs, options.qStopAfterAlignments, options.qStopAfterTrees))
    if q_do_orthologs:
        n_sp = len(speciesToUse)
        wd = files.FileHandler.GetWorkingDirectory_Write()
        wd_files_test = wd + "Files_test/"
        fh = []
        try:
            if not os.path.exists(wd_files_test):
                os.mkdir(wd_files_test)
            for i_sp in range(n_sp):
                di = wd_files_test + "Sp%d/" % i_sp
                if not os.path.exists(di):
                    os.mkdir(di)
                for j_sp in range(1):  # We only create a linear number of ortholog files now
                    fnij = di + "Sp%d.txt" % j_sp
                    fh.append(open(fnij, 'w'))
            # create a few extra files to be safe
            for i_extra in range(n_extra):
                fh.append(open(wd_files_test + "Extra%d.txt" % i_extra, 'w'))
            # close the files again and delete
            for fhh in fh:
                fhh.close()
            DeleteDirectoryTree(wd_files_test)
        except IOError as e:
            if str(e).startswith("[Errno 24] Too many open files"):
                util.number_open_files_exception_advice(len(speciesToUse), False)
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                util.Fail()
            else:
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                print("ERROR: Attempted to open required files for OrthoFinder run but an unexpected error occurred. \n\nStacktrace:")
                raise
    return options

def main(args=None):    
    try:
        if args is None:
            args = sys.argv[1:]
       # Create PTM right at start
        ptm_initialised = parallel_task_manager.ParallelTaskManager_singleton()
        print("")
        print(("OrthoFinder version %s Copyright (C) 2014 David Emms\n" % util.version))
        prog_caller = GetProgramCaller()
        
        options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault, user_specified_M = ProcessArgs(prog_caller, args)
        
        files.InitialiseFileHandler(options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault)     
        print("Results directory: %s" % files.FileHandler.GetResultsDirectory1())

        CheckDependencies(options, user_specified_M, prog_caller, files.FileHandler.GetWorkingDirectory1_Read()[0])
            
        # if using previous Trees etc., check these are all present - Job for orthologues
        if options.qStartFromBlast and options.qStartFromFasta:
            # 0. Check Files
            speciesInfoObj, speciesToUse_names = ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast)
            print("\nAdding new species in %s to existing analysis in %s" % (fastaDir, continuationDir))
            # 3. 
            speciesInfoObj = ProcessesNewFasta(fastaDir, options.dna, speciesInfoObj, speciesToUse_names)
            files.FileHandler.LogSpecies()
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4.
            seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
            # 5.
            speciesXML = GetXMLSpeciesInfo(speciesInfoObj, options) if options.speciesXMLInfoFN else None
            # 6.    
            util.PrintUnderline("Dividing up work for BLAST for parallel processing")
            CreateSearchDatabases(speciesInfoObj, options, prog_caller)
            # 7.  
            RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
            # 8.
            speciesNamesDict = SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
            gathering.DoOrthogroups(options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML)
            # 9.
            if not options.qStopAfterGroups:
                GetOrthologues(speciesInfoObj, options, prog_caller)   
        elif options.qStartFromFasta:
            # 3. 
            speciesInfoObj = None
            speciesInfoObj = ProcessesNewFasta(fastaDir, options.dna)
            files.FileHandler.LogSpecies()
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4
            seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
            # 5.
            speciesXML = GetXMLSpeciesInfo(speciesInfoObj, options) if options.speciesXMLInfoFN else None
            # 6.    
            util.PrintUnderline("Dividing up work for BLAST for parallel processing")
            CreateSearchDatabases(speciesInfoObj, options, prog_caller)
            # 7. 
            RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
            # 8.
            speciesNamesDict = SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
            gathering.DoOrthogroups(options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML)
            # 9. 
            if not options.qStopAfterGroups:
                GetOrthologues(speciesInfoObj, options, prog_caller)
        elif options.qStartFromBlast:
            # 0.
            speciesInfoObj, _ = ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast)
            files.FileHandler.LogSpecies()
            print("Using previously calculated BLAST results in %s" % (files.FileHandler.GetWorkingDirectory1_Read()[0]))
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4.
            seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
            # 5.
            speciesXML = GetXMLSpeciesInfo(speciesInfoObj, options) if options.speciesXMLInfoFN else None
            # 8
            speciesNamesDict = SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
            gathering.DoOrthogroups(options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML)
            # 9
            if not options.qStopAfterGroups:
                GetOrthologues(speciesInfoObj, options, prog_caller)
        elif options.qStartFromGroups:
            # 0.
            check_blast = not options.qMSATrees
            speciesInfoObj, _ = ProcessPreviousFiles(continuationDir, options.qDoubleBlast, check_blast=check_blast)
            files.FileHandler.LogSpecies()
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            # 9
            GetOrthologues(speciesInfoObj, options, prog_caller)
        elif options.qStartFromTrees:
            speciesInfoObj, _ = ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast, check_blast=False)
            files.FileHandler.LogSpecies()
            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            GetOrthologues_FromTrees(options)
        elif options.qFastAdd:
            # Prepare previous directory as database
            speciesInfoObj, speciesToUse_names = ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast, check_blast=False)
            # Check previous directory has been done with MSA trees
            if not acc.check_for_orthoxcelerate(continuationDir, speciesInfoObj):
                util.Fail()
            util.PrintUnderline("Creating orthogroup profiles")
            wd_list = files.FileHandler.GetWorkingDirectory1_Read()
            fn_diamond_db, q_hogs = acc.prepare_accelerate_database(continuationDir, wd_list, speciesInfoObj.nSpAll)
            print("\nAdding new species in %s to existing analysis in %s" % (fastaDir, continuationDir))
            speciesInfoObj = ProcessesNewFasta(fastaDir, options.dna, speciesInfoObj, speciesToUse_names)

            options = CheckOptions(options, speciesInfoObj.speciesToUse)
            seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
            # Add genes to orthogroups
            results_files = acc.RunSearch(options, speciesInfoObj, fn_diamond_db, prog_caller)

            # Clade-specific genes
            orphan_genes_version = 2
            speciesNamesDict = SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
            # if orphan_genes_version == 1:
            #     # v1 - This is unsuitable, it does an all-v-all search of all unassigned genes. Although these should have
            #     # been depleted of all genes that are not clade-specific, the resulting search still takes too long.
            #     clustersFilename_pairs, i_og_restart = clade_specific_orthogroups_v1(speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files, q_hogs)
            #     raise Exception("If q_hjogs then should be reading the N0.tsv file, not the original clusters")
            #     gathering.post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo, speciesNamesDict, options, speciesXML=None)
            if orphan_genes_version == 2:
                # v2 - Infer rooted species tree from new rooted gene trees, identify new species-clades & search within these
                clustersFilename_pairs, i_og_restart = BetweenCoreOrthogroupsWorkflow(continuationDir, speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files, q_hogs)

                # Infer clade-specific orthogroup gene trees
                gathering.post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo,
                                                      speciesNamesDict, options, speciesXML=None)
                if options.speciesTreeFN is None:
                    # No user species tree, use the one we've just inferred
                    options.speciesTreeFN = files.FileHandler.GetSpeciesTreeResultsFN(None, True)
            if not options.qStopAfterGroups:
                GetOrthologues(speciesInfoObj, options, prog_caller, i_og_restart)
        else:
            raise NotImplementedError
            ptm = parallel_task_manager.ParallelTaskManager_singleton()
            ptm.Stop()
        if not options.save_space and not options.qFastAdd:
            # split up the orthologs into one file per species-pair
            split_ortholog_files.split_ortholog_files(files.FileHandler.GetOrthologuesDirectory())
        d_results = os.path.normpath(files.FileHandler.GetResultsDirectory1()) + os.path.sep
        print("\nResults:\n    %s" % d_results)
        util.PrintCitation(d_results)
        files.FileHandler.WriteToLog("OrthoFinder run completed\n", True)
    except Exception as e:
        print(str(e))
        parallel_task_manager.print_traceback(e)
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()
        raise
    ptm = parallel_task_manager.ParallelTaskManager_singleton()
    ptm.Stop()
