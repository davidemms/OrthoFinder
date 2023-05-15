# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 15:20:28 2016

@author: david
"""

import os
import sys
import time
import subprocess
import fileinput
from collections import defaultdict, Counter

from . import util
from . import tree as tree_lib
from . import files, program_caller

def WriteGeneralOptions(filename, baseDir, qRunSingley, iogs4):
    x="""######## First, data files ########
BASEDIR=%s

RESULT=$(BASEDIR)phyldog/
PATH=$(RESULT)

genelist.file=$(RESULT)ListGenes.opt
init.species.tree=mrp
species.tree.file=$(BASEDIR)/Trees_ids/SpeciesTree_user_ids_rooted.txt
species.names.file=$(RESULT)ListSpecies.txt
starting.tree.file=$(RESULT)StartingTree.tree
output.tree.file=$(RESULT)OutputSpeciesTree.tree
output.duplications.tree.file=$(RESULT)OutputSpeciesTree_ConsensusDuplications.tree
output.losses.tree.file=$(RESULT)OutputSpeciesTree_ConsensusLosses.tree
output.numbered.tree.file=$(RESULT)OutputSpeciesTree_ConsensusNumbered.tree

######## Second, options ########
optimization.topology=no
branchProbabilities.optimization=average_then_branchwise
branch.expected.numbers.optimization=average_then_branchwise
spr.limit=5
time.limit=10000

### Added to remove warnings ###
reconciliation.model=DL
output.file.suffix=.txt
debug=0

# From specific file ... but required the variables

input.sequence.format=Fasta
output.reconciled.tree.file=$(RESULT)$(DATA).ReconciledTree
output.duplications.tree.file=$(RESULT)$(DATA).DuplicationTree
output.losses.tree.file=$(RESULT)$(DATA).LossTree
#output.numbered.tree.file=$(RESULT)OutputSpeciesTree_ConsensusNumbered.tree

use.quality.filters=false""" % baseDir
    if qRunSingley:
        for i in iogs4:
            base, ext = os.path.splitext(filename)
            og = "OG%07d" % i
            outFN = base + "_" + og + ext
            with open(outFN, 'w') as outfile:
                outfile.write(x.replace("ListGenes", "ListGenes_" + og).replace("mrp", "user"))
    else:
        with open(filename, 'w') as outfile: outfile.write(x)

def WriteOGOptions(phyldogDir, iogs4, exclude):
    basedir = phyldogDir + "../"
    x = """######## First, data files ########

BASEDIR=%s
RESULT=$(BASEDIR)phyldog/Results/
DATA=%s

taxaseq.file=$(BASEDIR)phyldog/$(DATA).map.txt
input.sequence.file=$(BASEDIR)Alignments_ids/$(DATA).fa

input.sequence.sites_to_use=all
input.sequence.max_gap_allowed=66%%
init.gene.tree=bionj

######## Second, model options ########
alphabet=Protein
model=LG08

######## Output options #########
gene.tree.file=$(RESULT)$(DATA).GeneTree
output.reconciled.tree.file=$(RESULT)$(DATA).ReconciledTree
output.duplications.tree.file=$(RESULT)$(DATA).DuplicationTree
output.losses.tree.file=$(RESULT)$(DATA).LossTree
output.numbered.tree.file=$(RESULT)$(DATA).NumberedTree

######## Finally, optimization options ########
optimization.topology=yes
optimization.topology.algorithm_nni.method=fast
optimization.tolerance=0.01
optimization.method_DB.nstep=0
optimization.topology.numfirst=false
optimization.topology.tolerance.before=100
optimization.topology.tolerance.during=100
optimization.max_number_f_eval=1000000
optimization.final=none
optimization.verbose=0
optimization.message_handler=none
optimization.profiler=none
optimization.reparametrization=no"""
    exclude = set(exclude)
    for i in iogs4:
        if i in exclude: continue
        ogName = "OG%07d" % i
        with open(phyldogDir + ogName + ".opt", 'w') as outfile: 
            outfile.write(x % (basedir, ogName))
    
def WriteListSpecies(filename, speciesToUse):
    with open(filename, 'w') as outfile:
        for i in speciesToUse:
            outfile.write("%d\n" % i)

def WriteGeneMaps(outputDir, iogs4, ogs, exclude):
    exclude = set(exclude)
    for i in iogs4:
        og = ogs[i]
        if i in exclude: continue
        genesForSpecies = defaultdict(list)
        for seq in og:
            name = seq.ToString()
            genesForSpecies[name.split("_")[0]].append(name)
        with open(outputDir + "OG%07d.map.txt" % i, 'w') as outfile:
            for species, genes in genesForSpecies.items():
                outfile.write("%s:%s\n" % (species, ";".join(genes)))

#def WriteGeneMaps(phyldogDir, ogs):
#    for i, og in enumerate(ogs):
#        with open(          

def CleanAlignmentsForPhyldog(phyldogDir, iogs4, ogs):
    """
    Remove * character
    Remove any orthogroups composed entierly of identical sequences
    Return alignments to be excluded
    """
    # 1. Remove * character
    for i in iogs4:
        for line in fileinput.FileInput(phyldogDir + "../Alignments_ids/OG%07d.fa" % i, inplace=True):
            if not line.startswith(">"): line=line.replace("*","-")
            sys.stdout.write(line)
    # 2. Remove any orthogroups composed entierly of identical sequences
    exclude = []
    for i in iogs4:
        og = ogs[i]
        with open(phyldogDir + "../Alignments_ids/OG%07d.fa" % i, 'r') as infile:
            seqs = []
            for line in infile:
                if line.startswith(">"):
                    seqs.append("")
                else:
                    seqs[-1] += line.rstrip()
        # 2a. check at least 4 sequences are different
        c = Counter(seqs)
        if len(c) < 4: exclude.append(i)
    print(("%d excluded alignments" % len(exclude)))
    print("Running PHYLDOG")
    return set(exclude)
 
def ProcessSpeciesTree(phyldogDir):
    species_tree_rooted_fn = phyldogDir + "OutputSpeciesTree_ConsensusNumbered.tree.txt"
    # Label nodes of species tree
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    species_tree_rooted.name = "N0"    
    leaf_node_labels = dict()
    for n in species_tree_rooted.traverse():
        if n.is_root(): continue
        elif n.is_leaf():
            sp, node_name = n.name.split("_") 
            leaf_node_labels[sp] = node_name
            n.name = sp
        else:
            n.name = "N%s" % n.name
    ret_species_tree_fn = phyldogDir + "Species_tree_labelled.tre"
    species_tree_rooted.write(outfile=ret_species_tree_fn)
    return ret_species_tree_fn
    
def WriteStandardFiles(phyldogDir, speciesToUse, qRunSingley, iogs4):
    WriteGeneralOptions(phyldogDir + "GeneralOptions.opt", phyldogDir + "../", qRunSingley, iogs4)
#    with open(phyldogDir + "listGenes_generic.txt", 'w') as outfile: outfile.write(phyldogDir + "OG_generic.opt:1")
    WriteListSpecies(phyldogDir + "ListSpecies.txt", speciesToUse)

def WriteListGenes(phyldogDir, iogs4, exclude, qRunSingley):
    if qRunSingley:
        for i in iogs4:
            if i in exclude: continue
            with open(phyldogDir + "ListGenes_OG%07d.opt" % i, 'w') as outfile:
                    outfile.write(phyldogDir + "OG%07d.opt:%s\n" % (i, str(os.stat( phyldogDir + "../Alignments_ids/OG%07d.fa" % i )[6])))   # phyldog prepareData.py method
    
    else:
        with open(phyldogDir + "ListGenes.opt", 'w') as outfile:
            for i in range(nOGs):
                if i in exclude: continue
                outfile.write(phyldogDir + "OG%07d.opt:%s\n" % (i, str(os.stat( phyldogDir + "../Alignments_ids/OG%07d.fa" % i )[6])))   # phyldog prepareData.py method
    

def Setup(phyldogDir, iogs4, ogs, speciesToUse, qRunSingley):
    if not os.path.exists(phyldogDir): os.mkdir(phyldogDir)
    if not os.path.exists(phyldogDir + "Results/"): os.mkdir(phyldogDir + "Results/")
    WriteStandardFiles(phyldogDir, speciesToUse, qRunSingley, iogs4)
    exclude = CleanAlignmentsForPhyldog(phyldogDir, iogs4, ogs)
    WriteOGOptions(phyldogDir, iogs4, exclude)
    WriteGeneMaps(phyldogDir, iogs4, ogs, exclude)
    WriteListGenes(phyldogDir, iogs4, exclude, qRunSingley)
    
def RunPhyldogAnalysis(phyldogDir, iogs4, ogs, speciesToUse, nParallel):
    qRunSingley = True
    Setup(phyldogDir, iogs4, ogs, speciesToUse, qRunSingley)
    start = time.time()
    if qRunSingley:
        cmds = [["mpirun -np 2 phyldog param=%s%s" % (phyldogDir, "GeneralOptions_OG%07d.opt" % i)] for i in iogs4]
        program_caller.RunParallelCommands(nParallel, cmds, qListOfList=False)
    else:
        popen = subprocess.Popen("mpirun -np %d phyldog param=GeneralOptions.opt" % nParallel, shell=True, cwd=phyldogDir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        popen.communicate()
    stop = time.time()
    print(("%f seconds" % (stop-start)))
    return ProcessSpeciesTree(phyldogDir)