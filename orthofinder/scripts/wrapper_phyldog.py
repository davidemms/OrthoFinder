# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 15:20:28 2016

@author: david
"""

import os
from collections import defaultdict

def WriteGeneralOptions(filename, baseDir):
    x="""######## First, data files ########
DATA=OG0000000
BASEDIR=%s

RESULT=$(BASEDIR)phyldog/
PATH=$(RESULT)

genelist.file=$(RESULT)listGenes_generic.txt
init.species.tree=user
#species.tree.file=$(BASEDIR)Trees_ids/SpeciesTree_ids_0_rooted.txt
species.tree.file=/home/david/projects/OrthoFinder/Development/phyldog/SpeciesTree_ids_0_rooted_from_elsewhere.txt
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
taxaseq.file=$(RESULT)/$(DATA).map.txt
input.sequence.file=$(BASEDIR)Alignments_ids/$(DATA).fa
input.sequence.format=Fasta
output.reconciled.tree.file=$(RESULT)$(DATA).ReconciledTree
output.duplications.tree.file=$(RESULT)$(DATA).DuplicationTree
output.losses.tree.file=$(RESULT)$(DATA).LossTree
#output.numbered.tree.file=$(RESULT)OutputSpeciesTree_ConsensusNumbered.tree

use.quality.filters=false""" % baseDir
    with open(filename, 'wb') as outfile: outfile.write(x)

def WriteOGOptions(filename):
    x = """######## First, data files ########

#DATA= defined on command line

input.sequence.sites_to_use=all
input.sequence.max_gap_allowed=66%
init.gene.tree=bionj

######## Second, model options ########
alphabet=Protein
model=LG08+F(initFreqs=observed )
rate_distribution=Gamma(n=4,alpha=1)

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
    with open(filename, 'wb') as outfile: outfile.write(x)
    
def WriteListSpecies(filename, speciesToUse):
    with open(filename, 'wb') as outfile:
        for i in speciesToUse:
            outfile.write("%d\n" % i)

def WriteGeneMaps(outputDir, ogs):
    for i, og in enumerate(ogs):
        genesForSpecies = defaultdict(list)
        for seq in og:
            name = seq.ToString()
            genesForSpecies[name.split("_")[0]].append(name)
        with open(outputDir + "OG%07d.map.txt" % i, 'wb') as outfile:
            for species, genes in genesForSpecies.items():
                outfile.write("%s:%s\n" % (species, ";".join(genes)))

#def WriteGeneMaps(phyldogDir, ogs):
#    for i, og in enumerate(ogs):
#        with open(            
    
def WriteStandardFiles(phyldogDir, speciesToUse):
    WriteGeneralOptions(phyldogDir + "GeneralOptions.opt", phyldogDir + "../")
    with open(phyldogDir + "listGenes_generic.txt", 'wb') as outfile: outfile.write(phyldogDir + "OG_generic.opt:1")
    WriteOGOptions(phyldogDir + "OG_generic.opt")
    WriteListSpecies(phyldogDir + "ListSpecies.txt", speciesToUse)
    
def Setup(phyldogDir, ogs, speciesToUse):
    if not os.path.exists(phyldogDir): os.mkdir(phyldogDir)
    WriteStandardFiles(phyldogDir, speciesToUse)
    WriteGeneMaps(phyldogDir, ogs)
    
def RunPhyldogAnalysis(phyldogDir, ogs, speciesToUse):
    Setup(phyldogDir, ogs, speciesToUse)