#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 16:24:50 2016

@author: david
"""

"""
Combinations to run:
- each main arg with
- different combinations of:
    -file arrangements
    - what order commands are called on a directory structure
    - which of the allowable ways of specify the directory/file is used
    (These are what I need to test, not different combinations of with/without -t argument etc)

1. Each of the individual commands

"""
import os
import shutil
import datetime
import subprocess

baseDir = os.path.dirname(os.path.realpath(__file__)) + os.sep
qBinary = False
orthofinder = baseDir + "../orthofinder/orthofinder.py"
my_env = os.environ.copy()
my_env["PATH"] = "/home/david/software/ncbi-blast-2.2.28+/bin:" + my_env["PATH"]

workspace = baseDir + "ArgumentCombinations/"
fileStore = workspace + "FileStore/"
fasta =  "Mycoplasma/"
fastaExtra = "Mycoplasma_extra/"
dirs = [fasta, fastaExtra]

d_fasta =  workspace + fasta
d_fastaExtra =  workspace + fastaExtra
#blastDir = ""
#groupsDir = ""
#orthologuesDir = ""

#def RunOrthoFinder(commands):
#    capture = subprocess.Popen("python %s %s" % (orthofinder, commands), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
#    stdout = "".join([x for x in capture.stdout])
#    stderr = "".join([x for x in capture.stderr])
#    return stdout, stderr

def Date():
    return datetime.date.today().strftime("%b%d")

class Workspace():
    def __init__(self):
        pass
    def __enter__(self):
        for d in dirs:
            shutil.copytree(fileStore + d, workspace + d)
    def __exit__(self, type, value, traceback):
        for d in dirs:
            shutil.rmtree(workspace + d)

def CreateCleanDirectories():
    for d in dirs:
        if os.path.exists(workspace + d):
            shutil.rmtree(workspace + d)
    for d in dirs:
        shutil.copytree(fileStore + d, workspace + d)
    
def PrintHeading(message):
    print("\n\n\n\n\n****** " + message + " ******\n")

def run_combinations():
    """
    Run through a set of OrthoFinder analyses
    """
    CreateCleanDirectories()
    
    PrintHeading("Get version")
    capture = subprocess.Popen(["python", orthofinder, "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
    stdout = "".join([x for x in capture.stdout])   
    iVersion = stdout.find("version") + len("version") + 1
    version = stdout[iVersion:].split(" ", 1)[0]
    
#    with Workspace():
    # Input: 4 fasta files
    # Output: ortholgoues with 4 possible roots for the species tree
    PrintHeading("Running from FASTA files")
    subprocess.call(["python", orthofinder, "-f", d_fasta], env=my_env)
    d_ogs = d_fasta + "Results_%s/" % Date() 
    d_orths = d_ogs + "Orthologues_%s/" % Date() 
    
    
    # Input: From trees from one of the 4 possible roots
    PrintHeading("Running from Trees (4 possible species tree roots, one specified)")
#    raise Exception("Running from Trees (4 possible species tree roots, one specified)")
    subprocess.call(["python", orthofinder, "-ft", d_orths, "-s", d_orths + "Orthologues_using_outgroup_2/SpeciesTree_rooted_at_outgroup_2.txt"], env=my_env)
    
    PrintHeading("Extra Species")
    subprocess.call(["python", orthofinder, "-b", d_ogs + "WorkingDirectory", "-f", d_fastaExtra], env=my_env)
    
    PrintHeading("From groups, only write sequences and exit")
    subprocess.call(["python", orthofinder, "-os", "-fg", d_ogs + "WorkingDirectory/clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt" % version], env=my_env)
    
    PrintHeading("Rerun from BLAST but stop at groups")
    subprocess.call(["python", orthofinder, "-b", d_ogs + "WorkingDirectory", "-og"], env=my_env)
    
    PrintHeading("From groups (specify file, specify species tree)")
    subprocess.call(["python", orthofinder, "-fg", d_ogs + "WorkingDirectory/clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt" % version, "-s", d_ogs + "WorkingDirectory/Orthologues_%s/Orthologues_using_outgroup_1/SpeciesTree_rooted_at_outgroup_1.txt" % Date()], env=my_env)

    PrintHeading("From groups, use MSAs (specify file, specify species tree)")
    subprocess.call(["python", orthofinder, "-M", "msa", "-fg", d_ogs + "WorkingDirectory/clusters_OrthoFinder_v%s_I1.5.txt_id_pairs.txt" % version, "-s", d_ogs + "WorkingDirectory/Orthologues_%s/Orthologues_using_outgroup_1/SpeciesTree_rooted_at_outgroup_1.txt" % Date()], env=my_env)
 
       
    
    
    
#    main_args = ["-f " + fastaDir,
#                 "-b " + blastDir,
#                 "-f " + fastaDir + " -b " + blastDir,
#                 "-fg " + groupsDir,
#                 "-ft " + orthologuesDir]
#                 #

if __name__ == "__main__":
    run_combinations()