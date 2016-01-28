# OrthoFinder - Accurate inference of orthologous gene groups made easy!
**Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy, Genome Biology 16:157**

Available from:

http://www.stevekellylab.com/software/orthofinder

https://github.com/davidemms/OrthoFinder

What's New
==========
**Jan. 2016**: Added the ability to **add and remove species**.

**Sept. 2015**: Added the **trees_for_orthogroups.py** utility to automatically calculate multiple sequence alignments and gene trees for the orthogroups calcualted using OrthoFinder.

Usage
=====
OrthoFinder runs as a single command that takes as input a directory of fasta files, one per species, and outputs a file containing the orthologous groups of genes from these species. 

**python orthofinder.py -f fasta_directory -t number_of_processes**

For example, if you want to run it using 16 processors in parallel on the example dataset move to the directory containing orthofinder.py and call:

**python orthofinder.py -f ExampleDataset -t 16**

Once complete your results will be in ExampleDataset/Results_\<date\>/

See below for details on:
- adding extra species to a previous analysis
- quickly running OrthoFinder on a subset of species from a previous analysis
- running OrthoFinder from pre-computed BLAST search results
- preparing files in the format required by OrthoFinder so you can run the BLAST searches yourself

Output File Format
==================
OrthoFinder generates three output files 

**1) OrthologousGroups.csv** is a comma separated text file. Each row comprises a single orthogroup and contains all the genes that belong to that orthogroup. The genes are organized into separate columns where each column corresponds to a single species.

**2) OrthologousGroups.txt** is a tab separated text file that is identical in format to the output file from OrthoMCL. This enables OrthoFinder to easily slot into existing bioinformatic pipelines.

**3) OrthologousGroups_UnassignedGenes.csv** is a comma separated text file that is identical in format to OrthologousGroups.csv but contains all of the genes that were not assigned to any orthogroup.


Installing Dependencies
=======================
OrthoFinder is written to run on linux and requires the following to be installed and in the system path:

1. Python 2.7 (version 3 isn't currently supported) together with the scipy libraries stack 

2. BLAST+ 

3. The MCL graph clustering algorithm 

To use the trees_for_orthogroups.py utility there are two additional dependencies which should be installed and in the system path:

1. MAFFT

2. FastTree 

Brief instructions are given below although users may wish to refer to the installation notes provided with these packages for more detailed instructions. BLAST+ and MCL must be in the system path so that the can be OrthoFinder can run them, instructions are provided below.

python and scipy
----------------
Up-to-date and clear instructions are provided here: http://www.scipy.org/install.html, be sure to chose a version using python 2.7. As websites can change, an alternative is to search online for "install scipy".

BLAST+
------
Executables are found here ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ (instructions are currently at http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/ and in more detail in the BLAST+ user manual). As websites can change, an alternative is to search online for "install BLAST+". 

1. Instructions are provided for installing BLAST+ on various flavours of Linux on the 'Standalone BLAST Setup for Unix' page of the BLAST+ Help manual currently at http://www.ncbi.nlm.nih.gov/books/NBK1762/. 
2. Follow the instructions under "Configuration" in the BLAST+ help manual to add BLAST+ to the PATH environment variable.

MCL
---
mcl is available in the repositories for some linux distributions and so can be installed in the same way as any other package. E.g. on Ubuntu "sudo apt-get install mcl". Alternatively it can be built from source which will likely require the build-essential or equivalent package on the Linux distribution being used. Instructions are provided on the MCL webpage.  


Setting up and running OrthoFinder
==================================
Once the required dependencies have been installed, OrthoFinder can be setup and run on the small example data included in the package as follows:

1. Save OrthoFinder-master.zip and unpack it 
2. Open a terminal and cd into the directory OrthoFinder-master
3. python orthofinder.py -f ExampleDataset/
4. If everything was successful the output generated will end with a line giving the location of the results file containing the orthologous groups.

The command for running OrthoFinder on any dataset is:

**python orthofinder.py -f directory_containing_fasta_files -t number_of_processes**

where:
directory_containing_fasta_files is a directory containing the fasta files (with filename extensions .fa or .fasta) for the species of interest, one species per fasta file. It is best to use a fasta file with only the longest transcript variant per gene.
number_of_processes is an optional argument that can be used to run the initial BLAST all-versus-all queries in parallel. As the BLAST queries are by far the time-consuming step it is best to use at least as many BLAST processes as there are CPUs on the machine. 

Adding extra species
====================
OrthoFinder allows you to add extra species to an analysis without re-running the time-consuming BLAST searches:

**python orthofinder.py -b previous_orthofinder_directory_containing_blast_results -f new_fasta_directory**

Will add each species from the new_fasta_directory to existing set of species, reuse all the previous BLAST results, perform only the new BLAST searches required for the new species and recalculate the orthogroups.

Removing Species
================
OrthoFinder allows you to remove species from a previous analysis. In the WorkingDirectory from a previous analysis there is a file called SpeciesIDs.txt. Comment out any species to be removed from the analysis using a '#' character and then run OrthoFinder using: 

**python orthofinder.py -b previous_orthofinder_directory_containing_blast_results**

Preparing files without running BLAST or calculating orthogroups
================================================================
The '-p' option will prepare the files in the format required by OrthoFinder and print the set of BLAST searches that need to be performed. This is useful if you want to manage the BLAST searches yourself. For example, you may want to distribute them across multiple machines. When the BLAST searches have been completed then orthogroups can be calcualted as described in the section, "Using pre-computed BLAST results".

Using pre-computed BLAST results
================================
It is possible to run OrthoFinder with pre-computed BLAST results provided they are in the correct format. The command is simply:

**python orthofinder.py -b directory_with_processed_fasta_and_blast_results**

The correctly formatted files can be prepared as described in the section, "Preparing files without running BLAST or calculating orthogroups". The files from a successful run of OrthoFinder are also in the required format and so can be used for follow-on analyses. 

These above two methods are strongly recommended for preparing files for use with the "-b" option. For completeness, the required files and their formats are described below.

The files that must be in directory_with_processed_fasta_and_blast_results are:
- a fasta file for each species
- BLAST results files for each species pair
- SequenceIDs.txt
- SpeciesIDs.txt

Examples of the format required for the files can be seen by running the supplied test dataset and looking in the working directory created. A description is given below.

Fasta files
-----------
Species0.fa  
Species1.fa  
etc.  

Within each fasta file the accessions for the sequences should be of the form x_y where x is the species ID, matching the number in the filename and y is the sequence ID starting from 0 for the sequences in each species. 

So the first few lines of start of Species0.fa would look like
```
>0_0
MFAPRGK...

>0_1
MFAVYAL...

>0_2
MTTIID...
```

And the first few lines of start of Species1.fa would look like
```
>1_0
MFAPRGK...

>1_1
MFAVYAL...

>1_2
MTTIID...
```

BLAST results files
-------------------
For each species pair x, y there should be a BLAST results file Blastx_y.txt where x is the index of the query fasta file and y is the index of the species used for the database. Similarly, there should be a BLAST results file Blasty_x.txt where y is the index of the query fasta file and x is the index of the species used for the database. The tabular BLAST output format 6 should be used. The query and hit IDs in the BLAST results files should correspond to the IDs in the fasta files.

**Aside, reducing BLAST computations:** Note that since the BLAST queries are by far the most computationally expensive step, considerable time could be saved by only performing n(n+1)/2 of the species versus species BLAST queries instead of n^2, where n is the number of species. This would be done by only searching Species<x>.fa against the BLAST database generated from Species<y>.fa if x <= y. The results would give the file Blastx_y.txt and then this file could be used to generate the Blasty_x.txt file by swapping the query and hit sequence on each line in the results file. This should have only a small effect on the generated orthologous groups.

SequenceIDs.txt
--------------- 
The SequenceIDs.txt give the translation from the IDs of the form x_y to the original accessions. An example line would be:
```
0_42: gi|290752309|emb|CBH40280.1| 
```
The IDs should be in order, i.e.
```
0_0: gi|290752267|emb|CBH40238.1|
0_1: gi|290752268|emb|CBH40239.1|
0_2: gi|290752269|emb|CBH40240.1|
...
...
1_0: gi|284811831|gb|AAP56351.2|
1_1: gi|284811832|gb|AAP56352.2|
...
```

SpeciesID.txt
-------------
The SpeciesIDs.txt file gives the translation from the IDs for the species to the original fasta file, e.g.:
```
0: Mycoplasma_agalactiae_5632_FP671138.faa
1: Mycoplasma_gallisepticum_uid409_AE015450.faa
2: Mycoplasma_genitalium_uid97_L43967.faa
3: Mycoplasma_hyopneumoniae_AE017243.faa
```

Orthobench with pre-computed BLAST results 
------------------------------------------
The BLAST pre-calculated BLAST results files etc. for the Orthobench dataset are available for download as are the original fasta files. 


Output orthologous groups using the orthoxml format
===================================================
Orthologous groups can be output using the orthoxml format. This is requested by adding '-x speciesInfoFilename' to the command used to call orthofinder, where speciesInfoFilename should be the filename (including the path if necessary) of a user prepared file providing the information about the species that is required by the orthoxml format. This file should contain one line per species and each line should contain the following 5 fields separated by tabs:

1. **fasta filename**: the filename (without path) of the fasta file for the species described on this line
2. **species name**: the name of the species
3. **NCBI Taxon ID**: the NCBI taxon ID for the species
4. **source database name**: the name of the database from which the fasta file was obtained (e.g. Ensembl)
5. **database fasta filename**: the name given to the fasta file by the database (e.g. Homo_sapiens.NCBI36.52.pep.all.fa)

so a single line of the file could look like this (where each field has been separated by a tab rather than just spaces):
```
HomSap.fa	Homo sapiens 36	Ensembl	Homo_sapiens.NCBI36.52.pep.all.fa
```

Information on the orthoxml format can be found here: http://orthoxml.org/0.3/orthoxml_doc_v0.3.html

Trees for Orthogroups
=====================
The 'trees_for_orthogroups.py' utility will automatically generate multiple sequence alignments and gene trees for each orthologous group generated by OrthoFinder. For example, once OrthoFinder has been run on the example dataset, trees_for_orthogroups can be run using:

**python trees_for_orthogroups.py ExampleDataset/Results_\<date\> -t 16**

where "ExampleDataset/Results_\<date\>" is the results directory from running OrthoFinder.

This will use MAFFT to generate the multiple sequence alignments and FastTree to generate the gene trees. Both of these programs need to be installed and in the system path.
