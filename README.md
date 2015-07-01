# OrthoFinder version 0.2.6 - Accurate inference of orthologous gene groups made easy!
Emms, D.M. and Kelly, S. (in press), OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthologous gene group inference accuracy, Genome Biology 

Available from:

http://www.stevekellylab.com/software/orthofinder

https://github.com/davidemms/OrthoFinder


Usage
=====
OrthoFinder runs as a single command that takes as input a directory of fasta files, one per species, and outputs a file containing the orthologous groups of genes from these species. 

python orthofinder.py -f fasta_directory -t number_of_processes

For example, if you want to run it using 16 processors in parallel on the example dataset move to the directory containing orthofinder.py and call:

python orthofinder.py -f ExampleDataset -t 16

Once complete your results will be in ExampleDataset/Results_\<date\>/OrthologoueGroups.txt

For details on running it from pre-computed BLAST search results see below.

Installing Dependencies
=======================
OrthoFinder is written to run on linux and requires the following to be installed and in the system path:
1. python 2.x (version 3 isn't currently supported) together with the scipy libraries stack 
2. BLAST+ 
3. The MCL graph clustering algorithm 

Brief instructions are given below although users may wish to refer to the installation notes provided with these packages for more detailed instructions. BLAST+ and MCL must be in the system path so that the can be OrthoFinder can run them, instructions are provided below.

python and scipy
----------------
Up-to-date and clear instructions are provided here: http://www.scipy.org/install.html.
As websites can change, an alternative is to search online for "install scipy".

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
1. Save OrthoFinder_v0.2.6.tar.gz 
2. Open a terminal and cd to the directory where OrthoFinder_v0.2.6.tar.gz was saved
3. tar xzf OrthoFinder_v0.2.6.tar.gz
4. cd OrthoFinder_v0.2.6
5. python orthofinder.py -f ExampleDataset/
6. If everything was successful the output generated will end with a line giving the location of the results file containing the orthologous groups.

The command for running OrthoFinder on any dataset is:

python orthofinder.py -f directory_containing_fasta_files -t number_of_processes

where:
directory_containing_fasta_files is a directory containing the fasta files (with filename extensions .fa or .fasta) for the species of interest, one species per fasta file. It is best to use a fasta file with only the longest transcript variant per gene.
number_of_processes is an optional argument that can be used to run the initial BLAST all-versus-all queries in parallel. As the BLAST queries are by far the time-consuming step it is best to use at least as many BLAST processes as there are CPUs on the machine. 

Using pre-computed BLAST results
================================
It is possible to run OrthoFinder with pre-computed BLAST results provided they are in the format specified below. The command is simply:

python orthofinder.py -b directory_with_processed_fasta_and_blast_results

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
>0_0
MFAPRGK...

>0_1
MFAVYAL...

>0_2
MTTIID...

And the first few lines of start of Species1.fa would look like
>1_0
MFAPRGK...

>1_1
MFAVYAL...

>1_2
MTTIID...

BLAST results files
-------------------
There should be BLAST results files with names of the form Blastx_y.txt and Blasty_x.txt for species pair x, y in BLAST output format 6. The query and hit IDs in the BLAST results files should correspond to the IDs in the fasta files.

Aside, reducing BLAST computations
Note that since the BLAST queries are by far the most computationally expensive step, considerable time could be saved by only performing n(n+1)/2 of the species versus species BLAST queries instead of n^2, where n is the number of species. This would be done by only searching Species<x>.fa against the BLAST database generated from Species<y>.fa if x <= y. The results would give the file BLASTx_y.txt and then this file could be used to generate the BLASTy_x.txt file by swapping the qury and hit sequence on each line in the results file. This should have only a small effect on the generated orthologous groups.

SequenceIDs.txt
--------------- 
The SequenceIDs.txt give the translation from the IDs of the form x_y to the original accessions. An example line would be:
0_42: Magnaporthe_oryzae_MGG_08114T0  
The IDs should be in order, i.e.
0_0: ...
0_1: ...
0_2: ...
...
...
1_0: ...
1_1: ...

SpeciesID.txt
-------------
The SpeciesIDs.txt file gives the translation from the IDs for the species to the original fasta file, e.g.:
0: /home/david/FungiFastaSmall/maor.fasta
1: /home/david/FungiFastaSmall/sace.fasta
2: /home/david/FungiFastaSmall/scpo.fasta
3: /home/david/FungiFastaSmall/sppu.fasta
4: /home/david/FungiFastaSmall/usma.fasta

Orthobench with pre-computer BLAST results 
------------------------------------------
The BLAST pre-calculated BLAST results files etc. for the Orthobench dataset are available for download as are the original fasta files. 


Output orthologous groups using the orthoxml format
===================================================
Orthologous groups can be output using the orthoxml format. This is requested by adding '-x speciesInfoFilename' to the command used to call orthofinder, where speciesInfoFilename should be the filename (including the path if necessary) of a user prepared file providing the information about the species that is required by the orthoxml format. This file should contain one line per species and each line should contain the following 5 fields separated by tabs:

fasta filename: the filename (without path) of the fasta file for the species described on this line
species name: the name of the species
NCBI Taxon ID: the NCBI taxon ID for the species
source database name: the name of the database from which the fasta file was obtained (e.g. Ensembl)
database fasta filename: the name given to the fasta file by the database (e.g. Homo_sapiens.NCBI36.52.pep.all.fa)

so a single line of the file could look like this (where each field has been separated by a tab rather than just spaces):
HomSap.fa	Homo sapiens 36	Ensembl	Homo_sapiens.NCBI36.52.pep.all.fa

Information on the orthoxml format can be found here: http://orthoxml.org/0.3/orthoxml_doc_v0.3.html
