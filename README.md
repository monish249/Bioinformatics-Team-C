# Bioinformatics-Team-C 

contents:
assemble_nd_anno.py
assanno.py

A command line program that will take bacterial forward and reverse reads from illumina sequencer. The program will trim the reads, assemble the genome using spades and megahit.

The program also takes reference genome and annotation as input and compares the assembles and chooses the best assembly. 

The resulting contigs are annotated using prokka

Installation:

conda create -n prok_genom_assanno
conda activate prok_genom_assanno
conda install -c conda-forge -c bioconda -c defaults prokka
conda install -c conda-forge -c bioconda biopython spades megahit quast fastqc trimmomatic unzip java-jdk --yes

we recommend to follow the exact order of installation in case you get "Can't locate Bio/Root/Version.pm in @inc" error while running prokka

copy and paste assanno.py and assemble_nd_anno.py

python assanno.py forward_reads reverse_reads output_dir reference_genome reference_annotation

