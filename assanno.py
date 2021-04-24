
  
#!/usr/bin/env python
"""
assanno is a command line program that takes forward reads(in fastq format),reverse reads (in fastq format), reference genome and reference annotation . The reads are processed with trimmomatic and a quality reports are generated using FastQC.  The paired reads from trimmomatic are used for different modes of spades and megahit assembly. Then the assemblies are evaluated using quast and the best assembly is  selected . Then the assembled genome is annotated using prokka.
Author: TeamC (Monish Kumar,Rumana Mustafa,Sakina Mandviwala, Ethan R Mason) ,Jessie Arce
"""
import sys
import argparse
import subprocess
import os

import assemble_nd_anno


def main():
    parser = argparse.ArgumentParser(description="Assembles and annotates prokaryotic genome ")
    parser.add_argument("F",
                        type=str,
                        metavar="<forward reads>",
                        help="Forward reads in FastQ format."
                        )
    parser.add_argument("R",
                        type=str,
                        metavar="<reverse reads>",
                        help="Reverse reads in FastQ format."
                        )
    parser.add_argument("outdir",
                        type=str,
                        metavar="<output directory>",
                        help="ass&anno results"
                        )
    parser.add_argument("r_genome",
                        type=str,
                        metavar="<ref_genome>",
                        help="Name of the reference genome file for quast.")

    parser.add_argument("r_anno",
                        type=str,
                        metavar="<ref_anno>",
                        help="Name of the reference annotation file for quast.")

    args = parser.parse_args()


    assemble_nd_anno.assemble_nd_anno(
        args.F,
        args.R,
        args.outdir,
        args.r_genome,
        args.r_anno,
    )

    
if __name__ == "__main__":
    main()
