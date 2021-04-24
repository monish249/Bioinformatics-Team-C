import os
import subprocess
import csv
import re

def run_trimmomatic(forward_reads,reverse_reads,forward_paired_reads,forward_unpaired_reads,reverse_paired_reads,reverse_unpaired_reads):
    """
    Arguments:
     forward_reads: filename of forward reads
     reverse_reads: filename of reverse reads
     forward_paired_reads: trimmomatic result file for forward paired reads
     forward_unparied_reads: result file for forward unpaired reads
     reverse_paired_reads: result file for reverse paired reads
     reverse_paired_reads: result file for reverse unpaired reads
    Results:
     Four result files for forward and reverse reads both paired and unpaired.
    """

    command = [
        "trimmomatic","PE",
        forward_reads,reverse_reads,
        forward_paired_reads,forward_unpaired_reads,
        reverse_paired_reads,reverse_unpaired_reads,
        "LEADING:10","TRAILING:10","SLIDINGWINDOW:5:20"
    ]

    subprocess.run(command)
    
def run_spades_default(forward_reads,reverse_reads,outdir):
    """Runs spade.py
    Arguments:
     forward_reads: filename of forward reads
     reverse_reads: filename of reverse reads
    Results:
     spades.py assembly in spades_output directory.
    """
    command = ["spades.py","-1",forward_reads,"-2",reverse_reads,"--only-assembler","-o",outdir+"_spades_default_output"]
    subprocess.run(command)

def run_spades_careful(forward_reads,reverse_reads,outdir):
    """Runs spade.py in careful mode
    Arguments:
     forward_reads: filename of forward reads
     reverse_reads: filename of reverse reads
    Results:
     spades.py assembly in spades_output directory.
    """
    command = ["spades.py","-k","21,33,55,77,99",
               "-1",forward_reads,"-2",reverse_reads,
               "--careful","--only-assembler",
               "-o",outdir+"_spades_careful_output"]
    subprocess.run(command)
    
def megahit_default(forward_reads,reverse_reads,outdir):  
    """Runs megahit
    Arguments:
     forward_reads: filename of forward reads
     reverse_reads: filename of reverse reads
    Results:
     megahit assembly in megahit_default_output directory.
    """
    command = ["megahit","-1",forward_reads,"-2",reverse_reads,"-o",outdir+"_megahit_default_output"]
    subprocess.run(command)
    
def megahit_mincount3(forward_reads,reverse_reads,outdir):  
    """Runs megahit mincount 3
    Arguments:
     forward_reads: filename of forward reads
     reverse_reads: filename of reverse reads
    Results:
     megahit assembly in megahit_mincount_output directory.
    """
    command = ["megahit","-1",forward_reads,"-2",reverse_reads,"-o",outdir+"_megahit_mincount3_output","--min-count","3"]
    subprocess.run(command)
    
def run_fastqc(forward_reads,reverse_reads,outdir):
    """Runs FastQC on forward and reverse reads.
   Arguments:
     forward_reads: filename of forward reads
    reverse_reads: filename of reverse reads
     outdir: directory for FastQC quality reports.
   Result:
     FastQC quality reports.
    """
    resultdir = outdir+'_fastqc_report/'
    os.mkdir(resultdir)
    subprocess.run(["fastqc",forward_reads,reverse_reads,"-o",resultdir])
    result_file1 = forward_reads.replace('.fastq','_fastqc.zip')
    result_file2 = reverse_reads.replace('.fastq','_fastqc.zip')
    fpath1 = os.path.join(resultdir,result_file1)
    fpath2 = os.path.join(resultdir,result_file2)
    subprocess.run(["unzip",fpath1])
    subprocess.run(["unzip",fpath2])
    spath1 = result_file1.replace('.zip','')
    spath2 = result_file2.replace('.zip','')
    path_f = os.path.join(spath1,"summary.txt")
    print("Forward read warning\n")
    for line in open(path_f,'r'):
        if re.search("WARN",line):
            print(line)

    path_r = os.path.join(spath2,"summary.txt")
    print("Reverse read warning\n")
    for line in open(path_r,'r'):
        if re.search("WARN",line):
            print(line)

def run_prokka (path, outdir):
    """Runs Prokka 
    Arguments:  path - path to the scaffold i.e contigs.fasta
    """
    prodir = outdir + "_prokka"
    command = ["prokka","-outdir",prodir,path]
    subprocess.run(command)


        
def assemble_nd_anno(forward_reads,reverse_reads,outdir,ref_genome,ref_anno):

    try:
        os.mkdir(outdir)
    except FileExistsError as error:
        print("Folder {} already exists.".format(outdir))
    else:

        ##for trimmomatic results
        base_Ffile = os.path.basename(forward_reads)
        base_Rfile = os.path.basename(reverse_reads)

        forward_paired = outdir+'_'+base_Ffile.replace('.','_paried.')
        forward_unpaired = outdir+'_'+base_Ffile.replace('.','_unparied.')
        reverse_paired = outdir+'_'+base_Rfile.replace('.','_paried.')
        reverse_unpaired = outdir+'_'+base_Rfile.replace('.','_unparied.')

        #trimmomatic
        run_trimmomatic(
            forward_reads,reverse_reads,
            forward_paired,forward_unpaired,
            reverse_paired,reverse_unpaired
        )
        
        run_fastqc(forward_paired,reverse_paired,outdir) #runs fastqc
        run_spades_default(forward_paired,reverse_paired,outdir) #runs spades default
        run_spades_careful(forward_paired,reverse_paired,outdir) #runs spades careful
        megahit_default(forward_paired,reverse_paired,outdir) #runs megahit default
        megahit_mincount3(forward_paired,reverse_paired,outdir) #runs megahit mincount
        
        command = [ "quast", "-o", outdir+"_quast_result","-R",ref_genome,
                   "-g",ref_anno, "-l","Spades_default,Spades_corrected,Megahit_default,Megahit_mincount",
                    os.path.join(outdir+"_spades_default_output","contigs.fasta"),
                    os.path.join(outdir+"_spades_careful_output","contigs.fasta"),
                    os.path.join(outdir+"_megahit_default_output","final.contigs.fa"),
                    os.path.join(outdir+"_megahit_mincount3_output","final.contigs.fa")
            
        ]
        
        subprocess.run(command) #runs quast
        
        
        #finding best assembly
        quast_path = os.path.join(outdir+"_quast_result","report.tsv")
        
        num = []
        with open(quast_path,'r') as f:
            read_tsv = csv.reader(f, delimiter="\t")
            for row in read_tsv:
                if 'Total aligned length' in row:
                    for char in row[1:]:
                        num.append(int(char))
        f.close()
        m = (max(num))
        i = num.index(m)
        
        if i == 0:
            ret_path = os.path.join(outdir+"_spades_default_output","contigs.fasta")
        elif i==1:
            ret_path = os.path.join(outdir+"_spades_careful_output","contigs.fasta")
        elif i==2:
            ret_path =os.path.join(outdir+"_megahit_default_output","final.contigs.fa")
        else:
            ret_path =os.path.join(outdir+"_megahit_mincount3_output","final.contigs.fa")
            
        #ret_path is the path to best assembly contigs file 
            
        run_prokka(ret_path,outdir) # runs prokka
        
        command = "mv "+outdir+'_*  '+outdir+'/' 
        subprocess.call(command, shell=True)   #moves every output folder to output directory
        


