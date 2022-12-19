#This is a snakemake rule that is used for the qc of short read data.
#This rule takes as input paired end read files and outputs: Trimmed reads and quality reports.
#The quality control tool used in this case is fastp.
#To use: snakemake --cores ncores --snakefile read_qc.smk sample_dir configfile
#Samples should be written to configfile first (use make_config.py).
#By default reports are stored in the same dir as the input reads.

#Import statements:
from sys import argv

#Declare path to configfile:
configfile: "../../config/config.yaml"

#Main rule:
rule all:
    input:
        expand("{sample_dir}{sample}/{sample}_{Rnum}_trimmed.fq.gz",
        sample_dir=config["INPUT_PATH"],
        sample = config["SAMPLES"],
        Rnum=['1', '2'])

#Rule to run fastp quality control:
rule fastp:
    input:
        r1="{sample_dir}{sample}/{sample}_1.fastq.gz",
        r2="{sample_dir}{sample}/{sample}_2.fastq.gz"
    output:
        r1out="{sample_dir}{sample}/{sample}_1_trimmed.fq.gz",
        r2out="{sample_dir}{sample}/{sample}_2_trimmed.fq.gz",
        html="{sample_dir}qc/{sample}/{sample}.html",
        json="{sample_dir}qc/{sample}/{sample}.json"
    conda:
        "../envs/fastp_env.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1out} -O {output.r2out} -h {output.html} -j {output.json}"
