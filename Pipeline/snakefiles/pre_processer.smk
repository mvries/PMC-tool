#This snakefile performs the downloading and pre processing of input samples.

#Import statements:
from snakemake.utils import min_version
min_version("6.0")

#set config file:
configfile: "../config/config.yaml"

#Determine input and output paths for snakemake:
output_dir = config["OUTPUT_PATH"]
phix_reference= config["PHIX_REF"]
plant_reference= config["PLANT_REF"]

################################################################################
###Execute the workflow:###

###Main rule:
rule all:
    input:
        expand(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_1.fq.gz",
        sample=config["SAMPLES"])

#The following rule downloads samples from the sequence read archive:
include: "../rules/pre_processing/sample_download.smk"

#The following rule is used to perform quality control of input reads with fastp:
include: "../rules/pre_processing/fastp.smk"

#The following rule uses bowtie2 to remove phix contamination:
include: "../rules/pre_processing/remove_phix.smk"

#The following rule uses bowtie2 to remove host plant contamination:
include: "../rules/pre_processing/remove_plant.smk"

#The following rules uses bbnorm function from bbmap to normalize fastq kmer depth:
include: "../rules/pre_processing/bbnorm.smk"
