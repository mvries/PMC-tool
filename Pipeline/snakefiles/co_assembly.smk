#This snakefile performs the linkage of SSU genes to bins produced by:
#MEGAHIT assemblies on pooled reads of sample files.

#Import statements:
from snakemake.utils import min_version
min_version("6.0")

#set config file:
configfile: "../config/config.yaml"

#Determine input and output paths for snakemake:
output_dir = config["OUTPUT_PATH"]
phix_reference= config["PHIX_REF"]
plant_reference= config["PLANT_REF"]

###############################################################################
###Execute main workflow:###

###Main rule:
rule all:
    input:
        expand(f'{output_dir}' + "MEGAHIT/{pool}/",
        pool=config["POOLS"])

#The following rule is used to pool the reads:
include: "../rules/pool_reads.smk"

#The followiung rule is used to assemble the pooled reads with MEGAHIT:
include: "../rules/megahit.smk"
