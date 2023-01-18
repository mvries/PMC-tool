#This snakefile performs the linkage of SSU genes to bins produced by:
#individual metaspades assemblies for each sample

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
        expand(f'{output_dir}' + "Quast/" + "{sample}" + "/report.html",
        sample=config["SAMPLES"])

#The following rules performs metaspades assemblies on individual samples:
include: "../rules/spades_assembly/metaspades.smk"

#The following rule produces an assemblie quality report using metaQUAST:
include: "../rules/spades_assembly/Quast.smk"
