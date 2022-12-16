#This is a snakemake rule that is used for the qc of short read data.
#This rule takes as input paired end read files and outputs a quality report.
#The quality control tool used in this case is fastp.

rule fastp:
