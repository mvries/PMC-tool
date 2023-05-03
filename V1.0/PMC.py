#This is a python wrapper that controls the PMC workflow.
#Running this script with the required inputs automatically performs analysis.
#To use: Python3 PMC.py sample_file phix.fasta plant.fasta output_dir/ ncores
#This is a python wrapper that controls the entire workflow.
#Running this script with the required inputs automatically performs analysis.
#To use: Python3 PMC.py sample_file configfile phix.fasta plant.fasta output_dir/ ncores

#Import statements:
from sys import argv
import subprocess
from scripts import make_config
from scripts import sample_pooler

#State directories and files needed:
sample_file = argv[1]
phix = argv[2]
plant = argv[3]
output_dir = argv[4]
cores = argv[5]

print("Pipeline initiated with the following arguments:")
print("sample file: " + sample_file)
print("phix: " + phix)
print("plant: " + plant)
print("outdir: " + output_dir)
print("cores: " + cores)

#Establish path to the snakemake config file:
configfile = "config/config.yaml"

#Reset config file for next run:
cmd = "rm config/config.yaml"
subprocess.check_call(cmd, shell=True)
cmd = "cp config/config_backup ../config/config.yaml"
subprocess.check_call(cmd, shell=True)

#Use the make_config script to prepare the config file for execution:
make_config.make_config(sample_file, configfile, phix, plant, output_dir)

#Prepare conda:
cmd = "conda config --set channel_priority strict"
subprocess.check_call(cmd, shell=True)

#Run the snakemake workflow that downloads and pre processes the input samples:
cmd = "snakemake --rerun-incomplete --cores " + str(cores) + " --use-conda " + "--snakefile Snakefiles/pre_processer.smk"
subprocess.check_call(cmd, shell=True)

#Clear the cache directory for SRA download:
print("Clearing cache Dir...")

cmd = "rm -r /dev/shm/sra"
subprocess.check_call(cmd, shell=True)

print("Starting MEGAHIT assemblies:")
#Perform the pooling of the samples:
#By default pools of 1 billion read pairs are made.
sample_pooler.make_pools(output_dir, configfile, 1000000000)
#Make directory for megahit output:
cmd = "mkdir " + str(output_dir) + "MEGAHIT"
subprocess.check_call(cmd, shell=True)
cmd = "snakemake --cores " + str(cores) + " --use-conda " + "--snakefile Snakefiles/co_assembly.smk"
subprocess.check_call(cmd, shell=True)
