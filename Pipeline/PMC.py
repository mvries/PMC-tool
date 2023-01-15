#This is a python wrapper that controls the PMC workflow.
#Running this script with the required inputs automatically performs analysis.
#To use: Python3 PMC.py sample_file phix.fasta plant.fasta output_dir/ ncores Assembly_type
#Possible assembly types: single (runs metaspades on each sample) co (runs megahit on pooled samples)

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
type = argv[6]

#Establish path to the snakemake config file:
configfile = "../config/config.yaml"

#Reset config file for next run (This is a safety measure in case of crash):
cmd = "rm " + configfile
subprocess.check_call(cmd, shell=True)
cmd = "cp ../config/config_backup " + configfile
subprocess.check_call(cmd, shell=True)

#Use the make_config script to prepare the config file for execution:
make_config.make_config(sample_file, configfile, phix, plant, output_dir)

#Prepare conda:
cmd = "conda config --set channel_priority strict"
subprocess.check_call(cmd, shell=True)

#Run the snakemake workflow that downloads and pre processes the input samples:
cmd = "snakemake --cores " + str(cores) + " --use-conda" + "--snakefile snakefiles/pre_processer.smk"
subprocess.check_call(cmd, shell=True)

#Run the appropriate snekmake assembly workflow depending on the assembly type:
if type == "single":
    cmd = "snakemake --cores " + str(cores) + " --use-conda" + "--snakefile snakefiles/spades_assembly.smk"
    subprocess.check_call(cmd, shell=True)
else:
    #Perform the pooling of the samples:
    sample_pooler.make_pools(output_dir, configfile, 1000000000)
    cmd = "snakemake --cores " + str(cores) + " --use-conda" + "--snakefile snakefiles/co_assembly.smk"
    subprocess.check_call(cmd, shell=True)
