#This is a python wrapper that controlls the entire workflow.
#Running this script with the required inputs automatically performs analysis.
#To use: Python3 PMC.py sample_file configfile phix.fasta plant.fasta output_dir/

#Import statements:
from sys import argv
import subprocess
from scripts import make_config

#State directories and files needed:
sample_file = argv[1]
configfile = argv[2]
phix = argv[3]
plant = argv[4]
output_dir = argv[5]



#Reset config file for next run:
cmd = "rm ../config/config.yaml"
subprocess.check_call(cmd, shell=True)
cmd = "cp ../config/config_backup ../config/config.yaml"
subprocess.check_call(cmd, shell=True)

#Use the make_config script to prepare the config file for execution:
make_config.make_config(sample_file, configfile, phix, plant, output_dir)

#Prepare conda:
cmd = "conda config --set channel_priority strict"
subprocess.check_call(cmd, shell=True)

#Run the snekmake workflow:
cmd = "snakemake --cores 4 --use-conda"
subprocess.check_call(cmd, shell=True)
