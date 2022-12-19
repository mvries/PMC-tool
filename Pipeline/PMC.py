#This is a python wrapper that controlls the entire workflow.
#Running this script with the required inputs automatically performs analysis.
#To use: Python3 PMC.py sample_dir/ Output_dir/ configfile

#Import statements:
from sys import argv
import subprocess
from scripts import make_config

#State directories
sample_dir = argv[1]
output_dir = argv[2]
configfile = argv[3]

#Use the make_config script to prepare the config file for execution:
make_config.make_config(sample_dir, output_dir, configfile)

#Prepare conda:
cmd = "conda config --set channel_priority strict"
subprocess.check_call(cmd, shell=True)

#Run the snekmake workflow:
cmd = "snakemake --cores 4 --use-conda"
subprocess.check_call(cmd, shell=True)

#Collect output to desired directory:
cmd = f"mkdir {output_dir}"
subprocess.check_call(cmd, shell=True)
cmd = f"mv {sample_dir}qc/ {output_dir}qc"
subprocess.check_call(cmd, shell=True)

#Reset config file for next run:
cmd = "rm ../config/config.yaml"
subprocess.check_call(cmd, shell=True)
cmd = "cp ../config/config_backup ../config/config.yaml"
subprocess.check_call(cmd, shell=True)
