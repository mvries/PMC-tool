#This is a python wrapper that controls the PMC workflow.
#Running this script with the required inputs automatically performs analysis.
#To use: Python3 PMC.py sample_file phix.fasta plant.fasta output_dir/ ncores

#Import statements:
from sys import argv
import subprocess
from Scripts import make_config
from Scripts import sample_pooler
from Scripts import contig_file

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
cmd = "cp config/config_backup config/config.yaml"
subprocess.check_call(cmd, shell=True)

#Use the make_config script to prepare the config file for execution:
make_config.make_config(sample_file, configfile, phix, plant, output_dir)

#Prepare conda:
cmd = "conda config --set channel_priority true"
subprocess.check_call(cmd, shell=True)

#Make environment if needed:
try:
    cmd = "mamba env create --file Environments/Snakemake.yaml"
    subprocess.check_call(cmd, shell=True)
except:
    pass

#Run the snakemake workflow that downloads and pre processes the input samples:
cmd = "conda run --no-capture-output -n snakemake snakemake --rerun-incomplete --cores " + str(cores) + " --use-conda " + "--snakefile Snakefiles/pre_processer.smk"
subprocess.check_call(cmd, shell=True)

#Clear the cache directory for SRA download:
print("Clearing cache Dir...")

cmd = "rm -r " + output_dir + "/sra"
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

"""
#After the individual binning tools are done we have to run some scripts so we can do the binning with magscot:
#Note that MAGSCOT needs to be in your path (This tool is not conda compatible), all the dependencies are in the magscot conda environment
#First we need to make a contig to bin file
contig_file.make_contig_to_bin(output_dir/metabat/, metabat,output_dir/metabat_to_bin.tsv)
contig_file.make_contig_to_bin(output_dir/concoct/, concoct,output_dir/concoct_to_bin.tsv)
contig_file.make_contig_to_bin(output_dir/maxbin/, maxbin,output_dir/maxbin_to_bin.tsv)

cmd = "cat " + str(output_dir) + "*_to_bin.tsv >> " + str(output_dir) + "final_bin.tsv"
subprocess.check_call(cmd, shell=True)

cmd = "conda env create -f Environments/magscot.yaml"
subprocess.check_call(cmd, shell=True)

cmd = "MAGScoT.R -i " + str(output_dir) + "final_bin.tsv " + "-o " + str(output_dir)/magscot/
"""




