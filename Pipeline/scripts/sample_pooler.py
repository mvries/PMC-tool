#This script is used to pool the input samples to a desired maximum read count per pool.
#To use the script: python3 sample_pooler.py input_dir/ configfile read read_depth
#Read depth = total number of reads desired per pool.

#Import statements:
from sys import argv
import subprocess
import glob

#This function is used to define pools of samples so that they are aprox desired read depth
def make_pools(input_dir, configfile, read_depth):
    print("Pooling samples: For large datasets this takes some time...")
    #The first step is counting the reads in each of the samples after qc:
    #Make dict to store reads counts for each sample:
    size_dict = {}
    #Make list of samples:
    sample_list = glob.glob(input_dir + 'bowtie2/plant/*')
    #Loop over the samples (extract id and read count):
    for sample in sample_list:
        #Extract ID:
        sample_id = sample.split('/')[-1]
        #Count forward reads:
        cmd = "zcat " + sample + "/" + "*_1* | wc -l"
        fw_count_full = subprocess.check_output(cmd, shell=True)
        fw_count_reads = fw_count_full / 4
        #Count reverse reads:
        cmd = "zcat " + sample + "/" + "*_2* | wc -l"
        rv_count_full = subprocess.check_output(cmd, shell=True)
        rv_count_reads = fw_count_full / 4
        #Calculate total reads for the sample:
        Total_reads_sample = fw_count_reads + rv_count_reads
        #Add sample and read count to the dict:
        size_dict[sample_id] = Total_reads_sample

    #Loop over the dictionary and split the samples into pools:
    #Counter to keep track of pool size:
    pool_size = 0
    #List to add samples to pool
    Pool = []
    #List to pool ... the pools
    pool_list = []
    for sample, size in size_dict.items():
        Pool.append(sample)
        pool_size += size
        if pool_size => read_depth:
            pool_list.append(pool)
            Pool = []
            pool_size = 0
        else:
            Pool = Pool

    #The next step is looping over the pool list and making the directories:
    #Declare a pool number variable to assist in writing the directories needed:
    print("Pooling done! " + str(len(pool_list)) + " pool(s) were made" )
    print("Writing pools to directories, this may also take some time...")
    pool_number = 0
    for pool in pool_list:
        pool_number += 1
        cmd = "mkdir " + input_dir + " Pools/pool" + str(pool_number)
        subprocess.check_call(cmd, shell=True)
        for sample in pool:
        cmd = "mv " + input_dir + "bowtie2/plant/" + sample + ' ' + input_dir + 'Pools/pool' + str(pool_number) + '/'
        subprocess.check_call(cmd, shell=True)

    #The final step is updating the configfile so that snakemake can use them:
    print("Writing done, updating configfile...")
    #Reset pool_number variable:
    pool_number = 0
    with open(configfile, 'a') as config_file:
        config_file.write("POOLS:\n")
        for pool in pool_list:
            pool_number += 1
            config_file.write(" pool" + str(pool_number) + ":\n")

    print("Pooling is done, starting the assembly stage of the pipeline!")

#Main function of the script:
def main():
    input_dir = argv[1]
    configfile = argv[2]
    read_depth = argv[3]
    make_pools(input_dir, configfile, read_depth)

#Main statement:
if __name__ == "__main__":
    main()
