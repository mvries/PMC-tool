#This script is used to pool the input samples to a desired maximum read count per pool.
#Read depth = total number of reads desired per pool.

#Note: This scripts contains quite a few hardcoded directories; do not be afraid!
#These directories are created as a core part of the pipeline.

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
    sample_list = glob.glob(input_dir + 'bbmap/*')
    #Loop over the samples (extract id and read count):
    print("Counting reads:")
    for sample in sample_list:
        #Extract ID:
        sample_id = sample.split('/')[-1]
        #Count forward reads:
        cmd = "zcat " + sample + "/" + "*normalized_1* | wc -l"
        fw_count_full = int(subprocess.check_output(cmd, shell=True))
        fw_count_reads = int(fw_count_full) / 4
        #Count reverse reads:
        cmd = "zcat " + sample + "/" + "*normalized_2* | wc -l"
        rv_count_full = int(subprocess.check_output(cmd, shell=True))
        rv_count_reads = int(rv_count_full) / 4
        #Calculate total reads for the sample:
        Total_reads_sample = int(fw_count_reads + rv_count_reads)
        #Add sample and read count to the dict:
        size_dict[sample_id] = Total_reads_sample
        print(sample_id + " has " + str(Total_reads_sample) + " reads.")

    #Loop over the dictionary and split the samples into pools:
    #Counter to keep track of pool size:
    pool_size = 0
    #List to add samples to pool
    Pool = []
    #List to pool ... the pools
    pool_list = []
    for sample, size in size_dict.items():
        print("pooling " + str(sample) + " of size " + str(size))
        Pool.append(sample)
        pool_size += size
        print("total pool size " + str(pool_size))
        if pool_size >= read_depth:
            pool_list.append(Pool)
            Pool = []
            pool_size = 0
    if len(Pool) > 0:
        pool_list.append(Pool)


    #The next step is looping over the pool list and making the directories:
    #Declare a pool number variable to assist in writing the directories needed:
    print("Pooling done! " + str(len(pool_list)) + " pool(s) were made" )
    print("Writing pools to directories, this may also take some time...")
    pool_number = 0
    cmd = "mkdir " + input_dir + "Pools/"
    subprocess.check_call(cmd, shell=True)
    for pool in pool_list:
        pool_number += 1
        cmd = "mkdir " + input_dir + "Pools/pool" + str(pool_number)
        subprocess.check_call(cmd, shell=True)
        for sample in pool:
            cmd = "mkdir " + input_dir + "Pools/pool" + str(pool_number) + "/" + sample
            subprocess.check_call(cmd, shell=True)
            cmd = "mv " + input_dir + "bbmap/" + sample + '/*norm* ' + input_dir + 'Pools/pool' + str(pool_number) + '/' + sample + "/"
            subprocess.check_call(cmd, shell=True)
        print("Pool " + str(pool_number) + " written.")

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
