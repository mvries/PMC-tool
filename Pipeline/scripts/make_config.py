#This is a script that is used to write a config file.
#This config file is used during the execution of the snakemake workflow.
#The config file contains a list of input samples in .yaml format
#To run this script: python3 make_config.py sample_dir path_to_yaml

def process_samples(sample_dir, path_to_config):
    #This function is used to process the input samples:
    #Required input: Path to Dir where sample reads are stored
    #Output: List of samples and paths in yaml format

    #Import statements:
    import glob

    #Make list to store paths:
    path_list = []
    #Make list to store sample names:
    sample_list = []
    #Glob is used to generate a list of paths to each sample
    sample_dir = sample_dir + "*"
    path_list = glob.glob(sample_dir)
    print(path_list)

    #Loop over path list to get the sample names and append to samples list:
    for i in path_list:
        sample = i.split("/")[-1]
        sample_list.append(sample)

    #Open config file and write in the right format (.yaml):
    with open(path_to_config, 'a') as config_file:
        config_file.write('samples:\n')

        #Loop over both the sample and path list to write them together:
        for i in range(0, len(sample_list)):
            #Make empty string with tab to store entry for writing to file:
            entry = "\t"
            #add sample name to entry
            entry += (sample_list[i] + ': ')
            #add path to sample reads to entry
            entry += (path_list[i])
            #Write entries to config file:
            config_file.write(entry + "\n")

#Main function:
def main():
    #Import statement
    from sys import argv
    #Establish path to sample dir:
    sample_dir = argv[1]
    #Establish path to config_file
    path_to_config = argv[2]
    #Run process sample function to write the config file
    process_samples(sample_dir, path_to_config)

if __name__ == "__main__":
    main()
