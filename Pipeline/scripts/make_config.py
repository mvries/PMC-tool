#This is a script that is used to write a config file.
#This config file is used during the execution of the snakemake workflow.
#The config file contains samples, input and output paths in .yaml format
#To run this script: python3 make_config.py sample_dir/ outdir/
#Sample dir structure: sample_dir/{Sample}/{Sample}.fq

def make_config(sample_dir, output_dir, configfile):
    #This function is used to process the input samples:
    #Required input: Path to Dir where sample reads are stored
    #Output: List of samples and paths in yaml format

    #Import statements:
    import glob

    #Make list to store paths:
    path_list = []
    #Make list to store sample names:
    sample_list = []
    #Glob is used to generate a list of paths to each sample in the input dir.
    #The 'root' path is saved and written to the config file later.
    sample_dir_full = sample_dir + "*"
    path_list = glob.glob(sample_dir_full)

    #Loop over path list to get the sample names and append to samples list:
    for i in path_list:
        sample = i.split("/")[-1]
        sample_list.append(sample)

    #Open config file and write in the right format (.yaml):
    with open(configfile, 'a') as config_file:
        config_file.write('SAMPLES:\n')

        #Loop over both the sample and path list to write them together:
        for i in range(0, len(sample_list)):
            #Make empty string with space to store entry for writing to file:
            entry = " "
            #add sample name to entry
            entry += (sample_list[i] + ': ')
            #add path to sample reads to entry
            entry += (path_list[i])
            #Write entries to config file:
            config_file.write(entry + "\n")

        #Add empty line to keep the file readable:
        config_file.write("\n")

        #Add the input and output path to the config file so snakemake can use:
        config_file.write("INPUT_PATH:\n")
        config_file.write(f' {sample_dir}\n')
        config_file.write("OUTPUT_PATH:\n")
        config_file.write(f' {output_dir}\n')

#Main function:
def main():
    #Import statement
    from sys import argv
    #Establish path to sample dir:
    sample_dir = argv[1]
    #Establish path to desired output dir:
    output_dir = argv[2]
    #Establish path to config file:
    config_file = argv[3]
    #Run process sample function to write the config file
    make_config(sample_dir, output_dir, configfile)

if __name__ == "__main__":
    main()
