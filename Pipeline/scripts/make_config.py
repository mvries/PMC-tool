#This is a script that is used to write a config file.
#This config file is used during the execution of the snakemake workflow.
#The config file contains samples, input and output paths in .yaml format
#To run this script: python3 make_config.py sample_dir/ outdir/ phix_file.fa plant_genome.fa
#Sample dir structure: sample_dir/{Sample}/{Sample}.fq

def make_config(sample_dir, output_dir, configfile, phix, plant):
    #This function is used to process the input samples:
    #Required input: Path to Dir where sample reads are stored
    #Output: List of samples and paths in yaml format

    #Import statements:
    import glob

    #Make list to store sample paths:
    sample_path_list = []
    #Make list to store sample names:
    sample_list = []
    #Make list to store reference genome paths:


    #Glob is used to generate a list of paths to each sample in the input dir.
    #The 'root' path is saved and written to the config file later.
    sample_dir_full = sample_dir + "*"
    sample_path_list = glob.glob(sample_dir_full)

    #Loop over path list to get the sample names and append to samples list:
    for i in sample_path_list:
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
            entry += (sample_path_list[i])
            #Write entries to config file:
            config_file.write(entry + "\n")

        #Add empty line to keep the file readable:
        config_file.write("\n")

        #Add the sample references and output paths to the config file so snakemake can use:
        config_file.write("INPUT_PATH:\n")
        config_file.write(f' {sample_dir}\n')
        config_file.write("OUTPUT_PATH:\n")
        config_file.write(f' {output_dir}\n')
        config_file.write("PHIX_REF:\n")
        config_file.write(f' {phix}\n')
        config_file.write("PLANT_REF:\n")
        config_file.write(f' {plant}\n')
        #Add empty line to keep the file readable:
        config_file.write("\n")



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
    #Establish path to Phix ref:
    phix = argv[4]
    #Establish path to plant ref:
    plant = argv[5]
    #Run process sample function to write the config file
    make_config(sample_dir, output_dir, configfile, phix, plant)

if __name__ == "__main__":
    main()
