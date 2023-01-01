#This is a script that is used to write a config file.
#This config file is used during the execution of the snakemake workflow.
#The config file contains samples and output paths in .yaml format
#To run this script: python3 make_config.py samplefile configfile phix.fa plant.fa outdir/

def make_config(sample_file, configfile, phix, plant, output_dir):
    #This function is used to process the input samples:

    #Loop over sample file to make list of samples:
    sample_list = []
    with open(sample_file, "r") as sample_file:
        for sample in sample_file:
            sample_list.append(sample)

    #Open config file and write in the right format (.yaml):
    with open(configfile, 'a') as config_file:
        config_file.write('\n')
        config_file.write('SAMPLES:\n')

        for sample in sample_list:
            entry = " "
            entry += sample
            config_file.write(entry + "\n")


        #Add empty line to keep the file readable:
        config_file.write("\n")

        #Add the sample references and output paths to the config file so snakemake can use:
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
    #Establish path to sample file:
    sample_file = argv[1]
    #Establish path to configfiler:
    configfile = argv[2]
    #Establish path to phix reference genome:
    phix = argv[3]
    #Establish path to plant ref:
    plant = argv[4]
    #Establish path to output_dir:
    output_dir = argv[5]
    #Run function to write the config file
    make_config(sample_file, configfile, phix, plant, output_dir)

if __name__ == "__main__":
    main()
