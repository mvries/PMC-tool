#This script is used to extract the bins from the output of magscot.
#To run python3 extract_bins.py contigfile binfile outdir


#Import statements
from sys import argv 

#This function makes a dictionary of the contig file.
#Keys are contig ids, values are the respective sequences.
def make_contig_library(contigfile):
	print("Making contig library:")
	contig_lib = {}
	with open(contigfile, 'r') as reader:
		key = ""
		sequence = ""
		for line in reader:
			if line[0] == ">":
				if sequence == "":
					key = line.strip(">")
					key_1 = key.strip("\n")
				else:
					sequence_1 = sequence.strip("\n")
					contig_lib[key_1] = sequence_1
					key = line.strip(">")
					key_1 = key.strip("\n")
					sequence = ""
			else:
				sequence += line
		sequence_1 = sequence.strip("\n")
		contig_lib[key_1] = sequence
		return contig_lib

#The following functions extracts the bins and generates the output fastas.		
def extract_bins(contig_lib, binfile, outdir):
	print("Extracting bins:")
	with open(binfile, 'r') as reader:
		for line in reader:
			split = line.split("\t")
			binid = split[0]
			contig = split[1]
			if binid == "binnew":
				continue
			else:
				output_bin = str(outdir) + str(binid) + ".fa"
				with open(output_bin, "a") as appender:
					sequence = contig_lib[contig.strip("\n")]
					appender.write(">" + str(contig))
					appender.write(sequence + "\n")
					
#Main statement:
def main():
	contigfile = argv[1]
	binfile = argv[2]
	outdir = argv[3]
	contig_lib = make_contig_library(contigfile)
	extract_bins(contig_lib, binfile, outdir)

if __name__ == "__main__":
	main()
