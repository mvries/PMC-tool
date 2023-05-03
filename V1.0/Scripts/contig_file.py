#This script is used to produce a contig to bin file for the bins made.
#This file is needed as input for the magscot binning and refinement tool.
#To run: python3 contig_file.py bin_dir/ binning_tool contig_to_bin_file.txt

#Import statements:
from sys import argv
import glob

def make_contig_to_bin(bindir, bintool, outfile):
	#First a list is made of the bins in the bindir:
	binlist = glob.glob(bindir + "*.fa*")
	#The output file is created:
	with open(outfile, "w+") as writer:
		#Declare a variable that tracks the number of the bin:
		bin_num = 0
		#The script loops over the bins in the bindir:
		for i in binlist:
			bin_num += 1
			#Prepare the bin entry for writing:
			entry_1 = str(bintool) + "_bin_" + str(bin_num)
			#Now the contigs have to be extracted from the bins:
			with open(i, "r") as bin_reader:
				#Make a list to store the contigs in this bin:
				contig_list = []
				#Loop over each line in the bin fasta:
				#Contigs ids are stored in the list created above:
				for line in bin_reader:
					if line[0] == ">":
						contig_id = line.strip(">")
						contig_id_2 = contig_id.strip("\n")
						contig_list.append(contig_id_2)
				#Now the script loops over the list of contigs ids:
				#The contig id, bin and binning tool are written to out:
				for contig in contig_list:
					writer.write(entry_1 + "\t" + str(contig) + '\t' +str(bintool) + '\n')
					
#Main statement:
def main():
	bindir = argv[1]
	bintool = argv[2]
	outfile = argv[3]
	make_contig_to_bin(bindir, bintool, outfile)

#Main switch:
if __name__ == "__main__":
	main()
