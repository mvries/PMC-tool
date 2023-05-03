#This script is used to extract the magscot bins to a share fasta file so they can be used for blast search.


#Import statement:
from sys import argv



#This funtion is used to extract the individual bins to a joined fasta.
def extract_joined_fasta(bin_dir, good_bins_file, output_dir):
	with open(good_bins_file, "r") as reader:
		for i in reader:
			name = i.strip("\n")
			with open(bin_dir + name + ".fa", "r") as reader2:
				with open(output_dir + "joined_bins.fasta", "a") as writer:
					writer.write(">" + name + "\n")
					for line in reader2:
						if line[0] != ">":
							writer.write(line)
						else:
							continue
	
#Main statement:
def main():
	bin_dir = argv[1]
	good_bins_file = argv[2]
	output_dir = argv[3]
	extract_joined_fasta(bin_dir, good_bins_file, output_dir)

if __name__ == "__main__":
	main()
