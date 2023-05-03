#Import statments:
from sys import argv


def fix_contigs(input_file, output_file):
	contig_number = 0
	with open(input_file, 'r') as reader:
		with open(output_file, "w+") as writer:
			for line in reader:
				if line[0] == ">":
					contig_number += 1
					writer.write(">k_" + str(contig_number) + "\n")
				else:
					writer.write(line)

def main():
	input_file = argv[1]
	output_file = argv[2]
	fix_contigs(input_file, output_file)
					
if __name__ == "__main__":
	main()
