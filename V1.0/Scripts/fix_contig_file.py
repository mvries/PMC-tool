#Import statements
from sys import argv

def fix_contig_file(input_file, output_file):
	with open(input_file, 'r') as reader:
		with open(output_file, "w+") as writer:
			for line in reader:
				split = line.split("\t")
				writer.write(split[0] + "\t" + "k_" + split[1] + "\t" + split[2])
			
#Main statement:
def main():
	input_file = argv[1]
	output_file = argv[2]
	fix_contig_file(input_file, output_file)

if __name__ == "__main__":
	main()
