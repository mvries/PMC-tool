#This script is used to perform an initial quality screening of the megahit contigs.
#Contigs that are smaller then 500bp are discarded from the analysis by default.
#To use pyhton3 contig_filter_megahit.py contig_file output_file contig_length

#Import statements:
from sys import argv

def filter_contigs(contig_file, output_file, contig_length):
    with open(contig_file, "r") as contigs:
        with open(output_file, "a") as output:
            writer = 0
            for line in contigs:
                if line[0] == ">":
                    length = line.split(" ")[-1]
                    length_int = length.split("=")[-1]
                    if length_int >= contig_length:
                        writer = 1
                    else:
                        writer = 0
                else:
                    writer = writer
                if writer == 1:
                    output.write(line)
                else:
                    writer = writer

#Main statement:
def main():
    contig_file = argv[1]
    output_file = argv[2]
    contig_length = argv[3]
    filter_contigs(contig_file, output_file, contig_length)

if __name__ == "__main__":
    main()
