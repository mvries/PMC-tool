#To run python3 binplot.py raw_checkmfile norm_checkmfile spades_checkmfile


from sys import argv
import matplotlib.pyplot as plt

#The following function is used to parse the checkm file:
def parse_checkm(checkm_file):
        #The problem is that the output of checkm qa is space seperated uneven
        #This file is parsed so that a single space seperates the columns:
        with open(checkm_file, "r") as reader:
                fixed_lines = []
                good_bins = []
                headerline = 0
                for line in reader:
                        new_line = ""
                        prevcar = ""
                        #First the header lines are skipped:
                        if line[0] == "-":
                                headerline += 1
                        if headerline == 2:
                                if line[0] != "-":
                                        #Now we loop over the characters in the lines to fix them:
                                        for char in line:
                                                if char == " ":
                                                        if prevcar == "":
                                                                continue
                                                        if prevcar == " ":
                                                                continue
                                                        if prevcar != " ":
                                                                new_line += " "
                                                                prevcar = " "
                                                else:   
                                                        new_line += char
                                                        prevcar = char
                        #Some lines will be empty and are discared.
                        #The fixed lines are returned as a list of lists.
                        if new_line == "":
                                continue
                        else:
                                fixed_lines.append(new_line.strip("\n"))
                #Now we selecct bins with completeness - contamination > 70%            
                for i in fixed_lines:
                        split = i.split(" ")
                        contamination = float(split[-3])
                        completeness = float(split[-4])
                        if completeness - contamination > 0:
                                good_bins.append(i)
                #We store the lists of total bins and good bins
                #Some plots need all the bins some plots only need the good ones.
                output_list = []
                output_list.append(fixed_lines)
                output_list.append(good_bins)
                return output_list

#The following function is used to make a point plot of bin scores:
def make_point_plot(checkm_parsed_norm, checkm_parsed_raw, checkm_spades, output_dir):
        #First we loop over the parsed checkM file to extract relevant data.
	checkm_parsed_all = checkm_parsed_norm[0]
	Contamination = []
	Completeness = []
	for entrie in checkm_parsed_all:
		entrie_split = entrie.split(" ")
		Contamination.append(float(entrie_split[-3]))
		Completeness.append(float(entrie_split[-4]))
	plt.scatter(Completeness, Contamination, c="b", alpha=0.5)

	checkm_parsed_all = checkm_parsed_raw[0]
	Contamination = []
	Completeness = []
	for entrie in checkm_parsed_all:
		entrie_split = entrie.split(" ")
		Contamination.append(float(entrie_split[-3]))
		Completeness.append(float(entrie_split[-4]))
	plt.scatter(Completeness, Contamination, c="r", alpha=0.5)

	checkm_spades_all = checkm_spades[0]
	Contamination = []
	Completeness = []
	for entrie in checkm_spades_all:
		entrie_split = entrie.split(" ")
		Contamination.append(float(entrie_split[-3]))
		Completeness.append(float(entrie_split[-4]))
	plt.scatter(Completeness, Contamination, c="g", alpha=0.8)
	plt.gca().invert_yaxis()
	plt.axvline(x=70)
	plt.ylabel("Contamination %")
	plt.xlabel("Completeness %")
	plt.title("Contamination and completeness scores of MAG sets")
	plt.legend(["Normalised pooled assembly (n=186)", "Original pooled assembly (n=211)", "Individual assemblies (n=80)"])
	plt.savefig(output_dir + "bin_score_scatterplot.png")
	plt.close()






def main():
	raw_checkm_file = argv[1]
	norm_checkm_file = argv[2]
	spades_checkm_file = argv[3]
	raw_parsed = parse_checkm(raw_checkm_file)
	norm_parsed = parse_checkm(norm_checkm_file)
	spades_parsed = parse_checkm(spades_checkm_file)
	output_dir = argv[4]
	make_point_plot(norm_parsed, raw_parsed, spades_parsed, output_dir)



if __name__ == "__main__":
	main()
