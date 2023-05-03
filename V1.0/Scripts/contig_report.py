
#This script handles the data processing of the analysis pipeline.
#This script mainly makes output plots.

#To run: Python3 contig_report.py output_dir gene_abbundance_file checkmfile_raw checkm_file_norm barrnap_file magscot_contig2bin_file CAToutput phyloflash fasta BAT_output
#This requires the following input files:
#1: Gene_abbundance_file: (Predicted ORFs) in following format (Tabsep):
#Name	Assembly_orfs	Bin_orfs

#2: CheckM qa output file for the bins to be analyzed.

#3: Output file of barrnap SSU gene predictions

#4: Magscot contig to bin file

#CAT output with taxomy names (please use --only_official naming in CAT)

#Import required packages:
from sys import argv
import matplotlib.pyplot as plt

############################ INPUT PARSING #############################

#The following function is used to parse the gene abbundance file:
def parse_gene_abbundance(abbundance_file):
	#This file has a specified format.
	#Parsing in this case implies extracting the data to lists.
	outlist = []
	assemblie_names = []
	binned_orfs = []
	unbinned_orfs = []
	with open(abbundance_file, "r") as reader:
		for line in reader:
			entrie = line.split("\t")
			if len(entrie) == 1:
				raise Exception("The provided gene abbundance file is not tab seperated!")
			if len(entrie) > 3:
				raise Exception("The provided gene abbundance file has more than 2 columns!")
			assemblie_names.append(entrie[0])
			binned = entrie[2].strip("\n")
			unbinned = int(entrie[1]) - int(entrie[2])
			binned_orfs.append(int(binned))
			unbinned_orfs.append(int(unbinned))
	outlist.append(assemblie_names)
	outlist.append(binned_orfs)
	outlist.append(unbinned_orfs)
	return outlist


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

#The following function is used to produce a txt file of the high quality bin names:
def make_good_bin_file(checkm_parsed, output_dir):
	with open(output_dir + "good_bins.txt", "w") as writer:
		for i in checkm_parsed[1]:
			writer.write(i.split(" ")[0] + "\n")
	
#The following function is used to produce a txt file of the SSU gene names:		
def make_SSU_name_file(phyloflash_fasta, output_dir):
	with open(phyloflash_fasta, "r") as reader:
		with open(output_dir + "SSU_names.txt", "w") as writer:
			for line in reader:
				if line[0] == ">":
					writer.write(line.strip(">"))
				else:
					continue
			
	
#The following function is used to parse the barrnap SSU gene prediction:
#Only the 16S gene predictions are extracted since these are linked:
def parse_barrnap(barrnap_file):
	contigs = []
	with open(barrnap_file, "r") as reader:
		for line in reader:
			#We only look at the identifier lines:
			if line[0] == ">":
				line_strip = line.strip(">")
				line_split_1 = line_strip.split("_")
				line_split_2 = line_strip.split(":")
				if line_split_1[0] == "16S":
					contigs.append(line_split_2[-2])
				else:
					continue
			else:
				continue
	return contigs

#The following function is used to parse the contig_to_bin file from magscot:
def parse_contig_to_bin(contig_to_bin_file):
	#We loop over the contig to bin file and assign values to lib..
	#Lib has configuration contig = key | bin = value
	contig_to_bin_lib = {}
	with open(contig_to_bin_file, "r") as reader:
		for line in reader:
			line_split_1 = line.split("\t")
			bin_name = str(line_split_1[0]).split("/")[-1]
			contig = str(line_split_1[-1]).strip("\n")
			contig_to_bin_lib[contig] = bin_name
	return contig_to_bin_lib
	
#The following function is used to parse the CAT output of conig taxonomy.	
def parse_cat_output(CAT_output):
	#We open the input files and read trough the file.
	with open(CAT_output) as reader:
		#Declare variables were we can store the taxonmy classification:
		Kingdom = 0
		Phyla = 0
		Class = 0
		Order = 0
		Family = 0
		Genus = 0
		Species = 0
		Total_contigs = 0
		headerline = 0
		output_list = []
		#We loop over the lines and add the assignments to the variables above:
		for line in reader:
			if line[0] == "#":
				headerline += 1
				if headerline == 1:
					Total_contigs = int(line.split(" ")[-4].replace(",", ""))
				else:
					continue
			else:
				line_split = line.split("\t")
				if line_split[0] == "superkingdom":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Kingdom += int(line_split[-3])
				if line_split[0] == "phylum":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Phyla += int(line_split[-3])
				if line_split[0] == "class":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Class += int(line_split[-3])
				if line_split[0] == "order":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Order += int(line_split[-3])
				if line_split[0] == "family":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Family += int(line_split[-3])
				if line_split[0] == "genus":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Genus += int(line_split[-3])
				if line_split[0] == "species":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Species += int(line_split[-3])
						
		output_list.append(Total_contigs)
		output_list.append(Kingdom)
		output_list.append(Phyla)
		output_list.append(Class)
		output_list.append(Order)
		output_list.append(Family)
		output_list.append(Genus)
		output_list.append(Species)
		return output_list

def parse_bat_output(BAT_output):
	#We open the input files and read trough the file.
	with open(BAT_output) as reader:
		#Declare variables were we can store the taxonmy classification:
		Kingdom = 0
		Phyla = 0
		Class = 0
		Order = 0
		Family = 0
		Genus = 0
		Species = 0
		Total_bins = 0
		headerline = 0
		output_list = []
		#We loop over the lines and add the assignments to the variables above:
		for line in reader:
			if line[0] == "#":
				headerline += 1
				if headerline == 1:
					Total_bins = int(line.split(" ")[-5])
				else:
					continue
			else:
				line_split = line.split("\t")
				if line_split[0] == "superkingdom":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Kingdom += int(line_split[-1])
				if line_split[0] == "phylum":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Phyla += int(line_split[-1])
				if line_split[0] == "class":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Class += int(line_split[-1])
				if line_split[0] == "order":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Order += int(line_split[-1])
				if line_split[0] == "family":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Family += int(line_split[-1])
				if line_split[0] == "genus":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Genus += int(line_split[-1])
				if line_split[0] == "species":
					if line_split[1] == "no support":
						continue
					if line_split[1] == "NA":
						continue
					else:
						Species += int(line_split[-1])
						
		output_list.append(Total_bins)
		output_list.append(Kingdom)
		output_list.append(Phyla)
		output_list.append(Class)
		output_list.append(Order)
		output_list.append(Family)
		output_list.append(Genus)
		output_list.append(Species)
		return output_list
	

######################## PLOT GENERATION ###############################

#The following function is used to create a bar chart of gene abbundance:
def make_bar_chart(abbundance_parsed, output_dir):
	assemblie_names = ["Normalised", "Original"]
	binned_orfs = abbundance_parsed[1]
	unbinned_orfs = abbundance_parsed[2]
	#Make the plot and save to png file:
	plt.bar(assemblie_names, binned_orfs, color="r")
	plt.bar(assemblie_names, unbinned_orfs, bottom=binned_orfs, color="b")
	plt.xlabel("Assembly")
	plt.ylabel("Ammount of predicted ORFs")
	plt.legend(["Binned", "Unbinned"])
	plt.title("Total ORFs / Binned ORFs")
	plt.savefig(output_dir + "gene_abbundance_hist.png")
	plt.close()
	
#The following function is used to make a point plot of bin scores:
def make_point_plot(checkm_parsed_norm, checkm_parsed_raw, output_dir):
	#First we loop over the parsed checkM file to extract relevant data.
	checkm_parsed_all = checkm_parsed_norm[0]
	Contamination = []
	Completeness = []
	for entrie in checkm_parsed_all:
		entrie_split = entrie.split(" ")
		Contamination.append(float(entrie_split[-3]))
		Completeness.append(float(entrie_split[-4]))
	plt.scatter(Completeness, Contamination, c="b")

	checkm_parsed_all = checkm_parsed_raw[0]
	Contamination = []
	Completeness = []
	for entrie in checkm_parsed_all:
		entrie_split = entrie.split(" ")
		Contamination.append(float(entrie_split[-3]))
		Completeness.append(float(entrie_split[-4]))
	plt.scatter(Completeness, Contamination, c="r")
	plt.gca().invert_yaxis()
	plt.ylabel("Contamination %")
	plt.xlabel("Completeness %")
	plt.title("Contamination and completeness scores of bins")
	plt.legend(["Normalised", "Original"])
	plt.savefig(output_dir + "bin_score_scatterplot.png")
	plt.close()

def make_taxonomy_barplot(checkm_parsed, Cat_output_parsed, Bat_output_parsed, output_dir):
	#Now we make the barplot for the bins
	#We have to give the new values to the variables above for the contigs:
	assignment = "Bins n=(" + str(Bat_output_parsed[0]) + ")"
	Total_assignments = sum(Bat_output_parsed[1:])
	K_percent = (Bat_output_parsed[1]/Total_assignments)*100
	P_percent = (Bat_output_parsed[2]/Total_assignments)*100
	C_percent = (Bat_output_parsed[3]/Total_assignments)*100
	O_percent = (Bat_output_parsed[4]/Total_assignments)*100
	F_percent = (Bat_output_parsed[5]/Total_assignments)*100
	G_percent = (Bat_output_parsed[6]/Total_assignments)*100
	S_percent = (Bat_output_parsed[7]/Total_assignments)*100
	#Calculate total hight of previous bars to stack the barplot:
	P_height = K_percent
	C_height = K_percent + P_percent
	O_height = K_percent + P_percent + C_percent
	F_height = K_percent + P_percent + C_percent + O_percent
	G_height = K_percent + P_percent + C_percent + O_percent + F_percent
	S_height = K_percent + P_percent + C_percent + O_percent + F_percent + G_percent
	#Next the barplot is made for the contigs:
	plt.bar(assignment, K_percent, color="k")
	plt.bar(assignment, P_percent, bottom=K_percent, color="b")
	plt.bar(assignment, C_percent, bottom=C_height, color="g")
	plt.bar(assignment, O_percent, bottom=O_height, color="r")
	plt.bar(assignment, F_percent, bottom=F_height, color="c")
	plt.bar(assignment, G_percent, bottom=G_height, color="m")
	plt.bar(assignment, S_percent, bottom=S_height, color="y")

	#Now we make the barplot for the contigs
	#We have to give the new values to the variables above for the contigs:
	assignment = "Contigs n=(" + str(Cat_output_parsed[0]) + ")"
	Total_assignments = sum(Cat_output_parsed[1:])
	K_percent = (Cat_output_parsed[1]/Total_assignments)*100
	P_percent = (Cat_output_parsed[2]/Total_assignments)*100
	C_percent = (Cat_output_parsed[3]/Total_assignments)*100
	O_percent = (Cat_output_parsed[4]/Total_assignments)*100
	F_percent = (Cat_output_parsed[5]/Total_assignments)*100
	G_percent = (Cat_output_parsed[6]/Total_assignments)*100
	S_percent = (Cat_output_parsed[7]/Total_assignments)*100
	#Calculate total hight of previous bars to stack the barplot:
	P_height = K_percent
	C_height = K_percent + P_percent
	O_height = K_percent + P_percent + C_percent
	F_height = K_percent + P_percent + C_percent + O_percent
	G_height = K_percent + P_percent + C_percent + O_percent + F_percent
	S_height = K_percent + P_percent + C_percent + O_percent + F_percent + G_percent
	#Next the barplot is made for the contigs:
	plt.bar(assignment, K_percent, color="k")
	plt.bar(assignment, P_percent, bottom=K_percent, color="b")
	plt.bar(assignment, C_percent, bottom=C_height, color="g")
	plt.bar(assignment, O_percent, bottom=O_height, color="r")
	plt.bar(assignment, F_percent, bottom=F_height, color="c")
	plt.bar(assignment, G_percent, bottom=G_height, color="m")
	plt.bar(assignment, S_percent, bottom=S_height, color="y")
	
	plt.title("Taxonomic classifications of bins and contigs")
	plt.xlabel("Classification level")
	plt.ylabel("Percentage of classifications")
	plt.legend(["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])
	plt.savefig(output_dir + "taxonomy_barplot.png")
	plt.close()



"""
#The following function is used to make stacked bar charts for taxonomic assignments.
def make_taxonomy_barplot(checkm_parsed, Cat_output_parsed, output_dir):
	checkm_parsed_all = checkm_parsed[1]
	#For the bin classification plot we loop over parsed checkm file and extract data:
	assignment = "Bins n=(" + str(len(checkm_parsed_all)) + ")"
	Kingdom = 0
	Phylum = 0
	Class = 0
	Order = 0
	Family = 0
	Genus = 0
	Species = 0
	#Now we loop over the taxonomic assignment and add to the variables above.
	for entrie in checkm_parsed_all:
		entrie_split = entrie.split(" ")
		entrie_assignment = entrie_split[1][0]
		if entrie_assignment == "k":
			Kingdom += 1
		if entrie_assignment == "p":
			Phylum += 1
		if entrie_assignment == "c":
			Class += 1
		if entrie_assignment == "o":
			Order += 1
		if entrie_assignment == "f":
			Family += 1
		if entrie_assignment == "g":
			Genus += 1
		if entrie_assignment == "s":
			Species += 1
	#Next we calculate the percentage of each classification:
	Total_assignments = Kingdom + Phylum + Class + Order + Family + Genus + Species
	K_percent = (Kingdom/Total_assignments)*100
	P_percent = (Phylum/Total_assignments)*100
	C_percent = (Class/Total_assignments)*100
	O_percent = (Order/Total_assignments)*100
	F_percent = (Family/Total_assignments)*100
	G_percent = (Genus/Total_assignments)*100
	S_percent = (Species/Total_assignments)*100
	#Calculate total hight of previous bars to stack the barplot:
	P_height = K_percent
	C_height = K_percent + P_percent
	O_height = K_percent + P_percent + C_percent
	F_height = K_percent + P_percent + C_percent + O_percent
	G_height = K_percent + P_percent + C_percent + O_percent + F_percent
	S_height = K_percent + P_percent + C_percent + O_percent + F_percent + G_percent
	#Next the barplot is made for the bins:
	plt.bar(assignment, K_percent, color="k")
	plt.bar(assignment, P_percent, bottom=K_percent, color="b")
	plt.bar(assignment, C_percent, bottom=C_height, color="g")
	plt.bar(assignment, O_percent, bottom=O_height, color="r")
	plt.bar(assignment, F_percent, bottom=F_height, color="c")
	plt.bar(assignment, G_percent, bottom=G_height, color="m")
	plt.bar(assignment, S_percent, bottom=S_height, color="y")
	
	#Now we make the barplot for the contigs
	#We have to give the new values to the variables above for the contigs:
	assignment = "Contigs n=(" + str(Cat_output_parsed[0]) + ")"
	Total_assignments = sum(Cat_output_parsed[1:])
	K_percent = (Cat_output_parsed[1]/Total_assignments)*100
	P_percent = (Cat_output_parsed[2]/Total_assignments)*100
	C_percent = (Cat_output_parsed[3]/Total_assignments)*100
	O_percent = (Cat_output_parsed[4]/Total_assignments)*100
	F_percent = (Cat_output_parsed[5]/Total_assignments)*100
	G_percent = (Cat_output_parsed[6]/Total_assignments)*100
	S_percent = (Cat_output_parsed[7]/Total_assignments)*100
	#Calculate total hight of previous bars to stack the barplot:
	P_height = K_percent
	C_height = K_percent + P_percent
	O_height = K_percent + P_percent + C_percent
	F_height = K_percent + P_percent + C_percent + O_percent
	G_height = K_percent + P_percent + C_percent + O_percent + F_percent
	S_height = K_percent + P_percent + C_percent + O_percent + F_percent + G_percent
	#Next the barplot is made for the contigs:
	plt.bar(assignment, K_percent, color="k")
	plt.bar(assignment, P_percent, bottom=K_percent, color="b")
	plt.bar(assignment, C_percent, bottom=C_height, color="g")
	plt.bar(assignment, O_percent, bottom=O_height, color="r")
	plt.bar(assignment, F_percent, bottom=F_height, color="c")
	plt.bar(assignment, G_percent, bottom=G_height, color="m")
	plt.bar(assignment, S_percent, bottom=S_height, color="y")
	
	plt.legend(["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])
	plt.savefig(output_dir + "taxonomy_barplot.png")
	plt.close()
"""

#The following function is used to make a plot of ssu genes to bins:
def make_ssu_distribution_plot(barrnap_parsed, contig_to_bin_file_parsed, checkm_parsed, output_dir):
	checkm_parsed_good = checkm_parsed[1]
	#Determine total number of bins so percentage can be calculated:
	total_bins = len(checkm_parsed_good)
	#Make list to save the bins with SSU links:
	Linked_bins = []
	#We extract the bins that have 16S genes in them (according to barrnap)
	for contig in barrnap_parsed:
		bin_name = contig_to_bin_file_parsed[contig]
		Linked_bins.append(bin_name)
	#Calculate ammount of linked bins:
	num_linked_bins = len(Linked_bins)
	#Now we calculate the percentage of linked bins:
	linked_bins_percent = (num_linked_bins/total_bins)*100
	#Next we see how many bins are linked uniquely:
	unique_list = []
	for bin_name in Linked_bins:
		if bin_name not in unique_list:
			unique_list.append(bin_name)
	#Calculate the number of unique bins with links:
	unique_bins = len(unique_list)
	#Calculate the percentage of uniquely linked bins:
	unique_bins_percent = (unique_bins/total_bins)*100
	Total_minus_unique = linked_bins_percent - unique_bins_percent
	#Calculate percentage of bins with no SSU link:
	unlinked = total_bins - num_linked_bins
	#Calculate percentage of bins with no SSU genes:
	Unlinked_percentage = (unlinked/total_bins)*100
	#Now we make the plot:
	labels = ["No SSU","Unique SSU","Multiple SSU"]
	data = []
	data.append(Unlinked_percentage)
	data.append(unique_bins_percent)
	data.append(Total_minus_unique)
	plt.pie(data, labels = labels, autopct='%1.0f%%')
	plt.title("Percentage of bins with predicted SSU genes n=(" + str(total_bins) + ")")
	plt.savefig(output_dir + "SSU_to_bin_piechart.png")
	plt.close()


############################## MAIN STATEMENT ##########################


#Main statement:
def main():
	#Declare the input files and output dir:
	output_dir = argv[1]
	abbundance_file = argv[2]
	checkm_file_raw = argv[3]
	checkm_file_norm = argv[4]
	barrnap_file = argv[5]
	contig_to_bin_file = argv[6]
	CAT_output = argv[7]
	phyloflash_fasta = argv[8]
	BAT_output = argv[9]
	#Run the functions of the script:
	
	abbundance_parsed = parse_gene_abbundance(abbundance_file)
	checkm_raw_parsed = parse_checkm(checkm_file_raw)
	checkm_norm_parsed = parse_checkm(checkm_file_norm)
	barrnap_parsed = parse_barrnap(barrnap_file)
	contig_to_bin_file_parsed = parse_contig_to_bin(contig_to_bin_file)
	Cat_output_parsed = parse_cat_output(CAT_output)
	Bat_output_parsed = parse_bat_output(BAT_output)
	
	make_SSU_name_file(phyloflash_fasta, output_dir)
	"""
	make_good_bin_file(checkm_parsed, output_dir)
	"""
	make_bar_chart(abbundance_parsed, output_dir)
	make_point_plot(checkm_raw_parsed, checkm_norm_parsed, output_dir)
	"""
	make_taxonomy_barplot(checkm_parsed, Cat_output_parsed, Bat_output_parsed, output_dir)
	"""
	make_ssu_distribution_plot(barrnap_parsed, contig_to_bin_file_parsed, checkm_norm_parsed, output_dir)
if __name__ == "__main__":
	main()
