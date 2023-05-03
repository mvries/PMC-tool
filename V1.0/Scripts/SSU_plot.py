

import matplotlib.pyplot as plt


from sys import argv



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
                        if completeness - contamination > float(70):
                                good_bins.append(i)
                #We store the lists of total bins and good bins
                #Some plots need all the bins some plots only need the good ones.
                output_list = []
                output_list.append(fixed_lines)
                output_list.append(good_bins)
                return output_list




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



#Main statement:
def main():
        #Declare the input files and output dir:
        output_dir = argv[1]
        checkm_file = argv[2]
        barrnap_file = argv[3]
        contig_to_bin_file = argv[4]
        #Run the functions of the script:

        checkm_parsed = parse_checkm(checkm_file)
        barrnap_parsed = parse_barrnap(barrnap_file)
        contig_to_bin_file_parsed = parse_contig_to_bin(contig_to_bin_file)
        make_ssu_distribution_plot(barrnap_parsed, contig_to_bin_file_parsed, checkm_parsed, output_dir)

if __name__ == "__main__":
        main()

