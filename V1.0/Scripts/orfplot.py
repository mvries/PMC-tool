import matplotlib.pyplot as plt
from sys import argv

#The following function is used to parse the gene abbundance file:
def parse_gene_abbundance(abbundance_file):
        #This file has a specified format.
        #Parsing in this case implies extracting the data to lists.
        outlist = []
        assemblie_names = []
        binned_orfs = []
        unbinned_orfs = []
        percentages = []
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

        for i in range(0, len(assemblie_names)):
            binned = binned_orfs[i]
            unbinned = unbinned_orfs[i]
            percentage = binned / (unbinned + binned) * 100
            percentages.append(percentage)
        outlist.append(percentages)
        return outlist

#The following function is used to create a bar chart of gene abbundance:
def make_bar_chart(abbundance_parsed, output_dir):
        assemblie_names = ["Normalised assembly", "Original assembly", "Individual assemblies"]
        binned_orfs = abbundance_parsed[1]
        unbinned_orfs = abbundance_parsed[2]
        percentages = abbundance_parsed[3]
        #Make the plot and save to png file:
        graph = plt.bar(assemblie_names, binned_orfs, color="r")
        plt.bar(assemblie_names, unbinned_orfs, bottom=binned_orfs, color="b")
        plt.xlabel("Assembly type")
        plt.ylabel("Amount of predicted ORFs")
        plt.legend(["Binned", "Unbinned"])
        plt.title("Total ORFs / Binned ORFs")

        i = 0
        for p in graph:
            width = p.get_width()
            height = p.get_height()
            x, y = p.get_xy()
            plt.text(x+width/2, y+height*1.01, str(round(percentages[i], 2))+'%', ha='center', weight='bold')
            i+=1

        plt.savefig(output_dir + "gene_abbundance_hist.png")
        plt.close

def main():
	abbundance_file = argv[1]
	output_dir = argv[2]
	abbundance_parsed = parse_gene_abbundance(abbundance_file)
	make_bar_chart(abbundance_parsed, output_dir)



if __name__ == "__main__":
	main()
