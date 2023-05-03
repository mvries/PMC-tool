Linkage of SSU genes to bins requires the following preperation steps:


First a blast search is performed (Use the blast conda environment in the Environments folder):
1. Create blast database:
makeblastdb -type nucl -in "SSU_sequences.fasta"

2. Run Blast for all bins ($BinFastaFile)
blastn \
-query $BinFastaFile \
-db "SSU_gene blastdb" \
-max_target_seqs 50 \
-outfmt 6 \
-evalue 0.000001 \
-num_threads 20 \
-out blastn-bins_Finals-full-nr99.50-hits-fasta.outfmt6

3. Filter Blast results (min ident 97% and min coverage 100bp)
cat blastn-bins_Finals-full-nr99.50-hits-fasta.outfmt6 | \
awk '$3 > 97' | awk '$4 > 100' | cut -f1,2,12 > filtered-BlastOut.tab

4. Create Blast results matrix
Needed Input files:
    filtered-BlastOut.tab
    list for all 16SrRNA names in text-format: "16SrRNAnames.txt"
    list for all bin names in text-format: "BinNames.txt"

Run reformat script (reformat-blast-results.sh)

reformat-blast-results.sh filtered-BlastOut.tab 16SrRNAnames.txt BinNames.txt

#final output: "matrix-Blast.tab"

