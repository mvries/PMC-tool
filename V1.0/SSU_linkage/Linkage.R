### Define input files ###

args1 <- commandArgs()

bins_file = as.character(args1[6])
bins_metadata_file = as.character(args1[7])
SSU_file = as.character(args1[8])
map_unam_file = as.character(args1[9])
blast_mat_file = as.character(args1[10])
SSU_lengths = as.character(args1[11])
res_folder = as.character(args1[12])

### Load input files ###

bin_abbundance_raw = read.csv(bins_file, sep = "\t", header = T)
SSU_abbundance_raw = read.csv(SSU_file, sep = "\t", header = T)

unambiguous_raw = read.csv(map_unam_file, sep = "\t", header = T)

blast_mat_raw = read.csv(blast_mat_file, sep = "\t", header = T)
blast_mat_2 = head(blast_mat_raw, - 1)

bins_metadata = read.csv(bins_metadata_file, sep = "\t", header = T)

SSU_lengths_file = read.csv(SSU_lengths, sep = "\t", header = T)


### Formatting of the input data ###

bin_abbundance_rownames <- bin_abbundance_raw[,1]
rownames(bin_abbundance_raw) <- bin_abbundance_rownames
bin_abbundance_matrix = as.matrix(bin_abbundance_raw[,-1])

SSU_abbundance_rownames <- SSU_abbundance_raw[,1]
rownames(SSU_abbundance_raw) <- SSU_abbundance_rownames
SSU_abbundance_matrix = as.matrix(SSU_abbundance_raw[,-1])

unambiguous_rownames <- unambiguous_raw[,1]
rownames(unambiguous_raw) <- unambiguous_rownames
unambiguous_matrix = as.matrix(unambiguous_raw[,-1])
colnames(unambiguous_matrix) <- sub('^out_', '', colnames(unambiguous_matrix))

blast_mat_rownames <- blast_mat_2[,1]
rownames(blast_mat_2) <- blast_mat_rownames
blast_mat = as.matrix(blast_mat_2[,-1])

metadata_rownames <- bins_metadata[,1]
metadata_matrix = as.matrix(bins_metadata[,-1])
rownames(metadata_matrix) <- metadata_rownames

SSU_lengths_rownames <- SSU_lengths_file[,1]
SSU_lengths_matrix = as.matrix(SSU_lengths_file[,-1])
rownames(SSU_lengths_matrix) <- SSU_lengths_rownames


### Normalisations ###
# First we normalise the bin abbundance (raw counts to TPM).
bin_names = rownames(bin_abbundance_matrix)

for(name in bin_names) {
        length = metadata_matrix[name,]
        counts = bin_abbundance_matrix[name,]
        newcount = counts / length
        bin_abbundance_matrix[name,] = newcount
}

bin_tpm_matrix <- t( t(bin_abbundance_matrix) * 1e6 / colSums(bin_abbundance_matrix))

#We also normalise the SSU gene abbundance (raw counts to TPM).
SSU_names = rownames(SSU_lengths_matrix)

for(name in SSU_names) {
        length = SSU_lengths_matrix[name,]
        counts = SSU_abbundance_matrix[name,]
        newcount = counts / length
        SSU_abbundance_matrix[name,] = newcount
}

SSU_tpm_matrix <- t( t(SSU_abbundance_matrix) * 1e6 / colSums(SSU_abbundance_matrix))

#In addition, we also need to normalize the reads after splitmapping (unbiquely mapped / total mapped) for each bin SSU pair.

norm_unambiguous_matrix <- t( t(unambiguous_matrix) / colSums(unambiguous_matrix))

#Then finally, we normalise the results of the blasting (Total length allignment / Total allignment of SSU genes)
norm_blast_mat <- t( t(blast_mat) / colSums(blast_mat) )

### Cleaning up steps ###

# We have to remove NaNs from the matrices after normalisation.

norm_unambiguous_matrix[is.nan(norm_unambiguous_matrix)] <- 0
norm_blast_mat[is.nan(norm_blast_mat)] <- 0

### Calculations  ###

# First we calculate the integrated scores for the correlations:
#For the correlation part we have to calculate pearsons and spearman correlations.

pearson_matrix = cor(t(bin_tpm_matrix), t(SSU_abbundance_matrix), method = "pearson")
spearmans_matrix = cor(t(bin_tpm_matrix), t(SSU_abbundance_matrix), method = "spearman")

#We take the absolute values of the spearman and pearson matrices:

absolute_pearson_mat = abs(pearson_matrix)
absolute_spearman_mat = abs(spearmans_matrix)

#Then we multiply the values in both matrices to obtain the correlation score:

correlation_score_matrix = absolute_pearson_mat * absolute_spearman_mat

#Next we calculate the linkage score of bin SSU pairs:

correlation_minus1 = 1 - correlation_score_matrix
unambiguous_minus1 = 1 - norm_unambiguous_matrix
blast_minus1 = 1 - norm_blast_mat

blast_minus1[blast_minus1==0] <- 0.1
unambiguous_minus1[unambiguous_minus1==0] <- 0.1

correlation_root = correlation_minus1^(1/3)
unambiguous_root = unambiguous_minus1^(1/3)
blast_root = blast_minus1^(1/3)

correlation_root_t = t(correlation_root)

### Loop over and extract results ###

bins = colnames(blast_mat)
genes = rownames(blast_minus1)


scoremat = blast_root

for(bin in bins) {
        for(gene in genes) {
                cor = correlation_root_t[gene, bin]
                blast = blast_root[gene, bin]
                count = unambiguous_root[gene, bin]
		score = cor * blast * count
		scoremat[gene, bin] = 1 - score
        }
}


bin_high_score = max(colSums(scoremat))
SSU_high_score = max(rowSums(scoremat))

scoremat_norm1 = scoremat / bin_high_score
scoremat_norm2 = scoremat / SSU_high_score

final_mat = scoremat_norm1 * scoremat_norm2

top_hits_all=lapply(1:ncol(final_mat),function(x){sort(final_mat[,x],decreasing = T)[1:2]})
names(top_hits_all)=colnames(final_mat)

sink("bin-to-SSU.tophits")
print(top_hits_all)
