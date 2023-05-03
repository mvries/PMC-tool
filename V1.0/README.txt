This github repo was originally made for my master thesis project "Building an integrated plant microbiome MAG catalogue".
A copy of my thesis is included in this git repo that contains aditional information on the tool and results that have been obtained.

This is in part a snakemake workflow.
The snakemake workflow handles the sample downloading, quality control, normalisation, assembly and subsequent binning with MAGSCoT.
The snakemake part of this workflow is handled by the PMC.py wrapper, the instructions can be found in the documentation of that script.

After the binning, Taxonomy is assigned to the MAG set using an SSU linkage approach first reported by Lesker et al.
The scripts and the instructions for the linkage procedure are in the linkage directory.

In addition, the scripts that were used to make the plots for my thesis with the data from this pipeline can also be found in the scripts diorectory.
