rule metabat:
    threads:
        config["PARAMS"]["METABAT"]["P"]
    params:
        P=config["PARAMS"]["METABAT"]["P"]
    input:
        assembly=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa",
        bam=f'{output_dir}' + "minimap/assemblies/" + "{pool}.bam"
    output:
        {output_dir} + "metabat/{pool}/"
    conda:
        "../Environments/metabat.yaml"
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth {output}metabat_depth {input.bam} && metabat2 -i {input.assembly} -a {output}metabat_depth -o {output}"
