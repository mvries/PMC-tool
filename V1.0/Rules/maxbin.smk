rule maxbin:
    threads:
        config["PARAMS"]["METABAT"]["P"]
    params:
        P=config["PARAMS"]["METABAT"]["P"]
    input:
        assembly=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa",
        depth={output_dir} + "metabat/{pool}/metabat_depth"
    output:
        {output_dir} + "maxbin/{pool}/"
    conda:
        "../../envs/concoct.yaml"
    shell:
        "run_maxbin.pl -contig {input.assembly} -abund {input.depth} o {output}"
