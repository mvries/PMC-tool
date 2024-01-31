rule concoct:
    threads:
        config["PARAMS"]["METABAT"]["P"]
    params:
        P=config["PARAMS"]["METABAT"]["P"]
    input:
        assembly=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa",
        depth={output_dir} + "metabat/{pool}/metabat_depth"
    output:
        {output_dir} + "concoct/{pool}/"
    conda:
        "../Environments/concoct.yaml"
    shell:
        "concoct -t {params.P} -i {input.assembly} --coverage-file {input.depth} -o {output}"
