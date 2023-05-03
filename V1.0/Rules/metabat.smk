rule metabat:
    threads:
        config["PARAMS"]["QUAST"]["P"]
    params:
        P=config["PARAMS"]["QUAST"]["P"]
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa"
    output:
        f'{output_dir}' + "Quast/" + "{pool}" + "/report.html"
    conda:
        "../../envs/metabat.yaml"
    shell:
        "metabat2 -i {input} -t {params.P} -o {output_dir}Quast/{wildcards.pool} {input}"
