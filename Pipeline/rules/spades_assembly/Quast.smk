rule Quast:
    threads:
        config["PARAMS"]["QUAST"]["P"]
    params:
        P=config["PARAMS"]["QUAST"]["P"]
    input:
        f'{output_dir}' + "spades/" + "{sample}" + "/contigs.fasta"
    output:
        f'{output_dir}' + "Quast/" + "{sample}" + "/report.html"
    conda:
        "../envs/Quast.yaml"
    shell:
        "metaquast --max-ref-number 0 -t {params.P} -o {output_dir}Quast/{wildcards.sample} {input}"
