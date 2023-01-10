rule Quast:
    threads:
        config["PARAMS"]["QUAST"]["P"]
    params:
        P=["PARAMS"]["QUAST"]["P"]
    input:
        f'{output_dir}' + "spades/" + "{sample}" + "/contigs.fasta"
    output:
        f'{output_dir}' + "Quast/" + "{sample}" + "/report.tsv"
    conda:
        "../envs/Quast.yaml"
    shell:
        "metaquast -t {params.P} -o {output_dir}Quast/{sample} {input}"
