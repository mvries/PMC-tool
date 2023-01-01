rule metaspades:
    threads:
        config["PARAMS"]["SPADES"]["P"]
    params:
        P=config["PARAMS"]["SPADES"]["P"],
    input:
        r1=f'{output_dir}' + "sra_download/" + "{sample}/{sample}_plant_removed_1.fq.gz",
        r2=f'{output_dir}' + "sra_download/" + "{sample}/{sample}_plant_removed_2.fq.gz"
    output:
        outdir=directory(f'{output_dir}' + "spades/" + "{sample}/"),
    conda:
        "../envs/metaspades.yaml"
    shell:
        "spades.py --meta -t {params.P} -1 {input.r1} -2 {input.r2} -o {output.outdir}"
