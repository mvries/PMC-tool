rule megahit:
    threads:
        config["PARAMS"]["MEGAHIT"]["P"]
    params:
        P=config["PARAMS"]["MEGAHIT"]["P"]
    input:
        fw=f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq",
        rv=f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq"
    output:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa"
    conda:
        "../../envs/megahit.yaml"
    shell:
        "megahit --presets meta-large -f -m 0.9 -t {params.P} -1 {input.fw} -2 {input.rv} -o {output_dir}MEGAHIT/{wildcards.pool} --out-prefix {wildcards.pool}"
