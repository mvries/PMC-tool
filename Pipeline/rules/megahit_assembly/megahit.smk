rule megahit:
    threads:
        1
    params:
        P=config["PARAMS"]["MEGAHIT"]["P"]
    input:
        fw=f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq.gz",
        rv=f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq.gz"
    output:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa"
    conda:
        "../../envs/megahit.yaml"
    shell:
        "heaptrack -o {output_dir}heaptrack megahit -f -m 0.9 -t {params.P} -1 {input.fw} -2 {input.rv} -o {output_dir}MEGAHIT/{wildcards.pool} --out-prefix {wildcards.pool}"
