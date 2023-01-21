rule megahit:
    threads:
        1
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
        "heaptrack -o {output_dir}heaptrack megahit -f -m 0.9 -t {params.P} -1 {input.fw} -2 {input.rv} -o {output_dir}MEGAHIT/{wildcards.pool} --out-prefix {wildcards.pool}"

rule filter_contigs:
    threads:
        1
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa"
    output:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs_filter.fa"
    shell:
        "python3 ../../scripts/contig_filter_megahit.py {input} {output} 500"
