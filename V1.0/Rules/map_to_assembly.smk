rule map_to_assembly:
    threads:
        config["PARAMS"]["MINIMAP"]["P"]
    params:
        P=config["PARAMS"]["MINIMAP"]["P"]
    input:
        assembly=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa",
        fw=f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq",
        rv=f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq"
    output:
        f'{output_dir}' + "minimap/assemblies/" + "{pool}.bam"
    conda:
        "../../envs/Quast.yaml"
    shell:
        "minimap2 -ax sr {input.assembly} -t {params.P} {input.fw} {input.rv} | samtools view -b -@ {params.p} | samtools sort -@ {params.p} -o {output}"
