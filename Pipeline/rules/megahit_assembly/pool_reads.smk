#This snakemake rule is used to pool reads from a set of input fastq files:
rule pool_reads:
    input:
        f'{output_dir}' + "Pools/{pool}"
    threads:
        8
    output:
        fw_unzip=temporary(f'{output_dir}' + "Pools/{pool}/{sample}/{sample}_plant_removed_1.fq"),
        rv_unzip=temporary(f'{output_dir}' + "Pools/{pool}/{sample}/{sample}_plant_removed_2.fq"),
        fw=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq"),
        rv=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq")
    shell:
        "unpigz -k -p 8 {input}/*/*_1*.f*q.* && cat {output.fw_unzip} > {output.fw} && unpigz -k -p 8 {input}/*/*_2*.f*q.* && cat {output.rv_unzip} > {output.rv}"
