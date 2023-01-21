#This snakemake rule is used to pool reads from a set of input fastq files:
rule pool_reads:
    input:
        f'{output_dir}' + "Pools/{pool}"
    threads:
        8
    output:
        fw_unzip=temporary(f'{output_dir}' + "Pools/{pool}/unzipped_reads_1"),
        rv_unzip=temporary(f'{output_dir}' + "Pools/{pool}/unzipped_reads_2"),
        fw=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq"),
        rv=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq")
    shell:
        "unpigz -c -f -k -p 8 {input}/*/*_1*.f*q.* > {output.fw_unzip} && cat {output.fw_unzip} > {output.fw} && unpigz -c -f -k -p 8 {input}/*/*_2*.f*q.* > {output.rv_unzip} && cat {output.rv_unzip} > {output.rv}"
