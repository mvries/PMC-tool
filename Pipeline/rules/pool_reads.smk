#This snakemake rule is used to pool reads from a set of input fastq files:
rule pool_reads:
    input:
        f'{output_dir}' + "Pools/{pool}"
    output:
        fw=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq.gz"),
        rv=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq.gz")
    shell:
        "cat {input}/*/*_1*.f*q.* > {output.fw} && cat {input}/*/*_2*.f*q.* > {output.rv}"
