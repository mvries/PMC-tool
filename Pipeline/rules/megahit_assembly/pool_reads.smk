#This snakemake rule is used to pool reads from a set of input fastq files:
rule pool_reads:
    input:
        f'{output_dir}' + "Pools/{pool}"
    threads:
        1
    output:
        fw=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq"),
        rv=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq"),
        fw_zip=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq.gz"),
        rv_zip=temporary(f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq.gz")
    shell:
        "zcat {input}/*/*_1*.f*q.* > {output.fw} && gzip -k {output.fw} && zcat {input}/*/*_2*.f*q.* > {output.rv} && gzip -k {output.rv}"
