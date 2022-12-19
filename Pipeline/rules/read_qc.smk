#This is a snakemake rule that is used for the qc of short read data.
#This rule takes as input paired end read files and outputs: Trimmed reads and quality reports.
#The quality control tool used in this case is fastp.
#By default reports are stored in the same dir as the input reads.

#Rule to run fastp quality control:
rule fastp:
    input:
        r1="{sample_dir}{sample}",
        r2="{sample_dir}{sample}"
    output:
        r1out="{sample_dir}{sample}/{sample}_1_trimmed.fq.gz",
        r2out="{sample_dir}{sample}/{sample}_2_trimmed.fq.gz",
        html="{sample_dir}qc/{sample}/{sample}.html",
        json="{sample_dir}qc/{sample}/{sample}.json"
    conda:
        "../envs/fastp_env.yaml"
    shell:
        "fastp -i {input.r1}/*1*f*q.gz -I {input.r2}/*2*f*q.gz -o {output.r1out} -O {output.r2out} -h {output.html} -j {output.json}"
