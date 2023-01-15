re	#Rule to build the reference index:
rule bowtie_build_phix:
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    input:
        f'{phix_reference}'
    output:
        f'{phix_reference}' + '.rev.1.bt2'
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build {input} {input}"

#Rule to map reads to the reference:
rule bowtie_map_phix:
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"],
        S=config["PARAMS"]["BOWTIE2"]["S"]
    input:
        r1=f'{output_dir}' + 'fastp_qc/' + '{sample}/{sample}_1_trimmed.fq.gz',
        r2=f'{output_dir}' + 'fastp_qc/' + '{sample}/{sample}_2_trimmed.fq.gz',
        index=f'{phix_reference}' + '.rev.1.bt2'
    output:
        sam=temporary(f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_phix_bt2.sam")
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2 {params.S} -p {params.P} -x {phix_reference} -1 {input.r1} -2 {input.r2} -S {output.sam}"

rule remove_phix:
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
    input:
        f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_phix_bt2.sam"
    output:
        bam1=temporary(f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_1.bam"),
        bam2=temporary(f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_2.bam"),
        sorted=temporary(f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_sorted.bam"),
        out1=temporary(f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_phix_removed_1.fq.gz"),
        out2=temporary(f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_phix_removed_2.fq.gz")
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "samtools view -@ {params.P} -bS {input} -o {output.bam1} && samtools view {output.bam1} -@ {params.P} -b -f 12 -F 256 -o {output.bam2} && samtools sort {output.bam2} -n -@ {params.P} -o {output.sorted} && samtools fastq -N -@ {params.P} -1 {output.out1} -2 {output.out2} {output.sorted}"
