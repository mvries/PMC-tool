#Rule to build the reference index:
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
        r1=f'{output_dir}' + 'sra_download/' + '{sample}/{sample}_1_trimmed.fq.gz',
        r2=f'{output_dir}' + 'sra_download/' + '{sample}/{sample}_2_trimmed.fq.gz',
        index=f'{phix_reference}' + '.rev.1.bt2'
    output:
        sam=temporary(f'{output_dir}' + "sra_download/" + "{sample}/{sample}_phix_bt2.sam")
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
        f'{output_dir}' + "sra_download/" + "{sample}/{sample}_phix_bt2.sam"
    output:
        out1=temporary(f'{output_dir}' + "sra_download/" + "{sample}/{sample}_phix_removed_1.fq.gz"),
        out2=temporary(f'{output_dir}' + "sra_download/" + "{sample}/{sample}_phix_removed_2.fq.gz")
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "samtools view -@ {params.P} -bS {input} | samtools view -@ {params.P} -b -f 12 -F 256 | samtools sort -n -@ {params.P} | samtools fastq -N -@ {params.P} -1 {output.out1} -2 {output.out2}"
