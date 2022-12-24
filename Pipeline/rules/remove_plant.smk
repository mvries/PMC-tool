#Rule to build the reference index:
rule bowtie_build_plant:
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    input:
        f'{plant_reference}'
    output:
        f'{plant_reference}' + '.rev.1.bt2'
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build {input} {input}"

#Rule to map reads to the reference:
rule bowtie_map_plant:
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"],
        S=config["PARAMS"]["BOWTIE2"]["S"]
    input:
        r1=f'{sample_dir}' + "{sample}" + "/{sample}_phix_removed_1.fq.gz",
        r2=f'{sample_dir}' + "{sample}" + "/{sample}_phix_removed_2.fq.gz",
        index=f'{plant_reference}' + '.rev.1.bt2'
    output:
        sam=temporary(f'{sample_dir}' + "{sample}" + "/{sample}_plant_bt2.sam"),
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2 {params.S} -p {params.P} -x {plant_reference} -1 {input.r1} -2 {input.r2} -S {output.sam}"

rule remove_plant:
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
    input:
        f'{sample_dir}' + "{sample}" + "/{sample}_plant_bt2.sam"
    output:
        out1=f'{sample_dir}' + "{sample}" + "/{sample}_plant_removed_1.fq.gz",
        out2=f'{sample_dir}' + "{sample}" + "/{sample}_plant_removed_2.fq.gz"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "samtools view -@ {params.P} -bS {input} | samtools view -@ {params.P} -b -f 12 -F 256 | samtools sort -n -@ {params.P} | samtools fastq -@ {params.P} -1 {output.out1} -2 {output.out2}"
