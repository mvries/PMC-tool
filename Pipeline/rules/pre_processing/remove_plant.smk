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
        "../../envs/bowtie2.yaml"
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
        r1=f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_phix_removed_1.fq.gz",
        r2=f'{output_dir}' + "bowtie2/phix/" + "{sample}/{sample}_phix_removed_2.fq.gz",
        index=f'{plant_reference}' + '.rev.1.bt2'
    output:
        sam=temporary(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_bt2.sam")
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "bowtie2 {params.S} -p {params.P} -x {plant_reference} -1 {input.r1} -2 {input.r2} -S {output.sam}"

rule remove_plant:
	threads:
		config["PARAMS"]["BOWTIE2"]["P"]
	params:
		P=config["PARAMS"]["BOWTIE2"]["P"]
	input:
		f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_bt2.sam"
	output:
		bamplant1=temporary(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_1_plant.bam"),
		bamplant2=temporary(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_2_plant.bam"),
		sortedplant=temporary(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_sorted_plant.bam"),
		out1=temporary(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_1.fq.gz"),
		out2=temporary(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_2.fq.gz")
	conda:
		"../../envs/bowtie2.yaml"
	shell:
		"samtools view -@ {params.P} -bS {input} -o {output.bamplant1} && samtools view {output.bamplant1} -@ {params.P} -b -f 12 -F 256 -o {output.bamplant2} && samtools sort {output.bamplant2} -n -@ {params.P} -o {output.sortedplant} && samtools fastq -N -@ {params.P} -1 {output.out1} -2 {output.out2} {output.sortedplant}"
