rule metaspades:
    threads:
        config["PARAMS"]["SPADES"]["P"]
    params:
        P=config["PARAMS"]["SPADES"]["P"],
        M=config["PARAMS"]["SPADES"]["M"]
    input:
		r1=(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_1.fq.gz"),
		r2=(f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_2.fq.gz")
    output:
        contigs=f'{output_dir}' + "spades/" + "{sample}" + "/contigs.fasta"
    conda:
        "../../envs/metaspades.yaml"
    shell:
        "spades.py --meta -m {params.M} -t {params.P} -1 {input.r1} -2 {input.r2} -o {output_dir}spades/{wildcards.sample}/"
