rule bbnorm:
    input:
    	r1=f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_1.fq.gz",
    	r2=f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_2.fq.gz"
    output:
        out1=f'{output_dir}' + "bbmap/" + "{sample}/{sample}_normalized_1.fq.gz",
        out2=f'{output_dir}' + "bbmap/" + "{sample}/{sample}_normalized_2.fq.gz",
        toss1=f'{output_dir}' + "bbmap/" + "{sample}/{sample}_tossed_1.fq.gz",
        toss2=f'{output_dir}' + "bbmap/" + "{sample}/{sample}_tossed_2.fq.gz",
        hist=f'{output_dir}' + "bbmap/" + "{sample}/{sample}.hist"
    threads:
        config["PARAMS"]["BBMAP"]["P"]
    params:
        depth=config["PARAMS"]["BBMAP"]["D"]
    conda:
        "../../envs/bbmap.yaml"
    shell:
        "bbnorm.sh in={input.r1} in2={input.r2} out={output.out1} out2={output.out2} outt={output.toss1} outt2={output.toss2} hist={output.hist} target={params.depth}"
