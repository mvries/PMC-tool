rule bbnorm:
    input:
    	r1=f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_1.fq.gz",
    	r2=f'{output_dir}' + "bowtie2/plant/" + "{sample}/{sample}_plant_removed_2.fq.gz"
    output:
        out1=f'{output_dir}' + "bbmap/" + "{sample}/{sample}_normalized_1.fq.gz",
        out2=f'{output_dir}' + "bbmap/" + "{sample}/{sample}_normalized_2.fq.gz",
        hist1=f'{output_dir}' + "bbmap/" + "{sample}/{sample}.hist1",
        hist2=f'{output_dir}' + "bbmap/" + "{sample}/{sample}.hist2"
    threads:
        config["PARAMS"]["BBMAP"]["P"]
    params:
        depth=config["PARAMS"]["BBMAP"]["D"],
        P=config["PARAMS"]["BBMAP"]["P"]
    conda:
        "../../envs/bbmap.yaml"
    log:
        f'{output_dir}' + "bbmap/" + "{sample}/{sample}_bbnorm.log"
    shell:
        "(bbnorm.sh target={params.depth} minprob=0.6 prefiltersize=0.50 prefilter=True in={input.r1} in2={input.r2} threads={params.P} out={output.out1} outt={output.out2} hist={output.hist1} histout={output.hist2} -Xmx100g) 2>{log}"
