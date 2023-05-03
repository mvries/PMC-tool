rule map_to_assembly:
    threads:
        config["PARAMS"]["MINIMAP"]["P"]
    params:
        P=config["PARAMS"]["MINIMAP"]["P"]
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa"
    output:
        f'{output_dir}' + "Quast/" + "{pool}" + "/report.html"
    conda:
        "../../envs/Quast.yaml"
    shell:
        "minimap2 -ax sr {input.assembly} -t {params.P} -o {output_dir}Quast/{wildcards.pool} {input}"
        
        
        
MINIMAP:
  P: 32
  I: 500
  K: 1500M

    ./minimap2 -ax sr ref.fa read1.fa read2.fa > aln.sam
