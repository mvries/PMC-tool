rule make_catalogue:
    threads:
        1
    input:
        contigs=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa",
        report=f'{output_dir}' + "Quast/" + "{pool}" + "/report.html"
    output:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna.gz"
    conda:
        "../../envs/vamb.yaml"
    shell:
        "concatenate.py {output} {input.contigs}"

rule make_assembly_index:
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna.gz"
    output:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna.gz" + '.rev.1.bt2'
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "bowtie2-build {input} {input}"
