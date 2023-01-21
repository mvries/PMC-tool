rule make_catalogue:
    threads:
        1
    input:
        contigs=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa",
        report=f'{output_dir}' + "Quast/" + "{pool}" + "/report.html"
    output:
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna.gz")
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
        zipped=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna.gz",
        unzipped=temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna")
    output:
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.1.bt2'),
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.2.bt2'),
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.3.bt2'),
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.4.bt2'),
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.rev.1.bt2'),
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.rev.2.bt2')
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "gunzip -k -f {input.zip} > {input.unzipped} && bowtie2-build {input.unzipped} {input.unzipped}"
