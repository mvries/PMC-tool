rule filter_contigs:
    threads:
        1
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa"
    output:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs_filter.fa"
    shell:
        "python3 ../../scripts/contig_filter_megahit.py {input} {output} 500"

rule make_catalogue:
    threads:
        1
    input:
        contigs=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs_filter.fa",
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
    output:
        unzipped=temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna"),
        i1=temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.1.bt2'),
        i2=temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.2.bt2'),
        i3=temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.3.bt2'),
        i4=temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.4.bt2'),
        rev1=temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.rev.1.bt2'),
        rev2=temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.rev.2.bt2')
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "gunzip -k -f {input.zipped} > {output.unzipped} && bowtie2-build --threads {params.P} {output.unzipped} {output.unzipped}"

rule map_to_assembly:
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"],
        S=config["PARAMS"]["BOWTIE2"]["S"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    input:
        catalogue=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna",
        r1=f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq",
        r2=f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq",
        i1=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.1.bt2',
        i2=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.2.bt2',
        i3=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.3.bt2',
        i4=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.4.bt2',
        rev1=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.rev.1.bt2',
        rev2=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna" + '.rev.2.bt2',
        contigs=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs_filter.fa"
    output:
        sam=temporary(f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.sam")
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "bowtie2 {params.S} -p {params.P} -x {input.catalogue} -1 {input.r1} -2 {input.r2} -S {output.sam}"

rule make_bam:
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    input:
        sam=f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.sam"
    output:
        bam=f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.bam"
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "samtools view -F 3584 -b --threads {params.P} {input.sam} > {output.bam}"
