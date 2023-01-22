rule filter_contigs:
    threads:
        1
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa"
    output:
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs_filter.fa")
    shell:
        "python3 scripts/contig_filter_megahit.py {input} {output} 500"

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
        bam=temporary(f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.bam")
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "samtools view -F 3584 -b --threads {params.P} {input.sam} > {output.bam}"

rule sort_bam:
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    input:
        f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.bam"
    output:
        temporary(f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2_sorted.bam")
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "samtools sort {input} --threads {params.P} -m 10G -o {output}"

rule jgi:
    threads:
        1
    input:
        f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2_sorted.bam"
    output:
        temporary(f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.raw.jgi")
    conda:
        "../../envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output} {input}"

rule cut_column1to3:
    threads:
        1
    input:
        f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.raw.jgi"
    output:
        temporary(f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.raw_13.jgi")
    shell:
        "cut -f1-3 {input} > {output}"

rule cut_column4to5:
    threads:
        1
    input:
        f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.raw.jgi"
    output:
        temporary(f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.raw_45.jgi")
    shell:
        "cut -f1-3 --complement {input} > {output}"

rule paste_abundances:
    threads:
        1
    input:
        cut13=f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.raw_13.jgi",
        cut45=f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}_bt2.raw_45.jgi"
    output:
        temporary(f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}.jgi.abundance.dat")
    shell:
        "paste {input.cut13} {input.cut45} > {output}"

rule vamb:
    input:
        jgi=f'{output_dir}' + "bowtie2/assembly/" + "{pool}/{pool}.jgi.abundance.dat",
        contigs=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna"
    output:
        f'{output_dir}' + "vamb/" + "{pool}/clusters.tsv",
        f'{output_dir}' + "vamb/" + "{pool}/latent.npz",
        f'{output_dir}' + "vamb/" + "{pool}/lengths.npz",
        f'{output_dir}' + "vamb/" + "{pool}/log.txt",
        f'{output_dir}' + "vamb/" + "{pool}/model.pt",
        f'{output_dir}' + "vamb/" + "{pool}/mask.npz",
        f'{output_dir}' + "vamb/" + "{pool}/tnf.npz"
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
        O=config["PARAMS"]["BOWTIE2"]["O"]
        M=config["PARAMS"]["BOWTIE2"]["M"]
        MIN=config["PARAMS"]["BOWTIE2"]["MIN"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    conda:
        "../../envs/vamb.yaml"
    shell:
        "vamb --outdir {output_dir}vamb/{pool}/ --fasta {input.contigs} --jgi {input.jgi} -o {params.O} -m {params.M} --minfasta {params.MIN}"
