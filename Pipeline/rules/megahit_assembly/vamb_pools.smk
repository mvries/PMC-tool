rule filter_contigs:
    threads:
        1
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa"
    output:
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs_filter.fa")
    shell:
        "python3 scripts/contig_filter_megahit.py {input} {output} 1000"

rule make_catalogue:
    threads:
        1
    input:
        contigs=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs_filter.fa",
        report=f'{output_dir}' + "Quast/" + "{pool}" + "/report.html"
    output:
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna")
    conda:
        "../../envs/vamb.yaml"
    shell:
        "concatenate.py {output} {input.contigs} -m 1000"

rule index:
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna"
    output:
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.flt.mmi")
    params:
        I=config["PARAMS"]["MINIMAP"]["I"]
    threads:
        1
    conda:
        "../../envs/minimap2.yaml"
    shell:
        "minimap2 -I {params.I} -d {output} {input}"

rule dict:
    input:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna"
    output:
        temporary(f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.flt.dict")
    threads:
        1
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "samtools dict {input} | cut -f1-3 > {output}"

rule map_to_assembly:
    params:
        P=config["PARAMS"]["MINIMAP"]["P"],
        K=config["PARAMS"]["MINIMAP"]["K"]
    threads:
        config["PARAMS"]["MINIMAP"]["P"]
    input:
        dict=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.flt.dict",
        r1=f'{output_dir}' + "Pools/{pool}/pooled_reads_1.fq",
        r2=f'{output_dir}' + "Pools/{pool}/pooled_reads_2.fq",
        mmi=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.flt.mmi"
    output:
        f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2.bam"
    log:
        f'{output_dir}' + "minimap2/" + "{pool}/mm2.log"
    conda:
        "../../envs/minimap2.yaml"
    shell:
        '''minimap2 -t {params.P} -ax sr {input.mmi} {input.r1} {input.r2} -N 5 -K {params.K} | grep -v "^@" | cat {input.dict} - | samtools view -F 3584 -b - > {output} 2>{log}'''

rule sort_bam:
    params:
        P=config["PARAMS"]["BOWTIE2"]["P"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    input:
        f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2.bam"
    output:
        f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2_sorted.bam"
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        "samtools sort {input} --threads {params.P} -m 10G -o {output}"

rule jgi:
    threads:
        1
    input:
        f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2_sorted.bam"
    output:
        temporary(f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2.raw.jgi")
    conda:
        "../../envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output} {input}"

rule cut_column1to3:
    threads:
        1
    input:
        f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2.raw.jgi"
    output:
        temporary(f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2_13.raw.jgi")
    shell:
        "cut -f1-3 {input} > {output}"

rule cut_column4to5:
    threads:
        1
    input:
        f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2.raw.jgi"
    output:
        temporary(f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2_45.raw.jgi")
    shell:
        "cut -f1-3 --complement {input} > {output}"

rule paste_abundances:
    threads:
        1
    input:
        cut13=f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2_13.raw.jgi",
        cut45=f'{output_dir}' + "minimap2/" + "{pool}/{pool}_mm2_45.raw.jgi"
    output:
        temporary(f'{output_dir}' + "minimap2/" + "{pool}/{pool}.jgi.abundance.dat")
    shell:
        "paste {input.cut13} {input.cut45} > {output}"

rule vamb:
    input:
        jgi=f'{output_dir}' + "minimap2/" + "{pool}/{pool}.jgi.abundance.dat",
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
        P=config["PARAMS"]["VAMB"]["P"],
        O=config["PARAMS"]["VAMB"]["O"],
        M=config["PARAMS"]["VAMB"]["M"],
        MIN=config["PARAMS"]["VAMB"]["MIN"]
    threads:
        config["PARAMS"]["BOWTIE2"]["P"]
    conda:
        "../../envs/vamb.yaml"
    shell:
        "rm -r {output_dir}vamb/{wildcards.pool}/;"
        "vamb --outdir {output_dir}vamb/{wildcards.pool}/ --fasta {input.contigs} --jgi {input.jgi} -o {params.O} -m {params.M} --minfasta {params.MIN}"

rule checkm:
    input:
        f'{output_dir}' + "vamb/" + "{pool}/clusters.tsv"
    output:
        f'{output_dir}' + "vamb/" + "{pool}/checkm.results"
    params:
        bins=f'{output_dir}' + "vamb/" + "{pool}/bins/",
        P=config["PARAMS"]["CHECKM"]["P"]
    threads:
        config["PARAMS"]["CHECKM"]["P"]
    conda:
        "../../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -f {output} -t {params.P} -x fna {params.bins} {output_dir}vamb/{wildcards.pool}"
