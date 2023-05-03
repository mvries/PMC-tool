rule download_sample:
    threads:
        config["PARAMS"]["SRA"]["P"]
    params:
        P=config["PARAMS"]["SRA"]["P"],
        M=config["PARAMS"]["SRA"]["M"],
        B=config["PARAMS"]["SRA"]["B"],
        C=config["PARAMS"]["SRA"]["C"]
    conda:
        "../../envs/sra.yaml"
    output:
        r1=temporary(f'{output_dir}' + "sra_download/" + "{sample}/{sample}_1.fastq.gz"),
        r2=temporary(f'{output_dir}' + "sra_download/" + "{sample}/{sample}_2.fastq.gz")
    shell:
        "fasterq-dump --split-files {wildcards.sample}  -t /dev/shm/sra -c {params.C} -b {params.B} -e {params.P} -m {params.M} -O {output_dir}sra_download/{wildcards.sample}/ -p && pigz --fast {output_dir}sra_download/{wildcards.sample}/*"
