rule fastp:
    threads:
        config["PARAMS"]["FASTP"]["P"]
    params:
        Q=config["PARAMS"]["FASTP"]["Q"],
        N=config["PARAMS"]["FASTP"]["N"],
        E=config["PARAMS"]["FASTP"]["E"],
        P=config["PARAMS"]["FASTP"]["P"]
    input:
        r1=f'{output_dir}' + "sra_download/" + "{sample}/{sample}_1.fastq.gz",
        r2=f'{output_dir}' + "sra_download/" + "{sample}/{sample}_2.fastq.gz"
    output:
        r1out=temporary(f'{output_dir}' + 'fastp_qc/' + '{sample}/{sample}_1_trimmed.fq.gz'),
        r2out=temporary(f'{output_dir}' + 'fastp_qc/' + '{sample}/{sample}_2_trimmed.fq.gz'),
        html=f'{output_dir}' + "fastp_qc/{sample}/" + '{sample}.html',
        json=f'{output_dir}' + "fastp_qc/{sample}/" + '{sample}.json'
    conda:
        "../../envs/fastp_env.yaml"
    shell:
        "fastp -w {params.P} -q {params.Q} -n {params.N} -e {params.E} -i {input.r1} -I {input.r2} -o {output.r1out} -O {output.r2out} -h {output.html} -j {output.json}"
