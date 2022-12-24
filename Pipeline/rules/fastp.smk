#This is a snakemake rule that is used for the qc of short read data.
#This rule takes as input paired end read files and outputs: Trimmed reads and quality reports.
#The quality control tool used in this case is fastp.
#By default reports are stored in the same dir as the input reads.

#Rule to run fastp quality control:
rule fastp:
    threads:
        config["PARAMS"]["FASTP"]["P"]
    params:
        Q=config["PARAMS"]["FASTP"]["Q"],
        N=config["PARAMS"]["FASTP"]["N"],
        E=config["PARAMS"]["FASTP"]["E"],
        P=config["PARAMS"]["FASTP"]["P"]
    input:
        sample=f'{sample_dir}' + "{sample}"
    output:
        r1out=temporary(f'{sample_dir}' + '{sample}' + '/{sample}_1_trimmed.fq.gz'),
        r2out=temporary(f'{sample_dir}' + '{sample}' + '/{sample}_2_trimmed.fq.gz'),
        html=f'{output_dir}' + "fastp_qc/{sample}/" + '{sample}.html',
        json=f'{output_dir}' + "fastp_qc/{sample}/" + '{sample}.json'
    conda:
        "../envs/fastp_env.yaml"
    shell:
        "fastp -w {params.P} -q {params.Q} -n {params.N} -e {params.E} -i {input}/*1*f*q.gz -I {input}/*2*f*q.gz -o {output.r1out} -O {output.r2out} -h {output.html} -j {output.json}"
