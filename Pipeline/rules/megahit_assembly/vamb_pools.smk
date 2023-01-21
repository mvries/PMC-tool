rule make_catalogue:
    input:
        contigs=f'{output_dir}' + "MEGAHIT/{pool}/{pool}.contigs.fa",
        report=f'{output_dir}' + "Quast/" + "{pool}" + "/report.html"
    output:
        f'{output_dir}' + "MEGAHIT/{pool}/{pool}.catalogue.fna.gz"
    conda:
        "../../envs/vamb.yaml"
    shell:
        "concatenate.py {output} {input.contigs}"
