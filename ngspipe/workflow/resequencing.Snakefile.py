configfile: "config.yaml"


rule all:
    input:
        report="report.html"


rule bwa_map:
    input:
        fa="data/chr19.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        bwa_bam=temp("mapped_reads/{sample}.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output} 2> {log}"


rule samtools_sort:
    input:
        bam="mapped_reads/{sample}.bam"
    output:
        bam_sort="sorted_reads/{sample}.bam"
    log:
        "logs/samtools_sort/{sample}.sorted.bam.log"
    threads: 4
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output} 2>{log}"


rule samtools_index:
    input:
        bam_sort="sorted_reads/{sample}.bam"
    output:
        bam_index="sorted_reads/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.samtools_index.log"
    threads: 4
    shell:
        "samtools index {input} 2>{log}"


rule bcftools_call:
    input:
        fa="data/chr16.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        bcftools_call="calls/all.vcf"
    log:
        "logs/bcftools_call/{sample}.bcftools_call.log"
    threads: 8
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output} 2>{log}"


rule report:
    input:
        bcftools_call="calls/all.vcf"
    output:
        report="report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        An example variant calling workflow
        ===================================

        Reads were mapped to the mouse
        reference chr19 genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])

