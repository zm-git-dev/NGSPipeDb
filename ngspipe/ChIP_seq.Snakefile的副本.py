# -*- coding: utf-8 -*-
import os
import pandas as pd
from os.path import join

# ----------------------------------------------------------------------- #
# sample information #
#
smpList = pd.read_csv(config["samplesList"], index_col=0, header=None)
SAMPLES = list(smpList.index)[0:]
# ----------------------------------------------------------------------- #


rule all:    #全局关键词，作用整个snakemake文件中的流程
    input:    #input用于定义整个流程需要输出的文件
        clean_R1 = expand(join(config['resultsDir'], "chipseq", "clean_fastq", "{sample}"+config["read1Suffix"]),sample=SAMPLES),    #cutadapt输出的R1.fastq
        clean_R2 = expand(join(config['resultsDir'], "chipseq", "clean_fastq", "{sample}"+config["read2Suffix"]),sample=SAMPLES),     #cutadapt输出的R2.fastq
        rmdup_bam = expand(join(config['resultsDir'], "chipseq", "picard", "{sample}.sorted_rmdup.bam"),sample=SAMPLES),
        #expand("bam/{sample}_bwa_mm10_sorted_rmdup.bam",sample=SAMPLES),    #输出比对后bam文件
        #expand("macs2_result/10-me_peaks.broadPeak")    #输出peaks

#expand作用是用于遍历所有SAMPLE中的样本，省略了循环的写法，整个snakemake file里的变量需要在rule all里指出，这样变量才会正常起作用,变量引用方式 {变量name}


rule bwa_index:
    input:
        genomeFa = config['genomeFasta']
    output:
        indexOk = touch(join(config['resultsDir'], "chipseq","bwa_index", "index.ok"))
    log:
        join(config['resultsDir'], "chipseq","bwa_index", "index.log")
    threads: 8
    params:
        prefix = join(config['resultsDir'], "chipseq","bwa_index", "bwa_genome")
    shell:
        "bwa index -p {params.prefix} {input.genomeFa} 1>{log} 2>&1"   #{input[0]}调用input关键字下的第一行，即"clean_fastq/{sample}_R1.fastq.gz"。

rule bwa_mapping:
    input:
        clean_R1 = join(config['resultsDir'], "chipseq", "clean_fastq", "{sample}"+config["read1Suffix"]),    #cutadapt输出的R1.fastq
        clean_R2 = join(config['resultsDir'], "chipseq", "clean_fastq", "{sample}"+config["read2Suffix"]),    #cutadapt输出的R2.fastq
        indexOk = touch(join(config['resultsDir'], "chipseq","bwa_index", "index.ok"))
    output:
        bwa_bam = join(config['resultsDir'], "chipseq", "mapping", "{sample}.sam")
    log:
        join(config['resultsDir'], "chipseq","mapping", "{sample}.log")
    params:
        prefix = join(config['resultsDir'], "chipseq","bwa_index", "bwa_genome")
    threads: 8
    shell:
        "bwa mem -t 10 -M {params.prefix} \
        {input.clean_R1} {input.clean_R2} 1>{output.bwa_bam}  2>{log}"   #{input[0]}调用input关键字下的第一行，即"clean_fastq/{sample}_R1.fastq.gz"。

rule bam_sort:
    input:
        bwa_bam = join(config['resultsDir'], "chipseq", "mapping", "{sample}.sam")
    output:
        bwa_bam_sort = join(config['resultsDir'], "chipseq", "sorted_bam", "{sample}.sorted.bam")
    log:
        join(config['resultsDir'], "chipseq", "sorted_bam", "{sample}.log")
    threads: 4
    shell:
        "samtools sort -O BAM -@ 4 {input.bwa_bam} \
        -o {output.bwa_bam_sort} \
        -T {output}.tem "    #将output文件以临时文件形式保存，流程结束自动删除以节约服务器空间

rule picard_remove_duplication:
    input:
        bwa_bam_sort = join(config['resultsDir'], "chipseq", "sorted_bam", "{sample}.sorted.bam")
    output:
        rmdup_bam = join(config['resultsDir'], "chipseq", "picard", "{sample}.sorted_rmdup.bam"),
        rmdup_matrix = join(config['resultsDir'], "chipseq", "picard", "{sample}.sorted_rmdup.matrix")
    log:
        join(config['resultsDir'], "chipseq", "picard", "{sample}.log")
    threads: 4
    shell:
        "picard MarkDuplicates\
        REMOVE_DUPLICATES=true \
        I= {input.bwa_bam_sort} \
        O={output.rmdup_bam} M={output.rmdup_matrix}"

rule macs2_call_peak:    #有control的call peaks
    input:
        treat="bam/10-me_bwa_mm10_sorted_rmdup.bam",    #实验组bam
        control="bam/10-Input_bwa_mm10_sorted_rmdup.bam"    #配对的对照bam
    output:
        "macs2_result/10-me_peaks.broadPeak"
    params:    #params关键词作用，创建局部变量，引用方式 {params.name}
        outdir="macs2_result",
        head_outfile="10-me"
    log:
        "macs2_result/10-me_peaks.log"
    threads:4
    shell:
        "macs2 callpeak -t {input.treat} -c {input.control} -f BAM  -g hs -n {params.head_outfile}  \
         --outdir {params.outdir} \
         --nomodel --extsize 200 \
         -B --broad --broad-cutoff  0.01 "


