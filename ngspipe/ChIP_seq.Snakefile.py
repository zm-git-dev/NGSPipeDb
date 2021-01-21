SAMPLES = {"10-me","10-Input"}    #区分不同样本的名字
rule all:    #全局关键词，作用整个snakemake文件中的流程
    input:    #input用于定义整个流程需要输出的文件
        expand("clean_fastq/{sample}_R1.fastq.gz",sample=SAMPLES),    #输出去接头后的R1 fastq
        expand("clean_fastq/{sample}_R2.fastq.gz",sample=SAMPLES),    #输出去接头后的R2 fastq
        expand("bam/{sample}_bwa_mm10_sorted_rmdup.bam",sample=SAMPLES),    #输出比对后bam文件
        expand("macs2_result/10-me_peaks.broadPeak")    #输出peaks
#expand作用是用于遍历所有SAMPLE中的样本，省略了循环的写法，整个snakemake file里的变量需要在rule all里指出，这样变量才会正常起作用,变量引用方式 {变量name}
rule cutadapt:    #cutadapt关键词，定义cutdapt的流程，主要包括input，output，shell三个部分
    input:
        raw_R1="rawdata/{sample}_R1.fastq.gz",    #rawdata的R1.fastq
        raw_R2="rawdata/{sample}_R2.fastq.gz"    #rawdata的R2.fastq
    output:
        clean_R1="clean_fastq/{sample}_R1.fastq.gz",    #cutadapt输出的R1.fastq
        clean_R2="clean_fastq/{sample}_R2.fastq.gz"    #cutadapt输出的R2.fastq
    log:
        "clean_fastq/{sample}_cutadapt.log"
    threads: 4    #所需的核心数，与-j保持一致
    shell:    #cutadapt的代码
        "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -u 5 -u -0 -U 5 -U -0 -m 30 -j 4 -o {output.clean_R1} -p {output.clean_R2} {input.raw_R1} {input.raw_R2} 1>{log} 2>&1"

rule bwa_mapping:   
    input:
        clean_R1="clean_fastq/{sample}_R1.fastq.gz",    #每行以逗号分隔为标志
        clean_R2="clean_fastq/{sample}_R2.fastq.gz",
        index="chr19.fa"
    output:
        bwa_bam="bam/{sample}_bwa_mm10.sam"
    log:
        "bam/{sample}_bwa_mm10.log"
    threads: 8
    shell:
        "bwa mem -t 10 -M {input.index} \
        {input[0]} {input[1]} 1>{output}  2>{log}"   #{input[0]}调用input关键字下的第一行，即"clean_fastq/{sample}_R1.fastq.gz"。
rule bam_sort:
    input:
        bwa_bam="bam/{sample}_bwa_mm10.sam"
    output:
        bwa_bam_sort="bam/{sample}_bwa_mm10_sorted.bam"
    log:
        "bam/{sample}_bwa_mm10_sorted.bam"
    threads: 4
    shell:
        "samtools sort -O BAM -@ 4 {input} \
        -o {output} \
        -T {output}.tem "    #将output文件以临时文件形式保存，流程结束自动删除以节约服务器空间
rule picard_remove_duplication:
    input:
        bwa_bam_sort="bam/{sample}_bwa_mm10_sorted.bam"
    output:
        rmdup_bam="bam/{sample}_bwa_mm10_sorted_rmdup.bam",
        rmdup_matrix="bam/{sample}_bwa_mm10_sorted_rmdup.matrix"
    log:
        "bam/{sample}_bwa_mm10_sorted_rmdup.log"
    threads: 4
    shell:
        "picard MarkDuplicates\
        REMOVE_DUPLICATES=true \
        I= {input} \
        O={output[0]}  M={output[1]}"
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


