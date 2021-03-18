# -*- coding: utf-8 -*-
import os
from os.path import join
import sys
import pandas as pd

# relative path
snake_dir = workflow.basedir # all configfile, scripts, restructuretext, ens are relative to snakefile (this file)
working_dir = os.getcwd() # input and output path are relative to current working directory

# get notice
receiver_email = 'nobody'

# ----------------------------------------------------------------------- #
# sample information #
#
smpList = pd.read_csv(config["samplesList"], index_col=0, header=None)
SAMPLES = list(smpList.index)[0:]
# ----------------------------------------------------------------------- #


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# detail parameters in pipe #
#
# 1. sampling data
# for test the pipe, you can choose to the part of the input file, can be whole,head:40000,tail:40000,random:0.5,random:40000
sampling_method = 'links' # tail, seqkit_number, seqkit_proportion, head, tail
sampling_data_outdir = join(config["resultsDir"], "sampling_data", "sampling_data_by_{}".format(sampling_method))

# 2. raw reads qc
qc_method = 'trim-galore' # trimomatic
qc_outdir = join(config["resultsDir"], "rawReads_qc", "rawReads_qc_by_{}".format(qc_method))

# 3. mapping to genome
mapping_method = 'bwa' # bowtie, bwa, 
mapping_outdir = join(config["resultsDir"], "mapping", "mapping_by_{}".format(mapping_method))
genome_index_prefix = 'genome'

# 4. call variation
call_variation_method = 'bcftools' # samtools, bcftools or gatk
variation_outdir = join(config["resultsDir"], "call_variation", "call_by_{}".format(call_variation_method))

# 5. variation annotation
variation_annotation_method = 'snpeff' # snpeff
annotation_outdir = join(config["resultsDir"], "variation_annotation", "annotation_by_{}".format('variation_annotation_method'))
annoFormat = 'gtf22' # gff3, gff2, gtf22
dbname = 'ngspipe-resequencing',
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# ------------------------------
#report: "report/workflow.rst"
report_outdir = join(config["reportsDir"], "report.html")
# ------------------------------

rule all:
    input:
        #
        # 1. sampling data #
        sampling_data_result                = expand(join(sampling_data_outdir, "{sample}"+config["read1Suffix"]), sample=SAMPLES),
        #
        # 2. raw reads qc #
        rawReads_qc_result                 = expand(join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"), sample=SAMPLES),
        #
        # 3. mapping #
        mapping_result                     = expand(join(mapping_outdir, "{sample}", "{sample}.sorted.bam"), sample=SAMPLES),
        #
        # 4. call variation #
        variation_result                   = expand(join(variation_outdir, "{sample}", "{sample}.vcf"), sample=SAMPLES),
        #
        # 5. annotation #
        #
        annotation_result                  = expand(join(annotation_outdir, "{sample}", "{sample}.snpeff.vcf"), sample=SAMPLES),
        #
        # 6. report #
        #report_result    = join(config['reportsDir'], "report.ok"),
        #


onsuccess:
    print("""
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Workflow finished, no error <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 ▄▄▄▄▄▄▄▄▄▄▄  ▄         ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄ 
▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
▐░█▀▀▀▀▀▀▀▀▀ ▐░▌       ▐░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ 
▐░▌          ▐░▌       ▐░▌▐░▌          ▐░▌          ▐░▌          ▐░▌          ▐░▌          
▐░█▄▄▄▄▄▄▄▄▄ ▐░▌       ▐░▌▐░▌          ▐░▌          ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ 
▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░▌          ▐░▌          ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
 ▀▀▀▀▀▀▀▀▀█░▌▐░▌       ▐░▌▐░▌          ▐░▌          ▐░█▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀█░▌ ▀▀▀▀▀▀▀▀▀█░▌
          ▐░▌▐░▌       ▐░▌▐░▌          ▐░▌          ▐░▌                    ▐░▌          ▐░▌
 ▄▄▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄█░▌ ▄▄▄▄▄▄▄▄▄█░▌
▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
 ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀ 
                                                                                           
    """)
    #shell("python ngspipe/scripts/sendmail.py {}".format(receiver_email))
    shell("python ngspipe/scripts/sendmail0129.py -r {} -t {} -d {}".format(receiver_email, "success", join(working_dir, ".snakemake/log/")))
    # NGSPipeDB_source_code/.snakemake/log/


onerror:
    print("An error occurred")
    print("""
>>>>>>>>>>>>>>>>> Workflow finished, no error <<<<<<<<<<<<<<<<<<<<
 ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄ 
▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌
▐░▌          ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌       ▐░▌▐░▌       ▐░▌
▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░▌       ▐░▌▐░█▄▄▄▄▄▄▄█░▌
▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌
▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀█░█▀▀ ▐░█▀▀▀▀█░█▀▀ ▐░▌       ▐░▌▐░█▀▀▀▀█░█▀▀ 
▐░▌          ▐░▌     ▐░▌  ▐░▌     ▐░▌  ▐░▌       ▐░▌▐░▌     ▐░▌  
▐░█▄▄▄▄▄▄▄▄▄ ▐░▌      ▐░▌ ▐░▌      ▐░▌ ▐░█▄▄▄▄▄▄▄█░▌▐░▌      ▐░▌ 
▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌
 ▀▀▀▀▀▀▀▀▀▀▀  ▀         ▀  ▀         ▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀         ▀ 
                                                                 
    """)
    shell("python ngspipe/scripts/sendmail0129.py -r {} -t {} -d {}".format(receiver_email, "error", join(working_dir, ".snakemake/log/")))

include: join("rules", "1.sampling_data_by_{}.Snakefile.py".format(sampling_method))
include: join("rules", "2.rawReads_qc_by_{}.Snakefile.py".format(qc_method))
include: join("rules", "3.genome_align_by_{}.Snakefile.py".format(mapping_method))
include: join("rules", "10.call_variation_by_{}.Snakefile.py".format(call_variation_method))
include: join("rules", "11.variation_annotation_by_{}.Snakefile.py".format(variation_annotation_method))
