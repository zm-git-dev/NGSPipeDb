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
sample_control = SAMPLES[0]
sample_treated = SAMPLES[1]
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
mapping_method = 'bwa' # bowtie2, bwa, 
mapping_outdir = join(config["resultsDir"], "mapping", "mapping_by_{}".format(mapping_method))
genome_index_prefix = 'genome'

# 4. remove_duplication
rmduplicate_method = 'samtools' # picard or macs2 or samtools, use picard will raise "EOFError: Compressed file ended before the end-of-stream marker was reached", macs2's output is bed format
remove_duplication_outdir = join(config["resultsDir"], "remove_duplication", "rmduplicate_by_{}".format(rmduplicate_method))

# 5. macs2_call_peak
call_peak_method = 'macs2' # macs2
call_peak_outdir = join(config["resultsDir"], "call_peak", "call_peak_by_{}".format(call_peak_method))

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
        # 4. remove duplication #
        rmdup_bam                          = expand(join(remove_duplication_outdir, "{sample}", "{sample}.sorted_rmdup.bam"), sample=SAMPLES),
        #
        # 5. macs2_call_peak #
        #
        callok                             = join(call_peak_outdir, "callpeak.ok"),
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
include: join("rules", "remove_duplicate_by_{}.Snakefile.py".format(rmduplicate_method))
include: join("rules", "call_peak_by_{}.Snakefile.py".format(call_peak_method))
