import os
import sys
import pandas as pd
from os.path import join

# ----------------------------------------------------------------------------
# sample to use
smpList = pd.read_csv(config["samplesList"], index_col=0, header=None)
SAMPLES = list(smpList.index)
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# outputdir
#
sampling_data_outdir = join(config["resultsDir"], "sampling_data")
# ----------------------------------------------------------------------------

rule all:
    input:
        sampling_data_by_head_result        = expand(join(sampling_data_outdir, "sampling_data_by_head", "{sample}"+config["read1Suffix"]), sample=SAMPLES),
        

onsuccess:
    print("Workflow finished, no error")
    #shell("mail -s 'Workflow finished, no error' 296373256@qq.com")

onerror:
    print("An error occurred")
    #shell("mail -s 'an error occurred' 296373256@qq.com ")

include: join("modules", "sampling_data.seqkit.Snakefile.py")
#include: join("modules", "align.hisat2_stringtie.Snakefile.py")
#include: join("modules", "stat.rseqc.Snakefile.py")
#include: join("modules", "quant.HtSeqCount.Snakefile.py")
#include: join("modules", "quant.stringtie.Snakefile.py")
#include: join("modules", "quant.featureCounts.Snakefile.py")
#include: join("modules", "report.markdown.Snakefile.py")
#include: join("modules", "qc.trim_galore.Snakefile.py")