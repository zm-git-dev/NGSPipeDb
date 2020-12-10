rule rawReads_data_stat:
    input:
        read1 = join(config["samplesDir"], "{sample}"+config["read1Suffix"]),
        read2 = join(config["samplesDir"], "{sample}"+config["read2Suffix"])
    output:
        read1 = join(rawReads_outdir, "{sample}", "{sample}_1.seqkit.summary"),
        read2 = join(rawReads_outdir, "{sample}", "{sample}_2.seqkit.summary"),
        read1_gc_qual = join(rawReads_outdir, "{sample}", "{sample}_1.seqkit.gc_qual"),
        read2_gc_qual = join(rawReads_outdir, "{sample}", "{sample}_2.seqkit.gc_qual")
    threads: 1
    log:
        join(rawReads_outdir, "{sample}", "run.log"),
    shell:
        '''
        seqkit stat -j {threads} -a {input.read1} >{output.read1} 2>{log};
        seqkit stat -j {threads} -a {input.read2} >{output.read2} 2>>{log};
        seqkit fx2tab -n -g -q -j {threads} {input.read1}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.read1_gc_qual} 2>>{log};
        seqkit fx2tab -n -g -q -j {threads} {input.read2}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.read2_gc_qual} 2>>{log};
        '''

rule rawReads_data_stat_ok:
    input:
        read1 = expand(join(rawReads_outdir, "{sample}", "{sample}_1.seqkit.summary"), sample=SAMPLES)
    output:
        rawReads_ok = touch(join(rawReads_outdir, 'statistic.completed'))
    shell:
        ""