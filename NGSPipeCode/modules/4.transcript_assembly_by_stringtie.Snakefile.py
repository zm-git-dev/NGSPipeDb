# ----------------------------------------------------------------------------
# assembly each sample
# ----------------------------------------------------------------------------
rule transcript_assembly_by_stringtie:
    message:
        '''
        stringtie
        '''
    input:
        sorted_bam = join(junction_align_outdir, "{sample}", "{sample}.sorted.bam"),
        genomeAnno = config["genomeAnno"]
    output:
        gtf = join(transcript_assembly_outdir, "{sample}", "{sample}.gtf")
    benchmark:
        join(transcript_assembly_outdir, "{sample}", "benchmark.txt")
    log:
        join(transcript_assembly_outdir, "{sample}", "stringtie.log")
    threads: 5
    shell:
        '''
        stringtie -p {threads} -G {input.genomeAnno} -o {output.gtf} {input.sorted_bam} 1>{log} 2>&1;
        '''

# ----------------------------------------------------------------------------
# merged_gtf
# ----------------------------------------------------------------------------
rule merge_stringtieResult_by_stringtieMerge:
    message:
        '''
        stringtieMerge
        '''
    input: 
        gtfs = expand(join(transcript_assembly_outdir, "{sample}", "{sample}.gtf"), sample = SAMPLES),
        genomeAnno = config["genomeAnno"]
    output:
        mergedGtf = join(transcript_assembly_outdir, "merged.gtf")
    benchmark:
        join(transcript_assembly_outdir, "benchmark.txt")
    log:
        join(transcript_assembly_outdir, "merged.log")
    shell:
        '''
        stringtie --merge -G {input.genomeAnno} -o {output.mergedGtf} {input.gtfs} 1>{log} 2>&1;
        '''
