annotation_gbrowse_outdir = join(gbrowse_outdir, 'annotation')

rule merge_gbrowse:
    message:
        '''
        ------------------------------
        8. merge all gbrowse database
        ------------------------------
        '''
    input:
        gtfgzip = join(annotation_gbrowse_outdir, 'trasncript.gtf.gz')
    output:
        merged_gbrowse = join(gbrowse_outdir, 'merged_gbrowse.ok')
    shell:
        '''
        '''

rule igv_annotation_gtf:
    message:
        '''
        ------------------------------
        8. merge all blast database
        ------------------------------
        '''
    input:
        refgtf = config['genomeAnno'],
        transcript_assembly = join(transcript_assembly_outdir, "merged.gff"),
    output:
        gtfgzip = join(annotation_gbrowse_outdir, 'trasncript.gff.gz')
    benchmark:
        join(annotation_gbrowse_outdir, "benchmark.txt")
    log:
        join(annotation_gbrowse_outdir, "run.log")
    conda:
        '../envs/agat.yaml'
    shell:
        '''
        gzip -c {input.transcript_assembly} > {output.gtfgzip} 2>{log};
        bgzip -c merged.sorted.gff >merged.sorted.gff.bgzip
        bedtools sort -i merged.gff merged.sorted.gff
        tabix -p gff merged.sorted.gff.bgzip;
        tabix -p gff trasncript.gtf.bgzip
        bgzip -c merged.sorted.gff >merged.sorted.gff.bgzip
        '''

"""
rule igv_refrence_fasta:
    input:
    output:
    shell:
        '''
        '''

rule igv_mapping_bam:
    input:
    output:
    shell:
        '''
        '''



rule genomebrowse:
    input:
        genome_gtf = join(config['result_dir'], 'database', 'data', 'genome.gff3'),
        genome_fasta = config['genome_fasta'],
        bigGenePred = config['bigGenePred']
    output:
        genome_bed = join(config['result_dir'], 'database', 'data', 'genome.bed'),
        genome_sortedbed = join(config['result_dir'], 'database', 'data', 'genome.sortedbed'),
        genome_bigbed = join(config['result_dir'], 'database', 'data', 'genome.sorted.bb'),
        genome_bigwig = join(config['result_dir'], 'database', 'data', 'genome.2bit'),
        genome_size = join(config['result_dir'], 'database', 'data', 'genome.size'),
        genome_pred = join(config['result_dir'], 'database', 'data', 'genome.pred'),
    conda: 'envs/genomebrowse.yaml'
    log: join(config['result_dir'], 'database', 'data', 'genomebrowse.log'),
    shell:
        '''
        faSize -detailed {input.genome_fasta} >{output.genome_size} 2>{log};
        gff3ToGenePred {input.genome_gtf} {output.genome_pred} 1>>{log} 2>&1;
        genePredToBed {output.genome_pred} {output.genome_bed} 1>>{log} 2>&1;
        bedSort {output.genome_bed} {output.genome_sortedbed} 1>>{log} 2>&1;
        #python bed12anno.py Function_annotation.xls mRNA_old2new.csv PVDB.genome.sorted.bed >PVDB.genome.sorted.bed15 2>>{log};
        #bedToBigBed -type=bed12+3 -as=bigGenePred.as {output.genome_sortedbed} {output.genome_size} {output.genome_bigbed} 1>>{log} 2>&1;
        bedToBigBed -type=bed12 -as={input.bigGenePred} {output.genome_sortedbed} {output.genome_size} {output.genome_bigbed} 1>>{log} 2>&1;
        faToTwoBit {input.genome_fasta} {output.genome_bigwig} 1>>{log} 2>&1;
        #bam2wig
        '''

# ----------------------------------------------------------------------------
#step2 mapping
step2_outdir = join(config["result_dir"], "align.hisat2_stringtie", \
    "step2_mapping_genome_by_hisat2")
# ----------------------------------------------------------------------------
rule bam2bigwig:
    input:
        bam = join(step2_outdir, "{sample}", "{sample}.sorted.bam"),
    output:
        sorted_bam = temp(join(config['result_dir'], 'database', 'sorted_bam', '{sample}.sorted.bam')),
        bigwig = join(config['result_dir'], 'database', 'bigwig', '{sample}.bigWig'),
    conda: 'envs/bam2bigwig.yaml'
    log: join(config['result_dir'], 'database', 'bigwig', '{sample}.bigWig.log'),
    shell:
        '''
        samtools sort {input.bam} -o {output.sorted_bam} 1>{log} 2>&1;
        bam2wig.sh {output.sorted_bam} > {output.bigwig} 2>>{log};
        # bedtools genomecov -ibam MyChIP.bam -bg > MyChIP.bedGraph;
        # bedGraphToBigWig MyChIP.bedGraph mm10.chrom.sizes MyChIP.bigWig;
        '''

rule genomebrowse:
    input:
        genome_gtf = join(config['result_dir'], 'database', 'data', 'genome.gff3'),
        genome_fasta = config['genome_fasta'],
        bigGenePred = config['bigGenePred']
    output:
        genome_bed = join(config['result_dir'], 'database', 'data', 'genome.bed'),
        genome_sortedbed = join(config['result_dir'], 'database', 'data', 'genome.sortedbed'),
        genome_bigbed = join(config['result_dir'], 'database', 'data', 'genome.sorted.bb'),
        genome_bigwig = join(config['result_dir'], 'database', 'data', 'genome.2bit'),
        genome_size = join(config['result_dir'], 'database', 'data', 'genome.size'),
        genome_pred = join(config['result_dir'], 'database', 'data', 'genome.pred'),
    conda: 'envs/genomebrowse.yaml'
    log: join(config['result_dir'], 'database', 'data', 'genomebrowse.log'),
    shell:
        '''
        faSize -detailed {input.genome_fasta} >{output.genome_size} 2>{log};
        gff3ToGenePred {input.genome_gtf} {output.genome_pred} 1>>{log} 2>&1;
        genePredToBed {output.genome_pred} {output.genome_bed} 1>>{log} 2>&1;
        bedSort {output.genome_bed} {output.genome_sortedbed} 1>>{log} 2>&1;
        #python bed12anno.py Function_annotation.xls mRNA_old2new.csv PVDB.genome.sorted.bed >PVDB.genome.sorted.bed15 2>>{log};
        #bedToBigBed -type=bed12+3 -as=bigGenePred.as {output.genome_sortedbed} {output.genome_size} {output.genome_bigbed} 1>>{log} 2>&1;
        bedToBigBed -type=bed12 -as={input.bigGenePred} {output.genome_sortedbed} {output.genome_size} {output.genome_bigbed} 1>>{log} 2>&1;
        faToTwoBit {input.genome_fasta} {output.genome_bigwig} 1>>{log} 2>&1;
        #bam2wig
        '''

# ----------------------------------------------------------------------------
#step2 mapping
step2_outdir = join(config["result_dir"], "align.hisat2_stringtie", \
    "step2_mapping_genome_by_hisat2")
# ----------------------------------------------------------------------------
rule bam2bigwig:
    input:
        bam = join(step2_outdir, "{sample}", "{sample}.sorted.bam"),
    output:
        sorted_bam = temp(join(config['result_dir'], 'database', 'sorted_bam', '{sample}.sorted.bam')),
        bigwig = join(config['result_dir'], 'database', 'bigwig', '{sample}.bigWig'),
    conda: 'envs/bam2bigwig.yaml'
    log: join(config['result_dir'], 'database', 'bigwig', '{sample}.bigWig.log'),
    shell:
        '''
        samtools sort {input.bam} -o {output.sorted_bam} 1>{log} 2>&1;
        bam2wig.sh {output.sorted_bam} > {output.bigwig} 2>>{log};
        # bedtools genomecov -ibam MyChIP.bam -bg > MyChIP.bedGraph;
        # bedGraphToBigWig MyChIP.bedGraph mm10.chrom.sizes MyChIP.bigWig;
        '''
"""