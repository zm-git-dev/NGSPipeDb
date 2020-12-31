genomedb_outdir = join(blastdb_outdir, 'nucl_genomedb')
transcriptdb_outdir = join(blastdb_outdir, 'nucl_transcriptdb')
proteindb_outdir = join(blastdb_outdir, 'prot_proteindb')

rule merge_blastdb:
    message:
        '''
        ------------------------------
        8. merge all blast database
        ------------------------------
        '''
    input:
        genomeFastaDb = touch(join(genomedb_outdir, "genomefasta_blastdb.ok"))
    output:
        merged_db = touch(join(blastdb_outdir, "merged_blastdb.ok"))

rule makeblastdb_genomefasta_by_blast:
    message:
        '''
        ------------------------------
        8. generate genome blast database
        ------------------------------
        '''
    input:
        genomeFasta = config['genomeFasta']
    output:
        genomeFastaDb = touch(join(genomedb_outdir, "genomefasta_blastdb.ok"))
    log:
        join(genomedb_outdir, "run.log")
    benchmark:
        join(genomedb_outdir, "benchmark.txt")
    conda:
        '../envs/blast.yaml'
    params:
        genomedb_prefix = 'genome'
    shell:
        '''
        makeblastdb -in {input.genomeFasta} -dbtype nucl -out {genomedb_outdir}/{params.genomedb_prefix} 1>{log} 2>&1;
        '''

"""
rule makeblastdb_transcriptfasta_by_blast:
    message:
        '''
        ------------------------------
        8. generate genome blast database
        ------------------------------
        '''
    input:
        genomeFasta = config['genomeFasta']
    output:
        genomeFastaDb = touch(join(blastdb_outdir, "nucl", "transcriptfasta_blastdb.ok"))
    log:
        join(blastdb_outdir, "run.log")
    benchmark:
        join(blastdb_outdir, "benchmark.txt")
    conda:
        '../envs/blast.yaml'
    params:
        genomedb_prefix = 'genome'
    shell:
        '''
        makeblastdb -in {input.genomeFasta} -dbtype nucl -o {blastdb_outdir}/{params.genomedb_prefix} 1>{log} 2>&1;
        '''

rule makeblastdb_proteinfasta_by_blast:
    message:
        '''
        ------------------------------
        8. generate genome blast database
        ------------------------------
        '''
    input:
        genomeFasta = config['genomeFasta']
    output:
        genomeFastaDb = touch(join(blastdb_outdir, "prot", "proteinfasta_blastdb.ok"))
    log:
        join(blastdb_outdir, "run.log")
    benchmark:
        join(blastdb_outdir, "benchmark.txt")
    conda:
        '../envs/blast.yaml'
    params:
        genomedb_prefix = 'genome'
    shell:
        '''
        makeblastdb -in {input.genomeFasta} -dbtype nucl -o {blastdb_outdir}/{params.genomedb_prefix} 1>{log} 2>&1;
        '''



rule blastdb:
    input:
        genome_gtf = config['genome_gtf'],
        genome_fasta = config['genome_fasta'],
        #genome_protein = 'if you have, give a path',
        #genome_cds = 'if you have, give a path',
    output:
        blastdb = directory(join(config['result_dir'], 'database', 'blastdb'))
    conda: 'envs/blast.yaml'
    log: join(config['result_dir'], 'database', 'blastdb', 'makeblasdb.log')
    shell:
        '''
        gffread -g {input.genome_fasta} -G {input.genome_gtf} -x {output}/cds.fa 1>{log} 2>&1;
        gffread -g {input.genome_fasta} -G {input.genome_gtf} -y {output}/protein.fa 1>>{log} 2>&1;
        makeblastdb -in {input.genome_fasta} -dbtype nucl -out {output}/genome 1>>{log} 2>&1;
        makeblastdb -in {output}/protein.fa -dbtype prot -out {output}/protein 1>>{log} 2>&1;
        makeblastdb -in {output}/cds.fa -dbtype nucl -out {output}/cds 1>>{log} 2>&1;
        '''

"""