rule trinity:
    input:
        unmappedfq1 = expand(config["resultfolder"] + "/unmapped/{sample}/{sample}_1.fq", sample=SAMPLES),
        unmappedfq2 = expand(config["resultfolder"] + "/unmapped/{sample}/{sample}_2.fq", sample=SAMPLES),
    output:
        trinitysample = config["resultfolder"] + "/trinity/trinitysample.txt",
        trinityout = touch(config["resultfolder"] + "/trinity/trinity_out.ok")
    params:
        unmappfolder = config["resultfolder"] + "/unmapped/",
        trinityfolder = config["resultfolder"] + "/trinity/trinity_out"
    conda: "env/trinity.yaml"
    threads: 40
    log: trinitysample = config["resultfolder"] + "/trinity/run_trinity.log",
    shell:
        "python script/sample4trinity.py {params.unmappfolder} >{output.trinitysample};"
        "Trinity --seqType fq --max_memory 110G --output {params.trinityfolder} --samples_file {output.trinitysample} --CPU {threads} 1>{log} 2>&1;"
        
rule QCsamples:
    input:
        gene_count = config["resultfolder"] + "/assembly_final/gene.counts.norm.tsv",
    output:
        ok = touch(config["resultfolder"] + "/sampleQC/sampleQC.ok")
    params:
        countFolder = config["resultfolder"] + "/sampleQC"
    conda: "env/trinity.yaml"
    log: config["resultfolder"] + "/sampleQC/QCsamples.log"
    shell:
        "cd {params.countFolder};"
        # generate a barplot that sums frag counts per replicate across all samples.
        "PtR --matrix ../../result/assembly_final/gene.tsv --samples ../../samples.txt --log2 --barplot_sum_counts >../../{log} 2>&1;"
        # generate a sample correlation matrix plot
        "PtR --matrix ../../result/assembly_final/gene.tsv --samples ../../samples.txt --log2 --CPM --sample_cor_matrix >../../{log}1 2>&1;"
        # PCA
        "PtR --matrix ../../result/assembly_final/gene.tsv --samples ../../samples.txt --log2 --CPM --prin_comp 2 >../../{log}2 2>&1"

rule DEanalysis:
    input:
        gene_count = config["resultfolder"] + "/assembly_final/gene.tsv",
    output:
        ok = touch(config["resultfolder"] + "/diff/diff.ok")
    conda: "env/trinity.yaml"
    log: config["resultfolder"] + "/diff/diff.log"
    params:
        diffFolder = config["resultfolder"] + "/diff"
    shell:
        "run_DE_analysis.pl --matrix {input.gene_count} --method DESeq2 --samples_file samples.txt --output {params.diffFolder} >{log} 2>&1;"
        # Extracting and clustering differentially expressed transcripts
        "cd {params.diffFolder};"
        "analyze_diff_expr.pl --matrix ../../result/assembly_final/gene.tsv --samples ../../samples.txt >../../{log}1 2>&1;"
        # Automatically Partitioning Genes into Expression Clusters
        "define_clusters_by_cutting_tree.pl -R diffExpr.P0.001_C2.matrix.RData --Ptree 60 >../../{log}2 2>&1;"

rule unmappedReads:
    input:
        config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        unmappedbam = temp(config["resultfolder"] + "/unmapped/{sample}/{sample}.unmapped.bam"),
        unmappedfq1 = config["resultfolder"] + "/unmapped/{sample}/{sample}_1.fq",
        unmappedfq2 = config["resultfolder"] + "/unmapped/{sample}/{sample}_2.fq"
    log:
        unmappedbamlog = config["resultfolder"] + "/unmapped/{sample}/{sample}.unmappedbamlog.log",
        bedtools = config["resultfolder"] + "/unmapped/{sample}/{sample}.bedtools.log"
    params:

    conda: "env/mapping.yaml"
    threads: 5
    shell:
        "samtools view -b -f 4 -@ {threads} {input}|samtools sort -@ {threads} -n - 2>{log.unmappedbamlog}|samtools fixmate - {output.unmappedbam} 1>>{log.unmappedbamlog} 2>&1;"
        "bedtools bamtofastq -i {output.unmappedbam} -fq {output.unmappedfq1} -fq2 {output.unmappedfq2} 1>{log.bedtools} 2>&1;"

rule trinity:
    input:
        unmappedfq1 = expand(config["resultfolder"] + "/unmapped/{sample}/{sample}_1.fq", sample=SAMPLES),
        unmappedfq2 = expand(config["resultfolder"] + "/unmapped/{sample}/{sample}_2.fq", sample=SAMPLES),
    output:
        trinitysample = config["resultfolder"] + "/trinity/trinitysample.txt",
        trinityout = touch(config["resultfolder"] + "/trinity/trinity_out.ok")
    params:
        unmappfolder = config["resultfolder"] + "/unmapped/",
        trinityfolder = config["resultfolder"] + "/trinity/trinity_out"
    conda: "env/trinity.yaml"
    threads: 40
    log: trinitysample = config["resultfolder"] + "/trinity/run_trinity.log",
    shell:
        "python script/sample4trinity.py {params.unmappfolder} >{output.trinitysample};"
        "Trinity --seqType fq --max_memory 110G --output {params.trinityfolder} --samples_file {output.trinitysample} --CPU {threads} 1>{log} 2>&1;"

rule unigenePlusGenome:
    input: 
        "unigene.fa",
        "genome.fa",
    output:
        "Ref.fa"
    shell:
        "cat genome.fa unigene.fa >Ref.fa"

rule refIndex:
    input: "Ref.fa"
    output:
        touch(config["resultfolder"] + "/refIndex/index.ok")
    conda: "env/mapping.yaml"
    log: config["resultfolder"] + "/refIndex/index.log"
    params:
        outfolder = config["resultfolder"] + "/genomeIndex"
    shell:
        "hisat2-build {input} {params.outfolder}/ref 1>{log} 2>&1"