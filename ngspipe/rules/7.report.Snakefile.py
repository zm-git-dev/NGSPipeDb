# global description
report: join(snake_dir, "reports/workflow.rst")

"""
rule dag:
    input:
    output:
        "dag.svg"
    shell:
        '''
        snakemake --snakefile NGSPipeCode/Snakefile --configfile NGSPipeCode/config.yaml --dag | dot -Tsvg > dag.svg
        '''

rule pipeline:
    input:
        'dag.svg'
    output:
        report(join(config['reportsDir'], 'report', "fig1.svg"), caption="reports/fig1.rst", category="Step 1: 设计自己的pipeline")
    shell:
        '''
        cp {input} {output};
        '''
"""

rule report_all:
    message:
        '''
        ------------------------------
        merge report result
        ------------------------------
        '''
    input:
        workflow_report = join(config['reportsDir'], 'workflow', "workflow.png"),
        rawdata_report = join(config['reportsDir'], '1.rawreads_stat', "rawreads_product.csv"),
        cleandata_report = join(config['reportsDir'], '2.cleanreads_stat', "cleanreads_product.csv"),
        bam_report = expand(join(config['reportsDir'], '3.mapping_stat', "{sample}.bam_stat.txt"), sample=SAMPLES),
        gene_boxplot_pdf = join(config['reportsDir'], '4.exp_stat', "gene_fpkm_boxplot.pdf"),
    output:
        touch(join(config['reportsDir'], "report.ok"))
    shell:
        '''
        '''

rule workflow_report:
    message:
        '''
        ------------------------------
        workflow_report
        ------------------------------
        '''
    input:
        workflow_report = join(snake_dir, "imgs", "workflow.png")
    output:
        workflow_report = join(config['reportsDir'], 'workflow', "workflow.png")
    shell:
        '''
        cp {input.workflow_report} {output.workflow_report};
        '''

rule rawreads_stat_report:
    message:
        '''
        ------------------------------
        rawreads_stat_report
        ------------------------------
        '''
    input:
        rawdata_report = join(rawReads_outdir, "reads_product.csv")
    output:
        rawdata_report = report(join(config['reportsDir'], '1.rawreads_stat', "rawreads_product.csv"), caption=join(snake_dir, "reports/rawreads_stat.rst"), category="Step 1: rawreads")
    shell:
        '''
        cp {input.rawdata_report} {output.rawdata_report};
        '''


rule cleanreads_stat_report:
    message:
        '''
        ------------------------------
        cleanreads_stat_report
        ------------------------------
        '''
    input:
        cleandata_report = join(cleanReads_outdir, "reads_product.csv")
    output:
        cleandata_report = report(join(config['reportsDir'], '2.cleanreads_stat', "cleanreads_product.csv"), caption=join(snake_dir, "reports/cleanreads_stat.rst"), category="Step 2: clean reads")
    shell:
        '''
        cp {input.cleandata_report} {output.cleandata_report};
        '''

rule cleanreads_mapping_report:
    message:
        '''
        ------------------------------
        cleanreads_mapping_report
        ------------------------------
        '''
    input:
        bam_stat = join(bam_outdir, "{sample}", "{sample}.bam_stat.txt")
    output:
        bam_report = report(join(config['reportsDir'], '3.mapping_stat', "{sample}.bam_stat.txt"), caption=join(snake_dir, "reports/cleanreads_stat.rst"), category="Step 3: mapping")
    shell:
        '''
        cp {input.bam_stat} {output.bam_report};
        '''

rule expression_report:
    message:
        '''
        ------------------------------
        7. expression_report
        ------------------------------
        '''
    input:
        gene_fpkm = join(quantify_outdir, "gene_fpkm_all_samples.tsv")
    output:
        gene_boxplot_pdf = report(join(config['reportsDir'], '4.exp_stat', "gene_fpkm_boxplot.pdf"), caption=join(snake_dir, "reports/expression.rst"), category="Step 4: show table")
    shell:
        '''
        python {snake_dir}/scripts/gene_exp_boxplot.py {input.gene_fpkm} {output.gene_boxplot_pdf};
        '''


"""
rule assembled_transcript_report:
    message:
        '''
        ------------------------------
        assembled_transcript_report
        ------------------------------
        '''
    input:
        'NGSPipeCode/reports/fig2.png'
    output:
        report(join(config['reportsDir'], 'report', "fig2.png"), caption=join(snake_dir, "reports/fig2.rst"), category="Step 2: 统计")
    shell:
        "cp {input} {output}"


"""