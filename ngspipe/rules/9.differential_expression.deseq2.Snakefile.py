# ----------------------------------------------------------------------------
# input are raw count, replicates must >=3
# ----------------------------------------------------------------------------
rule differential_expression_analysis_by_deseq2:
    message:
        '''
        ------------------------------
        differential_expression_analysis_by_deseq2
        ------------------------------
        '''
    input:
        gene_count = join(quantify_outdir, "gene.csv"),
        condition = config["conditionPath"]
    output:
        diff_outputok = touch(join(diff_outdir, "diff.ok"))
    conda:
        '../envs/requirements_exp_r_env.yaml'
    benchmark:
        join(diff_outdir, "benchmark.txt")
    log:
        join(diff_outdir, "diff.log")
    shell:
        '''
        Rscript ngspipe/scripts/deseq2_call_diff_with_counts_matrix.R --countsMatrix {input.gene_count} --conditionFile {input.condition} --resultDir {diff_outdir} 1>{log} 2>&1;
        '''

rule diff_gene_add_note:
    message:
        '''
        ------------------------------
        add annotation to differential gene
        ------------------------------
        '''
    input:
        diff_outputok = join(diff_outdir, "diff.ok"),
        gtf = config["genomeAnno"]
    output:
        diff_outputok = touch(join(diff_anno_outdir, "diff.ok"))
    conda:
        '../envs/requirements_exp_r_env.yaml'
    benchmark:
        join(diff_anno_outdir, "benchmark.txt")
    log:
        join(diff_anno_outdir, "diff.log")
    shell:
        '''
        python ngspipe/scripts/add_anno2diff_result.py {input.gtf} {diff_outdir} {diff_anno_outdir} 1>{log} 2>&1;
        '''
