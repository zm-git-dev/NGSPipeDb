report: "reports/workflow.rst"

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

rule b:
    input:
        'reports/fig2.png'
    output:
        report(join(config['reportsDir'], 'report', "fig2.png"), caption="reports/fig2.rst", category="Step 2: 统计")
    shell:
        "cp {input} {output}"


rule d:
    input:
        'reports/test.csv'
    output:
        report(join(config['reportsDir'], 'report', "test.csv"), caption="reports/table.rst", category="Step 3: show table")
    shell:
        "cp {input} {output}"
