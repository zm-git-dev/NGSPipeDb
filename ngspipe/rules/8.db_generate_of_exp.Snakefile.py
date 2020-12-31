rule sqlite3_auto_by_expression_data:
    message:
        '''
        ------------------------------
        8. generate expression data sqlite3
        ------------------------------
        '''
    input:
        gene_fpkm = join(quantify_outdir, "gene_fpkm_all_samples.tsv")
    output:
        exp_sqlite3 = join(exp_db_outdir, "exp.sqlite3")
    log:
        join(exp_db_outdir, "run.log")
    benchmark:
        join(exp_db_outdir, "benchmark.txt")
    shell:
        '''
        pip install simplesqlite 1>{log} 2>&1;
        python {snake_dir}/scripts/exp2sqlite3.py {input.gene_fpkm} {output.exp_sqlite3} 1>>{log} 2>&1;
        '''

rule sqlite3ForModel_auto_by_django:
    message:
        '''
        ------------------------------
        8. generate expression data sqlite3
        ------------------------------
        '''
    input:
        exp_sqlite3 = join(exp_db_outdir, "exp.sqlite3"),
    output:
        exp_django_model = join(config["djangoCode"], "geneExpAtlas", "models.py"),
    params:
        db_name = "expDb",
    log:
        join(exp_db_outdir, "change_model_py.log")
    shell:
        '''
        # inspectdb数据模型反向生成
        # make sure the right db name in ngsdb/ngsdb/setting.py
        python {config[djangoCode]}/manage.py inspectdb --database expDb|perl -ne 'if(/primary_key/){{s/null=True/null=False/}};print "$_"' > {output.exp_django_model} 2>{log};
        # 使用makemigrations创建迁移
        python {config[djangoCode]}/manage.py makemigrations geneExpAtlas 1>{log} 2>&1;
        # 使用migrate执行迁移
        python {config[djangoCode]}/manage.py migrate 1>{log} 2>&1;
        '''

rule expression_blast:
    input:
    output:
    shell:
        '''
        '''

rule expression_heatmap:
    input:
    output:
    shell:
        '''
        '''