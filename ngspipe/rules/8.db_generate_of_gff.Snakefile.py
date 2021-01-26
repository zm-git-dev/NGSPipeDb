
rule gtfTosqlite3:
    message:
        '''
        ------------------------------
        8. generate expression data sqlite3
        ------------------------------
        '''
    input:
        transcript_assembly = config['genomeAnno']
    output:
        gffdb = join(gffdb_outdir, "gtf.sqlite3")
    log:
        join(gffdb_outdir, 'run.log')
    conda:
        '../envs/gffutils.yaml'
    shell:
        '''
        # 这个数据库中的表是否没有index?
        gffutils-cli create --output {output.gffdb} {input.transcript_assembly} >{log} 2>&1;
        '''

rule gtfTosqlite3ForModel_auto_by_django:
    message:
        '''
        ------------------------------
        8. generate expression data sqlite3
        ------------------------------
        '''
    input:
        gffdb = join(gffdb_outdir, "gtf.sqlite3")
    output:
        gffdjango_model = join(config["djangoCode"], "geneAnno", "models.py"),
    params:
        db_name = "gffDb",
        db_app = 'geneAnno',
    log:
        join(gffdb_outdir, "change_model_py.log")
    shell:
        '''
        # inspectdb数据模型反向生成
        # make sure the right db name in ngsdb/ngsdb/setting.py
        python {config[djangoCode]}/manage.py inspectdb --database {params.db_name}|perl -ne 'if(/\s+id = /){{s/null=True/null=False, primary_key=True/}}print $_' > {output.gffdjango_model} 2>{log};
        # 使用makemigrations创建迁移
        python {config[djangoCode]}/manage.py makemigrations {params.db_app} 1>>{log} 2>&1;
        # 使用migrate执行迁移
        python {config[djangoCode]}/manage.py migrate 1>>{log} 2>&1;
        '''
