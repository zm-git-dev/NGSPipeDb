rule migration_model:
    message:
        '''
        ------------------------------
        8. after model generated, make migration
        ------------------------------
        '''
    input:
        gffdjango_model = join(config["djangoCode"], "geneAnno", "models.py"),
        exp_django_model = join(config["djangoCode"], "geneExpAtlas", "models.py"),
    output:
        migrationok = touch(join(migration_outdir, "migration.ok")),
    params:
        db_name = "gffDb",
        db_app = 'geneAnno',
    log:
        join(migration_outdir, "change_model_py.log")
    shell:
        '''
        pip install wooey clustergrammer sklearn pandas==0.25.3 1>{log} 2>&1;
        python {config[djangoCode]}/manage.py makemigrations 1>>{log} 2>&1;
        # django makemigrate
        python {config[djangoCode]}/manage.py migrate 1>>{log} 2>&1;
        '''
