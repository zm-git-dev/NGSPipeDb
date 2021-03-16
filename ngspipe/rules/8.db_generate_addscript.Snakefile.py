rule add_script2wooey:
    message:
        '''
        ------------------------------
        8. add_script2wooey
        ------------------------------
        '''
    input:
        makemigration                   = join(migration_outdir, "migration.ok"),
    output:
        addscript                       = touch(join(addscript_outdir, "addscript.ok")),
    log:
        join(addscript_outdir, "run.log")
    benchmark:
        join(addscript_outdir, "benchmark.txt")
    params:
        dna2rna_path = join('ngsdb', 'wooey', 'wooey_scripts', 'dna2rna.py'),
        dna2rna_name = 'dna2rna_convert',
    shell:
        '''
        # remove scripts in wooey database
        python ngspipe/scripts/remove_script_in_wooey_db.py 1>{log} 2>&1
        # add scripts in current machine
        python ngsdb/manage.py addscript {params.dna2rna_path} --name {params.dna2rna_name} 1>>{log} 2>&1;
        '''