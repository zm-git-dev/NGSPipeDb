import os
from os.path import join
import sys
import pandas as pd

# relative path
snake_dir = workflow.basedir # all configfile, scripts, restructuretext, ens are relative to snakefile (this file)
working_dir = os.getcwd() # input and output path are relative to current working directory

# configfile
configfile: join(snake_dir, "config", "rnaseq.config.yaml")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# detail parameters in pipe #

# 5. quantification
quantify_method = 'stringtie' # htseqcounts or featurecounts
quantify_outdir = join(config["resultsDir"], "quantify", "quantify_by_{}".format(quantify_method))

# 4. transcript assembly
transcript_assembly_method = 'stringtie' # star
transcript_assembly_outdir = join(config["resultsDir"], "transcript_assembly", "transcript_assembly_by_{}".format(transcript_assembly_method))


# blast
# 1. expression matrix database create
exp_db_outdir = join(config["dbDir"], "exp")
anno_db_outdir = join(config["dbDir"], "anno")

# sqlite3
# 2. blastdb
blastdb_outdir = join(config["dbDir"], "blastdb")

# 3. gffutils
gffdb_outdir = join(config["dbDir"], "gff_sqlite3")

# 4. genomebrowse
gbrowse_outdir = join(config["dbDir"], "gbrowse")

rule all:
    input:
        # 1. exp database create #
        sqlite3_exp                     = join(exp_db_outdir, "exp.sqlite3"),
        exp_django_model        = join(config["djangoCode"], "geneExpAtlas", "models.py"),

        # 2. makeblastdb
        merged_db               = join(blastdb_outdir, "merged_blastdb.ok"),

        # 3. gff database create
        gffdb                   = join(gffdb_outdir, "gtf.sqlite3"),
        gffdjango_model         = join(config["djangoCode"], "geneAnno", "models.py"),

        # 4. genome browse
        merged_gbrowse = join(gbrowse_outdir, 'merged_gbrowse.ok'),

        

onsuccess:
    print("""
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Workflow finished, no error <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 ▄▄▄▄▄▄▄▄▄▄▄  ▄         ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄ 
▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
▐░█▀▀▀▀▀▀▀▀▀ ▐░▌       ▐░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ 
▐░▌          ▐░▌       ▐░▌▐░▌          ▐░▌          ▐░▌          ▐░▌          ▐░▌          
▐░█▄▄▄▄▄▄▄▄▄ ▐░▌       ▐░▌▐░▌          ▐░▌          ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ 
▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░▌          ▐░▌          ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
 ▀▀▀▀▀▀▀▀▀█░▌▐░▌       ▐░▌▐░▌          ▐░▌          ▐░█▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀█░▌ ▀▀▀▀▀▀▀▀▀█░▌
          ▐░▌▐░▌       ▐░▌▐░▌          ▐░▌          ▐░▌                    ▐░▌          ▐░▌
 ▄▄▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄█░▌ ▄▄▄▄▄▄▄▄▄█░▌
▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
 ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀ 

    """)
    #shell("python NGSPipe/scripts/sendmail.py {}".format(receiver_email))
    # NGSPipeDB_source_code/.snakemake/log/
    
    

onerror:
    print("An error occurred")
    #shell("mail -s 'an error occurred' 296373256@qq.com ")

include: join("rules", "8.db_generate_of_exp.Snakefile.py")
include: join("rules", "8.db_generate_of_gff.Snakefile.py")
include: join("rules", "8.db_generate_of_blastdb.Snakefile.py")
include: join("rules", "8.db_generate_of_genomebrowser.Snakefile.py")
