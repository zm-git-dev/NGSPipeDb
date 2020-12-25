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
#

rule all:
    input:
        # 1. database create #
        sqlite3_exp                     = join(exp_db_outdir, "exp.sqlite3"),
        exp_django_model = join(config["djangoCode"], "geneExpAtlas", "models.py"),
        

onsuccess:
    print("""
    Workflow finished, no error
                ............                zhangxuan@T640P 
         .';;;;;.       .,;,.            --------------- 
      .,;;;;;;;.       ';;;;;;;.         OS: Deepin 20 x86_64 
    .;::::::::'     .,::;;,''''',.       Host: PowerEdge T640 
   ,'.::::::::    .;;'.          ';      Kernel: 5.4.50-amd64-desktop 
  ;'  'cccccc,   ,' :: '..        .:     Uptime: 1 hour, 47 mins 
 ,,    :ccccc.  ;: .c, '' :.       ,;    Packages: 2097 (dpkg) 
.l.     cllll' ., .lc  :; .l'       l.   Shell: bash 5.0.3 
.c       :lllc  ;cl:  .l' .ll.      :'   Resolution: 1920x1058 
.l        'looc. .   ,o:  'oo'      c,   WM: _NET_SUPPORTING_WM_CHECK: window id # 0x400001 
.o.         .:ool::coc'  .ooo'      o.   Icons: bloom [GTK2/3] 
 ::            .....   .;dddo      ;c    Terminal: /dev/pts/0 
  l:...            .';lddddo.     ,o     CPU: Intel Xeon Gold 5218R (80) @ 803MHz 
   lxxxxxdoolllodxxxxxxxxxc      :l      GPU: NVIDIA Quadro P620 
    ,dxxxxxxxxxxxxxxxxxxl.     'o,       Memory: 2601MiB / 128539MiB 
      ,dkkkkkkkkkkkkko;.    .;o;
        .;okkkkkdl;.    .,cl:.                                   
            .,:cccccccc:,.

    """)
    #shell("python NGSPipeCode/script/sendmail.py {}".format(receiver_email))
    # NGSPipeDB_source_code/.snakemake/log/
    
    

onerror:
    print("An error occurred")
    #shell("mail -s 'an error occurred' 296373256@qq.com ")

include: join("rules", "8.db_generate_of_exp.Snakefile.py")
