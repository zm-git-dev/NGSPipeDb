## nature protocol

```shell
conda create -n exp_r_env mamba -c conda-forge
conda activate exp_r_env
mamba env update --file ngspipe/envs/requirements_exp_r_env.yaml --prune

```

## Step-by-step database generate <a name="Step-by-Step-Database"></a>

We use django project to constructed our NGSDb. We have pareparied many apps for you. 嵌入. Please have a look:

django apps | description | package
---- | ----------- | -------
home | home page   | 
geneExpAtlas | table
network      |
igv          | genome browse | [IGV](https://github.com/igvteam/igv.js/wiki)
blastplus    | ncbi blast +  | NCBI
ngstools     | wooey         |
efp          | efp browse    |
clustergrammer-js https://clustergrammer.readthedocs.io/clustergrammer_js.html

AQUAMIS will provide you with an interactive, browser-based report, showing the most important measures of your data on the first sight. 
All tables in the report can be sorted and filtered.
The table on the first tab shows the key values for a quick estimation of the success of your sequencing experiment and the assembly. 
On the second tab, there is a more detailed table, giving many additional measures.
 Additionally to the tables, many measures are provided as graphical feedback. On the third tab, you see plots which are generated for one complete sequencing experiment. On the fourth tab, there are plots which each show measures on one specific dataset. 

### 1. 安装环境 <a name="DatabaseRequirement"></a>

```shell
django
wooey
clustergrammer
sklearn
pandas=0.25.3
```

### 2. convert analysis result to sqlite3 file <a name="Table2Sqlite3"></a>

modify model.py
```
snakemake -s ngspipe/db_generate.Snakefile.py -p -j1 
```

### 3. config <a name="DatabaseConfig"></a>

1. 修改 `mysite/mysite/settings.py`

```shell
# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    # add custom app
    'geneAnno', # gene annotation from nr/nt/pfam/go/kegg
    'geneExpAtlas', # gene expression matrix
    'blast', # blast tool
]
```

### 4. start server <a name="RunServer"></a>

Starting server by run `python ngsdb/manage.py runserver`, visit sebsite on http://127.0.0.1:8000.

### 5. add your custome data

1. table data in a page
2. custome script in wooey

python manage.py makemigrations
python manage.py migrate
python manage.py addscript ../ngsdb/tools/test1.py

If you want use wooey tools on a task model `celery -A ProjectName worker -c 1 --beat -l info` optional.
