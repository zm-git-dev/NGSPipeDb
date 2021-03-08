
# NGSDb

__Table of Contents:__

1. [Pre-prepare](#pre)
3. [Install the NGSDb conda environments](#NGSPipeDbEnv)
4. [Downloading the NGSDb source code](#NGSPipeDbSource)
5. [Downloading the NGSPipeDb test files](#Testdata)
6. [Run test data](#RunTest)
7. [Generate report](#Report)
8. [Run your custome data](#RunRawdata)


## one step to view the NGSPipe RNA-Seq analysis results in a browser

```shell
bash ngspipe/scripts/one_step_view_database.sh
```
test
Now you can viste your website on http://127.0.0.1:8000. All result are stored in `results`.
- Example of report <sub>[![html](https://img.icons8.com/ios/20/000000/html-filetype.png)](http://www.liu-lab.com)</sub>.
- Example of database <sub>[![html](https://img.icons8.com/dotty/25/000000/copy-link.png)](http://www.liu-lab.com)</sub>.

## Step-by-step to generate database with your own data <a name="Step-by-Step-Database"></a>

### 1. pre-prepare

Let's assume that you have the data you need to generate the database, such as gene annotation files, gene expression matrix and analysis reports. Before further steps, please [install CONDA](../NGSPipe-RNA-seq/#Miniconda) and [download the latest source code of NGSPipeDb](../NGSPipe-RNA-seq/#NGSPipeDbSource). Modify the project name and enter the project directory.

```shell
mv NGSPipeDb species_sample_transcript_analysis_by_NGSPipeDb
cd species_sample_transcript_analysis_by_NGSPipeDb
```

### 2. Install the NGSDb conda environments <a name="DatabaseRequirement"></a>

```shell
#conda create -n ngsdb python=3.8 -c conda-forge
#conda activate ngsdb
#conda env update -n ngsdb --file ngspipe/envs/requirements_ngsdb.yaml --prune

mamba create -n ngsdb python=3.8 -c conda-forge -y
mamba env update -n ngsdb --file ngspipe/envs/requirements_ngsdb.yaml --prune
conda activate ngsdb
```
mac 系统 source ~/.bash_profile

### 3. download test data

```shell
wget expression
wget differential
wget annotation
```

### 4. convert NGSPipe analysis results to sqlite3 format <a name="Table2Sqlite3"></a>

First, edit `ngspipe/config/ngsdb.config.yaml` file:

```yaml
# path relative to where you run snakemake

# gene structure annotation
genomeAnno: "testdata/GRCm38.83.chr19.gtf" # gene annotation file, can be gtf or gff
genomeFasta: "testdata/chr19.fa" # genome sequence

# result directory of NGSPipe
#testdata_resultsDir: "results/result"
#estdata_reportsDir: "results/report"
dbDir: "results/sqlite3"

djangoCode: "ngsdb"
exp_data: "results/result/quantify/quantify_by_stringtie/gene_fpkm_all_samples.tsv"
report_path: ""
gbrowse_data: ""

# sample list file
samplesList: "testdata/samples.xls"
# sample file directory
samplesDir: "testdata"
# fastq suffix, read1
read1Suffix: "_R1.fq.gz"
read2Suffix: "_R2.fq.gz"
# replict can by 1,2,3
replict_num: 3
# condition for differential expression by deseq2
```

Second, edit `ngspipe/workflow/db_generate.Snakefile.py` file:

```python
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
annotation_gbrowse_outdir = join(gbrowse_outdir, 'annotation')

# include modules
include: join("rules", "8.db_generate_of_exp.Snakefile.py")
include: join("rules", "8.db_generate_of_gff.Snakefile.py")
include: join("rules", "8.db_generate_of_blastdb.Snakefile.py")
include: join("rules", "8.db_generate_of_genomebrowser.Snakefile.py")
```

Then, run db_generate workflow to generate database files.

```shell
snakemake -s ngspipe/db_generate.Snakefile.py --configfile ngspipe/config/ngsdb.config.yaml -p -j 1
```

### 3. config <a name="DatabaseConfig"></a>

edit `mysite/mysite/settings.py` file

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

## overview the NGSDb apps

We use django project to constructed our NGSDb. We have pareparied many apps for you. embeed. Please have a look:

django apps | description | package
---- | ----------- | -------
report | ngs analysis report   | snakemake
geneExpAtlas | table | datatable js
network      | table | datatable js
igv          | genome browse | [IGV](https://github.com/igvteam/igv.js/wiki)
blastplus    | ncbi blast +  | NCBI
ngstools     | python script         | wooey
efp          | efp browse    |
heatmap      | expression    | clustergrammer-js https://clustergrammer.readthedocs.io/clustergrammer_js.html
