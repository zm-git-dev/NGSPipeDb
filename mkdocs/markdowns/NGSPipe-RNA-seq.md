# RNA-seq analysis

__Table of Contents:__

1. [Installing wget and git](#BasicLinux)
2. [Installing Miniconda3](#Miniconda)
3. [Installing the NGSPipeDb conda environments](#NGSPipeDbEnv)
4. [Downloading the NGSPipeDb source code](#NGSPipeDbSource)
5. [Downloading the NGSPipeDb test files](#Testdata)
6. [Run test data](#RunTest)
7. [Generate report](#Report)
8. [Run your custome data](#RunRawdata)

## Quick Start - One time installation of components necessary for an individual user <a name="QuickStarted"></a>

Three commands to start analysing [test data]():

```shell
# download ngspipedb to anywhere you want
git clone git://github.com/xuanblo/NGSPipeDb.git && mv NGSPipeDb mouse_transcriptome_analysis && cd mouse_transcriptome_analysis
# download test data and create environment
bash ngspipe/scripts/one_step_ranseq_test.sh
```

Now you can viste your website on http://127.0.0.1:8000. All result are stored in `results`.
- Example of report <sub>[![html](https://img.icons8.com/ios/20/000000/html-filetype.png)](http://www.liu-lab.com)</sub>.
- Example of database <sub>[![html](https://img.icons8.com/dotty/25/000000/copy-link.png)](http://www.liu-lab.com)</sub>.


If you have more time, then we recommend you configure ngspipedb according to your needs. For more details, please see [step by step](#step-by-step) bellow.



## Step-by-step RNA-seq workflow <a href="Step-by-Step-RNASeq"></a>

Although included in this section are step-by-step instructions, it is assumed that the user has a basic understanding of the [nix command line interface](https://en.wikipedia.org/wiki/Command-line_interface). Also, best practice RNA-seq analysis is plus.

### 1. Installing wget and git <a name="BasicLinux"></a>

To get some of the required software packages, we will use the command line tools called [wget](http://www.gnu.org/software/wget/) and [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).  *wget* is a popular tool for downloading things off of the internet.  *git* is a distributed version control system which we will use to checkout the NGSPipeDb code.

__Note__: These tools are already pre-installed in most systems, but if you are unsure whether or not you have *wget* enter `wget` and if the return is `wget: command not found`, then you will have to install *wget*.  Do likewise for *git*.

### 2. Installing Miniconda3 <a name="Miniconda"></a>

NGSPipeDb relies on the conda package manager for installation and dependency resolution, so you will need to install [conda](https://conda.io/docs/user-guide/install/index.html) first.

We will be using the [Miniconda3](http://conda.pydata.org/miniconda.html) package management system (aka __CONDA__) to manage all of the software packages that __NGSPipeDb__ is dependent on. 

Use following commands to retrieve and then run the Minicoda3 installation script:
1. `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`  
2. `bash Miniconda3-latest-Linux-x86_64.sh`  
    - While running the installation script, follow the commands listed on screen, and press the _enter_ key to scroll. __Make sure to answer `yes` when asked if you want to prepend Miniconda3 to PATH.__
    - Close your terminal, open a new one and you should now have Conda working! You could, alternatively, run `source ~/.bashrc` to initiate conda.
    - Test by entering: `conda update conda`
        - Press `y` to confirm the conda updates
3. Mamba is a reimplementation of the conda package manager in C++, the fast conda-alternative. If you already have conda: 
    - `conda install mamba -c conda-forge`
    
__Note__: you will only have to install Minicoda3 once.

### 3. Downloading the NGSPipeDb project <a name="NGSPipeDbSource"></a>

To install the latest stable version of NGSPipeDb, please clone the git repository to your system.

    cd /path/to/where_you_want
    git clone https://www.github.com/xuanblo/NGSPipeDb 
   or 

    cd /path/to/where_you_want
    git clone git://www.github.com/xuanblo/NGSPipeDb

If you want to use specific version, please cleckout the [release](https://github.com/xuanblo/NGSPipeDb/releases) and issue the following commands:

    cd /path/to/where_you_want
    wget https://github.com/xuanblo/NGSPipeDb/archive/NGSPipeDb_xxx.tar.gz
    tar -xf NGSPipeDb_xxx.tar.gz

### 4. Installing the NGSPipeDb conda environments <a name="NGSPipeDbEnv"></a>

We are now ready to use conda to install the software packages which NGSPipeDb is dependent on.

First, you will need to create the conda environments needed by the various workflows. The following command will create a new conda environment containing Snakemake and python called "ngspipedb" into which conda is installed, the default environment path is `~/miniconda/conda/env/ngspipedb`.

    mamba create -c conda-forge -c bioconda --name ngspipedb snakemake python=3.8

Next, to analysis NGS data some bioinformatics tools need to be installed.

    mamba env update -n ngspipedb --file ngspipe/envs/requirement.yaml --prune

__Note__: By default, all dependencies and tools should automatically be installed on the first execution. In case of download problems, please see how to [install conda env local](), [share conda env]() for more detail.

<font color='red'>由于差异表达需要安装的包与前面的环境有不兼容，所以暂时独立了出来。</font>
```shell
conda create -n exp_r_env mamba -c conda-forge
conda activate exp_r_env
mamba env update --file ngspipe/envs/requirements_exp_r_env.yaml --prune
```

### 5. Downloading the test files <a name="Testdata"></a>

__NGSPipeDb__ is dependent on reference files which can be found for the supported species listed below: [download link](http://www.liu-lab.com/ngspipedb)

To download RNA-seq test data into `./testdata`, simply run the bash script `download_testdata.sh`:

```shell
cd NGSPipeDb
bash ngspipe/scripts/download_testdata.sh testdata
```

or you download the reference files that you need and then untarring then in a directory called `testdata`.

```shell
wget http://www.liu-lab.com/ngspipedb/rnaseq_testdata.tar.gz
tar -zxvf rnaseq_testdata.tar.gz
```

Make sure you have the following directory structure by `tree testdata`:

    testdata/
    ├── GRCm38.83.chr19.gtf
    ├── chr19.fa
    ├── control_R1.fq.gz
    ├── control_R2.fq.gz
    ├── samples.xls
    ├── treated_R1.fq.gz
    └── treated_R2.fq.gz

### 6. run RNA-seq analysis on test data <a name="RunTest"></a>

We provied a simple RNA-seq workflow for you to take a glance of NGSPipe. In RNA-seq analysis part, it contains 7 step analysis:

```python
1. sampling data
2. raw reads qc
3. junction alignmnet
4. transcript assembly
5. quantification
6. statistic
7. diff
```

First, we need to check the `ngspipe/file`. This will not execute anything, but display what would be done.

    snakemake -s ngspipe/workflow/rnaseq_analysis.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -np

If nothing goes wrong, you can generate a dag plot:

    snakemake -s ngspipe/workflow/rnaseq_analysis.Snakefile.py --dag|dot -Tpng > dag.png

Now you can do RNA-seq analysis by just one simply command.

    snakemake -s ngspipe/workflow/rnaseq_analysis.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -p -j 10

The final data files are put in the folder `results`. Please check you result file `tree -d -L 2 results/`, it may like this:

    results/
    ├── report
    │   ├── 1.rawreads_stat
    │   ├── 2.cleanreads_stat
    │   ├── 3.mapping_stat
    │   ├── 4.exp_stat
    │   └── workflow
    ├── result
    │   ├── junction_align
    │   ├── quantify
    │   ├── rawReads_qc
    │   ├── sampling_data
    │   ├── statistic
    │   └── transcript_assembly
    └── sqlite3
        ├── blastdb
        ├── exp
        ├── gbrowse
        └── gff_sqlite3

__Note__: If you encounter any problem in this step, please turn to `TroubleShooting` for help.

### 7. rgenerate report <a name="Report"></a>

If all goes well, the proper analysis will be followed by the making of the html report using Snakemake to a html report file with pictures and tables.

    snakemake --snakefile ngspipe/workflow/rnaseq_analysis.Snakefile.py --report results/Report/report.html

The final report should appear as `results/Report/report.html`. This report is a single html file with all in it and can be sent to customers/colleagues as a final report. It is nicer than a PDF version because of large tables and figures which would suffer from page breaks and it can be viewed on any device supporting html include smartphones :-).

__Note__: Internet Explorer is not supported. If you get connected error in this step, you can solve this problem by edit file `ngspipedb_py38_conda_env/lib/python3.8/site-packages/snakemake/report/report.html.jinja2` to change `https://raw.githubusercontent.com/eligrey/FileSaver.js/2.0.0/src/FileSaver.js` to `https://cdnjs.cloudflare.com/ajax/libs/FileSaver.js/2.0.0/FileSaver.js`


### 8. run your custome data <a name="RunRawdata"></a>

NGSPipe is built to be used routinely. To ensure a maximum comparability of the results, you can copy default `ngspipe/config/rnaseq.config.yaml` and `ngspipe/workflow/rnaseq_analysis.Snakefile.py` file from the ngspipe source directory:

```shell
cp ngspipe/config/rnaseq.config.yaml ngspipe/workflow/rnaseq_analysis.Snakefile.py ./
```

We will explain how to edit and configure these files shortly below.

#### Rawdata sequence data

Raw data files can either be fastq, fastq.gz formated files. If your raw data are located in __somewhere else__, you can copy them to `rawdata`, or create soft links like `ln -s ../yoursamplepath/*.fq.gz rawdata/`.

    rawdata/
    ├── lung-rep1_R1.fq.gz
    ├── lung-rep1_R2.fq.gz
    ├── lung-rep2_R1.fq.gz
    ├── lung-rep2_R2.fq.gz
    ├── lung-rep3_R1.fq.gz
    ├── lung-rep3_R2.fq.gz
    ├── liver-rep1_R1.fq.gz
    ├── liver-rep1_R2.fq.gz
    ├── liver-rep2_R1.fq.gz
    ├── liver-rep2_R2.fq.gz
    ├── liver-rep3_R1.fq.gz
    ├── liver-rep3_R2.fq.gz

As recommended above, if all of your raw data are located in __rawdata__, then create a `rawdata/samples.xls` file like:

    lung-rep1
    lung-rep2
    lung-rep3
    liver.rep1
    liver.rep2
    liver.rep3

__Note__: You cannot mix Paired-end and Single-end samples within the same NGSPipe run as this will cause an ERROR. NGSPipe only support Paired-end samples.

#### reference data

You can download reference data from NCBI, Ensembl, or anywhere else to `Genomedata`. The most import file is genome in Fasta formant and gene annotation in GFF/GTF format.

    testdata/
    ├── GRCm38.83.chr19.gtf
    └── chr19.fa

#### edit config file

In this section, you will need to specify every term to match your own machine, reference genome, and data.

```yaml
# path relative to where you run snakemake

## reference ##
genomeAnno: "testdata/GRCm38.83.chr19.gtf" # gene annotation file, can be gtf or gff
genomeFasta: "testdata/chr19.fa" # genome sequence

## output directory ##
resultsDir: "results/result"
reportsDir: "results/report"
dbDir: "results/sqlite3"

## input ##
samplesList: "testdata/samples.xls" # sample file
samplesDir: "testdata" # sample file directory
read1Suffix: "_R1.fq.gz" # fastq suffix, read1
read2Suffix: "_R2.fq.gz"
replict_num: 3 # replict can by 1,2,3
# condition for differential expression by deseq2
```

#### condition for sample compare

To compare each sample, a group information is needed.

    sample_id,Sample,Tissue
    control-0,control,normal
    control-1,control,normal
    control-2,control,normal
    treated-0,treated,tumor
    treated-1,treated,tumor
    treated-2,treated,tumor

#### edit snakefile

edit file `ngspipe/workflow/rnaseq_analysis.Snakefile.py` for advance setting, such as sampling data method, mapping tool, and email address to receive run log.

```python
# relative path
snake_dir = workflow.basedir # all configfile, scripts, restructuretext, ens are relative to snakefile (this file)
working_dir = os.getcwd() # input and output path are relative to current working directory

# get notice
receiver_email = 'youemailaddress'

# ----------------------------------------------------------------------- #
# sample information #
#
smpList = pd.read_csv(config["samplesList"], index_col=0, header=None)
SAMPLES = list(smpList.index)[0:] # how many sample you want to run.
# ----------------------------------------------------------------------- #


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# detail parameters in pipe #
#
# 1. sampling data
# for test the pipe, you can choose to the part of the input file, can be whole,head:40000,tail:40000,random:0.5,random:40000
sampling_method = 'links' # tail, seqkit_number, seqkit_proportion, head, tail
sampling_data_outdir = join(config["resultsDir"], "sampling_data", "sampling_data_by_{}".format(sampling_method))

# 2. raw reads qc
qc_method = 'trim-galore' # trimomatic
qc_outdir = join(config["resultsDir"], "rawReads_qc", "rawReads_qc_by_{}".format(qc_method))

# 3. junction alignmnet
junction_align_method = 'hisat2' # star
junction_align_outdir = join(config["resultsDir"], "junction_align", "junction_align_by_{}".format(junction_align_method))
genome_index_prefix = "genome"
# rna-seq sequencing type, can be fr-firststrand, none, fr-secondstrand
rna_library = "" # "--rna-strandness RF"(fr-firststrand) or "--rna-strandness FR"(fr-secondstrand)

# 4. transcript assembly
transcript_assembly_method = 'stringtie' # star
transcript_assembly_outdir = join(config["resultsDir"], "transcript_assembly", "transcript_assembly_by_{}".format(transcript_assembly_method))

# 5. quantification
quantify_method = 'stringtie' # htseqcounts or featurecounts
quantify_outdir = join(config["resultsDir"], "quantify", "quantify_by_{}".format(quantify_method))

# 6. statistic
statistic_data_all = [
                  '0.genomeFa', 
                  '0.genomeAnno', 
                  '1.rawReads', 
                  '2.cleanReads', 
                  '2.multiqc', 
                  '3.bam', 
                  '4.mergedGtf', 
                  '5.exp',
                  ]
genomeFa_outdir, genomeAnno_outdir, rawReads_outdir, cleanReads_outdir, multiqc_outdir, bam_outdir, mergedGtf_outdir, exp_outdir \
  = [join(config["resultsDir"], "statistic", "statistic_data_of_{}".format(i)) for i in statistic_data_all]

statistic_data_choose = [
                  #'0.genomeFa', 
                  #'0.genomeAnno', 
                  '1.rawReads', 
                  '2.cleanReads', 
                  #'2.multiqc', 
                  '3.bam', 
                  #'4.mergedGtf', 
                  #'5.exp',
                  ]
stat_outdir = join(config["resultsDir"], "statistic")

# 7. diff gene discovery
diff_outdir = join(config["resultsDir"], "diff", "diff_by_{}".format('deseq2'))
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# ------------------------------
#report: "report/workflow.rst"
report_outdir = join(config["reportsDir"], "report.html")
# ------------------------------

# include modules
include: join("rules", "1.sampling_data_by_{}.Snakefile.py".format(sampling_method))
include: join("rules", "2.rawReads_qc_by_{}.Snakefile.py".format(qc_method))
include: join("rules", "3.junction_align_by_{}.Snakefile.py".format(junction_align_method))
include: join("rules", "4.transcript_assembly_by_{}.Snakefile.py".format(transcript_assembly_method))
include: join("rules", "5.quant_by_{}.Snakefile.py".format(quantify_method))
include: join("rules", "6.statistic_data_of_bam.Snakefile.py")
include: join("rules", "6.statistic_data_of_rawReads.Snakefile.py")
include: join("rules", "6.statistic_data_of_cleanReads.Snakefile.py")
include: join("rules", "7.report.Snakefile.py")
include: join("rules", "9.differential_expression.deseq2.Snakefile.py")

```

__Note__: The script/env/include path is always relative to the Snakefile containing the directive (in contrast to the input, output and log file paths, which are relative to the working directory). 

#### report

edit `ngspipe/report/*.rst` that will be added at the end of the report. For example, edit the `ngspipe/report/rawreads_stat.rst` file to include a text describing the 'Statistic' of the experiment. This text will be added to the report as static section and is one of the two report sections that can be edited by the end-user.

```matlab

summary of reads produced
=========================

The details of the data quality are as follows:

#. Sample name: Sample name
#. Raw reads: Count raw sequence data, with four rows in a unit, count the number of sequencing sequences in each file
#. Clean reads: The calculation method is the same as that of Raw reads, except that the statistical files are filtered sequencing data, and subsequent bioinformatics analysis is based on Clean reads
#. Clean bases: Multiply the number of sequencing sequences by the length of the sequencing sequence and convert it to G as the unit
#. Error rate: base sequencing error rate
#. Q20, Q30: respectively represent the percentage of bases whose Phred value is greater than 20 and 30 to the total bases
#. GC content: the percentage of the total number of bases G and C to the total number of bases
```

#### run snakemake


```shell
# dry run, use -n parameter only print task plan, -p print commands
snakemake -s ngspipe/workflow/rnaseq_analysis.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -np

# run pipe
snakemake -s ngspipe/workflow/rnaseq_analysis.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -p -j <thread number>

# after run pipe successfully, report can be generated.
snakemake --snakefile ngspipe/workflow/rnaseq_analysis.Snakefile.py --report results/Report/report.html
```