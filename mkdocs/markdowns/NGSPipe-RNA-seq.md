---
title: 
date: 2021-03-07 09:51:19
tags:
  - Markdown
  - rnaseq
categories: module

---

# RNA-seq analysis

## Quick Start - One time installation of components necessary for RNA-Seq analysis <a name="QuickStarted"></a>

=== "Linux & WSL"

    1.download ngspipedb to anywhere you want
    ```shell
    git clone git://github.com/xuanblo/NGSPipeDb.git
    mv NGSPipeDb mouse_transcriptome_analysis
    cd mouse_transcriptome_analysis
    ```

    2.create conda environment, download test data, and run RNA-Seq analysis
    ```shell
    bash ngspipe/scripts/one_step_ranseq_test.sh linux
    ```

=== "MacOSX"

    1.download ngspipedb to anywhere you want
    ```shell
    git clone git://github.com/xuanblo/NGSPipeDb.git
    mv NGSPipeDb mouse_transcriptome_analysis
    cd mouse_transcriptome_analysis
    ```

    2.create conda environment, download test data, and run RNA-Seq analysis
    ```shell
    bash ngspipe/scripts/one_step_ranseq_test.sh macos
    ```

!!! info
    If you have more time, then we recommend you configure ngspipedb according to your needs. For more details, please see [step by step](#step-by-step) bellow.

## Step-by-step RNA-seq workflow on testdata <a href="Step-by-Step-RNASeq"></a>

Although included in this section are step-by-step instructions, it is assumed that the user has a basic understanding of the [nix command line interface](https://en.wikipedia.org/wiki/Command-line_interface). Also, it would be better if the user has basic knowledge about [snakemake](), [conda]() and [best practice RNA sequence analysis](), but it is not required. You can find some easy-to-learn matierals in [linux & shell](../linux) and [RNASeq background](../ngs#rnaseq).

### 1. Install wget and git <a name="BasicLinux"></a>

To get some of the required software packages, we will use the command line tools called [wget](http://www.gnu.org/software/wget/) and [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).  *wget* is a popular tool for downloading things off of the internet.  *git* is a distributed version control system which we will use to checkout the NGSPipeDb code.

!!! note
    These tools are already pre-installed in most systems, but if you are unsure whether or not you have *wget* enter `wget` and if the return is `wget: command not found`, then you will have to install *wget*.  Do likewise for *git*.

### 2. Install Miniconda3 <a name="Miniconda"></a>

NGSPipeDb relies on the conda package manager for installation and dependency resolution, so you will need to install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) first.

We will be using the [Miniconda3](http://conda.pydata.org/miniconda.html) package management system (aka __CONDA__) to manage all of the software packages that __NGSPipeDb__ is dependent on. 

Use following commands to retrieve and then run the Minicoda3 installation script:

1.download miniconda3
```shell
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

2.install miniconda3
```shell
bash Miniconda3-latest-Linux-x86_64.sh
```
!!! important
    While running the installation script, follow the commands listed on screen, and press the _enter_ key to scroll. Make sure to answer `yes` when asked if you want to prepend Miniconda3 to PATH. After that, close your terminal, open a new one and you should now have Conda working! You could, alternatively, run `#!shell source ~/.bashrc` to initiate conda.  

3.Test by entering: `conda update conda`. Press `y` to confirm the conda updates.  

4.Finally, `conda install mamba -c conda-forge`. Mamba is a reimplementation of the conda package manager in C++, the fast conda-alternative.

!!! Note
    You will only have to install Minicoda3 once.

### 3. Download NGSPipeDb source code <a name="NGSPipeDbSource"></a>

We offer two download options.

=== "latest stable version"

    To install the latest stable version of NGSPipeDb, please clone the git repository to your system. First, `cd /path/to/where_you_want`, and then
    ```shell
    git clone https://www.github.com/xuanblo/NGSPipeDb
    ```

    !!! error
        If you meet problem with network, please try to download the "specific version" or modify the command to:  
        `git clone git://www.github.com/xuanblo/NGSPipeDb`.

=== "specific version"

    If you want to use specific version, please have a look at the [change log](../changelog) and different [release version](https://github.com/xuanblo/NGSPipeDb/releases) first.
    
    Then, `#!shell cd /path/to/where_you_want`
    ```shell
    wget https://github.com/xuanblo/NGSPipeDb/archive/NGSPipeDb_v0.0.1.tar.gz
    tar -xf NGSPipeDb_v0.0.1.tar.gz
    ```

### 4. Choose a good project name

Because all analysis will be done directly in the source code directory. So you can change the directory name in the following way.

=== "download by git"

    ```shell
    mv NGSPipeDb species_sample_transcript_analysis_by_NGSPipeDb
    ```

=== "download by wget"

    ```shell
    mv NGSPipeDb-NGSPipeDb_v0.0.1 species_sample_transcript_analysis_by_NGSPipeDb
    ```

Enter directory `#!shell cd species_sample_transcript_analysis_by_NGSPipeDb` for further setps.

You can view the __NGSPipe__ code structure with `#!shell tree ./ngspipe/ -d -L 2`

    ./ngspipe/
    ├── config
    ├── envs
    ├── imgs
    ├── metadata
    ├── notebooks
    ├── reports
    ├── rules
    └── scripts

!!! note
    Files other than `ngspipe` are not required for RNA-Seq analysis.

### 5. Installing the NGSPipe RNA-Seq conda environments <a name="NGSPipeDbEnv"></a>

First, you will need to create a conda environments. For details, see [manage envirement](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) with conda.

```shell
conda create -c conda-forge -c bioconda --name ngspipe-rnaseq snakemake=5.30.2 python=3.8 seqkit=0.14.0
```

Next, to producess and analysis the RNA-Seq sequenceding data, some bioinformatics tools need to be installed. Chick to [see software we used](https://github.com/xuanblo/NGSPipeDb/blob/master/ngspipe/envs/requirements_rnaseq.yaml).

```shell
conda env update -n ngspipe-rnaseq --file ngspipe/envs/requirements_rnaseq.yaml --prune
```

Activate conda environment:

```shell
conda activate ngspipe-rnaseq
```

!!! Note
    By default, all softwares and packages are stored under `~/miniconda/conda/env/ngspipe-rnaseq`.

!!! error
    If the conda downloading encounters any problems, you can refer to how to [download pre-build ngspipe-ranseq env](../conda/#pre-build).

### 6. Download test files <a name="Testdata"></a>

__NGSPipe__ is dependent on reference files and raw sequence reads which can be found in [http://www.liu-lab.com/ngspipedb/testdata](http://www.liu-lab.com/ngspipedb/testdata)

To download the mouse RNA-seq test data into `./testdata`, you can choose auto download with a script or Manually:

=== "auto download by a shell script"

    ```shell
    bash ngspipe/scripts/download_testdata.sh testdata
    ```

=== "Manually download by wget"

    1.download the test files (160M) that you need and then untarring then in a directory called `testdata`.
    ```shell
    wget http://www.liu-lab.com/ngspipedb/rnaseq_testdata.tar.gz
    tar -zxvf rnaseq_testdata.tar.gz
    ```

    2. generate replicated samples:
    ```shell
    cd testdata
    python ../ngspipe/scripts/generate_replicat.py control_R1.fq.gz control_R2.fq.gz treated_R1.fq.gz treated_R2.fq.gz
    gunzip chr19.fa.gz
    gunzip GRCm38.83.chr19.gtf.gz
    cd ..
    ```

Make sure you have the following directory structure by command `tree testdata`:

    testdata
    ├── chr19.fa
    ├── condition.xls
    ├── control-0_R1.fq.gz
    ├── control-0_R2.fq.gz
    ├── control-1_R1.fq.gz
    ├── control-1_R2.fq.gz
    ├── control-2_R1.fq.gz
    ├── control-2_R2.fq.gz
    ├── control_R1.fq.gz
    ├── control_R2.fq.gz
    ├── GRCm38.83.chr19.gtf
    ├── samples.xls
    ├── treated-0_R1.fq.gz
    ├── treated-0_R2.fq.gz
    ├── treated-1_R1.fq.gz
    ├── treated-1_R2.fq.gz
    ├── treated-2_R1.fq.gz
    ├── treated-2_R2.fq.gz
    ├── treated_R1.fq.gz
    └── treated_R2.fq.gz

!!! warning
    The test data is only used to verify that the analytical process is working properly and the results do not have r biological significance.

### 7. run RNA-seq analysis on test data <a name="RunTest"></a>

We provied a simple RNA-seq workflow for you to take a glance of NGSPipe. In RNA-seq analysis part, it contains 7 step analysis:

```python
1. sampling data (choose part of your data)
2. raw reads qc
3. junction align to genome
4. transcript assembly
5. gene quantification
6. statistic
7. differential gene analysis
```

First, we need to check the `ngspipe/file`. This will not execute anything, but display what would be done.

```shell
# dry run, use -n parameter only print task plan, -p print commands
snakemake -s ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -np
```

If nothing goes wrong, you can generate a dag plot:

```shell
snakemake -s ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml --dag|dot -Tpng > dag.png
```

Now you can do RNA-seq analysis by just one simply command.

```shell
snakemake -s ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -p -j 10
```

The final data files are put in the folder `results`. Please check you result file `tree -d -L 2 results/`, it may like this:

    results/rnaseq
    ├── report
    │   ├── 1.rawreads_stat
    │   ├── 2.cleanreads_stat
    │   ├── 3.mapping_stat
    │   ├── 4.exp_stat
    │   └── workflow
    └── result
        ├── diff
        ├── junction_align
        ├── quantify
        ├── rawReads_qc
        ├── sampling_data
        └── statistic

!!! Note
    If you encounter any problem in this step, please turn to `TroubleShooting` for help.

### 8. rgenerate report <a name="Report"></a>

If all goes well, the proper analysis will be followed by the making of the html report using Snakemake to a html report file with pictures and tables.

```shell
snakemake --snakefile ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml --report results/report/report.html
```

The final report should appear as `results/report/report.html`. This report is a single html file with all in it and can be sent to customers/colleagues as a final report. It is nicer than a PDF version because of large tables and figures which would suffer from page breaks and it can be viewed on any device supporting html include smartphones :-).

!!! Note
    Internet Explorer is not supported.

## run your custome data <a name="RunRawdata"></a>

NGSPipe is built to be used routinely. To ensure a maximum comparability of the results, you can copy default `ngspipe/config/rnaseq.config.yaml` and `ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py` file to the same directory:

```shell
cp ngspipe/config/rnaseq.config.yaml ngspipe/config/my_own_rnaseq.config.yaml 
cp ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py ngspipe/my_own_rnaseq_analysis.Snakefile.py
```
We will explain how to edit and configure these files shortly below.

### 1. Rawdata sequence data <a name="download_raw_data"></a>

Raw data files can either be fastq, fastq.gz formated files. `mkdir rawdata` and upload your own data to this directory. If your raw data are located in __somewhere else__, you can copy them to `rawdata`, or create soft links like `ln -s ../yoursamplepath/*.fq.gz rawdata/`.

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
    liver-rep1
    liver-rep2
    liver-rep3

!!! Note
    You cannot mix Paired-end and Single-end samples within the same NGSPipe run as this will cause an ERROR. NGSPipe only support Paired-end samples.

### 2. Reference data <a name="download_ref_data"></a>

You can download reference data from NCBI, Ensembl, or anywhere else to `genomedata`. The most import file is genome in Fasta formant and gene annotation in GFF/GTF format. Use the same method as rawdata does:

    genomedata/
    ├── GRCm38.83.chr19.gtf
    └── chr19.fa

### 3. Edit config file <a name="edit_config"></a>

In this section, you will need to specify every term to match your own machine, reference genome, and data.

```yaml
# path relative to where you run snakemake

## reference ##
genomeAnno: "testdata/GRCm38.83.chr19.gtf" # gene annotation file, can be gtf or gff
genomeFasta: "testdata/chr19.fa" # genome sequence

## output directory ##
resultsDir: "results/rnaseq"
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

!!! Note
    The input, output file paths are relative to the working directory. 

### 4. Condition for sample compare <a name="edit_condition"></a>

To perform gene differential analysis, please create a `rawdata/condition.xls` file.

    sample_id,Sample,Tissue
    lung-rep1,lung,normal
    lung-rep2,lung,normal
    lung-rep3,lung,normal
    liver-rep1,liver,tumor
    liver-rep2,liver,tumor
    liver-rep3,liver,tumor

### 5. edit snakefile <a name="edit_snakefile"></a>

Edit file `ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py` for advance setting, such as sampling data method, mapping tool, and email address to receive run log.

``` python
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

!!! Note
    The `script/env/include` path is always relative to the Snakefile containing the directive (in contrast to the input, output and log file paths, which are relative to the working directory). 

### 6. Custome report <a name="custom_report"></a>

edit `ngspipe/report/*.rst` that will be added at the end of the report. For example, edit the `ngspipe/report/rawreads_stat.rst` file to include a text describing the 'Statistic' of the experiment. This text will be added to the report as static section and is one of the two report sections that can be edited by the end-user.

```yaml

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

NGSPipe will provide you with an interactive, browser-based report, showing the most important measures of your data on the first sight. All tables in the report can be sorted and filtered. The table on the first tab shows the key values for a quick estimation of the success of your sequencing experiment and the assembly. On the second tab, there is a more detailed table, giving many additional measures. Additionally to the tables, many measures are provided as graphical feedback. On the third tab, you see plots which are generated for one complete sequencing experiment. On the fourth tab, there are plots which each show measures on one specific dataset.

### 7. Run snakemake <a name="custom_run"></a>

Run snakemake in a prebuild environment:

```shell
snakemake -s ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -p -j <thread number>
```

Or install dependencies while snakemake run. By default, all dependencies and tools should automatically be installed on the first execution. 

```shell
snakemake -s ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -p -j <thread number> --use-conda
```

See Specifying a location for an environment or run conda create --help for information on specifying a different path.

```shell
snakemake -s ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -p -j <thread number> -p ./ngspipe_rnaseq_macox64
```