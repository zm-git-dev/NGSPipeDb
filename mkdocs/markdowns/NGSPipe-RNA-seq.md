## RNA-seq analysis

    1. [Installing wget and git](#BasicLinux)
    2. [Installing Miniconda3](#Miniconda)
    3. [Installing the NGSPipeDb conda environments](#NGSPipeDbEnv)
    4. [Downloading the NGSPipeDb source code](#NGSPipeDbSource)
    5. [Downloading the NGSPipeDb test files](#Testdata)
    6. [Run test data](#RunTest)
    7. [Generate report](#Report)
    8. [Run your custome data](#RunRawdata)
    1. [Installing requirement](#DatabaseRequirement)
    2. [Convert table file to sqlite3](#Table2Sqlite3)
    3. [Database config](#DatabaseConfig)
    4. [Start server](#RunServer)

## Quick Start - One time installation of components necessary for an individual user <a name="QuickStarted"></a>

Three commands to start analysing test data:
```shell
# download ngspipedb to anywhere you want
git clone https://github.com/xuanblo/NGSPipeDb.git && mv NGSPipeDb mouse_transcriptome_analysis && cd mouse_transcriptome_analysis
# download test data and create environment
bash ngspipe/scripts/one_step_runtest.sh
```

Now you can viste your website on http://127.0.0.1:8000. All result are stored in `results`.
- Example of report <sub>[![html](https://img.icons8.com/ios/20/000000/html-filetype.png)](http://www.liu-lab.com)</sub>.
- Example of database <sub>[![html](https://img.icons8.com/dotty/25/000000/copy-link.png)](http://www.liu-lab.com)</sub>.


If you have more time, then we recommend you configure ngspipedb according to your needs. For more details, please see [step by step](#step-by-step) bellow.



## Step-by-step RNA-seq workflow <a href="Step-by-Step-RNASeq"></a>

Although included in this section are step-by-step instructions, it is assumed that the user has a basic understanding of the [nix command line interface](https://en.wikipedia.org/wiki/Command-line_interface). Also, best practice RNA-seq analysis is plus

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
__Note__: the xxx refers to the latest changeset of NGSPipeDb, so it will differ.

### 4. Installing the NGSPipeDb conda environments <a name="NGSPipeDbEnv"></a>

We are now ready to use conda to install the software packages which NGSPipeDb is dependent on.

First, you will need to create the conda environments needed by the various workflows. The following command will create a new conda environment containing Snakemake and python called "ngspipedb" into which conda is installed, the default environment path is `~/miniconda/conda/env/ngspipedb`.

    mamba create -c conda-forge -c bioconda --name ngspipedb snakemake python=3.8

Next, to analysis NGS data some bioinformatics tools need to be installed.

    mamba env update -n ngspipedb --file ngspipe/envs/requirement.yaml --prune

__Note__: By default, all dependencies and tools should automatically be installed on the first execution. In case of download problems, please see how to [install conda env local](), [share conda env]() for more detail.

```shell
conda create -n exp_r_env mamba -c conda-forge
conda activate exp_r_env
mamba env update --file ngspipe/envs/requirements_exp_r_env.yaml --prune
```

### 5. Downloading the test files <a name="Testdata"></a>

__NGSPipeDb__ is dependent on reference files which can be found for the supported species listed below:  [download link](http://www.liu-lab.com/ngspipedb)

To download RNA-seq test data into `./testdata`, simply run the bash script `download_testdata.sh`:

    cd NGSPipeDb
    bash ngspipe/scripts/download_testdata.sh testdata`

or you download the reference files that you need and then untarring then in a directory called `testdata`.

    wget http://www.liu-lab.com/ngspipedb/rnaseq_testdata.tar.gz
    tar -zxvf rnaseq_testdata.tar.gz

Make sure you have the following directory structure by `tree testdata`:

```
testdata/
├── GRCm38.83.chr19.gtf
├── chr19.fa
├── control_R1.fq.gz
├── control_R2.fq.gz
├── samples.xls
├── treated_R1.fq.gz
└── treated_R2.fq.gz
```

### 6. run RNA-seq analysis on test data <a name="RunTest"></a>

We provied a simple RNA-seq workflow for you to take a glance of NGSPipe. In RNA-seq analysis part, it contains 7 step analysis:

```python
1. sampling data
2. raw reads qc
3. junction alignmnet
4. transcript assembly
5. quantification
6. statistic
```

First, we need to check the workflow/path/file. This will not execute anything, but display what would be done.

    snakemake -s ngspipe/rnaseq_analysis.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -np

If nothing goes wrong, you can generate a dag plot:

    snakemake -s ngspipe/rnaseq_analysis.Snakefile.py --dag|dot -Tpng > dag.png

Now you can do RNA-seq analysis by juest one simply command.

    snakemake -s ngspipe/rnaseq_analysis.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -p -j 10

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

    snakemake --snakefile ngspipe/rnaseq_analysis.Snakefile.py --report results/Report/report.html

- the final report should appear as `results/Report/report.html`. This report is a single html file with all in it and can be sent to customers/colleagues as a final report. It is nicer than a PDF version because of large tables and figures which would suffer from page breaks and it can be viewed on any device supporting html include smartphones :-).

__Note__: Internet Explorer is not supported. If you get connected error in this step, you can solve this problem by edit file `ngspipedb_py38_conda_env/lib/python3.8/site-packages/snakemake/report/report.html.jinja2` to change `https://raw.githubusercontent.com/eligrey/FileSaver.js/2.0.0/src/FileSaver.js` to `https://cdnjs.cloudflare.com/ajax/libs/FileSaver.js/2.0.0/FileSaver.js`


### 8. run your custome data <a name="RunRawdata"></a>

NGSPipeDb is built to be used routinely. To ensure a maximum comparability of the results, a default config.yaml file is generated when calling the AQUAMIS.py wrapper script. The wrapper itself only allows configuring basic functionalities. The config.yaml can be initialized by starting AQUAMIS with the dry-run flag -n . Then, you can alter it to configure NGSPipeDb in more detail.

Running this code requires:

- installing **miniconda** and all required dependencies from the provided **environment.yaml** (```conda env create --name pinfish_3.6 --file=environment.yaml```)
- editing the **config.yaml** file to match your own machine, reference genome, and data
- adding the required data files in due locations (matching the yaml)
- edit the **Preamble.md** file to include a text describing the 'Aim' of the experiment. This text will be added to the report as first section and is one of the two report sections that can be edited by the end-user.
- edit `ngspipe/report/*.rst` that will be added at the end of the report.

- run Snakemake with ```snakemake --use-conda -j <thread number>```


[view the report hosted here](http://htmlpreview.github.io/?https://github.com/Nucleomics-VIB/Nanopore_Pinfish_Analysis/blob/master/Nanopore_Pinfish_Analysis.html)

rem: when something breaks the snake, or if you add more text/comments in the initial Rmd report, you can regenerate the report manually with ```R --slave -e 'rmarkdown::render("Nanopore_Pinfish_Analysis.Rmd", "html_document")'``` within the base project folder.

edit file `NGSPipeCode/config.yaml` for general data path or something. edit file `snakefile` for general data path or something.
Configuring the META files: config.yaml <a name="config"></a>
The config.yaml file has three main sections. __PATHS__, __PARAMS__, __SAMPLES__:

edit file `NGSPipeCode/Snakefile` for advance setting, such as sampling method, mapping tool, email address to receive run log.



3. create samples.xls, for example, if you have two samples named "control" and "treated", just create a text file (maybe named sample.xls) with one column and two rows.

resources/testdata/sample.info.xls:
```
control
treated
```

Raw data files can either be fastq, fastq.gz, or bam formated files. If your raw data are located in __somewhere else__, you can copy them to `rawdata`, or create soft links like `cd rawdata && ln -s yoursamplepath/*.fq.gz ./`.

As recommended above, if all of your raw data are located in __rawdata__, then create a `samples.xls` file like:

```
lung.rep1
lung.rep2
lung.rep3
liver.rep1
liver.rep2
liver.rep3
```

__If you did not follow the recommended best practice__ then you will have to specify the full paths here.

Each sample should be given a __NAME__ (arbitrary text) and a __PATH__


__IMPORTANT__: __You cannot mix Paired-end and Single-end samples within the same VIPER run as this will cause an ERROR__. If necessary, run all of one type of data first, followed by the net type of data after.

__*2. Copying over the META files:*__

The __META__ files (*config.yaml* and *metasheet.csv*) allow you to configure each run.  They are explained in much more detail below.  For now, we simply copy them from the viper source directory:
```
    cd PROJECT
    cp viper/config.yaml .
    cp viper/metasheet.csv .
```
__We will explain how to edit and configure these files shortly below__


In this section, you will need to specify the location of the following static reference files.

__The script path is always relative to the Snakefile containing the directive (in contrast to the input and output file paths, which are relative to the working directory).__

All paths in the snakefile are interpreted relative to the directory snakemake is executed in. 

__*3. custom your snakefile*__

__*4. run*__

```shell
# dry run, use -n parameter only print task plan, -p print commands
snakemake -np --snakefile NGSPipeCode/Snakefile --configfile NGSPipeCode/config.yaml

# run pipe
snakemake -p --snakefile NGSPipeCode/Snakefile --configfile NGSPipeCode/config.yaml -j1


```

__Note__: Input, output and log are relative to your execution directory. Other paths, such as Env and include, are relative to the snapfile path
