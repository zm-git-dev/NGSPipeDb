# NGSPipeDb

[toc]

## NGSPipeDb介绍

#### 为什么我们开发NGSPipeDb

不同类型组学数据的比较和整合研究有助于揭示不同类型的分子和生化反应复杂的调控机制。然而目前整合多组学的研究工作还很缺乏,尤其是基于不同类型组学数据关联分析的调控机制研究。

#### NGSPipeDb有什么特点

1. 可重复的数据分析。
2. 基于工作流的数据自动化分析流程，解决组学数据分析的难题。
3. 结构化数据存储和分析平台。

## NGSPipe使用

#### 1. prepare working directory

首先在命令行下执行如下命令：

```shell
# creating directory
mkdir jatropha_transcriptome_analysis

# enter directory
cd jatropha_transcriptome_analysis

# download NGSPipeDb source code
git clone http://xuanblo.github.com/NGSPipeDb
```

然后执行`tree -L 2 -d`，你将会看到如下的目录结构：

```shell
.
├── NGSDjangoCode # 主要存放Django网站代码
├── NGSDjangoDB # 数据库需要的标准数据，如表达矩阵，注释信息等等
├── NGSPipeCode # 主要存放转录组分析的snakemake流程代码
│   ├── modules # 可以组合的rule，环境，report配置
│   └── script # 画图和数据处理代码
├── NGSPipeOut # 主要存放转录组分析的结果
│   ├── Report # 报告
│   └── Result # 最终结果
├── Rawdata # 转录组分析的原始数据
└── Testdata # 转录组分析的测试数据

10 directories
```

code size control：`du -sh ./*`

```shell
0	./NGSDjangoCode
0	./NGSDjangoDB
50M	./NGSPipeCode
8.0K	./NGSPipeOut
0	./Rawdata
8.0K	./README.md
16K	./Testdata
```

#### 2. install conda environment

- On linux

```shell
cd NGSPipeDb

# download latest miniconda3 and install
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/Miniconda3-latest-Linux-x86_64.sh && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3

# activate base and set miniconda path
conda init

# Conda Prompt Customization
conda config --set env_prompt '({name}) '

source ~/.bashrc

# update conda, (optional)
conda update conda

# create conda visual environment, python version, snakemake version, env directory,django version
conda create -p ngspipedb_py38_conda_env python=3.8

# activate conda env
conda activate ./ngspipedb_py38_conda_env
conda install mamba -c conda-forge
mamba install --file req.txt

# exit env
conda deactivate
```

2. On Mac

```shell
# mac用户请用：
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

3. On Win10

```shell
hello
```

#### 3. using test data (Optional)

1. download test data

```shell
cd Testdata
sh RunMe.sh
```

2. generate replicate data (Optional)

```shell
# move your data to Rawdata and 
```

3. create samples.xls

```
```

#### 4. Custom data

1. your own data

2. download rawdata from NCBI

#### 5. config NGSPipe

edit file `NGSPipeCode/config.yaml`

```yaml
# gene annotation file, can be gtf or gff
genomeAnno: "../test_data/GRCm38.83.chr19.gtf"
# genome sequence
genomeFasta: "../test_data/chr19.fa"
# result directory of NGSPipe
resultsDir: "../NGSPipeOut/Result/20201102-StringtieMaxIntron1000"
# sample description file
samplesPath: "../Testdata/samples.xls"
# fastq suffix
fastq: "gz"
# rna-seq sequencing type, can be fr-firststrand, none, fr-secondstrand
rna_library: "fr-firststrand"
# for test the pipe, you can choose to the part of the input file, can be whole,head:40000,tail:40000,random:0.5,random:40000
preRun: "random:"
```

#### 6. run NGSPipe

```shell
snakemake -np --snakefile NGSPipeCode/Snakefile --configfile NGSPipeCode/config.yaml
snakemake -p --snakefile NGSPipeCode/Snakefile --configfile NGSPipeCode/config.yaml -j1
snakemake -p --snakefile NGSPipeCode/Snakefile --configfile NGSPipeCode/config.yaml -j1 --use-conda --conda-frontend mamba
snakemake -np --snakefile Snakefile.test.py --configfile config.yaml 
snakemake -s rnaseq.snakefile.py --use-conda -n 20
snakemake --snakefile NGSPipeCode/Snakefile --configfile NGSPipeCode/config.yaml --report report.html
```
rsync保留权限，scp不保留权限

## NGSDb使用

#### 0. env

```shell
pip install django
cd NGSDjangoCode
django-admin.py startproject fresh
```

查看目录`tree fresh`

```shell
fresh/
├── fresh
│   ├── __init__.py
│   ├── asgi.py
│   ├── settings.py
│   ├── urls.py
│   └── wsgi.py
└── manage.py

```
django-admin.py startproject jatrophaDb
cd jatrophaDb
python manage.py startapp home
python manage.py startapp geneExpAtlas
python manage.py startapp blast
python manage.py startapp geneAnno
python manage.py runserver
```

1 directory, 6 files
```

#### config

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

#### 1. using test data

#### 2. run NGSPipe for custom data

#### 3. generate sqlite

#### 4. run server

```shell
python manage.py runserver 0.0.0.0:8000
```

## make clean清除数据重来

1. NGSPipeData
```shell
# Testdata
find Testdata/* | grep -v '\(RunMe.sh\|samples.xls\)' | xargs rm
# NGSPipeOut
find NGSPipeOut/* | grep -v '\(RunMe.sh\|samples.xls\)' | xargs rm
```

2. NGSDBData
```shell
```

## Troubleshooting

#### 1. conda install software error
    error1:
    ```
    Collecting package metadata (current_repodata.json): done
    Solving environment: failed with initial frozen solve. Retrying with flexible solve.
    Solving environment: failed with repodata from current_repodata.json, will retry with next repodata source.
    Collecting package metadata (repodata.json): done
    Solving environment: failed with initial frozen solve. Retrying with flexible solve.
    Solving environment: / 
    Found conflicts! Looking for incompatible packages.
    This can take several minutes.  Press CTRL-C to abort.
    failed                                                                                                                                                                                                              

    UnsatisfiableError: The following specifications were found to be incompatible with each other:

    Output in format: Requested package -> Available versions
    ```
    error2:
    ```
    Collecting package metadata (current_repodata.json): done
    Solving environment: failed with repodata from current_repodata.json, will retry with next repodata source.
    Collecting package metadata (repodata.json): done
    Solving environment: failed

    PackagesNotFoundError: The following packages are not available from current channels:

    - ngspipedb_py36_conda_env.yaml

    Current channels:

    - https://conda.anaconda.org/bioconda/osx-64
    - https://conda.anaconda.org/bioconda/noarch
    - https://conda.anaconda.org/conda-forge/osx-64
    - https://conda.anaconda.org/conda-forge/noarch
    - https://repo.anaconda.com/pkgs/main/osx-64
    - https://repo.anaconda.com/pkgs/main/noarch
    - https://repo.anaconda.com/pkgs/r/osx-64
    - https://repo.anaconda.com/pkgs/r/noarch

    To search for alternate channels that may provide the conda package you're
    looking for, navigate to

        https://anaconda.org

    and use the search bar at the top of the page.
    ```
    Solution2:
    ```
    conda update conda
    conda update anaconda
    ```

## Exporting an environment file across platforms/conda环境分享

```shell
cd NGSPipeDB_source_code
# on linux
conda env export --no-builds -p ./ngspipedb_py36_conda_env >ngspipedb_py36_conda_env.yaml
# on mac os
conda install conda-devenv -c conda-forge
# 如果上一步安装出错
# conda update --force conda
conda create -f ngspipedb_py36_conda_env.yaml -p ./ngspipedb_py36_conda_env
conda devenv -f ngspipedb_py36_conda_env.yaml -n ngs
```

#### package conda env

```shell
# on source 
conda install -c conda-forge conda-pack
source activate ./ngspipedb_py36_conda_env
conda pack -p ./ngspipedb_py36_conda_env -j 2 -o ngspipedb_py36_conda_env_osx64.tar.gz
# target
mkdir -p ngspipedb_py36_conda_env
tar -xzf ngspipedb_py36_conda_env.tar.gz -C ngspipedb_py36_conda_env
conda-unpack
```

conda Solving environment

```shell
conda install mamba -c conda-forge
mamba install snakemake=5.27.4
conda config --set channel_priority strict
```

conda update env

source activate ./ngspipedb_py38_conda_env/
mamba env update --prefix ./ngspipedb_py38_conda_env/ --file req.yaml  --prune
