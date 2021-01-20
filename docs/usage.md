# Run snake-pipe-rna


## 1. install miniconda

- First, download miniconda from url (https://docs.conda.io/en/latest/miniconda.html)
- Then, install follow the guide (https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

For example, if your system is linux, use the following command:
```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh && rm -f bash Miniconda3-latest-Linux-x86_64.sh
```

## 2. create virtual environment

We suggest user to use snake-pipe-rna in a virtual conda environment.
```shell
# create virtual environment named rnaseq
conda create -n rnaseq python=3.6
# activate your environmnet
source activate rnaseq
# if your want exit this environmnet
conda deactivate
```

## 3. install snakemake

```shell
# using conda install
conda install -c bioconda -c conda-forge snakemake=5.5.4
```

## 4. download snake-pipe-rna

```shell
git clone https://github.com/xuanblo/snake-pipe-rna.git
```

## 5. run snake-pipe-rna
```
#test Snakefile
snakemake -np
#show DAG
snakemake --dag | dot -Tpdf > dag.pdf
```

### 5.1 mRNA-seq

```shell
# Download the file of interest (here using a loop)
# Note that it can be interesting to store the STDERR of wget.
for i in control_R1 control_R2 treated_R1 treated_R2;
  do wget --no-clobber http://pedagogix-tagc.univ-mrs.fr/courses/data/ngs/abd/${i}.fq.gz 2>/dev/null;
done

# Download chr19 sequence (mm10 version)
wget --no-clobber http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr19.fa.gz 2> chr19.fa_wget.log
gunzip chr19.fa.gz # uncompress the file

curl ftp://ftp.ensembl.org/pub/release-83/gtf/mus_musculus/Mus_musculus.GRCm38.83.chr.gtf.gz | gunzip -c | grep "^19" | sed 's/^19/chr19/' > GRCm38.83.chr19.gtf
```

### 5.2 lncRNA-seq

```shell
snakemake --use-conda -j 40
```

### 5.3 miRNA-seq

## 6. run auto-rna-server
