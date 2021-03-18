# resequcing analysis

## step-by-step run resequencing on testdata

### 1. pre-prepare download source code

1. [Install wget and git](../NGSPipe-RNA-seq/#BasicLinux)
2. [Install Miniconda3](../NGSPipe-RNA-seq/#Miniconda)
3. [Download NGSPipeDb source code](../NGSPipe-RNA-seq/#NGSPipeDbSource)

Modify the project name and enter the project directory.

```shell
mv NGSPipeDb species_sample_transcript_analysis_by_NGSPipeDb
cd species_sample_transcript_analysis_by_NGSPipeDb
```

### 2. create conda envirenment

```shell
mamba create -n ngspipe-resequencing python=3.8 -c conda-forge -y
mamba env update -n ngspipe-resequencing --file ngspipe/envs/requirements_resequencing.yaml --prune
conda activate ngspipe-resequecing
```
  - snpeff
  - gatk


### 3. download testdata

```
mkdir -p testdata && cd testdata
wget http://www.liu-lab.com/ngspipedb/testdata/control_R1.fq.gz
wget http://www.liu-lab.com/ngspipedb/testdata/treated_R1.fq.gz
wget http://www.liu-lab.com/ngspipedb/testdata/control_R2.fq.gz
wget http://www.liu-lab.com/ngspipedb/testdata/treated_R2.fq.gz
wget http://www.liu-lab.com/ngspipedb/testdata/chr19.fa.gz
wget http://www.liu-lab.com/ngspipedb/testdata/GRCm38.83.chr19.gtf.gz
gunzip chr19.fa.gz
gunzip GRCm38.83.chr19.gtf.gz
wget http://www.liu-lab.com/ngspipedb/testdata/sample_resequecing.xls
cd ..
```
