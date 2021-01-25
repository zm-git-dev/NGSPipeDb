## 使用教程

conda create -n ngspipe_chipseq 
conda activate ngspipe_chipseq
conda install mamba -c conda-forge

mamba env update --file ngspipe/envs/requirements_chipseq.yaml --prune

snakemake -s ngspipe/ChIP_seq.Snakefile.py --configfile ngspipe/config/chipseq.config.yaml -np