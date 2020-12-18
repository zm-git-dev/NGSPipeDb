# Download the file of interest (here using a loop)
# Note that it can be interesting to store the STDERR of wget.
for i in control_R1 control_R2 treated_R1 treated_R2;
  do wget --no-clobber http://www.liu-lab.com/ngspipedb/Testdata/${i}.fq.gz 2>/dev/null;
done

# Download chr19 sequence (mm10 version)
wget --no-clobber http://www.liu-lab.com/ngspipedb/Testdata/chr19.fa.gz 2> chr19.fa_wget.log
gunzip chr19.fa.gz # uncompress the file

curl http://www.liu-lab.com/ngspipedb/Testdata/GRCm38.83.chr19.gtf.gz | gunzip -c > GRCm38.83.chr19.gtf


# generate replicate sample, need install seqkit first, conda install seqkit -c bioconda
# python generate_replicat.py control_R1.fq.gz treated_R1.fq.gz control_R2.fq.gz treated_R2.fq.gz