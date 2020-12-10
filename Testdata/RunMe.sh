# Download the file of interest (here using a loop)
# Note that it can be interesting to store the STDERR of wget.
for i in control_R1 control_R2 treated_R1 treated_R2;
  do wget --no-clobber http://pedagogix-tagc.univ-mrs.fr/courses/data/ngs/abd/${i}.fq.gz 2>/dev/null;
done

# Download chr19 sequence (mm10 version)
wget --no-clobber http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr19.fa.gz 2> chr19.fa_wget.log
gunzip chr19.fa.gz # uncompress the file

curl ftp://ftp.ensembl.org/pub/release-83/gtf/mus_musculus/Mus_musculus.GRCm38.83.chr.gtf.gz | gunzip -c | grep "^19" | sed 's/^19/chr19/' > GRCm38.83.chr19.gtf


# 访问我们组的网站下载 
python generate_replicat.py control_R1.fq.gz treated_R1.fq.gz control_R2.fq.gz treated_R2.fq.gz
