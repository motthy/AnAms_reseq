#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=4G,s_vmem=4G
#$ -l d_rt=192:00:00
#$ -l s_rt=192:00:00

## set bwa path

mkdir ref index
cd ref

echo Downloading reference files from Cats-I...
curl -LO https://cat.annotation.jp/download/AnAms1.0/AnAms1.0.genome.fa.gz
echo Done.

echo Uncompressing files...
gunzip *.gz
echo Done.

cd ../
echo Making index by bwa
/your/bwa/path/bwa index -p index/AnAms1.0.genome.fa ref/AnAms1.0.genome.fa
echo Done.
