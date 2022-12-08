#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=4G,s_vmem=4G
#$ -l d_rt=192:00:00
#$ -l s_rt=192:00:00

mkdir ref index
cd ref

echo Downloading reference files from UCSC...
curl -LO https://hgdownload-test.gi.ucsc.edu/goldenPath/felCat9/bigZips/felCat9.fa.gz
echo Done.

echo Uncompressing files...
gunzip *.gz
echo Done.

cd ../
echo Making index by bwa
/your/bwa/path/bwa index -p index/felCat9.fa ref/felCat9.fa
echo Done.
