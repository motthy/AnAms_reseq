#!/bin/bash
#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=8G,s_vmem=8G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00

## make index for bwa ##

# apptainer container
BWA=/usr/local/biotools/b/bwa:0.7.17--pl5.22.0_2
SAMTOOLS=/usr/local/biotools/s/samtools:0.1.19--hfb9b9cc_8
PICARD=/usr/local/biotools/p/picard:2.26.4--hdfd78af_0

# ref directry
mkdir ref
cd ref

ln -s https://cat.annotation.jp/download/AnAms1.0/AnAms1.0.genome.fa.gz

gunzip AnAms1.0.genome.fa.gz

## make bwa index
apptainer exec $BWA bwa index -p AnAms1.0 AnAms1.0.genome.fa

## make GATK dict
apptainer exec $PICARD java -jar build/libs/picard.jar \
                         CreateSequenceDictionary AnAms1.0.genome.fa
## make faidx
apptainer exec $SAMTOOLS samtools faidx AnAms1.0.genome.fa
