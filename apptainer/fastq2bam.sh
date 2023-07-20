#!/bin/bash
#$ -S bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=20G,s_vmem=20G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00
#$ -V

## mapping ##

set -eu

# isolate id
ISOLATE=$1

# reference
REF=$(pwd)/ref/AnAms1.0.genome.fa

# apptainer container
BWA=/usr/local/biotools/b/bwa:0.7.17--pl5.22.0_2
SAMTOOLS=/usr/local/biotools/s/samtools:0.1.19--hfb9b9cc_8
PICARD=/usr/local/biotools/p/picard:2.26.4--hdfd78af_0

# working directory
cd $ISOLATE

# mapping w/BWA-MEM
apptainer exec $BWA bwa mem -M \
             $REF \
             ${ISOLATE}_1.trim.fq.gz \
             ${ISOLATE}_2.trim.fq.gz > ${}.sam

# samをbamに変換
apptainer exec $SAMTOOLS \
             samtools view -b -@ 8 ${ISOLATE}.sam > ${ISOLATE}.bam

# bamのソート
apptainer exec $SAMTOOLS \
             samtools sort -@ 8 -m 10G ${ISOLATE}.bam > ${ISOLATE}.sort.bam

# unalignd bam
apptainer exec $PICARD \
             picard FastqToSam \
             FASTQ=${ISOLATE}_1.trim.fq.gz \
             FASTQ2=${ISOLATE}_2.trim.fq.gz \
             OUTPUT=${ISOLATE}.uBAM.bam

# mapped bamとunaligned bamのマージ (realignment)
apptainer exec $PICARD \
             picard MergeBamAlignment \
             ALIGNED_BAM=${ISOLATE}.sort.bam \
             UNMAPPED_BAM=${ISOLATE}.uBAM.bam \
             REFERENCE_SEQUENCE=$REF \
             OUTPUT=${ISOLATE}.merge.bam

# unmapped bam
apptainer exec $SAMTOOLS \
             samtools view -f 4 -@ 8 ${ISOLATE}.merge.bam \
             -o ${ISOLATE}_unmapped.bam

# unmapped readsの抽出
apptainer exec $PICARD \
             picard SamToFastq \
             VALIDATION_STRINGENCY=SILENT \
             INPUT=${ISOLATE}_unmapped.bam \
             FASTQ=${ISOLATE}_unmapped.r1.fq \
             SECOND_END_FASTQ=${ISOLATE}_unmapped.r2.fq

gzip -c ${ISOLATE}_unmapped.r1.fq > ${ISOLATE}_unmapped.r1.fq.gz
gzip -c ${ISOLATE}_unmapped.r2.fq > ${ISOLATE}_unmapped.r2.fq.gz

# PCR duplicateの除去
apptainer exec $PICARD \
             picard MarkDuplicates \
             INPUT=${ISOLATE}.merge.bam \
             OUTPUT=${ISOLATE}.rmdup.bam \
             METRICS_FILE=${ISOLATE}.metrics \
             REMOVE_DUPLICATES=true \
             MAX_RECORDS_IN_RAM=1000000 \
             TMP_DIR=./tmp

# @RGの追記
apptainer exec $PICARD \
　　　　　　　　picard AddOrReplaceReadGroups \
　　　　　　　　INPUT=${ISOLATE}.rmdup.bam \
    　　　　　 OUTPUT=${ISOLATE}.rmdup.addRG.bam \
    　　　　　 RGID=${ISOLATE} \
             RGLB=${ISOLATE} \
             RGPU=${ISOLATE} \
             RGPL=illumina \
             RGSM=${ISOLATE} \
             VALIDATION_STRINGENCY=LENIENT

# bamのstatisticsとindex作成
apptainer exec $SAMTOOLS samtools stats ${ISOLATE}.rmdup.addRG.bam
apptainer exec $SAMTOOLS samtools index -@ 8 ${ISOLATE}.rmdup.addRG.bam

# remove unsorted sam/bam and uncompressed files
#rm ${ISOLATE}.sam
#rm ${ISOLATE}.bam
#rm ${ISOLATE}_unmapped.r1.fq
#rm ${ISOLATE}_unmapped.r2.fq
#rm ${ISOLATE}.merge.bam
