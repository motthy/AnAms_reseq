#!/bin/bash
#$ -S bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=20G,s_vmem=20G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00
#$ -V

set -eu

# isolate id
ISOLATE=$1

# reference
REF=$(pwd)/ref/AnAms1.0.genome.fa

# apptainer container
GATK=/usr/local/biotools/g/gatk4:4.2.2.0--hdfd78af_1
PICARD=/usr/local/biotools/p/picard:2.26.4--hdfd78af_0

# working directory
cd $ISOLATE

# select SNP
apptainer exec $GATK gatk --java-options "-Xmx4G" SelectVariants \
                    --variant variants_${ISOLATE}.genotype.g.vcf.gz \
                    --output variants_${ISOLATE}.snp.vcf.gz \
                    --reference $REF \
                    --select-type-to-include "SNP"

apptainer exec $GATK gatk --java-options "-Xmx4G" IndexFeatureFile \
     --feature-file variants_${ISOLATE}.snp.vcf.gz \
     --output variants_${ISOLATE}.snp.vcf.gz.tbi

# GATK VariantFiltration
## https://gatk.broadinstitute.org/hc/en-us/articles/360035532412?id=11097
## for SNP
#singularity exec -B $WORKDIR $GATK gatk --java-options "-Xmx4G" VariantFiltration \
#     --variant variants_${SRR_ACC}.snp.vcf.gz \
#     --reference $REF \
#     --output variants_${SRR_ACC}.snp.filt.vcf.gz \
#     -filter "QD < 2.0" --filter-name "QD2" \
#     -filter "MQ < 40.0" --filter-name "MQ40" \
#     -filter "FS > 60.0" --filter-name "FS60" \
#     -filter "SOR > 3.0" --filter-name "SOR3" \
#     -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#     -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

apptainer exec $GATK gatk --java-options "-Xmx4G" VariantFiltration \
     --variant variants_${ISOLATE}.snp.vcf.gz \
     --reference $REF \
     --output variants_${ISOLATE}.snp.filt.vcf.gz \
     --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
     --filter-name "snv_hard_filtering"

## index
apptainer exec $GATK gatk --java-options "-Xmx4G" IndexFeatureFile \
     --feature-file variants_${ISOLATE}.snp.filt.vcf.gz \
     --output variants_${ISOLATE}.snp.filt.vcf.gz.tbi

# Allele Balance filtering
## heterozygous calls (ABHet=ref/(ref+alt)) ABHet < 0.2 or ABHet > 0.8 were removed
apptainer exec $PICARD \
             java -jar build/libs/picard.jar FilterVcf \
             INPUT=variants_${ISOLATE}.snp.filt.vcf.gz \
             OUTPUT=variants_${ISOLATE}.snp.ABHet.filt.vcf.gz \
             MIN_AB=0.2
## index
apptainer exec $GATK gatk --java-options "-Xmx4G" IndexFeatureFile \
     --feature-file variants_${ISOLATE}.snp.ABHet.filt.vcf.gz \
     --output variants_${ISOLATE}.snp.ABHet.filt.vcf.gz.tbi
