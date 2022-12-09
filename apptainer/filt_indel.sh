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

# select INDEL < 6bp
apptainer exec $GATK gatk \
                    --java-options "-Xmx4G" SelectVariants \
                    --variant variants_${ISOLATE}.genotype.g.vcf.gz \
                    --output variants_${ISOLATE}.indel.vcf.gz \
                    --reference $REF \
                    --select-type-to-include "INDEL" \
                    --max-indel-size 5

## index
apptainer exec $GATK gatk IndexFeatureFile \
     --feature-file variants_${ISOLATE}.indel.vcf.gz  \
     --output variants_${ISOLATE}.indel.vcf.gz.tbi


# GATK VariantFiltration
## https://gatk.broadinstitute.org/hc/en-us/articles/360035532412?id=11097
## for indel
#singularity exec -B $WORKDIR $GATK gatk \
#     --java-options "-Xmx4G" VariantFiltration \
#     --variant variants_${SRR_ACC}.indel.vcf.gz \
#     --reference $REF \
#     --output variants_${SRR_ACC}.indel.filt.vcf.gz \
#     -filter "QD < 2.0" --filter-name "QD2" \
#     -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum20" \
#     -filter "InbreedingCoeff < -0.8" --filter-name "InbreedingCoeff-0.8" \
#     -filter "FS > 200.0" --filter-name "FS200" \
#     -filter "SOR > 10.0" --filter-name "SOR10"

apptainer exec $GATK gatk \
     --java-options "-Xmx4G" VariantFiltration \
     --variant variants_${ISOLATE}.indel.vcf.gz \
     --reference $REF \
     --output variants_${ISOLATE}.indel.filt.vcf.gz \
     --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0 || SOR > 10.0 ' \
     --filter-name "indel_hard_filtering"

# index
apptainer exec $GATK gatk --java-options "-Xmx4G" IndexFeatureFile \
     --feature-file variants_${ISOLATE}.indel.filt.vcf.gz  \
     --output variants_${ISOLATE}.indel.filt.vcf.gz.tbi

# Allele Balance filtering w/GATK VariantAnnotator
## heterozygous calls (ABHet=ref/(ref+alt)) ABHet < 0.2 or ABHet > 0.8 were removed
apptainer exec $PICARD \
             java -jar build/libs/picard.jar FilterVcf \
             INPUT=variants_${ISOLATE}.indel.filt.vcf.gz \
             OUTPUT=variants_${ISOLATE}.indel.ABHet.filt.vcf.gz \
             MIN_AB=0.2

# index
apptainer exec $GATK gatk --java-options "-Xmx4G" IndexFeatureFile \
     --feature-file variants_${ISOLATE}.indel.ABHet.filt.vcf.gz \
     --output variants_${ISOLATE}.indel.ABHet.filt.vcf.gz.tbi

