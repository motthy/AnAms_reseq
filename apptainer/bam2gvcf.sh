#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=10G,s_vmem=10G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00

## variant calling ##

# isolate id
ISOLATE=$1

# reference
REF=$(pwd)/ref/AnAms1.0.genome.fa

# apptainer container
GATK=/usr/local/biotools/g/gatk4:4.2.2.0--hdfd78af_1

# working directory
cd $ISOLATE

# GATK Haplotypecaller
apptainer exec --no-mount tmp $GATK gatk --java-options "-Xmx4G" HaplotypeCaller \
                    --input ${ISOLATE}.rmdup.bam \
                    --output ${ISOLATE}.g.vcf.gz \
                    --reference $REF \
                    --emit-ref-confidence GVCF\
                    -max-alternate-alleles 2

apptainer exec --no-mount tmp $GATK gatk --java-options "-Xmx4G" IndexFeatureFile \
     --feature-file ${ISOLATE}.g.vcf.gz \
     --output ${ISORATE}.g.vcf.gz.tbi

# GATK GenotypeGVCFs
apptainer exec --no-mount tmp $GATK gatk --java-options "-Xmx4G" GenotypeGVCFs \
                    --variant variants_${ISOLATE}.g.vcf.gz \
                    --output variants_${ISOLATE}.genotype.g.vcf.gz \
                    --reference $REF

apptainer exec --no-mount tmp $GATK gatk --java-options "-Xmx4G" IndexFeatureFile \
     --feature-file variants_${ISOLATE}.genotype.g.vcf.gz  \
     --output variants_${ISOLATE}.genotype.g.vcf.gz.tbi

