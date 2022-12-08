#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=10G,s_vmem=10G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00

ISOLATE=<set sample id>
FQ1=<set R1_fastq>
FQ2=<set R2_fastq>

# 1. readのQCとトリミング
bash read_preprocessing_fastp.sh $ISOLATE $FQ1 $FQ2

# 2. mapping w/BWA-MEM
# realignment
# marking PCR duplicate
bash fastq2bam.sh $ISOLATE

# 3. SNV calling w/GATK Haplotypecaller
bash bam2gvcf.sh $ISOLATE

# 4. Hard filltering suggested by GATK developers
# Allele Balance filtering w/GATK VariantAnnotator
# remove indel > 6bp
bash filt_snp.sh $ISOLATE
bash filt_indel.sh $ISOLATE

# 5. Annotation w/snpeff (optional)
## Would you use Variant Reflector@Cats-I ?
## https://cat.annotation.jp/tools/variant_reflector/

