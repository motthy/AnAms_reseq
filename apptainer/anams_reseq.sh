#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=10G,s_vmem=10G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00

# 1. readのQCとトリミング
bash read_preprocessing.sh

# 2. mapping w/BWA-MEM
# realignment
# marking PCR duplicate
bash fastq2bam.sh

# 3. SNV calling w/GATK Haplotypecaller
bash bam2gvcf.sh

# 4. Hard filltering suggested by GATK developers
# Allele Balance filtering w/GATK VariantAnnotator
# remove indel > 6bp
bash filt_snp.sh
bash filt_indel.sh

# 5. Annotation w/snpeff (optional)
## Would you use Variant Reflector@Cats-I ?
## https://cat.annotation.jp/tools/variant_reflector/

