# Cat genome AnAms1.0 variant call pipeline

This workflow is forked from  **rice_reseq** https://github.com/nigyta/rice_reseq

## Workflow

1. Read preprocessing (read_preprocessing.cwl)  
    1.1 FASTQ stats for raw reads (seqkit stats) __[stats report]__  
    1.2 Quality check for raw reads (FastQC) __[Report HTML]__  
    1.3 Adapter trimming and read QC (Trimmomatic) __[summary]__  
    1.4 FASTQ stats for preprocessed reads (seqkit stats) __[stats report]__  
    1.5 Read quality check (FastQC) __[Report HTML]__  
2. Read mapping and BAM conversion (fastq2bam.cwl)  
    2.1 Read mapping (BWA)  
    2.2 Convert SAM to BAM (samtools)  
    2.3 Sort BAM (samtools)  
    2.4 Create unmapped BAM (Picard FastqToSam)  
    2.5 Merge mapped and unmapped BAM (Picard Merge)  
    2.6 Remove duplicated reads (Picard MarkDuplicate) __[BAM, metrics file]__  
    2.7 Create BAM index (samtools index) __[BAM index (.bai)]__  
    2.8 Statistics for de-duplicated BAM (samtools stats) __[stats.txt]__  
    2.9 Extract unmapped reads from BAM (samtools, picard) __[FASTQ]__  
3. Variant calling, genotyping, filtering (bam2vcf.cwl)  
    3.1 Variant calling (GATK HaplotypeCaller) __[gVCF, tbi]__  
    3.2 Genotyping (GATK GenotypeGVCFs)  
    3.3 Filtering (GATK VariantFiltration)   
        Filtering condition: "QD < 5.0 || FS > 50.0 || SOR > 3.0 || MQ < 50.0 || MQRankSum < -2.5 || ReadPosRankSum < -1.0 || ReadPosRankSum > 3.5"  
    3.4 Variant selection for SNP and INDEL (GATK SelectVariants) __[VCF, tbi]__  
4. Variant Annotation (snpeff_all.cwl)  
    4.1 Build SnpEff database from reference files (SnpEff build)  
    4.2 Annotate variants (SnpEff) __[VCF, tbi, summary, effected genes]__  

## Output
|  File name  |  Description  | Step |
| ---- | ---- | --- |
| {outprefix}_read-stats-raw.tsv | Stats file for raw FASTQ files | 1.1 |
| {fastq1/2}_fastqc.html | FastQC report for raw reads | 1.2 |
| {outprefix}.trimmomatic.summary.txt | Summary of Trimmomatic result | 1.3 |
| {outprefix}_read-stats.tsv | Stats file for preprocessed FASTQ files | 1.4 |
| {fastq1/2}.trimmomatic-pe_fastqc.html | FastQC report for reads processed using Trimmomatoc | 1.5 |
| {outprefix}.rmdup.bam | Alignment result in BAM format, after de-duplication | 2.6 |
| {outprefix}.rmdup.bam.bai | BAM index for the file above | 2.7 |
| {outprefix}.rmdup.metrics | Metrics file for BAM de-duplication | 2.6 |
| {outprefix}.rmdup.bam.stats.txt | Statistics for de-duplicated BAM (generated using samtools-depth, in-house script) | 2.8 |
| {outprefix}.unmapped-read.{r1/r2}.fastq.gz | Reads not mapped to the reference genome | 2.9 |
| variants_{outprefix}.g.vcf.gz | gVCF file generated from GATK-HapolotypeCaller | 3.1 |
| variants_{outprefix}.g.vcf.gz.tbi | Tab index for the file above | 3.1 |
| variants_{outprefix}.varonly.vcf.gz | VCF file genotyped, filtered, and selected for variants | 3.4 |
| variants_{outprefix}.varonly.vcf.gz.tbi | Tab index for the file above | 3.4 |
| variants_{outprefix}.snpEff.vcf.gz | VCF file annotated using snpEFff | 4.2 |
| variants_{outprefix}.snpEff.vcf.gz.tbi | Tab index for the file above | 4.2 |
| snpEff_genes_{outprefix}.txt | snpEff result for effected genes| 4.2 |
| snpEff_summary_{outprefix}.html | Summary of snpEff | 4.2 |


