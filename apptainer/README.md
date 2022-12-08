# Cat genome AnAms1.0 variant call pipeline

This workflow is forked from  **rice_reseq** https://github.com/nigyta/rice_reseq

Modified to allow use of Apptainer

## usage
1. Set prefix (e.g., sample id) the valuable ISOLATE in each script
2. Prepair reference via qsub
```
qsub prepair_reference.sh
```
3. Execute wrapper script via qsub
```
qsub anams_reseq.sh
```

## Workflow

1. Read preprocessing (read_preprocessing.sh)
    1.1 FASTQ stats for raw reads (seqkit stats) __[stats report]__
    1.2 Quality check for raw reads (FastQC) __[Report HTML]__
    1.3 Adapter trimming and read QC (Trimmomatic) __[summary]__
    1.4 FASTQ stats for preprocessed reads (seqkit stats) __[stats report]__
    1.5 Read quality check (FastQC) __[Report HTML]__
2. Read mapping and BAM conversion (fastq2bam.sh)
    2.1 Read mapping (BWA)
    2.2 Convert SAM to BAM (samtools)
    2.3 Sort BAM (samtools)
    2.4 Create unmapped BAM (Picard FastqToSam)
    2.5 Merge mapped and unmapped BAM (Picard Merge)
    2.6 Remove duplicated reads (Picard MarkDuplicate) __[BAM, metrics file]__
    2.7 Create BAM index (samtools index) __[BAM index (.bai)]__
    2.8 Statistics for de-duplicated BAM (samtools stats) __[stats.txt]__
    2.9 Extract unmapped reads from BAM (samtools, picard) __[FASTQ]__
3. Variant calling, genotyping (bam2gvcf.sh)
    3.1 Variant calling (GATK HaplotypeCaller) __[gVCF, tbi]__
    3.2 Genotyping (GATK GenotypeGVCFs)
4. Filtering (filt_indel.sh, filt_snp.sh)
    3.3 Filtering for indel (GATK VariantFiltration)__[VCF, tbi]__
        Filtering condition: "QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0 || SOR > 10.0"
    3.4 Filtering for SNP (GATK VariantFiltration)__[VCF, tbi]__
        Filtering condition: "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    3.5 Allele Balance filtering for SNP and indel (Picard FilterVcf) __[VCF, tbi]__
        Heterozygous calls (ABHet=ref/(ref+alt)) ABHet < 0.2 or ABHet > 0.8 were removed
5. Variant Annotation
    Would you use [Variant Reflector@Cats-I](https://cat.annotation.jp/tools/variant_reflector/)?

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
| variants_{outprefix}.indel.filt.vcf.gz | VCF file genotyped, filtered, and selected for variants | 3.3 |
| variants_{outprefix}.indel.filt.vcf.gz.tbi | Tab index for the file above | 3.3 |
| variants_{outprefix}.snp.filt.vcf.gz | VCF file genotyped, filtered, and selected for variants | 3.4 |
| variants_{outprefix}.snp.filt.vcf.gz.tbi | Tab index for the file above | 3.4 |
| variants_{outprefix}.indel.ABHet.filt.vcf.gz | VCF file Allele Balance filtering for variants | 3.5 |
| variants_{outprefix}.indel.ABHet.filt.vcf.gz.tbi | Tab index for the file above | 3.5 |
| variants_{outprefix}.snp.ABHet.filt.vcf.gz | VCF file Allele Balance filtering for variants | 3.5 |
| variants_{outprefix}.snp.ABHet.filt.vcf.gz.tbi | Tab index for the file above | 3.5 |


