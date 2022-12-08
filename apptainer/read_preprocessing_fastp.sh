#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=8G,s_vmem=8G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00
#$ -V

## DDBJの/usr/local/resources/dra/fastqから目的のfastqを取得してtrimmingとQCを行う ###

# isolate and SRR acccession id
ISOLATE=24898
SRR_ACC=SRR1927133

# working directory
mkdir $ISOLATE
WORKDIR=/home/mhiromi/rhesus_macaque/ms/${ISOLATE}

# fastq
FQ1=/ddbj_database/dra/fastq/SRA248/SRA248488/SRX966727/SRR1927133_1.fastq.bz2
FQ2=/ddbj_database/dra/fastq/SRA248/SRA248488/SRX966727/SRR1927133_2.fastq.bz2
FQDIR=$(dirname $FQ1)

# singularity container
FASTQC=/usr/local/biotools/f/fastqc:0.11.9--hdfd78af_1
SEQKIT=/usr/local/biotools/s/seqkit:2.0.0--h9ee0642_0
FASTP=/usr/local/biotools/f/fastp:0.23.1--h79da9fb_0

# working directory
cd $ISOLATE

# 1.1 raw readsのstats
singularity exec -B $FQDIR -B $WORKDIR $SEQKIT seqkit stats $FQ1 $FQ2 -o read-stats-raw.tsv

# 1.2 raw readのqc
singularity exec -B $FQDIR -B $WORKDIR $FASTQC fastqc $FQ1 --nogroup -o .
singularity exec -B $FQDIR -B $WORKDIR $FASTQC fastqc $FQ2 --nogroup -o .

# 1.3 readsのtrimmingとQC
bunzip2 -c $FQ1 > ${SRR_ACC}_1.fq
bunzip2 -c $FQ2 > ${SRR_ACC}_2.fq

singularity exec -B $WORKDIR $FASTP fastp \
    -i ${SRR_ACC}_1.fq -o ${SRR_ACC}_1.trim_fq.gz \
    -I ${SRR_ACC}_2.fq -O ${SRR_ACC}_2.trim_fq.gz \
    -h ${SRR_ACC}_report.html \
    -j ${SRR_ACC}_report.json \
	-w 8
               
# 1.4 trimmed readsのstats
singularity exec -B $WORKDIR $SEQKIT seqkit stats ${SRR_ACC}_1.trim.fq.gz ${SRR_ACC}_2.trim.fq.gz -o read-stats.tsv

# 1.5 trimmed readのqc
singularity exec -B $WORKDIR $FASTQC fastqc ${SRR_ACC}_1.trim.fq.gz --nogroup -o .
singularity exec -B $WORKDIR $FASTQC fastqc ${SRR_ACC}_2.trim.fq.gz --nogroup -o .

# remove uncompressed fastq
#rm ${SRR_ACC}_1.fq
#rm ${SRR_ACC}_2.fq
