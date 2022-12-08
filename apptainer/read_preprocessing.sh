#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=8G,s_vmem=8G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00
#$ -V

## DDBJの/usr/local/resources/dra/fastqから目的のfastqを取得してtrimmingとQCを行う ###

# isolate id
ISOLATE=<set sample name>

# fastq
FQ1=<set R1_fastq>
FQ2=<set R1_fastq>
FQDIR=$(dirname $FQ1)

# working directory
mkdir $ISOLATE
WORKDIR=$(pwd)/${ISOLATE}



# singularity container
FASTQC=/usr/local/biotools/f/fastqc:0.11.9--hdfd78af_1
SEQKIT=/usr/local/biotools/s/seqkit:2.0.0--h9ee0642_0
TRIMMOMATIC=/usr/local/biotools/t/trimmomatic:0.39--hdfd78af_2

# working directory
cd $ISOLATE

# 1.1 raw readsのstats
apptainer exec --no-mount tmp -B $FQDIR -B $WORKDIR $SEQKIT seqkit stats $FQ1 $FQ2 -o read-stats-raw.tsv

# 1.2 raw readのqc
apptainer exec --no-mount tmp -B $FQDIR -B $WORKDIR $FASTQC fastqc $FQ1 --nogroup -o .
apptainer exec --no-mount tmp -B $FQDIR -B $WORKDIR $FASTQC fastqc $FQ2 --nogroup -o .

# 1.3 readsのtrimmingとQC
# bunzip2 -c $FQ1 > ${SRR_ACC}_1.fq
# bunzip2 -c $FQ2 > ${SRR_ACC}_2.fq

apptainer exec --no-mount tmp -B $WORKDIR $TRIMMOMATIC trimmomatic PE \
  -threads 8 \
  -phred33 \
  -trimlog out.trimmomatic.log.txt \
  -summary out.trimmomatic.summary.txt \
  $FQ1 \
  $FQ2 \
  ${ISOLATE}_1.trim.fq.gz \
  ${ISOLATE}_1.unpaired.fq.gz \
  ${ISOLATE}_2.trim.fq.gz \
  ${ISOLATE}_2.unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:20 \
  TRAILING:20 \
  SLIDINGWINDOW:10:20 \
  MINLEN:30

# 1.4 trimmed readsのstats
apptainer exec --no-mount tmp -B $WORKDIR $SEQKIT seqkit stats ${ISOLATE}_1.trim.fq.gz ${ISOLATE}_2.trim.fq.gz -o read-stats.tsv

# 1.5 trimmed readのqc
apptainer exec --no-mount tmp -B $WORKDIR $FASTQC fastqc ${ISOLATE}_1.trim.fq.gz --nogroup -o .
apptainer exec --no-mount tmp -B $WORKDIR $FASTQC fastqc ${ISOLATE}_2.trim.fq.gz --nogroup -o .

