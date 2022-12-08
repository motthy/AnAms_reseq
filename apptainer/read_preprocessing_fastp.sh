#$ -S /bin/bash
#$ -pe def_slot 8
#$ -cwd
#$ -l mem_req=8G,s_vmem=8G
#$ -l d_rt=720:00:00
#$ -l s_rt=720:00:00
#$ -V

## QC and trimming fastq ##

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
FASTP=/usr/local/biotools/f/fastp:0.23.1--h79da9fb_0

# working directory
cd $ISOLATE

# 1.1 raw readsのstats
apptainer exec --no-mount tmp -B $FQDIR -B $WORKDIR $SEQKIT seqkit stats $FQ1 $FQ2 -o read-stats-raw.tsv

# 1.2 raw readのqc
apptainer exec --no-mount tmp -B $FQDIR -B $WORKDIR $FASTQC fastqc $FQ1 --nogroup -o .
apptainer exec --no-mount tmp -B $FQDIR -B $WORKDIR $FASTQC fastqc $FQ2 --nogroup -o .

# 1.3 readsのtrimmingとQC
# bunzip2 -c $FQ1 > ${ISOLATE}_1.fq
# bunzip2 -c $FQ2 > ${ISOLATE}_2.fq

apptainer exec --no-mount tmp -B $WORKDIR $FASTP fastp \
    -i $FQ1 -o ${ISOLATE}_1.trim_fq.gz \
    -I $FQ2 -O ${ISOLATE}_2.trim_fq.gz \
    -h ${ISOLATE}_report.html \
    -j ${ISOLATE}_report.json \
	-w 8

# 1.4 trimmed readsのstats
apptainer exec --no-mount tmp -B $WORKDIR $SEQKIT seqkit stats ${ISOLATE}_1.trim.fq.gz ${ISOLATE}_2.trim.fq.gz -o read-stats.tsv

# 1.5 trimmed readのqc
apptainer exec --no-mount tmp -B $WORKDIR $FASTQC fastqc ${ISOLATE}_1.trim.fq.gz --nogroup -o .
apptainer exec --no-mount tmp -B $WORKDIR $FASTQC fastqc ${ISOLATE}_2.trim.fq.gz --nogroup -o .

