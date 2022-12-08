#$ -S /bin/bash                                                                                                                                   
#$ -pe def_slot 4                                                                                                                                 
#$ -cwd                                                                                                                                           
#$ -l s_vmem=8G,mem_req=8G


cwltool --singularity --outdir result ../workflows/fastq2bam.cwl fastq2bam.job.yaml
