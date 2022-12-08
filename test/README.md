# How to run the cat_reseq pipeline

## Download reference files and make index for bwa

```bash
% ./download_reference_and_makeindex.sh

```

Genome Fasta is downloaded to `ref` directory.
BWA index are into `index` directory.

## Run

- Help  

```bash
% cwltool ../workflows/cat_reseq-pipeline.cwl -h
  ```

- Run

  All the output files will be generated into the directory specified with `--outdir`. If not specified, files will be generated into the current directory.specifying a `job.yaml` file

```bash
% cwltool --outdir test_out ../workflows/cat_reseq-pipeline.cwl cat_reseq-pipeline.test.job.yaml
```

- Run with singularity (NIG-SC)  
  Use `--singularity`

```bash
% cwltool --singularity --outdir test_out ../workflows/cat_reseq-pipeline.cwl cat_reseq-pipeline.test.job.yaml
```

- For debugging (Use `cachedir` to keep cache files and resume the job)

```bash
%  cwltool --cachedir test_cache --outdir test_out ../workflows/cat_reseq-pipeline.cwl cat_reseq-pipeline.test.job.yaml  
```
  
