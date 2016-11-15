RNAseq-pipeline
===============

1. Make a directory to deposit subfolders and processed files. `$DIR`

2. Often, `$DIR/FASTQ` contains the fastq files for this experiment.

2. In $DIR, make text file called `SAMPLE_IDs.txt` with FASTQ file prefixes delimited with newline.
  
  Sample file possibilites for prefix `SAMPLE`: 
  
  Paired-end files: `SAMPLE_R1_001.fastq.gz` and `SAMPLE_R2_001.fastq.gz`
 
  Singe-end files : `SAMPLE.fastq.gz`
  
  *fq.gz files also supported*

3. Run pipeline with following call in `$DIR`:

  `/DIR/TO/rnaseq_run_all.sh ExprName SomeIdentifier Genome PE_Bool Stranded_Bool ERCC_Bool $DIR/FASTQ`

  ExprName: main identifier, experiment name
  SomeIdentifier: spurious additional identifier, date is a good thing use here
  Genome: supported genomes are `hg19, hg38, mm10, rn6`
  PE_Bool: `TRUE` If samples paired-ended
  Stranded_Bool: `TRUE` If samples are reverse-stranded (forward-stranded not compatible yet)
  ERCC_Bool: `TRUE` If ERCC mix 1 was added

4. Hidden run files will be created with all calls needed to run the pipeline. Each step is submitted to SGE cluster and queued to run sequentially. Steps that can be parallelized are submitted as array jobs.
