RNAseq-pipeline
===============

This pipeline requires several R packages. You can install them by running:

```
qsub pipeline_setup.sh
```

1. Make a directory to deposit subfolders and processed files. `$DIR`

2. Often, `$DIR/FASTQ` contains the fastq files for this experiment.

2. In $DIR, make text file called `SAMPLE_IDs.txt` with FASTQ file prefixes delimited with newline.
  
  Sample file possibilities for prefix `SAMPLE`: 
  
  Paired-end files: `SAMPLE_R1_001.fastq.gz` and `SAMPLE_R2_001.fastq.gz`
 
  Single-end files : `SAMPLE.fastq.gz`
  
  *fq.gz files also supported*
  
  If you reads are split in multiple files and you want to merge them specify in `SAMPLE_IDs.txt` a second column with the group identifiers and the boolean ${MERGE} in the next section. An example of such a `SAMPLE_IDs.txt` file would be
  
```
R10126_C1BP4ACXX_GAGATTCC_L005  1
R10126_C1BP4ACXX_GAGATTCC_L006  1
R10126_C1BP4ACXX_GAGATTCC_L007  2
R10126_C1BP4ACXX_GAGATTCC_L008  3
```

3. Run pipeline with following call in `$DIR`:

```
sh /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq_run_all.sh $ExprName $SomeIdentifier $Genome $PE_Bool $Stranded_Bool $ERCC_Bool $DIR/FASTQ $MERGE
```

  ExprName: main identifier, experiment name
  
  SomeIdentifier: spurious additional identifier, date is a good thing use here
  
  Genome: supported genomes are `hg19, hg38, mm10, rn6`
  
  PE_Bool: `TRUE` If samples paired-ended
  
  Stranded_Bool: `TRUE` If samples are reverse-stranded (forward-stranded not compatible yet)
  
  ERCC_Bool: `TRUE` If ERCC mix 1 was added
  
  MERGE: `TRUE` if you want to merge the files. The will be saved in a directory called `merged_fastq`. `FALSE` by default and doesn't need to be specified.

4. Hidden run files will be created with all calls needed to run the pipeline. Each step is submitted to SGE cluster and queued to run sequentially. Steps that can be parallelized are submitted as array jobs.
