RNAseq-pipeline
===============

This pipeline requires several R packages. You can install them by running:

```
qsub pipeline_setup.sh
```

1. Make a directory to deposit subfolders and processed files: `$DIR`.

  Often, `$DIR/FASTQ` contains the fastq files for this experiment but the fastq files don't necessarily have to be in this location.

2. In `$DIR`, make text file called `SAMPLE_IDs.txt` with FASTQ file prefixes delimited with newline.
  
  Sample file possibilities for prefix `SAMPLE`: 
  
  Paired-end files: __SAMPLE__`_R1_001.fastq.gz` and __SAMPLE__`_R2_001.fastq.gz`
 
  Single-end files : __SAMPLE__`.fastq.gz`
  
  *fq.gz files also supported*
  
  For example, a `SAMPLE_IDs.txt` file can look like this (these are paired-end samples):

  ```
  R10126_C1BP4ACXX_GAGATTCC_L005
  R10126_C1BP4ACXX_GAGATTCC_L006
  R10126_C1BP4ACXX_GAGATTCC_L007
  R10126_C1BP4ACXX_GAGATTCC_L008
  ```

  You can alternatively add the path to the files, which can be useful if you have the FASTQ files in different directories:

  ```
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L005
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L006
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10145_C1BM1ACXX/R10145_C1BM1ACXX_AGCGATAG_L005
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10145_C1BM1ACXX/R10145_C1BM1ACXX_AGCGATAG_L006
  ```
  
  If you reads are split in multiple files and you want to merge them specify in `SAMPLE_IDs.txt` a second column with the group identifiers and the boolean ${MERGE} in the next section. An example of such a `SAMPLE_IDs.txt` file would be
  
  ```
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L005   1
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L006   1
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10145_C1BM1ACXX/R10145_C1BM1ACXX_AGCGATAG_L005   2
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10145_C1BM1ACXX/R10145_C1BM1ACXX_AGCGATAG_L006   2
  ```

3. Run pipeline with following call in `$DIR`:

  ```
  sh /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq_run_all.sh $ExprName $SomeIdentifier $Genome $PE $Stranded $ERCC $FASTQ_DIR $MERGE $LARGE
  ```

  * __ExprName__: main identifier, experiment name
  * __SomeIdentifier__: spurious additional identifier, date is a good thing use here
  * __Genome__: supported genomes are `hg19, hg38, mm10, rn6`
  * __PE__: `TRUE` If samples paired-ended
  * __Stranded__: `TRUE` If samples are reverse-stranded (forward-stranded not compatible yet)
  * __ERCC__: `TRUE` If ERCC mix 1 was added
  * __FASTQ_DIR__: The path of the directory containing the FASTQ files. Use `""` if `SAMPLE_IDs.txt` already contains full paths (default: `""`).
  * __MERGE__: `TRUE` if you want to merge the files. The will be saved in a directory called `merged_fastq`. `FALSE` by default and doesn't need to be specified.
  * __LARGE__: `TRUE` if you want to use double the default memory settings. Useful for large projects (many samples and/or many reads). `FALSE` by default and doesn't need to be specified.

4. Hidden run files will be created with all calls needed to run the pipeline. Each step is submitted to SGE cluster and queued to run sequentially. Steps that can be parallelized are submitted as array jobs.

5. Completion emails: by default you will only get an email if a job failed. If you want to get completion emails add the empty file ${DIR}/.send_emails
