RNAseq-pipeline
===============

This pipeline requires several R packages. You can install them by running:

```
qsub pipeline_R_setup.sh
```

1. Make a directory to deposit subfolders and processed files: `$DIR`.

  Often, `$DIR/FASTQ` contains the fastq files for this experiment but the fastq files don't necessarily have to be in this location.

2. In `$DIR`, make text file called `SAMPLE_IDs.txt` with FASTQ file prefixes delimited with newline.
  
  Sample file possibilities for prefix _SAMPLE_: 
  
  Paired-end files: _SAMPLE_`_R1_001.fastq.gz` and _SAMPLE_`_R2_001.fastq.gz`
 
  Single-end files : _SAMPLE_`.fastq.gz`
  
  The following extensions are supported: `fastq.gz`, `fq.gz`, `fastq`, `fq`.
  
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
  
  If you reads are split in multiple files and you want to merge them specify in `SAMPLE_IDs.txt` a second column with the group identifiers and the boolean `${MERGE}` in the next section. An example of such a `SAMPLE_IDs.txt` file would be
  
  ```
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L005   1
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L006   1
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10145_C1BM1ACXX/R10145_C1BM1ACXX_AGCGATAG_L005   2
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10145_C1BM1ACXX/R10145_C1BM1ACXX_AGCGATAG_L006   2
  ```

3. Run pipeline with following call in `$DIR`:

  ```
  ## You need a compute node to run this script!
  qrsh
  bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq_run_all.sh $ExprName $SomeIdentifier $Genome $PE $Stranded $ERCC $FASTQ_DIR $MERGE $LARGE
  ```
  
  Note that you have to use __bash__ and not __sh__, otherwise you'll get an error like this:
  
  ```
  module: command not found
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
  * __FULLCOV__: `TRUE` if you want to create the fullCoverage object. Set to `FALSE` by default. Note that the fullCoverage object is not needed (by default) since we create the normalized mean BigWig already and can use it to define the ERs with `derfinder::findRegions()`, then use the resulting GRanges object and the paths to the BigWig files in `derfinder::getRegionCoverage(fullCov = NULL, files = bigWigs, regions = outputFrom_findRegions)` or alternatively write the regions to a BED file with rtracklayer, create the counts with [bwtool](https://github.com/CRG-Barcelona/bwtool) and then read them into R manually (similar to what we did in [recount-website](https://github.com/leekgroup/recount-website)).

4. Hidden run files will be created with all calls needed to run the pipeline. Each step is submitted to SGE cluster and queued to run sequentially. Steps that can be parallelized are submitted as array jobs.

5. Completion emails: by default you will only get an email if a job failed. If you want to get completion emails add the empty file `${DIR}/.send_emails`. You can create it with:

  ```
  touch ${DIR}/.send_emails
  ```
