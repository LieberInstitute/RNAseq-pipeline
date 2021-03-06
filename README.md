**Outdated**

Note that this pipeline is no longer under active development. You can find the latest version at https://github.com/LieberInstitute/SPEAQeasy which is described in the following pre-print: https://www.biorxiv.org/content/10.1101/2020.12.11.386789v1.



RNAseq-pipeline
===============

This pipeline requires several R packages. You can install them by running:

```
qsub pipeline_R_setup.sh
```

1. Make a directory to deposit subfolders and processed files: `$DIR`.

  Often, `$DIR/FASTQ` contains the fastq files for this experiment but the fastq files don't necessarily have to be in this location.

2. In `$DIR`, make text file called `samples.manifest` which has to be a manifest file just like the ones used in [Rail-RNA](http://rail.bio) and [Myrna](http://bowtie-bio.sourceforge.net/myrna/). The exception is that the paths have to be local as no URLs are supported. It has the following structure:


    1. (for a set of unpaired input reads) `<FASTQ FILE>(tab)<optional MD5>(tab)<sample label>`
    2. (for a set of paired input reads) `<FASTQ FILE 1>(tab)<optional MD5 1>(tab)<FASTQ FILE 2>(tab)<optional MD5 2>(tab)<sample label>`
  
  The following extensions are supported: `fastq.gz`, `fq.gz`, `fastq`, `fq`.
  
  For example, a `samples.manifest` file can look like this (these are paired-end samples):

  ```
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L005_R1_001.fastq.gz 0   /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L005_R2_001.fastq.gz   0   sample1
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L006_R1_001.fastq.gz 0   /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L006_R2_001.fastq.gz   0   sample2
  ```
  
  If you reads are split in multiple files and you want to merge them, simply repeat a sample id. The merged files will be saved at in a directory called `merged_fastq`.
  
  ```
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L005_R2_001.fastq.gz 0   /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L005_R2_001.fastq.gz   0   sample1
  /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L006_R2_001.fastq.gz 0   /dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L006_R2_001.fastq.gz   0   sample1
  ```

3. Run pipeline with following call in `$DIR`:

  ```
  ## You need a compute node to run this script!
  qrsh
  ## Change the path if you cloned this repo elsewhere!
  bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh --help
  ```
  
  Note that you have to use __bash__ and not __sh__, otherwise you'll get an error like.

  1. __experiment__: main identifier, experiment name.
  1. __prefix__: additional identifier (could be date or reference genome used).
  1. __reference__: supported genomes are `hg19, hg38, mm10, rn6`
  1. __stranded__: `FALSE` by default. The other valid options are `forward` and `reverse`.
  1. __ercc__: `FALSE` by default and will be turned to `TRUE` if specified. Specify if ERCC mix 1 was added.
  1. __cores__: defaults to 8. Specifies how many cores to use per job for the jobs that are parallelized.
  1. __large__: `TRUE` if you want to use double the default memory settings. Useful for large projects (many samples and/or many reads). `FALSE` by default and doesn't need to be specified.
  1. __fullcov__: `TRUE` if you want to create the fullCoverage object. Set to `FALSE` by default. Note that the fullCoverage object is not needed (by default) since we create the normalized mean BigWig already and can use it to define the ERs with `derfinder::findRegions()`, then use the resulting GRanges object and the paths to the BigWig files in `derfinder::getRegionCoverage(fullCov = NULL, files = bigWigs, regions = outputFrom_findRegions)` or alternatively write the regions to a BED file with rtracklayer, create the counts with [bwtool](https://github.com/CRG-Barcelona/bwtool) and then read them into R manually (similar to what we did in [recount-website](https://github.com/leekgroup/recount-website)).
  1. __bashfolder__: defaults to `/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh`. It's the directory where the shell files are located at. You only need to specify it if you cloned this repository somewhere else.
  1. __annofolder__: defaults to `/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation`. We currently don't have scripts for you to automatically reproduce that directory, but you can copy it to another location if you want to have more control on it and change some files.
  

4. Hidden run files will be created with all calls needed to run the pipeline. Each step is submitted to SGE cluster and queued to run sequentially. Steps that can be parallelized are submitted as array jobs. If you with to manually run a given step, check the default options by using the `--help` option. For example:

  ```
  ## Change the path if you cloned this repo elsewhere!
  bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step3-hisat2.sh --help
  ```

5. Completion emails: by default you will only get an email if a job failed. If you want to get completion emails add the empty file `${DIR}/.send_emails`. You can create it with:

  ```
  touch ${DIR}/.send_emails
  ```

6. Cluster queue: by default the shared queue is used, but if you want to specify one create the hidden file `${DIR}/.queue` with the name of the queue inside of it (no spaces, no new lines). For example:

  ```
  bluejay
  ```
  
  Do not save `shared` into `${DIR}/.queue`. Otherwise your jobs won't run.

For reproducibility purposes, the information about the version of the pipeline and contents of `ANNO_FOLDER` will be stored in `logs/pipeline_information.txt`.


## First time user

Note that you will need the [jaffelab](https://github.com/LieberInstitute/jaffelab) package installed among other packages. Most of them can be installed automatically by scripts in the pipeline, but `jaffelab` requires you to create the `GITHUB_PAT` environment variable following the instructions by `devtools::install_github()`.


## Read and alignment metrics output

  __FastQC Flags__ (pass/warn/fail)

Flag | Description
--- | ---
FQCbasicStats | Composite statistic: overall pass/warn/fail
perBaseQual | Quality values across all bases at each position in the FastQ file
perTileQual | Quality scores from each tile across all bases to see if there was loss in quality associated with only one part of the flowcell
perSeqQual | Report if a subset of sequences have universally low quality values
perBaseContent | G, A, T, and C should have roughly equal proportions at each base position
GCcontent | GC content across the length of each sequence
Ncontent | Percentage of base calls at each position for which an N was called
SeqLengthDist | Distribution of fragment sizes (should all be uniform length)
SeqDuplication | Degree of duplication for every sequence in the library
OverrepSeqs | Identical sequences that appear at a higher rate than expected
AdapterContent | Proportion of your library which has adapter sequences at each position
KmerContent | Measures the number of each kmer at each position in your library and then looks for significant deviations from an even coverage at all positions
  
  __FastQC metrics__ (There will be two sets of metrics for each sample if library is paired-end.)

Metric | Description
--- | ---
SeqLength | Length of the shortest and longest sequence in the set
percentGC | Overall %GC of all bases in all sequences
phredXX-XX | Median Phred score at read positions XX (taken at different positions depending on sequence length)
phredGT30 | Proportion of reads above Phred score of 30
phredGT35 | Proportion of reads above Phred score of 35
AdapterXX-XX | Adapter content for Illumina standard adapter at positions XX (taken at different positions depending on sequence length)
 
  __HISAT2 metrics__

Metric | Description
--- | ---
numReads | Total number of reads in fastq file
numMapped | Number of reads that were aligned by HISAT2 to the reference genome
numUnmapped | Number of reads that did not align to the reference genome
overallMapRate | Proportion of total reads mapped (numMapped / numReads)
concordMapRate | If library is paired-end, percentage of read pairs that mapped concordantly

