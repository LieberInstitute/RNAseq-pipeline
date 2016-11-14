RNAseq-pipeline
===============

1. Make a directory to deposit subfolders and processed files. $DIR
2. Make text file with FASTQ file name headers delimited \n and called `SAMPLE_IDs.txt'
  Sample name prefix possibilites : SAMPLE_R1_001.fastq.gz/fq.gz & SAMPLE_R2_001.fastq.gz/fq.gz
  Or : SAMPLE.fastq.gz/fq.gz
3. Run pipeline with /DIR/TO/run_all.sh ExprName SomeIdentifier Genome PE_Bool Stranded_Bool ERCC_Bool $DIR/FASTQ_FolderName 
