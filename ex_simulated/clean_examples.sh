#!/bin/bash

echo "**** Cleaning job starts ****"
date

rm */.*.sh */.*.R */.paired_end */.send_emails */*.txt */*.rda */*.Rdata */*.csv */.queue
rm -fr */Coverage */logs */FastQC */HISAT2_out */Counts merge/merged_fastq */Salmon_tx */ERs */Genotypes

echo "**** Cleaning job ends ****"
date
