#!/bin/bash

echo "**** Cleaning job starts ****"
date

## Reset the manifest file
mv merge/.samples_unmerged.manifest merge/samples.manifest

rm */.*.sh */.*.R */.paired_end */.send_emails */*.txt */*.rda */*.Rdata */*.csv
rm -fr */Coverage */logs */FastQC */HISAT2_out */Counts merge/merged_fastq */Salmon_tx

echo "**** Cleaning job ends ****"
date
