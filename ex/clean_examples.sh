#!/bin/bash

echo "**** Cleaning job starts ****"
date

## Reset the manifest file
mv merge/.samples_unmerged.manifest merge/samples.manifest

rm */.*.sh */.*.R */.paired_end */.send_emails */*.txt */*.rda */*.Rdata */*.csv */.queue
rm -fr */Coverage */logs */FastQC */HISAT2_out */Counts merge/merged_fastq */Salmon_tx */ERs

echo "**** Cleaning job ends ****"
date
