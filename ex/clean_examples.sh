#!/bin/bash

echo "**** Cleaning job starts ****"
date

## Reset the manifest file
mv merge/.samples_unmerged.manifest merge/samples.manifest

rm */.*.sh */.*.R */.paired_end */.send_emails */inferred_strandness_pattern.txt
rm -fr */Coverage */logs */FastQC */HISAT2_out */Counts merge/example

echo "**** Cleaning job ends ****"
date
