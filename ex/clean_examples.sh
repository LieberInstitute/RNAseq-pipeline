#!/bin/bash

echo "**** Cleaning job starts ****"
date

rm */.*.sh */.*.R */.paired_end */.send_emails */inferred_strandness_pattern.txt
rm -fr */Coverage */logs */FastQC */HISAT2_out */Counts 

echo "**** Cleaning job ends ****"
date
