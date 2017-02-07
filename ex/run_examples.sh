#!/bin/bash
#$ -cwd
#$ -N pipeline_ex
#$ -o ./pipeline_ex.txt
#$ -e ./pipeline_ex.txt
#$ -m e
echo "**** Job starts ****"
date

cat > .gitignore <<EOF
*/.*.R
*/.*.sh
*/.paired_end
*/Coverage
*/logs
*/FastQC
*/HISAT2_out
*/Counts
merge/merged_fastq
merge/.send_emails
*/*.txt
*/*.rda
*/*.csv
*/Salmon_tx
EOF

## Clean examples
bash clean_examples.sh

cd single_end
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "single" --reference "hg19" --cores 1 --fullcov "TRUE"

cd ../paired_end
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "paired" --reference "mm10" --cores 1

cd ../merge
touch .send_emails
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "merge" --reference "hg38" --stranded "TRUE" --cores 2

echo "**** Job ends ****"
date
