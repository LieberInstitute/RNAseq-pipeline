#!/bin/bash
#$ -cwd
#$ -N pipeline_ex
#$ -o ./pipeline_ex.txt
#$ -e ./pipeline_ex.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

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

## Add queue files
#echo "bluejay" > single_end/.queue
#echo "bluejay" > paired_end/.queue
#echo "bluejay" > merge/.queue

echo "**************************"
echo "Running single end example"
echo "**************************"
cd single_end
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "single" --reference "hg19" --cores 1 --fullcov "TRUE"

echo "**************************"
echo "Running paired end example"
echo "**************************"
cd ../paired_end
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "paired" --reference "mm10" --cores 1

echo "**************************"
echo "Running merging example"
echo "**************************"
cd ../merge
touch .send_emails
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "merge" --reference "hg38" --stranded "reverse" --cores 2

echo "**** Job ends ****"
date
