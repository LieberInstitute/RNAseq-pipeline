#!/bin/bash
#$ -cwd
#$ -N pipeline_ex_simulated
#$ -o ./pipeline_ex_simulated.txt
#$ -e ./pipeline_ex_simulated.txt
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
echo "bluejay" > paired_end_stranded/.queue
echo "bluejay" > paired_end_unstranded/.queue
echo "bluejay" > single_end_stranded/.queue
echo "bluejay" > single_end_unstranded/.queue

echo "**************************"
echo "Running single end unstranded example"
echo "**************************"
cd single_end_unstranded
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "single_unstranded" --reference "hg38" --cores 1

echo "**************************"
echo "Running paired end unstranded example"
echo "**************************"
cd ../paired_end_unstranded
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "paired_unstranded" --reference "hg38" --cores 1

echo "**************************"
echo "Running single end stranded example"
echo "**************************"
cd ../single_end_stranded
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "single_stranded" --reference "hg38" --cores 1 --stranded "reverse"

echo "**************************"
echo "Running paired end stranded example"
echo "**************************"
cd ../paired_end_stranded
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "paired_stranded" --reference "hg38" --cores 1  --stranded "reverse"

echo "**** Job ends ****"
date
