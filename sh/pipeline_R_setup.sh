#!/bin/bash
#$ -cwd
#$ -N pipeline_setup
#$ -o ./logs/pipeline_setup.txt
#$ -e ./logs/pipeline_setup.txt
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

mkdir -p logs

Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite(c('Biostrings', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db', 'biomaRt', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38', 'org.Mm.eg.db', 'BSgenome.Mmusculus.UCSC.mm10', 'org.Rn.eg.db', 'BSgenome.Rnorvegicus.UCSC.rn6', 'derfinder', 'bumphunter', 'LieberInstitute/jaffelab', 'devtools', 'getopt', 'BiocParallel', 'rafalib')); options(width = 120); devtools::session_info()"

echo "**** Job ends ****"
date
