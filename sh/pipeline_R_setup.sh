#!/bin/bash
#$ -cwd
#$ -N pipeline_setup
#$ -o ./logs/pipeline_setup.txt
#$ -e ./logs/pipeline_setup.txt
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

## Check that the user has the GITHUB_PAT variable which is needed for
## installing the private package LieberInstitute/jaffelab
if [[ "${GITHUB_PAT}" == "" ]]; then
    echo "Please follow the instructions in the R help page of devtools::install_github() to add the GITHUB_PAT environment variable. This will then allow you to install the private package https://github.com/LieberInstitute/jaffelab."
    exit 1
fi

Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite(c('Biostrings', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db', 'biomaRt', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38', 'org.Mm.eg.db', 'BSgenome.Mmusculus.UCSC.mm10', 'org.Rn.eg.db', 'BSgenome.Rnorvegicus.UCSC.rn6', 'derfinder', 'bumphunter', 'LieberInstitute/jaffelab', 'devtools', 'getopt', 'BiocParallel', 'rafalib', 'SummarizedExperiment')); options(width = 120); devtools::session_info()"

echo "**** Job ends ****"
date
