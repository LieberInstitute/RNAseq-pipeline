#!/bin/bash
#$ -cwd
#$ -N pipeline_setup
#$ -pe local 8
#$ -o ./logs/pipeline_setup.o.txt
#$ -e ./logs/pipeline_setup.e.txt
echo "**** Job starts ****"
date

mkdir -p logs

echo "**** Pipeline version: latest GitHub sha ****"
Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite(c('Biostrings', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db', 'biomaRt', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38', 'org.Mm.eg.db', 'BSgenome.Mmusculus.UCSC.mm10', 'org.Rn.eg.db', 'BSgenome.Rnorvegicus.UCSC.rn6', 'derfinder', 'bumphunter', 'LieberInstitute/jaffelab', 'devtools', 'getopt', 'BiocParallel', 'rafalib')); options(width = 120); devtools::session_info()"

echo "**** Job ends ****"
date
