#!/bin/bash
#$ -cwd
#$ -N pipeline_setup
#$ -o ./logs/pipeline_setup.o.txt
#$ -e ./logs/pipeline_setup.e.txt
echo "**** Job starts ****"
date

mkdir -p logs

echo -e "**** Pipeline version: GitHub sha ****\n${pipelineversion}"
Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite(c('Biostrings', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db', 'biomaRt', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38', 'org.Mm.eg.db', 'BSgenome.Mmusculus.UCSC.mm10', 'org.Rn.eg.db', 'BSgenome.Rnorvegicus.UCSC.rn6', 'derfinder', 'bumphunter', 'LieberInstitute/jaffelab', 'devtools', 'getopt', 'BiocParallel', 'rafalib', 'qualV')); options(width = 120); devtools::session_info()"

echo "**** Job ends ****"
date
