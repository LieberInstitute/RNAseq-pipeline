#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=55G,h_fsize=100G
#$ -N salmon_build
#$ -pe local 1
#$ -o build_indexes_log
#$ -e build_indexes_log
#$ -m a
echo "**** Job starts ****"
date

MAINDIR=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts

fastaHG38=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/gencode.v25.transcripts.fa

### hg38

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.8.2_linux_x86_64/bin/salmon index \
-t ${fastaHG38} -i ${MAINDIR}/salmon_0.8.2_index_gencode.v25.transcripts -p 1 --type quasi -k 31

echo "**** Job ends ****"
date
