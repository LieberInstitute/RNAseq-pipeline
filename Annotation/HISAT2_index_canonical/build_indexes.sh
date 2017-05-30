#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=18G,h_vmem=20G,h_fsize=100G
#$ -N hisat2_build
#$ -pe local 6
#$ -o build_indexes.log
#$ -e build_indexes.log
#$ -m a
echo "**** Job starts ****"
date

module load python/2.7.9


DATADIR=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/HISAT2_index_test
MAINDIR=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/HISAT2_index_canonical


### hg19

echo "Subsetting hg19"
## Based on https://www.biostars.org/p/49773/
~/.local/bin/faidx --regex "^chr" ${DATADIR}/GRCh37.primary_assembly.genome.fa > ${MAINDIR}/GRCh37.primary_assembly.genome.canonical.fa

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2-build -p 6 \
${MAINDIR}/GRCh37.primary_assembly.genome.canonical.fa ${MAINDIR}/hisat2_GRCh37primary



### hg38
echo "Subsetting hg38"
~/.local/bin/faidx --regex "^chr" ${DATADIR}/GRCh38.primary_assembly.genome.fa > ${MAINDIR}/GRCh38.primary_assembly.genome.canonical.fa

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2-build -p 6 \
${MAINDIR}/GRCh38.primary_assembly.genome.canonical.fa ${MAINDIR}/hisat2_GRCh38primary


echo "**** Job ends ****"
date
