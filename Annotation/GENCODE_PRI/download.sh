#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=18G,h_vmem=20G,h_fsize=100G
#$ -N download_gtf
#$ -o download.log
#$ -e download.log
#$ -m a

## Download Gencodev25 hg38 PRI gtf file and gunzip it while keeping the original .gz file
echo "**** Job starts ****"
date
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.primary_assembly.annotation.gtf.gz

gunzip < gencode.v25.primary_assembly.annotation.gtf.gz > gencode.v25.primary_assembly.annotation.gtf

echo "**** Job ends ****"
date
