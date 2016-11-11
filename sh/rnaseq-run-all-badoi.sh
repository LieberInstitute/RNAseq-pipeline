#!/bin/sh

## Usage
# sh /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh testrun run1 hg38 TRUE TRUE FALSE /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/test_runthroughAZ/fq
# sh /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh bs run1 hg38 FALSE FALSE FALSE /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/brainspan
# sh /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh fulltest sep23 hg38 TRUE TRUE TRUE /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/test_runthroughAZ/fq

# Define variables
EXPERIMENT=$1
PREFIX=$2
hgXX=$3
PE=$4
STRANDED=$5
ERCC=$6
FQ_FOLDER=$7

SH_FOLDER=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh
ANNO_FOLDER=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation

# Set variables for desired genome version
if [ $hgXX == "hg38" ] ; then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/hisat2_GRCh38primary
	CHRSIZES=${ANNO_FOLDER}/hg38.chrom.sizes.gencode
elif [ $hgXX == "hg19" ] ; then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/hisat2_GRCh37primary
	CHRSIZES=${ANNO_FOLDER}/hg19.chrom.sizes.gencode
elif [ $hgXX == "mm10" ] ; then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/gencode.vM11.annotation.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/hisat2_GRCm38primary
	CHRSIZES=${ANNO_FOLDER}/mm10.chrom.sizes.gencode
elif [ $hgXX == "rn6" ] ; then 
	GTF=${ANNO_FOLDER}/ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf
	HISATIDX=${ANNO_FOLDER}/ensembl/Rnor_6.0/hisat2_Rnor6.0toplevel
	CHRSIZES=${ANNO_FOLDER}/rn6.chrom.sizes.ensembl
else 
	echo "Enter hg19 or hg38, mm10 for mouse, or rn6 for rat." >&2
	exit 1
fi

## Find extension of fastq file names
FILENAME=$(ls $FQ_FOLDER | head -1)
if [ $(printf $FILENAME | tail -c 8) == "fastq.gz" ] ; then EXT="fastq.gz"
elif [ $(printf $FILENAME | tail -c 5) == "fq.gz" ] ; then EXT="fq.gz"
elif [ $(printf $FILENAME | tail -c 5) == "fastq" ] ; then EXT="fastq"
elif [ $(printf $FILENAME | tail -c 2) == "fq" ] ; then EXT="fq"
else 
	echo "Unrecognized fastq filename extension." >&2
	exit 1
fi


## create and submit all scripts

#if [ $ERCC == "TRUE" ] ; then
#sh ${SH_FOLDER}/step0-ercc.sh ${EXPERIMENT} ${PREFIX} ${PE} ${FQ_FOLDER} ${EXT}
#fi

#sh ${SH_FOLDER}/step1-fastqc.sh ${EXPERIMENT} ${PREFIX} ${PE} ${FQ_FOLDER} ${EXT}
#sh ${SH_FOLDER}/step2-trim.sh ${EXPERIMENT} ${PREFIX} ${PE} ${FQ_FOLDER} ${EXT}
#sh ${SH_FOLDER}/step3-hisat2.sh ${EXPERIMENT} ${PREFIX} ${PE} ${FQ_FOLDER} ${EXT} ${HISATIDX}
#sh ${SH_FOLDER}/step4-featureCounts.sh ${EXPERIMENT} ${PREFIX} ${STRANDED} ${GTF} ${hgXX} ${PE}
sh ${SH_FOLDER}/step5-coverage.sh ${EXPERIMENT} ${PREFIX} ${CHRSIZES}
#sh ${SH_FOLDER}/step6-makeRobjects.sh ${EXPERIMENT} ${PREFIX} ${hgXX} ${SH_FOLDER} ${PE} ${ERCC}


