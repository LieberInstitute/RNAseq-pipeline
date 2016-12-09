#!/bin/sh

## Usage
# qrsh
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
FQ_FOLDER=${7-""}
MERGE=${8-"FALSE"}
LARGE=${9-"FALSE"}
FULLCOV=${10-"FALSE"}


echo "**** Pipeline version: latest GitHub sha ****"
git --git-dir=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/.git rev-parse origin/master

## Try running R. If it fails it means that the user is on the login node.
Rscript -e "Sys.time()" &> .try_load_R
LOGNODE=$(grep force-quitting .try_load_R | wc -l)
if [ ${LOGNODE} != "0" ]
then
    echo "**** You are on the login node. Use qrsh to run this script ****"
    exit 1
fi
rm .try_load_R

SH_FOLDER=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh
ANNO_FOLDER=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation

## Check dependencies
echo "**** checking that R packages are present ****"
if [ -e ".missing_R_packages" ]
then
    echo "**** Installing R packages since some of them are missing ****"
    qsub ${SH_FOLDER}/pipeline_R_setup.sh
    rm .missing_R_packages
fi

echo "**** checking that RSeQC is installed ****"
module load python/2.7.9
python -c 'from qcmodule import SAM' &> .check_python_rseqc
RSEQC=$(cat .check_python_rseqc | wc -l)
if [ ${RSEQC} != "0" ]
then
    echo "**** Installing RSeQC: will take less than 5 minutes ****"
    echo "**** You can test that it successfully installed by running: python -c 'from qcmodule import SAM'    ****"
    pip install --user RSeQC
fi
rm .check_python_rseqc

# Set variables for desired genome version
if [ $hgXX == "hg38" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/hisat2_GRCh38primary
	CHRSIZES=${ANNO_FOLDER}/hg38.chrom.sizes.gencode
elif [ $hgXX == "hg19" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/hisat2_GRCh37primary
	CHRSIZES=${ANNO_FOLDER}/hg19.chrom.sizes.gencode
elif [ $hgXX == "mm10" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/gencode.vM11.annotation.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/hisat2_GRCm38primary
	CHRSIZES=${ANNO_FOLDER}/mm10.chrom.sizes.gencode
elif [ $hgXX == "rn6" ]
then 
	GTF=${ANNO_FOLDER}/ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf
	HISATIDX=${ANNO_FOLDER}/ensembl/Rnor_6.0/hisat2_Rnor6.0toplevel
	CHRSIZES=${ANNO_FOLDER}/rn6.chrom.sizes.ensembl
else 
	echo "Enter hg19 or hg38, mm10 for mouse, or rn6 for rat." >&2
	exit 1
fi

## Add full paths to SAMPLE_IDs.txt if necessary
if [ "${FQ_FOLDER}" != "" ]
then
    echo "Adding ${FQ_FOLDER} to SAMPLE_IDs.txt"
    mv SAMPLE_IDs.txt .SAMPLE_IDs_original.txt
    awk -v fold=${FQ_FOLDER} '{print fold"/" $0;}' .SAMPLE_IDs_original.txt > SAMPLE_IDs.txt
fi

## Find extension of fastq file names
FILEID=$(head -n 1 SAMPLE_IDs.txt | cut -f 1 -d " ")
Rscript ${SH_FOLDER}/find_extension.R -f ${FILEID}

if [ -e ".FILE_extension.txt" ]
then
    EXT=$(cat .FILE_extension.txt)
else
    exit 1
fi

## create and submit all scripts

if [ ${MERGE} == "TRUE" ]
then
    sh ${SH_FOLDER}/step00-merge.sh ${EXPERIMENT} ${PREFIX} ${PE} ${SH_FOLDER} ${LARGE}
fi

if [ ${ERCC} == "TRUE" ]
then
    sh ${SH_FOLDER}/step0-ercc.sh ${EXPERIMENT} ${PREFIX} ${PE} ${LARGE}
fi

sh ${SH_FOLDER}/step1-fastqc.sh ${EXPERIMENT} ${PREFIX} ${PE} ${LARGE}
sh ${SH_FOLDER}/step2-trim.sh ${EXPERIMENT} ${PREFIX} ${PE} ${LARGE}
sh ${SH_FOLDER}/step3-hisat2.sh ${EXPERIMENT} ${PREFIX} ${PE} ${HISATIDX} ${LARGE}
sh ${SH_FOLDER}/step4-featureCounts.sh ${EXPERIMENT} ${PREFIX} ${STRANDED} ${GTF} ${hgXX} ${PE} ${LARGE}
sh ${SH_FOLDER}/step5-coverage.sh ${EXPERIMENT} ${PREFIX} ${CHRSIZES} ${LARGE}
sh ${SH_FOLDER}/step5b-meanCoverage.sh ${EXPERIMENT} ${PREFIX} ${LARGE}
sh ${SH_FOLDER}/step6-makeRobjects.sh ${EXPERIMENT} ${PREFIX} ${hgXX} ${SH_FOLDER} ${PE} ${ERCC} ${LARGE} ${FULLCOV}
