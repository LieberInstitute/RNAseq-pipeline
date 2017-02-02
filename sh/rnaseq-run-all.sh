#!/bin/bash

## Usage
# qrsh
# bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh ${EXPERIMENT} ${PREFIX} ${hgXX} ${STRANDED} ${ERCC} ${FQ_FOLDER} ${CORES} ${LARGE} ${FULLCOV} ${BASH_FOLDER} ${ANNO_FOLDER}
# bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh testrun run1 hg38 TRUE TRUE FALSE /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/test_runthroughAZ/fq
# bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh bs run1 hg38 FALSE FALSE FALSE /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/brainspan
# bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh fulltest sep23 hg38 TRUE TRUE TRUE /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/test_runthroughAZ/fq

# Define variables
EXPERIMENT=$1
PREFIX=$2
hgXX=$3
STRANDED=$4
ERCC=$5
FQ_FOLDER=${6-""}
CORES=${7-8}
LARGE=${8-"FALSE"}
FULLCOV=${9-"FALSE"}
BASH_FOLDER=${10-"/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh"}
ANNO_FOLDER=${11-"/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation"}

## Try running R. If it fails it means that the user is on the login node.
Rscript -e "Sys.time()" &> .try_load_R
LOGNODE=$(grep force-quitting .try_load_R | wc -l)
if [ ${LOGNODE} != "0" ]
then
    echo "**** You are on the login node. Use qrsh to run this script ****"
    exit 1
fi
rm .try_load_R

## Create logs dir, otherwise scripts fail since they use the -o and -e
## options
mkdir -p logs

## Check dependencies
echo "**** checking that R packages are present ****"
if [ -f ".missing_R_packages" ]
then
    echo "**** Installing R packages since some of them are missing ****"
    qsub ${BASH_FOLDER}/pipeline_R_setup.sh
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

echo "**** checking that checksumdir is installed ****"
python -c "import checksumdir" &> .check_checksumdir
CHECKSUM=$(cat .check_checksumdir | wc -l)
if [ ${CHECKSUM} != "0" ]
then
    echo "**** Installing checksumdir ****"
    pip install --user checksumdir
fi
rm .check_checksumdir

## Save the information about the pipeline version and annotation folder
## for reproducibility purposes
REPODIR=$(dirname $BASH_FOLDER)
pipelineversion=$(git --git-dir=${REPODIR}/.git rev-parse origin/master)

echo "**** Computing the md5 for ${ANNO_FOLDER}, takes 2-3 minutes ****"
annofoldermd5=$(~/.local/bin/checksumdir -a md5 ${ANNO_FOLDER})

## Save the reproducibility info
echo -e "**** Pipeline version: GitHub sha ****\n${pipelineversion}\n**** BASH_FOLDER: ****\n${BASH_FOLDER}\n**** ANNO_FOLDER: ****\n${ANNO_FOLDER}\n**** md5sum for ANNO_FOLDER ****\n${annofoldermd5}\n**** ANNO_FOLDER contents ****" > logs/pipeline_information.txt
ls -lhtR ${ANNO_FOLDER} >> logs/pipeline_information.txt

# Set variables for desired genome version
if [ $hgXX == "hg38" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/hisat2_GRCh38primary
	CHRSIZES=${ANNO_FOLDER}/hg38.chrom.sizes.gencode
    BED=${ANNO_FOLDER}/RSeQC/hg38.bed
elif [ $hgXX == "hg19" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/hisat2_GRCh37primary
	CHRSIZES=${ANNO_FOLDER}/hg19.chrom.sizes.gencode
    BED=${ANNO_FOLDER}/RSeQC/hg19.bed
elif [ $hgXX == "mm10" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/gencode.vM11.annotation.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/hisat2_GRCm38primary
	CHRSIZES=${ANNO_FOLDER}/mm10.chrom.sizes.gencode
    BED=${ANNO_FOLDER}/RSeQC/mm10.bed
elif [ $hgXX == "rn6" ]
then 
	GTF=${ANNO_FOLDER}/ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf
	HISATIDX=${ANNO_FOLDER}/ensembl/Rnor_6.0/hisat2_Rnor6.0toplevel
	CHRSIZES=${ANNO_FOLDER}/rn6.chrom.sizes.ensembl
    BED=${ANNO_FOLDER}/RSeQC/rn6.bed
else 
	echo "Enter hg19 or hg38, mm10 for mouse, or rn6 for rat." >&2
	exit 1
fi

## Find extension of fastq file and whether to merge or not
## also add  full paths to SAMPLE_IDs.txt if necessary
Rscript ${BASH_FOLDER}/find_sample_info.R -s SAMPLE_IDs.txt -f "${FQ_FOLDER}"

if [ ! -f ".file_extensions.txt" ]
then
    echo "Error: could not find .file_extensions.txt"
    exit 1
fi

## create and submit all scripts

if [ -f ".requires_merging" ]
then
    rm .file_extensions.txt
    sh ${BASH_FOLDER}/step00-merge.sh ${EXPERIMENT} ${PREFIX} ${CORES} ${LARGE} ${BASH_FOLDER}
    rm .requires_merging
fi

if [ ${ERCC} == "TRUE" ]
then
    sh ${BASH_FOLDER}/step0-ercc.sh ${EXPERIMENT} ${PREFIX} ${CORES} ${LARGE}
fi

sh ${BASH_FOLDER}/step1-fastqc.sh ${EXPERIMENT} ${PREFIX} ${LARGE}
sh ${BASH_FOLDER}/step2-trim.sh ${EXPERIMENT} ${PREFIX} ${CORES} ${LARGE}
sh ${BASH_FOLDER}/step3-hisat2.sh ${EXPERIMENT} ${PREFIX} ${HISATIDX} ${BED} ${CORES} ${LARGE}
sh ${BASH_FOLDER}/step4-featureCounts.sh ${EXPERIMENT} ${PREFIX} ${STRANDED} ${GTF} ${hgXX} ${CORES} ${LARGE}
sh ${BASH_FOLDER}/step5-coverage.sh ${EXPERIMENT} ${PREFIX} ${CHRSIZES} ${LARGE}
sh ${BASH_FOLDER}/step5b-meanCoverage.sh ${EXPERIMENT} ${PREFIX} ${CHRSIZES} ${LARGE}
sh ${BASH_FOLDER}/step6-makeRobjects.sh ${EXPERIMENT} ${PREFIX} ${hgXX} ${ERCC} ${CORES} ${LARGE} ${FULLCOV} ${BASH_FOLDER}
