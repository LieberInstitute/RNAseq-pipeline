#!/bin/bash

## Usage information:
# bash rnaseq-run-all.sh --help

# Define variables
TEMP=$(getopt -o x:p:r:s:e:c:l:f:b:a:u:h --long experiment:,prefix:,reference:,stranded:,ercc:,cores:,large:,fullcov:,bashfolder:,annofolder:,unaligned:,help -n 'rnaseq-run-all' -- "$@")
eval set -- "$TEMP"

STRANDED="FALSE"
ERCC="FALSE"
LARGE="FALSE"
FULLCOV="FALSE"
CORES=8
BASH_FOLDER="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh"
ANNO_FOLDER="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation"
UNALIGNED="FALSE"


while true; do
    case "$1" in
        -x|--experiment)
            case "$2" in
                "") shift 2 ;;
                *) EXPERIMENT=$2 ; shift 2;;
            esac;;
        -p|--prefix)
            case "$2" in
                "") shift 2 ;;
                *) PREFIX=$2 ; shift 2;;
            esac;;
        -r|--reference)
            case "$2" in
                "") shift 2 ;;
                *) hgXX=$2 ; shift 2;;
            esac;;
        -s|--stranded)
            case "$2" in
                "") STRANDED="FALSE" ; shift 2;;
                *) STRANDED=$2; shift 2;;
            esac ;;
        -e|--ercc)
            case "$2" in
                "") ERCC="FALSE" ; shift 2;;
                *) ERCC=$2; shift 2;;
            esac ;;
        -c|--cores)
            case "$2" in
                "") CORES="8" ; shift 2;;
                *) CORES=$2; shift 2;;
            esac ;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -f|--fullcov)
            case "$2" in
                "") FULLCOV="FALSE" ; shift 2;;
                *) FULLCOV=$2; shift 2;;
            esac ;;
        -b|--bashfolder)
            case "$2" in
                "") BASH_FOLDER="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh"; shift 2;;
                *) BASH_FOLDER=$2; shift 2;;
            esac;;
        -a|--annofolder)
            case "$2" in
                "") ANNO_FOLDER="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation"; shift 2;;
                *) ANNO_FOLDER=$2; shift 2;;
            esac;;
        -u|--unaligned)
            case "$2" in
                "") UNALIGNED="FALSE" ; shift 2;;
                *) UNALIGNED=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\n  qrsh\nShort options:\n  bash rnaseq-run-all.sh -x -p -r (hg38, hg19, mm10, rn6) -s (default:FALSE) -e (default:FALSE) -c (default:8) -l (default:FALSE) -f (default:FALSE) -b (default:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh) -a (default:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation) -u (default:FALSE)\nLong options:\n  bash rnaseq-run-all.sh --experiment --prefix --reference (hg38, hg19, mm10, rn6) --stranded (default:FALSE) --ercc (default:FALSE) --cores (default:8) --large (default:FALSE) --fullcov (default:FALSE) --bashfolder (default:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh) --annofolder (default:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation) --unaligned (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

if [ ${STRANDED} == "FALSE" ]
then
    STRANDOPTION=""
elif [ ${STRANDED} == "forward" ]
then
    STRANDOPTION="forward"
elif [ ${STRANDED} == "reverse" ]
then
    STRANDOPTION="reverse"
else
    echo "The option --stranded has to either be 'FALSE', 'forward' or 'reverse'."
    exit 1
fi

## Try running R. If it fails it means that the user is on the login node.
Rscript -e "Sys.time()" &> .try_load_R
LOGNODE=$(grep force-quitting .try_load_R | wc -l)
if [ ${LOGNODE} != "0" ]
then
    echo "**** You are on the login node. Use qrsh to run this script ****"
    date
    exit 1
fi
rm .try_load_R

## Create logs dir, otherwise scripts fail since they use the -o and -e
## options
mkdir -p logs

## Check dependencies
echo "**** checking that R packages are present ****"
date
Rscript ${BASH_FOLDER}/check_R_packages.R
if [ -f ".missing_R_packages" ]
then
    echo "**** Installing R packages since some of them are missing ****"
    qsub ${BASH_FOLDER}/pipeline_R_setup.sh
    rm .missing_R_packages
fi

echo "**** checking that RSeQC is installed ****"
date
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
date
python -c "import checksumdir" &> .check_checksumdir
CHECKSUM=$(cat .check_checksumdir | wc -l)
if [ ${CHECKSUM} != "0" ]
then
    echo "**** Installing checksumdir ****"
    pip install --user checksumdir
fi
rm .check_checksumdir

# Set variables for desired genome version
if [ $hgXX == "hg38" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/hisat2_GRCh38primary
	CHRSIZES=${ANNO_FOLDER}/hg38.chrom.sizes.gencode
    BED=${ANNO_FOLDER}/RSeQC/hg38.bed
	SALMONIDX=${ANNO_FOLDER}/GENCODE/GRCh38_hg38/transcripts/salmon_index_gencode.v25.transcripts
	STEP6=TRUE
	STEP8=TRUE
elif [ $hgXX == "hg19" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/hisat2_GRCh37primary
	CHRSIZES=${ANNO_FOLDER}/hg19.chrom.sizes.gencode
    BED=${ANNO_FOLDER}/RSeQC/hg19.bed
	SALMONIDX=${ANNO_FOLDER}/GENCODE/GRCh37_hg19/transcripts/salmon_index_gencode.v25lift37.transcripts
	STEP6=TRUE
	STEP8=FALSE
elif [ $hgXX == "mm10" ]
then 
	GTF=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/gencode.vM11.annotation.gtf
	HISATIDX=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/hisat2_GRCm38primary
	CHRSIZES=${ANNO_FOLDER}/mm10.chrom.sizes.gencode
    BED=${ANNO_FOLDER}/RSeQC/mm10.bed
	SALMONIDX=${ANNO_FOLDER}/GENCODE/GRCm38_mm10/transcripts/salmon_index_gencode.vM11.transcripts
	STEP6=TRUE
	STEP8=FALSE
elif [ $hgXX == "rn6" ]
then 
	GTF=${ANNO_FOLDER}/ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf
	HISATIDX=${ANNO_FOLDER}/ensembl/Rnor_6.0/hisat2_Rnor6.0toplevel
	CHRSIZES=${ANNO_FOLDER}/rn6.chrom.sizes.ensembl
    BED=${ANNO_FOLDER}/RSeQC/rn6.bed
	STEP6=FALSE
	STEP8=FALSE
else
	echo "Error: enter hg19 or hg38, mm10 for mouse, or rn6 for rat as the reference." >&2
    date
	exit 1
fi

## Save the information about the pipeline version and annotation folder
## for reproducibility purposes
REPODIR=$(dirname $BASH_FOLDER)
pipelineversion=$(git --git-dir=${REPODIR}/.git rev-parse origin/master)

#echo "**** Computing the md5 for ${ANNO_FOLDER}, takes 2-3 minutes ****"
#date
#annofoldermd5=$(~/.local/bin/checksumdir -a md5 ${ANNO_FOLDER})
annofoldermd5="--skipped checksumdir--"

## Save the reproducibility info
echo "**** Saving the reproducibility information in logs/pipeline_information.txt ****"
date
echo -e "**** Pipeline version: GitHub sha ****\n${pipelineversion}\n**** BASH_FOLDER: ****\n${BASH_FOLDER}\n**** ANNO_FOLDER: ****\n${ANNO_FOLDER}\n**** md5sum for ANNO_FOLDER ****\n${annofoldermd5}\n**** ANNO_FOLDER contents ****" > logs/pipeline_information.txt
ls -lhtR ${ANNO_FOLDER} >> logs/pipeline_information.txt

## Find extension of fastq file and whether to merge or not
## also add  full paths to samples.manifest if necessary
echo "**** Finding the sample information ****"
date
Rscript ${BASH_FOLDER}/find_sample_info.R -s samples.manifest -o ${PWD}/merged_fastq


echo "**** Creating bash scripts for every step ****"

## Check that find_sample_info.R worked.
## It can fail in some situations that are checked by the R code.
if [ -f "find_sample_error" ]
then
    echo "Fatal error when running: ${BASH_FOLDER}/find_sample_info.R -s samples.manifest -o ${PWD}/merged_fastq"
    echo "Once you have fixed the error, manually remove the file find_sample_error before running this script again."
    exit 1
fi

## create and submit all scripts
if [ -f ".requires_merging" ]
then
    sh ${BASH_FOLDER}/step00-merge.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --cores ${CORES} --large ${LARGE} --bashfolder ${BASH_FOLDER}
    rm .requires_merging
fi

if [ ${ERCC} == "TRUE" ]
then
    sh ${BASH_FOLDER}/step0-ercc.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --cores ${CORES} --large ${LARGE} --stranded ${STRANDED}
fi

sh ${BASH_FOLDER}/step1-fastqc.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --large ${LARGE}
sh ${BASH_FOLDER}/step2-trim.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --cores ${CORES} --large ${LARGE}
sh ${BASH_FOLDER}/step3-hisat2.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --index ${HISATIDX} --bed ${BED} --cores ${CORES} --large ${LARGE} --stranded ${STRANDED} --unaligned ${UNALIGNED}
sh ${BASH_FOLDER}/step4-featureCounts.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --stranded ${STRANDED} --gtf ${GTF} --reference ${hgXX} --cores ${CORES} --large ${LARGE}
sh ${BASH_FOLDER}/step5-coverage.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --chrsizes ${CHRSIZES} --large ${LARGE}
sh ${BASH_FOLDER}/step5b-meanCoverage.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --chrsizes ${CHRSIZES} --large ${LARGE}

if [ ${STEP6} == "TRUE" ]
then
    sh ${BASH_FOLDER}/step6-txQuant.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --stranded ${STRANDED} --index ${SALMONIDX} --cores ${CORES} --large ${LARGE}
fi

sh ${BASH_FOLDER}/step7-makeRobjects.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --reference ${hgXX} --ercc ${ERCC} --cores ${CORES} --large ${LARGE} --fullcov ${FULLCOV} --bashfolder ${BASH_FOLDER} --stranded ${STRANDED}

if [ ${STEP8} == "TRUE" ]
then
    sh ${BASH_FOLDER}/step8-callVariants.sh --experiment ${EXPERIMENT} --prefix ${PREFIX} --reference ${hgXX} --cores ${CORES} --large ${LARGE}
	sh ${BASH_FOLDER}/step8b-mergeVariantCalls.sh --experiment ${EXPERIMENT} --prefix ${PREFIX}
fi



