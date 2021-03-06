#!/bin/bash

## Usage information:
# bash step7-makeRobjects.sh --help

# Define variables
TEMP=$(getopt -o x:p:r:e:l:f:b:a:s:h --long experiment:,prefix:,reference:,ercc:,large:,fullcov:,bashfolder:,annofolder:,stranded:,help -n 'step7-makeRobjects' -- "$@")
eval set -- "$TEMP"

ERCC="FALSE"
LARGE="FALSE"
FULLCOV="FALSE"
BASH_FOLDER="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh"
STRANDED="FALSE"

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
        -e|--ercc)
        case "$2" in
            "") ERCC="FALSE" ; shift 2;;
            *) ERCC=$2; shift 2;;
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
        -s|--stranded)
            case "$2" in
                "") STRANDED="FALSE" ; shift 2;;
                *) STRANDED=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step7-makeRobjects.sh -x -p -r (hg38, hg19, mm10, rn6) -e (default:FALSE) -l (default:FALSE) -f (default:FALSE) -b (default:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh) -s (default:FALSE)\nLong options:\n  bash step7-makeRobjects.sh --experiment --prefix --reference (hg38, hg19, mm10, rn6) --ercc (default:FALSE) --large (default:FALSE) --fullcov (default:FALSE) --bashfolder (default:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh) --stranded (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SHORT="Rcounts-${EXPERIMENT}"
sname="step7-${SHORT}.${PREFIX}"
SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
CORES=5

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=28G,h_vmem=30G,h_fsize=200G"
else
    MEM="mem_free=18G,h_vmem=20G,h_fsize=200G"
fi

if [ -f ".send_emails" ]
then
    EMAIL="e"
else
    EMAIL="a"
fi

if [ -f ".queue" ]
then
    SGEQUEUE="$(cat .queue),"
else
    SGEQUEUE=""
fi

if [ -f ".paired_end" ]
then
    PE="TRUE"
else
    PE="FALSE"
fi

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

if [ $hgXX == "mm10" ]; then SPEC="mouse";
elif [ $hgXX == "rn6" ]; then SPEC="rat";
else SPEC="human";
fi

cp ${BASH_FOLDER}/create_count_objects-${SPEC}.R ${MAINDIR}/.step7-create_count_objects-${SPEC}.R
cp ${BASH_FOLDER}/create_fullCov_object.R ${MAINDIR}/.step7-create_fullCov_object.R

# Construct shell files
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe local ${CORES}
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -hold_jid pipeline_setup,step4-featCounts-${EXPERIMENT}.${PREFIX},step6-txQuant-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "****"

Rscript ${MAINDIR}/.step7-create_count_objects-${SPEC}.R -o ${hgXX} -m ${MAINDIR} -e ${EXPERIMENT} -p ${PREFIX} -l ${PE} -c ${ERCC} -t ${CORES} -s ${STRANDED}

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call

if [[ ${FULLCOV} == "TRUE" ]]
then
    SHORT="fullCov-${EXPERIMENT}"
    sname="step7-${SHORT}.${PREFIX}"
    # Construct shell files
    echo "Creating script ${sname}"
    
    cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe local ${CORES}
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -hold_jid pipeline_setup,step5-coverage-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "****"


## Don't create the fullCoverage object by default:
## it's not needed since we create the mean bigwig already and can use it
## to define the ERs with derfinder::findRegions(), then use the resulting
## GRanges object and the paths to the BigWig files in
## derfinder::getRegionCoverage(fullCov = NULL, files = bigWigs, regions = outputFrom_findRegions)
## or alternatively write the regions to a BED file with rtracklayer,
## create the counts with bwtool and then read them into R manually
Rscript ${MAINDIR}/.step7-create_fullCov_object.R -o ${hgXX} -m ${MAINDIR} -e ${EXPERIMENT} -p ${PREFIX} -l ${PE} -f ${FULLCOV} -c ${CORES}


echo "**** Job ends ****"
date
EOF

    call="qsub .${sname}.sh"
    echo $call
    $call
fi
