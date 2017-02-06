#!/bin/bash

## Usage information:
# bash step00-merge.sh --help

# Define variables
TEMP=$(getopt -o x:p:c:lb:h --long experiment:,prefix:,cores:,large,bashfolder:,help -n 'step00-merge' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"
CORES=8
BASH_FOLDER="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh"

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
        -c|--cores)
            case "$2" in
                "") CORES="8" ; shift 2;;
                *) CORES=$2; shift 2;;
            esac ;;
        -l|--large) LARGE="TRUE"; shift ;;
        -b|--bashfolder)
            case "$2" in
                "") BASH_FOLDER="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh"; shift 2;;
                *) BASH_FOLDER=$2; shift 2;;
            esac;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step00-merge.sh -x -p -c (default:8) -l (default:FALSE) -b (default:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh)\nLong options:\n  bash step00-merge.sh --experiment --prefix --cores (default:8) --large (default:FALSE) --bashfolder (default:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="merge-${EXPERIMENT}"
sname="step00-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=6G,h_vmem=10G,h_fsize=150G"
else
    MEM="mem_free=3G,h_vmem=5G,h_fsize=150G"
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

# Construct shell files
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -hold_jid pipeline_setup
#$ -m ${EMAIL}

echo "**** Job starts ****"
date

Rscript ${BASH_FOLDER}/step00-merge.R -s ${MAINDIR}/samples.manifest -o ${MAINDIR}/${EXPERIMENT}/${PREFIX}/merged_fastq -c ${CORES}

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
