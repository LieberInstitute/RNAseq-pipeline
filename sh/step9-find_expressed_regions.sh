#!/bin/bash

## Usage information:
# bash step9-find_expressed_regions.sh --help

# Define variables
TEMP=$(getopt -o x:p:s:l:h --long experiment:,prefix:,chrsizes:,large:,help -n 'step9-find_expressed_regions' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"

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
        -s|--chrsizes)
            case "$2" in
                "") shift 2 ;;
                *) CHRSIZES=$2 ; shift 2;;
            esac;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step9-find_expressed_regions.sh -x -p -s -l (default:FALSE)\nLong options:\n  bash step9-find_expressed_regions.sh --experiment --prefix --chrsizes --large (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="findERs-${EXPERIMENT}"
sname="step9-${SHORT}.${PREFIX}"
CORES=10

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=25G,h_vmem=25G,h_fsize=200G"
else
    MEM="mem_free=15G,h_vmem=15G,h_fsize=200G"
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

# Directories
mkdir -p ${MAINDIR}/Coverage

cp ${BASH_FOLDER}/step9-find_expressed_regions.R ${MAINDIR}/.step9-find_expressed_regions.R

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe local ${CORES}
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -hold_jid pipeline_setup,step5b-meanCoverage-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "****"

## Locate mean files
meanFiles=Coverage/mean*.bw

for meanFile in meanFiles
do
    echo "Initializing script for \${meanFile}"
    Rscript ${MAINDIR}/.step9-find_expressed_regions.R -m \${meanFile} -o ${MAINDIR}/Coverage -i ${CHRSIZES} -c 10
done

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
