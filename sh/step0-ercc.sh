#!/bin/bash

## Usage information:
# bash step0-ercc.sh --help

# Define variables
TEMP=$(getopt -o x:p:c:l:s:h --long experiment:,prefix:,cores:,large:,stranded:,help -n 'step0-ercc' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"
CORES=8
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
        -s|--stranded)
            case "$2" in
                "") STRANDED="FALSE" ; shift 2;;
                *) STRANDED=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step0-ercc.sh -x -p -c (default:8) -l (default:FALSE) -s (default:FALSE)\nLong options:\n  bash step0-ercc.sh --experiment --prefix --cores (default:8) --large (default:FALSE) --stranded (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="ercc-${EXPERIMENT}"
sname="step0-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=6G,h_vmem=10G,h_fsize=100G"
else
    MEM="mem_free=3G,h_vmem=5G,h_fsize=100G"
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
    STRANDOPTION="--fr-stranded"
elif [ ${STRANDED} == "reverse" ]
then
    STRANDOPTION="--rf-stranded"
else
    echo "The option --stranded has to either be 'FALSE', 'forward' or 'reverse'."
    exit 1
fi


# Construct shell files
FILELIST=${MAINDIR}/samples.manifest
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 20
#$ -hold_jid pipeline_setup,step00-merge-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"
echo "****"
echo "Sample id: \$(cat ${MAINDIR}/samples.manifest | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
if [ $PE == "TRUE" ] 
then
    FILE2=\$(awk 'BEGIN {FS="\t"} {print \$3}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
fi
ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")
mkdir -p ${MAINDIR}/Ercc/\${ID}

if [ $PE == "TRUE" ]
then 
    ${SOFTWARE}/kallisto quant \
    -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx \
    -o ${MAINDIR}/Ercc/\${ID} -t ${CORES} ${STRANDOPTION} \
    \${FILE1} \${FILE2}
else
    ${SOFTWARE}/kallisto quant \
    -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx \
    -o ${MAINDIR}/Ercc/\${ID} -t ${CORES} --single ${STRANDOPTION} \${FILE1}
fi

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
