#!/bin/bash

## Usage information:
# bash step6-txQuant.sh --help

# Define variables
TEMP=$(getopt -o x:p:s:i:c:l:h --long experiment:,prefix:,stranded:,index:,cores:,large:,help -n 'step6-txQuant' -- "$@")
eval set -- "$TEMP"

STRANDED="FALSE"
LARGE="FALSE"
CORES=1


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
        -s|--stranded)
            case "$2" in
                "") STRANDED="FALSE" ; shift 2;;
                *) STRANDED=$2; shift 2;;
            esac ;;
        -i|--index)
            case "$2" in
                "") shift 2 ;;
                *) SALMONINDEX=$2 ; shift 2;;
            esac;;
        -c|--cores)
            case "$2" in
                "") CORES="1" ; shift 2;;
                *) echo "Ignoring --cores $2 and will use 1 core"; CORES="1"; shift 2;;
            esac ;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step6-txQuant.sh -x -p -s (default:FALSE) -i -c (default:4) -l (default:FALSE)\nLong options:\n  bash step6-txQuant.sh --experiment --prefix --stranded (default:FALSE) --index --cores (default:8) --large (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="txQuant-${EXPERIMENT}"
sname="step6-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=100G,h_vmem=120G,h_fsize=100G"
else
    MEM="mem_free=80G,h_vmem=90G,h_fsize=100G"
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


# Construct shell files
FILELIST=${MAINDIR}/samples.manifest
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
#$ -pe local 1
#$ -o ./logs/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 15
#$ -hold_jid pipeline_setup,step4-featCounts-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

FILE1=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
if [ $PE == "TRUE" ] 
then
    FILE2=\$(awk 'BEGIN {FS="\t"} {print \$3}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
fi
ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")


mkdir -p ${MAINDIR}/Salmon_tx/\${ID}

if [ $PE == "TRUE" ] ; then 
	${SOFTWARE}/Salmon-0.7.2_linux_x86_64/bin/salmon quant \
	-i $SALMONINDEX -p ${CORES} -l ISR \
	-1 \${FILE1} -2 \${FILE2} \
	-o ${MAINDIR}/Salmon_tx/\${ID}
else
	${SOFTWARE}/Salmon-0.7.2_linux_x86_64/bin/salmon quant \
	-i $SALMONINDEX -p ${CORES} -l U \
	-r \${FILE1} \
	-o ${MAINDIR}/Salmon_tx/\${ID}
fi


echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call

