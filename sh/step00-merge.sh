#!/bin/sh

## Usage
# ${SH_FOLDER}/step00-merge.sh ${EXPERIMENT} ${PREFIX} ${PE} ${SH_FOLDER} ${LARGE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
PE=$3
SH_FOLDER=$4
LARGE=${5-"FALSE"}

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="merge-${EXPERIMENT}"
sname="step00-${SHORT}.${PREFIX}"
pipelineversion=$(git --git-dir=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/.git rev-parse origin/master)

if [ -e "${MAINDIR}/.FILE_extension.txt" ]
then
    EXT=$(cat ${MAINDIR}/.FILE_extension.txt)
else
    echo "Error: could not find ${MAINDIR}/.FILE_extension.txt"
    exit 1
fi

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=6G,h_vmem=10G,h_fsize=150G"
else
    MEM="mem_free=3G,h_vmem=5G,h_fsize=150G"
fi

if [ -e ".send_emails" ]
then
    EMAIL="e"
else
    EMAIL="a"
fi

if [ -e ".queue" ]
then
    QUEUE=$(cat .queue)
else
    QUEUE="shared"
fi

# Construct shell files
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE},${MEM}
#$ -N ${sname}
#$ -pe local 8
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -hold_jid pipeline_setup
#$ -m ${EMAIL}

echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n${pipelineversion}"

Rscript ${SH_FOLDER}/step00-merge.R -s ${MAINDIR}/SAMPLE_IDs.txt -o ${MAINDIR}/${EXPERIMENT}/${PREFIX}/merged_fastq -p ${PE} -e ${EXT} -c 8

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
