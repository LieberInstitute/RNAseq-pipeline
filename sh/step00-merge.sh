#!/bin/sh

## Usage
# ${BASH_FOLDER}/step00-merge.sh ${EXPERIMENT} ${PREFIX} ${CORES} ${LARGE} ${BASH_FOLDER}

# Define variables
EXPERIMENT=$1
PREFIX=$2
CORES=${3-8}
LARGE=${4-"FALSE"}
BASH_FOLDER=${5-"/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh"}

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
    QUEUE=$(cat .queue)
else
    QUEUE="shared"
fi

if [ -f ".paired_end" ]
then
    PE="TRUE"
else
    PE="FALSE"
fi

# Construct shell files
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE},${MEM}
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -hold_jid pipeline_setup
#$ -m ${EMAIL}

echo "**** Job starts ****"
date

Rscript ${BASH_FOLDER}/step00-merge.R -s ${MAINDIR}/SAMPLE_IDs.txt -o ${MAINDIR}/${EXPERIMENT}/${PREFIX}/merged_fastq -c ${CORES}

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
