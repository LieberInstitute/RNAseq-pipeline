#!/bin/sh

## Usage
# ${SH_FOLDER}/step1-fastqc.sh ${EXPERIMENT} ${PREFIX} ${PE} ${LARGE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
PE=$3
LARGE=${4-"FALSE"}

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="fastqc-${EXPERIMENT}"
sname="step1-${SHORT}.${PREFIX}"

if [ -e "${MAINDIR}/.FILE_extension.txt" ]
then
    EXT=$(cat ${MAINDIR}/.FILE_extension.txt)
else
    echo "Error: could not find ${MAINDIR}/.FILE_extension.txt"
    exit 1
fi

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=10G,h_vmem=14G,h_fsize=100G"
else
    MEM="mem_free=5G,h_vmem=7G,h_fsize=100G"
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

# Directories
mkdir -p ${MAINDIR}/FastQC/Untrimmed
mkdir -p ${MAINDIR}/logs

# Construct shell files
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE},${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 100
#$ -hold_jid pipeline_setup,step00-merge-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date



FILEID=\$(awk "NR==\${SGE_TASK_ID}" $FILELIST )
ID=\$(basename "\${FILEID}")
mkdir -p ${MAINDIR}/FastQC/Untrimmed/\${ID}

if [ $PE == "TRUE" ]
then 
    ${SOFTWARE}/FastQC_v0.11.5/fastqc \
\${FILEID}_R1_001.${EXT} \${FILEID}_R2_001.${EXT} \
--outdir=${MAINDIR}/FastQC/Untrimmed/\${ID} --extract
else
    ${SOFTWARE}/FastQC_v0.11.5/fastqc \${FILEID}.${EXT} \
--outdir=${MAINDIR}/FastQC/Untrimmed/\${ID} --extract
fi

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
