#!/bin/sh

## Usage
# ${BASH_FOLDER}/step0-ercc.sh ${EXPERIMENT} ${PREFIX} ${CORES} ${LARGE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
CORES=${3-8}
LARGE=${4-"FALSE"}


SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="ercc-${EXPERIMENT}"
sname="step0-${SHORT}.${PREFIX}"

if [ ! -f "${MAINDIR}/.file_extensions.txt" ]
then
    echo "Error: could not find ${MAINDIR}/.file_extensions.txt"
    exit 1
fi

if [ ! -f "${MAINDIR}/.file_extensions.txt" ]
then
    echo "Error: could not find ${MAINDIR}/.file_extensions.txt"
    exit 1
fi

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
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE},${MEM}
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 20
#$ -hold_jid pipeline_setup,step00-merge-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date


FILEID=\$(awk "NR==\${SGE_TASK_ID}" $FILELIST )
ID=\$(basename "\${FILEID}")
mkdir -p ${MAINDIR}/Ercc/\${ID}
EXT=\$(awk "NR==\${SGE_TASK_ID}" ${MAINDIR}/.file_extensions.txt )

if [ $PE == "TRUE" ] ; then 
${SOFTWARE}/kallisto quant \
-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx \
-o ${MAINDIR}/Ercc/\${ID} -t ${CORES} --rf-stranded \
\${FILEID}_R1_001.${EXT} \${FILEID}_R2_001.${EXT}
else
${SOFTWARE}/kallisto quant \
-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx \
-o ${MAINDIR}/Ercc/\${ID} -t ${CORES} --single \${FILEID}.${EXT}

fi

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
