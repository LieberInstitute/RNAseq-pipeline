#!/bin/sh

## Usage
# ${SH_FOLDER}/step0-ercc.sh ${EXPERIMENT} ${PREFIX} ${PE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
PE=$3

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="ercc-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"

if [ -e "${MAINDIR}/.FILE_extension.txt" ]
then
    EXT=$(cat ${MAINDIR}/.FILE_extension.txt)
else
    echo "Error: could not find ${MAINDIR}/.FILE_extension.txt"
    exit 1
fi

# Construct shell files
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=5G,h_fsize=100G
#$ -N ${sname}
#$ -pe local 8
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 20
echo "**** Job starts ****"
date

echo "**** Pipeline version: latest GitHub sha ****"
git --git-dir=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/.git rev-parse origin/master


FILEID=\$(awk "NR==\$SGE_TASK_ID" $FILELIST )
ID=\$(basename "\${FILEID}")
mkdir -p ${MAINDIR}/Ercc/\${ID}

if [ $PE == "TRUE" ] ; then 
${SOFTWARE}/kallisto quant \
-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx \
-o ${MAINDIR}/Ercc/\${ID} -t 8 --rf-stranded \
\${FILEID}_R1_001.${EXT} \${FILEID}_R2_001.${EXT}
else
${SOFTWARE}/kallisto quant \
-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx \
-o ${MAINDIR}/Ercc/\${ID} -t 8 --single \${FILEID}.${EXT}

fi

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call

