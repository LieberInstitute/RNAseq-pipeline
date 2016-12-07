#!/bin/sh

## Usage
# ${SH_FOLDER}/step1-fastqc.sh ${EXPERIMENT} ${PREFIX} ${PE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
PE=$3

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="fastqc-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"

if [ -e "${MAINDIR}/.FILE_extension.txt" ]
then
    EXT=$(cat ${MAINDIR}/.FILE_extension.txt)
else
    echo "Error: could not find ${MAINDIR}/.FILE_extension.txt"
    exit 1
fi

# Directories
mkdir -p ${MAINDIR}/FastQC/Untrimmed
mkdir -p ${MAINDIR}/logs

# Construct shell files
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=7G,h_fsize=10G
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 100
echo "**** Job starts ****"
date

echo "**** Pipeline version: latest GitHub sha ****"
git --git-dir=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/.git rev-parse origin/master

FILEID=\$(awk "NR==\$SGE_TASK_ID" $FILELIST )
ID=\$(basename "\${FILEID}")
mkdir -p ${MAINDIR}/FastQC/Untrimmed/\${ID}

if [ $PE == "TRUE" ] ; then 
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
