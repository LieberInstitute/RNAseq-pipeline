#!/bin/sh

## Usage
# ${SH_FOLDER}/step00-merge.sh ${EXPERIMENT} ${PREFIX} ${PE} ${SH_FOLDER}

# Define variables
EXPERIMENT=$1
PREFIX=$2
PE=$3
SH_FOLDER=$4

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="merge-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"

if [ -e "${MAINDIR}/.FILE_extension.txt" ]
then
    EXT=$(cat ${MAINDIR}/.FILE_extension.txt)
else
    echo "Error: could not find ${MAINDIR}/.FILE_extension.txt"
    exit 1
fi

# Construct shell files
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=5G,h_fsize=150G
#$ -N ${sname}
#$ -pe local 8
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt

echo "**** Job starts ****"
date

echo "**** Pipeline version: latest GitHub sha ****"
git --git-dir=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/.git rev-parse origin/master

Rscript ${SH_FOLDER}/step00-merge.R -s ${MAINDIR}/SAMPLE_IDs.txt -m ${MAINDIR} -p ${PE} -e ${EXT} -c 8

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
