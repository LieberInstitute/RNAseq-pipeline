#!/bin/sh

## Usage
# ${SH_FOLDER}/step6-makeRobjects.sh ${EXPERIMENT} ${PREFIX} ${hgXX} ${SH_FOLDER} ${PE} ${ERCC} ${LARGE} ${FULLCOV}

# Define variables
EXPERIMENT=$1
PREFIX=$2
hgXX=$3
SH_FOLDER=$4
PE=$5
ERCC=$6
LARGE=${7-"FALSE"}
FULLCOV=${8-"FALSE"}

SHORT="Rcounts-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"
SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=20G,h_vmem=24G,h_fsize=200G"
else
    MEM="mem_free=5G,h_vmem=6G,h_fsize=200G"
fi

if [ -e ".send_emails" ]
then
    EMAIL="e"
else
    EMAIL="a"
fi

if [ $hgXX == "mm10" ]; then SPEC="mouse";
elif [ $hgXX == "rn6" ]; then SPEC="rat";
else SPEC="human";
fi

cp ${SH_FOLDER}/create_count_objects-${SPEC}.R ${MAINDIR}/.create_count_objects-${SPEC}.R
cp ${SH_FOLDER}/create_fullCov_object.R ${MAINDIR}/.create_fullCov_object.R

# Construct shell files
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l ${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -hold_jid pipeline_setup,coverage-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** Pipeline version: latest GitHub sha ****"
git --git-dir=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/.git rev-parse origin/master

Rscript ${MAINDIR}/.create_count_objects-${SPEC}.R -o ${hgXX} -m ${MAINDIR} -e ${EXPERIMENT} -p ${PREFIX} -l ${PE} -c ${ERCC}

## Don't create the fullcoverage object by default
Rscript ${MAINDIR}/.create_fullCov_object.R -o ${hgXX} -m ${MAINDIR} -e ${EXPERIMENT} -p ${PREFIX} -l ${PE} -f ${FULLCOV}


echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
