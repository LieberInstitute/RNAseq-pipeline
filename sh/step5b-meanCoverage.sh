#!/bin/sh

## Usage
# ${SH_FOLDER}/step5b-meanCoverage.sh ${EXPERIMENT} ${PREFIX} ${CHRSIZES} ${LARGE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
CHRSIZES=$3
LARGE=${4-"FALSE"}

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="meanCoverage-${EXPERIMENT}"
sname="step5b-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=200G,h_vmem=240G,h_fsize=100G"
else
    MEM="mem_free=100G,h_vmem=120G,h_fsize=100G"
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
mkdir -p ${MAINDIR}/Coverage

# Construct shell files
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE},${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -hold_jid pipeline_setup,step5-coverage-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date


## Locate normalized BigWig files and concatenate them in a space separated list
BIGWIGS=\$(while read line; do ID=\$(basename \${line}); echo "${MAINDIR}/Coverage/\${ID}.bw"; done < ${FILELIST} | paste -sd " ")

## Create mean of normalized bigwigs
module load wiggletools/default
module load ucsctools
wiggletools write ${MAINDIR}/Coverage/mean.wig mean \${BIGWIGS}
wigToBigWig ${MAINDIR}/Coverage/mean.wig ${CHRSIZES} ${MAINDIR}/Coverage/mean.bw

## Remove temp files
rm ${MAINDIR}/Coverage/mean.wig

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
