#!/bin/sh

## Usage
# ${SH_FOLDER}/step5-coverage.sh ${EXPERIMENT} ${PREFIX} ${CHRSIZES} ${LARGE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
CHRSIZES=$3
LARGE=${4-"FALSE"}

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="coverage-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=20G,h_vmem=40G,h_fsize=100G"
else
    MEM="mem_free=10G,h_vmem=20G,h_fsize=100G"
fi

if [ -e ".send_emails" ]
then
    EMAIL="e"
else
    EMAIL="a"
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
#$ -l ${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 40
#$ -hold_jid pipeline_setup,featCounts-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** Pipeline version: latest GitHub sha ****"
git --git-dir=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/.git rev-parse origin/master

FILEID=\$(awk "NR==\${SGE_TASK_ID}" $FILELIST )
ID=\$(basename "\${FILEID}")

## Normalizing bigwigs to 40 million 100 bp reads
module load python/2.7.9
module load ucsctools
python ~/.local/bin/bam2wig.py -s ${CHRSIZES} -i ${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam -t 4000000000 -o ${MAINDIR}/Coverage/\${ID}

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
