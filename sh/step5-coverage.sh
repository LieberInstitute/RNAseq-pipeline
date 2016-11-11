#!/bin/sh

## Usage
# ${SH_FOLDER}/step5-coverage.sh ${EXPERIMENT} ${PREFIX} ${CHRSIZES}

# Define variables
EXPERIMENT=$1
PREFIX=$2
CHRSIZES=$3

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="coverage-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"

# Directories
mkdir -p ${MAINDIR}/Coverage

# Construct shell files
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=20G,h_fsize=10G
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 40
#$ -hold_jid featCounts-${EXPERIMENT}.${PREFIX}
echo "**** Job starts ****"
date
ID=\$(awk "NR==\$SGE_TASK_ID" $FILELIST )

BAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam
BG=${MAINDIR}/Coverage/\${ID}.bedGraph
BGS=${MAINDIR}/Coverage/\${ID}.sorted.bedGraph
BW=${MAINDIR}/Coverage/\${ID}.bw

${SOFTWARE}/bedtools-2.26.0/bin/bedtools genomecov -ibam \$BAM -bga -split > \$BG
${SOFTWARE}/bedtools-2.26.0/bin/sortBed -i \$BG > \$BGS
${SOFTWARE}/bedGraphToBigWig \$BGS $CHRSIZES \$BW

rm \$BG \$BGS

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call

