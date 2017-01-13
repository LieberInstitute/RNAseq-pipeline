#!/bin/sh

## Usage
# ${SH_FOLDER}/step3-hisat2.sh ${EXPERIMENT} ${PREFIX} ${PE} ${HISATIDX} ${LARGE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
PE=$3
HISATIDX=$4
LARGE=${5-"FALSE"}

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="hisat2-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"
pipelineversion=$(git --git-dir=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/.git rev-parse origin/master)

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
    QUEUE=""
fi

# Directories
mkdir -p ${MAINDIR}/HISAT2_out/align_summaries

# Construct shell files
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE} ${MEM}
#$ -N ${sname}
#$ -pe local 8
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 15
#$ -hold_jid pipeline_setup,trim-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n${pipelineversion}"

FILEID=\$(awk "NR==\${SGE_TASK_ID}" $FILELIST )
ID=\$(basename "\${FILEID}")

if [ -e ${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_paired.fq.gz ] ; then
	## Trimmed, paired-end
	echo "HISAT2 alignment run on trimmed paired-end reads"
	FP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_paired.fq.gz
	FU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_unpaired.fq.gz
	RP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_paired.fq.gz
	RU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_unpaired.fq.gz
	
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p 8 \
	-x $HISATIDX -1 \$FP -2 \$RP -U \${FU},\${RU} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam --rna-strandness RF --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt
	
elif  [ -e ${MAINDIR}/trimmed_fq/\${ID}_trimmed.fq.gz ] ; then
	## Trimmed, single-end
	echo "HISAT2 alignment run on trimmed single-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p 8 \
	-x $HISATIDX -U ${MAINDIR}/trimmed_fq/\${ID}_trimmed.fq.gz \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt

elif [ $PE == "TRUE" ] ; then
	## Untrimmed, pair-end
	echo "HISAT2 alignment run on original untrimmed paired-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p 8 \
	-x $HISATIDX -1 \${FILEID}_R1_001.${EXT} -2 \${FILEID}_R2_001.${EXT} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam --rna-strandness RF --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt

else
	## Untrimmed, single-end
	echo "HISAT2 alignment run on original untrimmed single-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p 8 \
	-x $HISATIDX -U \${FILEID}.${EXT} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt
fi


###sam to bam
SAM=${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam
BAMACC=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.bam
BAMS=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted

#filter unmapped segments
${SOFTWARE}/samtools-1.2/samtools view -bh -F 4 \$SAM > \$BAMACC
${SOFTWARE}/samtools-1.2/samtools sort -@ 8 \$BAMACC \$BAMS
${SOFTWARE}/samtools-1.2/samtools index \${BAMS}.bam

rm \$SAM

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
