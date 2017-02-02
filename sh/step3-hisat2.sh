#!/bin/sh

## Usage
# ${BASH_FOLDER}/step3-hisat2.sh ${EXPERIMENT} ${PREFIX} ${HISATIDX} ${BED} ${CORES} ${LARGE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
HISATIDX=$3
BED=$4
CORES=${5-8}
LARGE=${6-"FALSE"}

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="hisat2-${EXPERIMENT}"
sname="step3-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=10G,h_vmem=14G,h_fsize=100G"
else
    MEM="mem_free=5G,h_vmem=7G,h_fsize=100G"
fi

if [ -f ".send_emails" ]
then
    EMAIL="e"
else
    EMAIL="a"
fi

if [ -f ".queue" ]
then
    QUEUE="$(cat .queue),"
fi

if [ -f ".paired_end" ]
then
    PE="TRUE"
else
    PE="FALSE"
fi

# Construct shell files
FILELIST=${MAINDIR}/samples.manifest
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE}${MEM}
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 15
#$ -hold_jid pipeline_setup,step2-trim-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

# Directories
mkdir -p ${MAINDIR}/HISAT2_out/align_summaries
mkdir -p ${MAINDIR}/HISAT2_out/infer_strandness

## Locate file and ids
FILE1=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
if [ $PE == "TRUE" ] 
then
    FILE2=\$(awk 'BEGIN {FS="\t"} {print \$3}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
fi
ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")

if [ -f ${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_paired.fastq ] ; then
	## Trimmed, paired-end
	echo "HISAT2 alignment run on trimmed paired-end reads"
	FP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_paired.fastq
	FU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_unpaired.fastq
	RP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_paired.fastq
	RU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_unpaired.fastq
	
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p ${CORES} \
	-x $HISATIDX -1 \$FP -2 \$RP -U \${FU},\${RU} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam --rna-strandness RF --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt
	
elif  [ -f ${MAINDIR}/trimmed_fq/\${ID}_trimmed.fastq ] ; then
	## Trimmed, single-end
	echo "HISAT2 alignment run on trimmed single-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p ${CORES} \
	-x $HISATIDX -U ${MAINDIR}/trimmed_fq/\${ID}_trimmed.fastq \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt

elif [ $PE == "TRUE" ] ; then
	## Untrimmed, pair-end
	echo "HISAT2 alignment run on original untrimmed paired-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p ${CORES} \
	-x $HISATIDX -1 \${FILE1} -2 \${FILE2} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam --rna-strandness RF --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt

else
	## Untrimmed, single-end
	echo "HISAT2 alignment run on original untrimmed single-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p ${CORES} \
	-x $HISATIDX -U \${FILE1} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt
fi


###sam to bam
SAM=${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam
ORIGINALBAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.bam
SORTEDBAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted

#filter unmapped segments
${SOFTWARE}/samtools-1.2/samtools view -bh -F 4 \${SAM} > \${ORIGINALBAM}
${SOFTWARE}/samtools-1.2/samtools sort -@ ${CORES} \${ORIGINALBAM} \${SORTEDBAM}
${SOFTWARE}/samtools-1.2/samtools index \${SORTEDBAM}.bam

## Clean up
rm \${SAM}
rm \${ORIGINALBAM}

## Run infer experiment
module load python/2.7.9
~/.local/bin/infer_experiment.py -i \${SORTEDBAM} -r ${BED} 1> ${MAINDIR}/HISAT2_out/infer_strandness/\${ID}.txt 2>&1

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call

## Process the output from infer_experiment.py for all samples
SHORT="infer-strandness-${EXPERIMENT}"
sname="step3b-${SHORT}.${PREFIX}"
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -hold_jid pipeline_setup,step3-hisat2-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

## Process the infer experiment info
Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step3b_infer_strandness.R

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
