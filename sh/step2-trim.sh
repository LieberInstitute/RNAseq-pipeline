#!/bin/sh

## Usage
# ${SH_FOLDER}/step2-trim.sh ${EXPERIMENT} ${PREFIX} ${PE} ${FQ_FOLDER} ${EXT}

# Define variables
EXPERIMENT=$1
PREFIX=$2
PE=$3
FQ_FOLDER=$4
EXT=$5

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="trim-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"

# Construct shell files
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=15G,h_fsize=10G
#$ -N ${sname}
#$ -pe local 8
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 8
#$ -hold_jid fastqc-${EXPERIMENT}.${PREFIX}
echo "**** Job starts ****"
date
ID=\$(awk "NR==\$SGE_TASK_ID" $FILELIST )

if [ $PE == "TRUE" ] ; then 
	REPORT1=${MAINDIR}/FastQC/Untrimmed/\${ID}/\${ID}_R1_001_fastqc/summary.txt
	REPORT2=${MAINDIR}/FastQC/Untrimmed/\${ID}/\${ID}_R2_001_fastqc/summary.txt
	RESULT1=\$(grep "Adapter Content" \$REPORT1 | cut -c1-4)
	RESULT2=\$(grep "Adapter Content" \$REPORT2 | cut -c1-4)

	if [[ \$RESULT1 == "FAIL" || \$RESULT2 == "FAIL" ]] ; then
		## trim, rerun fastQC
		echo "End 1 adapters: \$RESULT1"
		echo "End 2 adapters: \$RESULT2"
		echo "Trimming will occur."
		
		mkdir -p ${MAINDIR}/trimmed_fq
		FP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_paired.fq.gz
		FU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_unpaired.fq.gz
		RP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_paired.fq.gz
		RU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_unpaired.fq.gz
		
		## trim adapters
		java -jar ${SOFTWARE}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 \
		${FQ_FOLDER}/\${ID}_R1_001.${EXT} ${FQ_FOLDER}/\${ID}_R2_001.${EXT} \$FP \$FU \$RP \$RU \
		ILLUMINACLIP:${SOFTWARE}/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10:1 \
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

		## rerun fastqc
		mkdir -p ${MAINDIR}/FastQC/Trimmed/\${ID}
		${SOFTWARE}/FastQC_v0.11.5/fastqc \
		\$FP \$FU \$RP \$RU \
		--outdir=${MAINDIR}/FastQC/Trimmed/\${ID} --extract
	else
		echo "No trimming required!"
	fi

else
	## reads are single-end
	REPORT1=${MAINDIR}/FastQC/Untrimmed/\${ID}/\${ID}_fastqc/summary.txt
	RESULT1=\$(grep "Adapter Content" \$REPORT1 | cut -c1-4)

	if [ \$RESULT1 == "FAIL" ] ; then
		## trim, rerun fastQC
		echo "Adapters: \$RESULT1"
		echo "Trimming will occur."
		
		mkdir -p ${MAINDIR}/trimmed_fq
		OUT=${MAINDIR}/trimmed_fq/\${ID}_trimmed.${EXT}
		
		## trim adapters
		java -jar ${SOFTWARE}/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 \
		${FQ_FOLDER}/\${ID}.${EXT} \$OUT \
		ILLUMINACLIP:${SOFTWARE}/Trimmomatic-0.36/adapters/TruSeq2-SE.fa:2:30:10:1 \
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		## rerun fastqc
		mkdir -p ${MAINDIR}/FastQC/Trimmed/\${ID}
		${SOFTWARE}/FastQC_v0.11.5/fastqc \$OUT \
		--outdir=${MAINDIR}/FastQC/Trimmed/\${ID} --extract
	else
		echo "No trimming required!"
	fi
fi

	
echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
