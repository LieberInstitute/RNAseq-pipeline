#!/bin/bash

## Usage information:
# bash step2-trim.sh --help

# Define variables
TEMP=$(getopt -o x:p:c:l:h --long experiment:,prefix:,cores:,large:,help -n 'step2-trim' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"
CORES=8

while true; do
    case "$1" in
        -x|--experiment)
            case "$2" in
                "") shift 2 ;;
                *) EXPERIMENT=$2 ; shift 2;;
            esac;;
        -p|--prefix)
            case "$2" in
                "") shift 2 ;;
                *) PREFIX=$2 ; shift 2;;
            esac;;
        -c|--cores)
            case "$2" in
                "") CORES="8" ; shift 2;;
                *) CORES=$2; shift 2;;
            esac ;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step2-trim.sh -x -p -c (default:8) -l (default:FALSE)\nLong options:\n  bash step2-trim.sh --experiment --prefix --cores (default:8) --large (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="trim-${EXPERIMENT}"
sname="step2-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=30G,h_vmem=35G,h_fsize=100G"
else
    MEM="mem_free=20G,h_vmem=25G,h_fsize=100G"
fi

if [ -f ".send_emails" ]
then
    EMAIL="e"
else
    EMAIL="a"
fi

if [ -f ".queue" ]
then
    SGEQUEUE="$(cat .queue),"
else
    SGEQUEUE=""
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
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 5
#$ -hold_jid pipeline_setup,step1-fastqc-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"
echo "****"
echo "Sample id: \$(cat ${MAINDIR}/samples.manifest | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
FILEBASE1=\$(basename \${FILE1} | sed 's/.fq.gz//; s/.fq//; s/.fastq.gz//; s/.fastq//')
if [ $PE == "TRUE" ] 
then
    FILE2=\$(awk 'BEGIN {FS="\t"} {print \$3}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
    FILEBASE2=\$(basename \${FILE2} | sed 's/.fq.gz//; s/.fq//; s/.fastq.gz//; s/.fastq//')
fi
ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")

if [ $PE == "TRUE" ] ; then 
	REPORT1=${MAINDIR}/FastQC/Untrimmed/\${ID}/\${FILEBASE1}_fastqc/summary.txt
	REPORT2=${MAINDIR}/FastQC/Untrimmed/\${ID}/\${FILEBASE2}_fastqc/summary.txt
	RESULT1=\$(grep "Adapter Content" \$REPORT1 | cut -c1-4)
	RESULT2=\$(grep "Adapter Content" \$REPORT2 | cut -c1-4)

	if [[ \$RESULT1 == "FAIL" || \$RESULT2 == "FAIL" ]] ; then
		## trim, rerun fastQC
		echo "End 1 adapters: \$RESULT1"
		echo "End 2 adapters: \$RESULT2"
		echo "Trimming will occur."
		
		mkdir -p ${MAINDIR}/trimmed_fq
		FP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_paired.fastq
		FU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_unpaired.fastq
		RP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_paired.fastq
		RU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_unpaired.fastq
		
		## trim adapters
		java -Xmx512M -jar ${SOFTWARE}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads ${CORES} -phred33 \
		\${FILE1} \${FILE2} \$FP \$FU \$RP \$RU \
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
	REPORT1=${MAINDIR}/FastQC/Untrimmed/\${ID}/\${FILEBASE1}_fastqc/summary.txt
	RESULT1=\$(grep "Adapter Content" \$REPORT1 | cut -c1-4)

	if [[ \$RESULT1 == "FAIL" ]] ; then
		## trim, rerun fastQC
		echo "Adapters: \$RESULT1"
		echo "Trimming will occur."
		
		mkdir -p ${MAINDIR}/trimmed_fq
		OUT=${MAINDIR}/trimmed_fq/\${ID}_trimmed.fastq
		
		## trim adapters
		java -Xmx512M -jar ${SOFTWARE}/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads ${CORES} -phred33 \
		\${FILE1} \$OUT \
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
