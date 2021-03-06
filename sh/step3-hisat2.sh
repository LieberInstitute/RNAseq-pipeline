#!/bin/bash

## Usage information:
# bash step3-hisat2.sh --help

# Define variables
TEMP=$(getopt -o x:p:i:b:l:s:u:h --long experiment:,prefix:,index:,bed:,large:,stranded:,unaligned:,help -n 'step3-hisat2' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"
STRANDED="FALSE"
UNALIGNED="FALSE"

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
        -i|--index)
            case "$2" in
                "") shift 2 ;;
                *) HISATIDX=$2 ; shift 2;;
            esac;;
        -b|--bed)
            case "$2" in
                "") shift 2 ;;
                *) BED=$2 ; shift 2;;
            esac;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -s|--stranded)
            case "$2" in
                "") STRANDED="FALSE" ; shift 2;;
                *) STRANDED=$2; shift 2;;
            esac ;;
        -u|--unaligned)
            case "$2" in
                "") UNALIGNED="FALSE" ; shift 2;;
                *) UNALIGNED=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step3-hisat2.sh -x -p -i -b -l (default:FALSE) -s (default:FALSE) -u (default:FALSE)\nLong options:\n  bash step3-hisat2.sh --experiment --prefix --index --bed --large (default:FALSE) --stranded (default:FALSE) --unaligned (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="hisat2-${EXPERIMENT}"
sname="step3-${SHORT}.${PREFIX}"
CORES=8

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=12G,h_vmem=14G,h_fsize=200G"
else
    MEM="mem_free=5G,h_vmem=6G,h_fsize=200G"
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
    if [ ${STRANDED} == "FALSE" ]
    then
        STRANDOPTION=""
    elif [ ${STRANDED} == "forward" ]
    then
        STRANDOPTION="--rna-strandness FR"
    elif [ ${STRANDED} == "reverse" ]
    then
        STRANDOPTION="--rna-strandness RF"
    else
        echo "The option --stranded has to either be 'FALSE', 'forward' or 'reverse'."
        exit 1
    fi
else
    PE="FALSE"
    if [ ${STRANDED} == "FALSE" ]
    then
        STRANDOPTION=""
    elif [ ${STRANDED} == "forward" ]
    then
        STRANDOPTION="--rna-strandness F"
    elif [ ${STRANDED} == "reverse" ]
    then
        STRANDOPTION="--rna-strandness R"
    else
        echo "The option --stranded has to either be 'FALSE', 'forward' or 'reverse'."
        exit 1
    fi
fi

if [ ${UNALIGNED} == "FALSE" ]
then
    UNALIGNEDOPT=""
elif [ ${UNALIGNED} == "TRUE" ]
then
    mkdir -p ${MAINDIR}/HISAT2_out/unaligned
    UNALIGNEDOPT="--un-conc ${MAINDIR}/HISAT2_out/unaligned/\${ID}.fastq"
else
    echo "The option --unaligned has to either be 'FALSE' or 'TRUE'"
    exit 1
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
#$ -tc 15
#$ -hold_jid pipeline_setup,step2-trim-${EXPERIMENT}.${PREFIX}
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

# Directories
mkdir -p ${MAINDIR}/HISAT2_out/align_summaries
mkdir -p ${MAINDIR}/HISAT2_out/infer_strandness

if [ ${UNALIGNED} == "TRUE" ]
then
    mkdir -p ${MAINDIR}/HISAT2_out/unaligned
fi

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
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam ${STRANDOPTION} --phred33 \
    ${UNALIGNEDOPT} \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt
	
elif  [ -f ${MAINDIR}/trimmed_fq/\${ID}_trimmed.fastq ] ; then
	## Trimmed, single-end
	echo "HISAT2 alignment run on trimmed single-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p ${CORES} \
	-x $HISATIDX -U ${MAINDIR}/trimmed_fq/\${ID}_trimmed.fastq \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam ${STRANDOPTION} --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt

elif [ $PE == "TRUE" ] ; then
	## Untrimmed, pair-end
	echo "HISAT2 alignment run on original untrimmed paired-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p ${CORES} \
	-x $HISATIDX -1 \${FILE1} -2 \${FILE2} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam ${STRANDOPTION} --phred33 \
    ${UNALIGNEDOPT} \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt

else
	## Untrimmed, single-end
	echo "HISAT2 alignment run on original untrimmed single-end reads"
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p ${CORES} \
	-x $HISATIDX -U \${FILE1} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam ${STRANDOPTION} --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt
fi


###sam to bam
SAM=${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam
ORIGINALBAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.bam
SORTEDBAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted

#filter unmapped segments
echo "**** Filtering unmapped segments ****"
date
${SOFTWARE}/samtools-1.2/samtools view -bh -F 4 \${SAM} > \${ORIGINALBAM}
${SOFTWARE}/samtools-1.2/samtools sort -@ ${CORES} \${ORIGINALBAM} \${SORTEDBAM}
${SOFTWARE}/samtools-1.2/samtools index \${SORTEDBAM}.bam

## Clean up
rm \${SAM}
rm \${ORIGINALBAM}

## Run infer experiment
echo "**** Inferring strandedness with infer_experiment.py ****"
date
module load python/2.7.9
~/.local/bin/infer_experiment.py -i \${SORTEDBAM}.bam -r ${BED} 1> ${MAINDIR}/HISAT2_out/infer_strandness/\${ID}.txt 2>&1

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
#$ -N ${sname}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -hold_jid pipeline_setup,step3-hisat2-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "****"

## Process the infer experiment info
Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step3b_infer_strandness.R -o "HISAT2_out/infer_strandness" -p "inferred_strandness_pattern.txt"

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
