#!/bin/bash

## Usage information:
# bash step5-coverage.sh --help

# Define variables
TEMP=$(getopt -o x:p:s:l:h --long experiment:,prefix:,chrsizes:,large:,help -n 'step5-coverage' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"

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
        -s|--chrsizes)
            case "$2" in
                "") shift 2 ;;
                *) CHRSIZES=$2 ; shift 2;;
            esac;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step5-coverage.sh -x -p -s -l (default:FALSE)\nLong options:\n  bash step5-coverage.sh --experiment --prefix --chrsizes --large (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="coverage-${EXPERIMENT}"
sname="step5-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=18G,h_vmem=20G,h_fsize=100G"
else
    MEM="mem_free=12G,h_vmem=16G,h_fsize=100G"
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

# Directories
mkdir -p ${MAINDIR}/Coverage

# Construct shell files
FILELIST=${MAINDIR}/samples.manifest
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 40
#$ -hold_jid pipeline_setup,step3-hisat2-${EXPERIMENT}.${PREFIX},step3b-infer-strandness-${EXPERIMENT}.${PREFIX}
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

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

FILE1=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
if [[ $PE == "TRUE" ]]
then
    FILE2=\$(awk 'BEGIN {FS="\t"} {print \$3}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
fi
ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")
STRANDRULE=\$(cat inferred_strandness_pattern.txt)

## Normalizing bigwigs to 40 million 100 bp reads
module load python/2.7.9
module load ucsctools

## Can only use -d when the data is stranded
if [ "\${STRANDRULE}" == "none" ]
then
    python ~/.local/bin/bam2wig.py -s ${CHRSIZES} -i ${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam -t 4000000000 -o ${MAINDIR}/Coverage/\${ID}
elif [ "\${STRANDRULE}" == "1++,1--,2+-,2-+" ]
then
    python ~/.local/bin/bam2wig.py -s ${CHRSIZES} -i ${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam -t 4000000000 -o ${MAINDIR}/Coverage/\${ID} -d "1++,1--,2+-,2-+"
elif [ "\${STRANDRULE}" == "1+-,1-+,2++,2--"]
then
    python ~/.local/bin/bam2wig.py -s ${CHRSIZES} -i ${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam -t 4000000000 -o ${MAINDIR}/Coverage/\${ID} -d "1+-,1-+,2++,2–"
elif [ "\${STRANDRULE}" == "++,--" ]
then
    python ~/.local/bin/bam2wig.py -s ${CHRSIZES} -i ${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam -t 4000000000 -o ${MAINDIR}/Coverage/\${ID} -d "++,--"
elif [ "\${STRANDRULE}" == "+-,-+" ]
then
    python ~/.local/bin/bam2wig.py -s ${CHRSIZES} -i ${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam -t 4000000000 -o ${MAINDIR}/Coverage/\${ID} -d "+-,-+"
else
    echo "Found unexpected value in inferred_strandness_pattern.txt: \${STRANDRULE}"
fi

## Remove temp files
rm ${MAINDIR}/Coverage/\${ID}*.wig

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
