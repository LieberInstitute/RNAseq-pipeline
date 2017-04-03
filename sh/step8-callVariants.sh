#!/bin/bash

## Usage information:
# bash step8-callVariants.sh --help

# Define variables
TEMP=$(getopt -o x:p:r:c:l:h --long experiment:,prefix:,reference:,cores:,large:,help -n 'step8-callVariants' -- "$@")
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
		-r|--reference)
            case "$2" in
                "") shift 2 ;;
                *) hgXX=$2 ; shift 2;;
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
            echo -e "Usage:\nShort options:\n  bash step8-callVariants.sh -x -p -r (hg38, hg19, mm10, rn6) -c (default:8) -l (default:FALSE)\nLong options:\n  bash step8-callVariants.sh --experiment --prefix --reference (hg38, hg19, mm10, rn6) --cores (default:8) --large (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="callVariants-${EXPERIMENT}"
sname="step8-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=10G,h_vmem=12,h_fsize=100G"
else
    MEM="mem_free=5G,h_vmem=8G,h_fsize=100G"
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



# Reference genome files
if [ $hgXX == "hg38" ] ; then 
	BEDFILE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed
	FAFILE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/GRCh38.primary_assembly.genome.fa
elif [ $hgXX == "hg19" ] ; then 
	BEDFILE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg19.bed
	FAFILE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh37_hg19/GRCh37.primary_assembly.genome.fa
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
## See thread with Mark about -pe local at
## https://lists.johnshopkins.edu/sympa/arc/bithelp
# -pe local 1
#$ -o ./logs/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 50
#$ -hold_jid pipeline_setup,step3-hisat2-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

mkdir -p ${MAINDIR}/Genotypes/

ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")
BAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam

SNPBED=${BEDFILE}
GENOME=${FAFILE}
SNPTMP=${MAINDIR}/Genotypes/\${ID}_calledVariants.tmp
SNPOUT=${MAINDIR}/Genotypes/\${ID}_calledVariants.txt
module load samtools
samtools mpileup -l \${SNPBED} -AB -q0 -Q0 -d1000000 -t DP -f \${GENOME} \${BAM} | /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/pileVar.pl > \${SNPOUT}

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
