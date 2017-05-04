#!/bin/bash

## Usage information:
# bash step8-callVariants.sh --help

# Define variables
TEMP=$(getopt -o x:p:r:h --long experiment:,prefix:,reference:,help -n 'step8-callVariants' -- "$@")
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
		-r|--reference)
            case "$2" in
                "") shift 2 ;;
                *) hgXX=$2 ; shift 2;;
            esac;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step8-callVariants.sh -x -p -r (hg38, hg19, mm10, rn6) \nLong options:\n  bash step8-callVariants.sh --experiment --prefix --reference (hg38, hg19, mm10, rn6)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="callVariants-${EXPERIMENT}"
sname="step8-${SHORT}.${PREFIX}"


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
#$ -l ${SGEQUEUE}mem_free=2G,h_vmem=3G,h_fsize=100G
#$ -N ${sname}
#$ -o ./logs/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 100
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
echo "****"
echo "Sample id: \$(cat ${MAINDIR}/samples.manifest | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")"
echo "****"

mkdir -p ${MAINDIR}/Genotypes/

ID=\$(cat ${MAINDIR}/samples.manifest | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")
BAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam

SNPTMP=${MAINDIR}/Genotypes/\${ID}_tmp.vcf
SNPOUTGZ=${MAINDIR}/Genotypes/\${ID}.vcf.gz

module load bcftools
${SOFTWARE}/samtools-1.2/samtools mpileup -l ${BEDFILE} -AB -q0 -Q13 -d1000000 -uf ${FAFILE} \${BAM} -o \${SNPTMP}
bcftools call -mv -Oz \${SNPTMP} > \${SNPOUTGZ}

module load htslib
tabix -p vcf \${SNPOUTGZ}

rm \${SNPTMP}

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
