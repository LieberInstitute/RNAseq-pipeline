#!/bin/bash

## Usage information:
# bash step8b-callVariants.sh --help

# Define variables
TEMP=$(getopt -o x:p:h --long experiment:,prefix:,help -n 'step8b-mergeVariantCalls' -- "$@")
eval set -- "$TEMP"

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
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step8-callVariants.sh -x -p \nLong options:\n  bash step8-callVariants.sh --experiment --prefix"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="mergeVariantCalls-${EXPERIMENT}"
sname="step8b-${SHORT}.${PREFIX}"

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


# Construct shell files
FILELIST=${MAINDIR}/samples.manifest
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${SGEQUEUE}mem_free=8G,h_vmem=10G
#$ -N ${sname}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -hold_jid pipeline_setup,step8-callVariants-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "****"

VCFS=\$(cat ${MAINDIR}/samples.manifest | awk '{print "${MAINDIR}/Genotypes/"\$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
module load htslib
vcf-merge \${VCFS} | bgzip -c > ${MAINDIR}/Genotypes/mergedVariants.vcf.gz

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
