#!/bin/bash

## Usage information:
# bash step8b-callVariants.sh --help

# Define variables
TEMP=$(getopt -o x:p:r:c:l:h --long experiment:,prefix:,reference:,cores:,large:,help -n 'step8b-mergeVariantCalls' -- "$@")
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
#$ -l ${SGEQUEUE}${MEM}
#$ -N ${sname}
## See thread with Mark about -pe local at
## https://lists.johnshopkins.edu/sympa/arc/bithelp
#$ -o ./logs/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.\$TASK_ID.txt
#$ -hold_jid pipeline_setup,step8-callVariants-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

VCFS=$(cat ${MAINDIR}/samples.manifest | awk '{print "${MAINDIR}/Genotypes/"$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
vcf-merge \${VCFS} |bgzip -c > ${MAINDIR}/Genotypes/mergedVariants.vcf.gz


echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
