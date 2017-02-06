#!/bin/bash

## Usage information:
# bash step5b-meanCoverage.sh --help

# Define variables
TEMP=$(getopt -o x:p:s:l:h --long experiment:,prefix:,chrsizes:,large:,help -n 'step5b-meanCoverage' -- "$@")
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
            echo -e "Usage:\nShort options:\n  bash step5b-meanCoverage.sh -x -p -s -l (default:FALSE)\nLong options:\n  bash step5b-meanCoverage.sh --experiment --prefix --chrsizes --large (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="meanCoverage-${EXPERIMENT}"
sname="step5b-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=200G,h_vmem=240G,h_fsize=100G"
else
    MEM="mem_free=100G,h_vmem=120G,h_fsize=100G"
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
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -hold_jid pipeline_setup,step5-coverage-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

## Load required software
module load wiggletools/default
module load ucsctools

## Read the strandness information
STRANDRULE=\$(cat inferred_strandness_pattern.txt)

if [ \${STRANDRULE} == "none" ]
then
    ## Locate normalized BigWig files and concatenate them in a space separated list
    BIGWIGS=\$(while read line; do ID=\$(basename \${line}); echo "${MAINDIR}/Coverage/\${ID}.bw"; done < ${FILELIST} | paste -sd " ")
    
    ## Create mean of normalized bigwigs
    wiggletools write ${MAINDIR}/Coverage/mean.wig mean \${BIGWIGS}
    wigToBigWig ${MAINDIR}/Coverage/mean.wig ${CHRSIZES} ${MAINDIR}/Coverage/mean.bw
else
    for strand in Forward Reverse
    do
        echo "Processing strand \${strand}"
        ## Locate normalized BigWig files and concatenate them in a space separated list
        BIGWIGS=\$(while read line; do ID=\$(basename \${line}); echo "${MAINDIR}/Coverage/\${ID}.\${strand}.bw"; done < ${FILELIST} | paste -sd " ")
    
        ## Create mean of normalized bigwigs
        wiggletools write ${MAINDIR}/Coverage/mean.\${strand}.wig mean \${BIGWIGS}
        wigToBigWig ${MAINDIR}/Coverage/mean.\${strand}.wig ${CHRSIZES} ${MAINDIR}/Coverage/mean.\${strand}.bw
    done
fi

## Remove temp files
rm ${MAINDIR}/Coverage/mean*.wig

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
