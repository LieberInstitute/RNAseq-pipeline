#!/bin/bash

## Usage information:
# bash step4-featureCounts.sh --help

# Define variables
TEMP=$(getopt -o x:p:s:g:r:l:h --long experiment:,prefix:,stranded:,gtf:,reference:,large:,help -n 'step4-featureCounts' -- "$@")
eval set -- "$TEMP"

STRANDED="FALSE"
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
        -s|--stranded)
            case "$2" in
                "") STRANDED="FALSE" ; shift 2;;
                *) STRANDED=$2; shift 2;;
            esac ;;
        -g|--gtf)
            case "$2" in
                "") shift 2 ;;
                *) GTF=$2 ; shift 2;;
            esac;;
        -r|--reference)
            case "$2" in
                "") shift 2 ;;
                *) hgXX=$2 ; shift 2;;
            esac;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step4-featureCounts.sh -x -p -s (default:FALSE) -g -r (hg38, hg19, mm10, rn6) -l (default:FALSE)\nLong options:\n  bash step4-featureCounts.sh --experiment --prefix --stranded (default:FALSE) --gtf --reference (hg38, hg19, mm10, rn6) --large (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="featCounts-${EXPERIMENT}"
sname="step4-${SHORT}.${PREFIX}"
CORES=8

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=12G,h_vmem=12G,h_fsize=100G"
else
    MEM="mem_free=8G,h_vmem=8G,h_fsize=100G"
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

if [ ${STRANDED} == "FALSE" ]
then
    STRANDOPTION="0"
elif [ ${STRANDED} == "forward" ]
then
    STRANDOPTION="1"
elif [ ${STRANDED} == "reverse" ]
then
    STRANDOPTION="2"
else
    echo "The option --stranded has to either be 'FALSE', 'forward' or 'reverse'."
    exit 1
fi


# File name of featureCounts output
if [ $hgXX == "mm10" ] ; then 
	FCFILE="\${ID}_Gencode.M11.${hgXX}"
elif [ $hgXX == "rn6" ] ; then 
	FCFILE="\${ID}_Ensembl.rnor6.0.${hgXX}"
elif [ $hgXX == "hg38" ] ; then 
	FCFILE="\${ID}_Gencode.v25.${hgXX}"
else 
	FCFILE="\${ID}_Gencode.v25lift37.${hgXX}"
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
#$ -tc 30
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

## remove to fix issues
module unload conda_R 

# Directories
mkdir -p ${MAINDIR}/Counts/gene
mkdir -p ${MAINDIR}/Counts/exon
mkdir -p ${MAINDIR}/Counts/junction/tmpdir

FILE1=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
if [ $PE == "TRUE" ] 
then
    FILE2=\$(awk 'BEGIN {FS="\t"} {print \$3}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
fi
ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SGE_TASK_ID}")
BAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam

if [ $PE == "TRUE" ] ; then 
	# genes	
	${SOFTWARE}/subread-1.5.0-p3-source/bin/featureCounts \
	-s ${STRANDOPTION} -p -T ${CORES} -a $GTF \
	-o ${MAINDIR}/Counts/gene/${FCFILE}_Genes.counts \$BAM
	# exons	
	${SOFTWARE}/subread-1.5.0-p3-source/bin/featureCounts \
	-s ${STRANDOPTION} -p -O -f -T ${CORES} -a $GTF \
	-o ${MAINDIR}/Counts/exon/${FCFILE}_Exons.counts \$BAM
else
	# genes	
	${SOFTWARE}/subread-1.5.0-p3-source/bin/featureCounts \
	-s ${STRANDOPTION} -T ${CORES} -a $GTF \
	-o ${MAINDIR}/Counts/gene/${FCFILE}_Genes.counts \$BAM
	# exons	
	${SOFTWARE}/subread-1.5.0-p3-source/bin/featureCounts \
	-s ${STRANDOPTION} -O -f -T ${CORES} -a $GTF \
	-o ${MAINDIR}/Counts/exon/${FCFILE}_Exons.counts \$BAM
fi
	
# junctions	
OUTJXN=${MAINDIR}/Counts/junction/\${ID}_junctions_primaryOnly_regtools.bed
OUTCOUNT=${MAINDIR}/Counts/junction/\${ID}_junctions_primaryOnly_regtools.count
TMPDIR=${MAINDIR}/Counts/junction/tmpdir
TMPBAM=\${TMPDIR}/\${ID}.bam
#filter only primary alignments
${SOFTWARE}/samtools-1.2/samtools view -@ ${CORES} -bh -F 0x100 \$BAM > \${TMPBAM}
${SOFTWARE}/samtools-1.2/samtools index \${TMPBAM}

## Load python 2.7.9 since the default one cannot run:
# python
# import site
module load python/2.7.9
${SOFTWARE}/regtools/build/regtools junctions extract -i 9 -o \${OUTJXN} \${TMPBAM}
${SOFTWARE}/bed_to_juncs_withCount < \${OUTJXN} > \${OUTCOUNT}


echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call

cat > ${MAINDIR}/.${sname}_clean.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -N ${sname}_clean
#$ -o ./logs/${SHORT}_clean.txt
#$ -e ./logs/${SHORT}_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "****"

## Delete temporary files after they have been used
rm -rf ${MAINDIR}/Counts/junction/tmpdir

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}_clean.sh"
echo $call
$call
