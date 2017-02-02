#!/bin/sh

## Usage
# ${BASH_FOLDER}/step4-featureCounts.sh ${EXPERIMENT} ${PREFIX} ${STRANDED} ${GTF} ${hgXX} ${CORES} ${LARGE}

# Define variables
EXPERIMENT=$1
PREFIX=$2
STRANDED=$3
GTF=$4
hgXX=$5
CORES=${6-8}
LARGE=${7-"FALSE"}

SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}
SHORT="featCounts-${EXPERIMENT}"
sname="step4-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=20G,h_vmem=24G,h_fsize=100G"
else
    MEM="mem_free=10G,h_vmem=12G,h_fsize=100G"
fi

if [ -f ".send_emails" ]
then
    EMAIL="e"
else
    EMAIL="a"
fi

if [ -f ".queue" ]
then
    QUEUE="$(cat .queue),"
fi

if [ -f ".paired_end" ]
then
    PE="TRUE"
else
    PE="FALSE"
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
FILELIST=${MAINDIR}/SAMPLE_IDs.txt
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l ${QUEUE}${MEM}
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 10
#$ -hold_jid pipeline_setup,step3-hisat2-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date


# Directories
mkdir -p ${MAINDIR}/Counts/gene
mkdir -p ${MAINDIR}/Counts/exon
mkdir -p ${MAINDIR}/Counts/junction/tmpdir

FILEID=\$(awk "NR==\${SGE_TASK_ID}" $FILELIST )
ID=\$(basename "\${FILEID}")
BAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted.bam

if [ $PE == "TRUE" ] ; then 
	# genes	
	${SOFTWARE}/subread-1.5.0-p3-source/bin/featureCounts \
	-s 2 -p -T ${CORES} -a $GTF \
	-o ${MAINDIR}/Counts/gene/${FCFILE}_Genes.counts \$BAM
	# exons	
	${SOFTWARE}/subread-1.5.0-p3-source/bin/featureCounts \
	-s 2 -p -O -f -T ${CORES} -a $GTF \
	-o ${MAINDIR}/Counts/exon/${FCFILE}_Exons.counts \$BAM
else
	# genes	
	${SOFTWARE}/subread-1.5.0-p3-source/bin/featureCounts \
	-T ${CORES} -a $GTF \
	-o ${MAINDIR}/Counts/gene/${FCFILE}_Genes.counts \$BAM
	# exons	
	${SOFTWARE}/subread-1.5.0-p3-source/bin/featureCounts \
	-O -f -T ${CORES} -a $GTF \
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
#$ -o ./logs/${SHORT}_clean.o.txt
#$ -e ./logs/${SHORT}_clean.e.txt
#$ -hold_jid pipeline_setup,step4-featCounts-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date


## Delete temporary files after they have been used
rm -rf ${MAINDIR}/Counts/junction/tmpdir

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}_clean.sh"
echo $call
$call
