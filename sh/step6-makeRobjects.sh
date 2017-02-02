#!/bin/sh

## Usage
# ${BASH_FOLDER}/step6-makeRobjects.sh ${EXPERIMENT} ${PREFIX} ${hgXX} ${ERCC} ${CORES} ${LARGE} ${FULLCOV} ${BASH_FOLDER}

# Define variables
EXPERIMENT=$1
PREFIX=$2
hgXX=$3
ERCC=$4
CORES=${5-8}
LARGE=${6-"FALSE"}
FULLCOV=${7-"FALSE"}
BASH_FOLDER=${8-"/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh"}

SHORT="Rcounts-${EXPERIMENT}"
sname="step6-${SHORT}.${PREFIX}"
SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=20G,h_vmem=24G,h_fsize=200G"
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
    QUEUE="$(cat .queue),"
fi

if [ -f ".paired_end" ]
then
    PE="TRUE"
else
    PE="FALSE"
fi

if [ $hgXX == "mm10" ]; then SPEC="mouse";
elif [ $hgXX == "rn6" ]; then SPEC="rat";
else SPEC="human";
fi

cp ${BASH_FOLDER}/create_count_objects-${SPEC}.R ${MAINDIR}/.create_count_objects-${SPEC}.R
cp ${BASH_FOLDER}/create_fullCov_object.R ${MAINDIR}/.create_fullCov_object.R

# Construct shell files
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe local ${CORES}
#$ -l ${QUEUE}${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -hold_jid pipeline_setup,step4-featCounts-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date


Rscript ${MAINDIR}/.create_count_objects-${SPEC}.R -o ${hgXX} -m ${MAINDIR} -e ${EXPERIMENT} -p ${PREFIX} -l ${PE} -c ${ERCC} -t ${CORES}

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call

if [[ ${FULLCOV} == "TRUE" ]]
then
    SHORT="fullCov-${EXPERIMENT}"
    sname="step6-${SHORT}.${PREFIX}"
    # Construct shell files
    echo "Creating script ${sname}"
    
    cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe local ${CORES}
#$ -l ${QUEUE}${MEM}
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -hold_jid pipeline_setup,step5-coverage-${EXPERIMENT}.${PREFIX}
#$ -m ${EMAIL}
echo "**** Job starts ****"
date


## Don't create the fullCoverage object by default:
## it's not needed since we create the mean bigwig already and can use it
## to define the ERs with derfinder::findRegions(), then use the resulting
## GRanges object and the paths to the BigWig files in
## derfinder::getRegionCoverage(fullCov = NULL, files = bigWigs, regions = outputFrom_findRegions)
## or alternatively write the regions to a BED file with rtracklayer,
## create the counts with bwtool and then read them into R manually
Rscript ${MAINDIR}/.create_fullCov_object.R -o ${hgXX} -m ${MAINDIR} -e ${EXPERIMENT} -p ${PREFIX} -l ${PE} -f ${FULLCOV} -c ${CORES}


echo "**** Job ends ****"
date
EOF

    call="qsub .${sname}.sh"
    echo $call
    $call
fi
