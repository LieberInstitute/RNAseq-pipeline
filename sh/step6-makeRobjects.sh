#!/bin/sh

## Usage
# ${SH_FOLDER}/step6-makeRobjects.sh ${EXPERIMENT} ${PREFIX} ${hgXX} ${SH_FOLDER} ${PE} ${ERCC}

# Define variables
EXPERIMENT=$1
PREFIX=$2
hgXX=$3
SH_FOLDER=$4
PE=$5
ERCC=$6

SHORT="Rcounts-${EXPERIMENT}"
sname="${SHORT}.${PREFIX}"
SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=${PWD}

if [ $hgXX == "mm10" ]; then SPEC="mouse";
elif [ $hgXX == "rn6" ]; then SPEC="rat";
else SPEC="human";
fi

cp ${SH_FOLDER}/create_count_objects-${SPEC}.R ${MAINDIR}/.create_count_objects-${SPEC}.R

# Construct shell files
echo "Creating script ${sname}"

cat > ${MAINDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe local 12
#$ -l mem_free=3G,h_vmem=5G
#$ -N ${sname}
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -hold_jid coverage-${EXPERIMENT}.${PREFIX}
echo "**** Job starts ****"
date

rm -rf ${MAINDIR}/Counts/junction/tmpdir

Rscript ${MAINDIR}/.create_count_objects-${SPEC}.R $hgXX $MAINDIR $EXPERIMENT $PREFIX $PE $ERCC

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call

