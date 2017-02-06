#!/bin/bash
#$ -cwd
#$ -N pipeline_ex
#$ -o ./pipeline_ex.o.txt
#$ -e ./pipeline_ex.e.txt
#$ -m e
echo "**** Job starts ****"
date

cat > .gitignore <<EOF
*/.*.R
*/.*.sh
*/.paired_end
*/Coverage
*/logs
*/FastQC
EOF


cd single_end
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "single" --reference "hg19" --cores 1 --fullcov
rm .*.sh .*.R

cd ../paired_end
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "paired" --reference "mm10" --cores 1
rm .paired_end .*.sh .*.R

cd ../merge
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "merge" --reference "hg38" --stranded --cores 2
rm .paired_end .*.sh .*.R

echo "**** Job ends ****"
date
