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

## Clean examples
bash clean_examples.sh

cd single_end
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "single" --reference "hg19" --cores 1 --fullcov

cd ../paired_end
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "paired" --reference "mm10" --cores 1

cd ../merge
bash ../../sh/rnaseq-run-all.sh --experiment "example" --prefix "merge" --reference "hg38" --stranded --cores 2

echo "**** Job ends ****"
date
