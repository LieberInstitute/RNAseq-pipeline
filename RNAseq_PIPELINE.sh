#$ -l mf=6G,h_vmem=8G,h_stack=256M,h_fsize=100G
#$ -pe local 8
#$ -cwd
#$ -o ./sh/o.$TASK_ID.txt
#$ -e ./sh/e.$TASK_ID.txt
#$ -t 1-3
#
##############
FQ_FOLDER=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/test_runthroughAZ/fq
EXP_FOLDER=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/test_runthroughAZ
FILELIST=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/test_runthroughAZ/SAMPLE_IDs.txt
#"hg19" or "hg38"
hgXX="hg38"
#"yes" or "no"
ERCC="yes"
# also change option -t 1-# above based on num of samples
##############

ID=$(awk "NR==$SGE_TASK_ID" $FILELIST )
SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software

if [ $hgXX == "hg38" ] ; then 
GTF=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf
HISATIDX=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary
CHRALIAS=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/chrAliases_GRCh38_gencode2UCSC.csv
elif [ $hgXX == "hg19" ] ; then 
GTF=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf
HISATIDX=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh37_hg19/hisat2_GRCh37primary
CHRALIAS=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/chrAliases_GRCh37_gencode2UCSC.csv
else 
	echo "Enter hg19 or hg38 for hgXX variable." >&2
	exit 1
fi

############
###ercc spikein tpm quantification
if [ $ERCC == "yes" ] ; then
	mkdir -p ${EXP_FOLDER}/Ercc/${ID}
	${SOFTWARE}/kallisto quant \
	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx \
	-o ${EXP_FOLDER}/Ercc/${ID} -t 8 ${FQ_FOLDER}/${ID}_R1.fastq.gz ${FQ_FOLDER}/${ID}_R2.fastq.gz

	#R CMD BATCH ./sh/plot_ercc_spikein_results.R
fi


############
###fastQC on untrimmed
mkdir -p ${EXP_FOLDER}/FastQC/Untrimmed/${ID}
${SOFTWARE}/FastQC_v0.11.5/fastqc \
${FQ_FOLDER}/${ID}_R1.fastq.gz ${FQ_FOLDER}/${ID}_R2.fastq.gz \
--outdir=${EXP_FOLDER}/FastQC/Untrimmed/${ID} --extract


############
###check adapter content
###trim adapters if needed
###run HISAT2 on trimmed or untrimmed
REPORT1=${EXP_FOLDER}/FastQC/Untrimmed/${ID}/${ID}_R1_fastqc/summary.txt
REPORT2=${EXP_FOLDER}/FastQC/Untrimmed/${ID}/${ID}_R2_fastqc/summary.txt
RESULT1=$(grep "Adapter Content" $REPORT1 | cut -c1-4)
RESULT2=$(grep "Adapter Content" $REPORT2 | cut -c1-4)

if [[ $RESULT1 == "FAIL" || $RESULT2 == "FAIL" ]] ; then
	## trim, rerun fastQC, run hisat on trimmed
	
	mkdir -p ${EXP_FOLDER}/fastq/Trimmed
	FP=${EXP_FOLDER}/fastq/Trimmed/${ID}_trimmed_forward_paired.fq.gz
	FU=${EXP_FOLDER}/fastq/Trimmed/${ID}_trimmed_forward_unpaired.fq.g
	RP=${EXP_FOLDER}/fastq/Trimmed/${ID}_trimmed_reverse_paired.fq.gz
	RU=${EXP_FOLDER}/fastq/Trimmed/${ID}_trimmed_reverse_unpaired.fq.gz
	
	## trim adapters
	java -jar ${SOFTWARE}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 \
	${FQ_FOLDER}/${ID}_R1.fastq.gz ${FQ_FOLDER}/${ID}_R2.fastq.gz $FP $FU $RP $RU \
	ILLUMINACLIP:${SOFTWARE}/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10:1 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

	## rerun fastqc
	mkdir -p ${EXP_FOLDER}/FastQC/Trimmed/${ID}
	${SOFTWARE}/FastQC_v0.11.5/fastqc \
	${FQ_FOLDER}/${ID}_R1.fastq.gz ${FQ_FOLDER}/${ID}_R2.fastq.gz \
	--outdir=${EXP_FOLDER}/FastQC/Trimmed/${ID} --extract

	## run hisat on trimmed
	mkdir -p ${EXP_FOLDER}/HISAT2_out/align_summaries
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p 8 \
	-x $HISATIDX -1 $FP -2 $RP -U ${FU},${RU} \
	-S ${EXP_FOLDER}/HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF \
	2>${EXP_FOLDER}/HISAT2_out/align_summaries/${ID}_summary.txt
	
else 
	############
	###run hisat on untrimmed
	mkdir -p ${EXP_FOLDER}/HISAT2_out/align_summaries
	${SOFTWARE}/hisat2-2.0.4/hisat2 -p 8 \
	-x $HISATIDX -1 ${FQ_FOLDER}/${ID}_R1.fastq.gz -2 ${FQ_FOLDER}/${ID}_R2.fastq.gz \
	-S ${EXP_FOLDER}/HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF \
	2>${EXP_FOLDER}/HISAT2_out/align_summaries/${ID}_summary.txt
fi


############
###sam to bam
SAM=${EXP_FOLDER}/HISAT2_out/${ID}_hisat_out.sam
BAMACC=${EXP_FOLDER}/HISAT2_out/${ID}_accepted_hits.bam
BAMS=${EXP_FOLDER}/HISAT2_out/${ID}_accepted_hits_sorted

#filter unmapped segments
samtools view -bh -F 4 $SAM > $BAMACC
samtools sort -@ 8 $BAMACC $BAMS
samtools index ${BAMS}.bam

#rm $SAM $BAMACC


############
###featureCounts
mkdir -p ${EXP_FOLDER}/Counts/gene
mkdir -p ${EXP_FOLDER}/Counts/exon
mkdir -p ${EXP_FOLDER}/Counts/junction/tmpdir

BAM=${BAMS}.bam

# genes	
featureCounts -s 2 -A $CHRALIAS -a $GTF \
-o ${EXP_FOLDER}/Counts/gene/${ID}_Gencode.v25.${hgXX}_Genes.counts $BAM
# exons	
featureCounts -s 2 -T 8 -O -f -A $CHRALIAS -a $GTF \
-o ${EXP_FOLDER}/Counts/exon/${ID}_Gencode.v25.${hgXX}_Exons.counts $BAM

# junctions	
OUTJXN=${EXP_FOLDER}/Counts/junction/${ID}_junctions_primaryOnly_regtools.bed
OUTCOUNT=${EXP_FOLDER}/Counts/junction/${ID}_junctions_primaryOnly_regtools.count
TMPDIR=${EXP_FOLDER}/Counts/junction/tmpdir
TMPBAM=$TMPDIR/${ID}.bam
#filter only primary alignments
samtools view -bh -F 0x100 $BAM > $TMPBAM
samtools index $TMPBAM
regtools junctions extract -i 9 -o $OUTJXN $TMPBAM
${SOFTWARE}/bed_to_juncs_withCount < $OUTJXN > $OUTCOUNT

#rm -rf $TMPDIR


############
###coverage files from bam
mkdir -p ${EXP_FOLDER}/Coverage

BG=${EXP_FOLDER}/Coverage/${ID}.bedGraph
BGS=${EXP_FOLDER}/Coverage/${ID}.sorted.bedGraph
BW=${EXP_FOLDER}/Coverage/${ID}.bw
CHRSIZE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/${hgXX}.chrom.sizes.gencode

${SOFTWARE}/bedtools-2.26.0/bin/bedtools genomecov -ibam $BAM -bga -split > $BG
${SOFTWARE}/bedtools-2.26.0/bin/sortBed -i $BG > $BGS
${SOFTWARE}/bedGraphToBigWig $BGS $CHRSIZE $BW

#rm $BG $BGS




######## TO DO:

#star alignment?

#leo? mean coverage bigwig

###automated R script that creates count objects
#R CMD BATCH create_counts_objects.R




