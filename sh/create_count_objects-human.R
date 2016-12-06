## Required libraries
library('derfinder')
library('BiocParallel')
library('Biostrings')
library('GenomicRanges')
library('GenomicFeatures')
library('org.Hs.eg.db')
library('biomaRt')
library('stringr')
library('jaffelab')

##
args = commandArgs(TRUE)
hgXX = args[1]
MAINDIR = args[2]
EXPERIMENT = args[3]
PREFIX = args[4]
PE = args[5]
ERCC = args[6]

if (hgXX == "hg19") { 
	library('BSgenome.Hsapiens.UCSC.hg19')
} else if (hgXX == "hg38") { 
	library('BSgenome.Hsapiens.UCSC.hg38')
}

#hgXX="hg19"
#MAINDIR="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/tests/lake"
#EXPERIMENT="lake"
#PREFIX="ten"
#PE="FALSE"
#ERCC="FALSE"

source("/users/ajaffe/Lieber/lieber_functions_aj.R")
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")
RDIR="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/junction_txdb"
EXPNAME = paste0(EXPERIMENT,"_",PREFIX)


## read in pheno	
pd = data.frame(read.table(paste0(MAINDIR,"/SAMPLE_IDs.txt"), as.is=TRUE, header=F))
names(pd)[1] = "SAMPLE_ID"
N = length(pd$SAMPLE_ID)


############################################################ 
###### ercc plots
if (ERCC == TRUE ){
	sampIDs = as.vector(pd$SAMPLE_ID)

	##observed kallisto tpm
	erccTPM = sapply(sampIDs, function(x) {
	  read.table(paste0(MAINDIR,"/Ercc/",x,"/abundance.tsv"),header=T)$tpm
	})
	rownames(erccTPM) = read.table(paste0(MAINDIR,"/Ercc/",sampIDs[1],"/abundance.tsv"),
							header=T)$target_id
	#expected concentration
	spikeIns = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/ercc_actual_conc.txt",
								as.is=TRUE,row.names=2)
	##match row order
	spikeIns = spikeIns[match(rownames(erccTPM),rownames(spikeIns)),]

	pdf(paste0(MAINDIR,"/Ercc/ercc_spikein_check_mix1.pdf"),h=12,w=18)
	mypar(4,6)
	for(i in 1:ncol(erccTPM)) {
		plot(log2(spikeIns[,"concentration.in.Mix.1..attomoles.ul."]+1) ~ log2(erccTPM[,i]+1),
			xlab="Kallisto log2(TPM+1)", ylab="Mix 1: log2(Concentration+1)",
			main = colnames(erccTPM)[i],
			xlim = c(min(log2(erccTPM+1)),max(log2(erccTPM+1))))
		abline(0, 1, lty=2)
	}
	dev.off()

	mix1conc = matrix(rep(spikeIns[,"concentration.in.Mix.1..attomoles.ul."]), 
						nc = ncol(erccTPM), nr = nrow(erccTPM), byrow=FALSE)
	logErr = (log2(erccTPM+1) - log2(mix1conc+1))
	pd$ERCCsumLogErr = colSums(logErr)
		
	}
############################################################


### add bam file
pd$bamFile = paste0(MAINDIR, "/HISAT2_out/", pd$SAMPLE_ID, "_accepted_hits.sorted.bam")

### get alignment metrics
if (PE == TRUE) {
hisatStats = function(logFile) {
	y = scan(logFile, what = "character", sep= "\n", 
		quiet = TRUE, strip=TRUE)
		
	if (as.numeric(ss(ss(y[2], "\\(",2), "%")) == 100) {
	## 100% of reads paired
	reads = as.numeric(ss(y[1], " "))*2
	unaligned = as.numeric(ss(y[12], " "))
	o = c(trimmed="FALSE",
		numReads = reads,
		numMapped = reads - unaligned,
		numUnmapped = unaligned,
		overallMapRate = as.numeric(ss(y[15], "\\%"))/100,
		concordMapRate = (as.numeric(ss(ss(y[4], "\\(",2), "%"))+as.numeric(ss(ss(y[5], "\\(",2), "%")))/100)
	} else {
	## Combo of paired and unpaired (from trimming)
	reads = as.numeric(ss(y[2], " "))*2 + as.numeric(ss(y[15], " "))
	unaligned = as.numeric(ss(y[12], " ")) + as.numeric(ss(y[16], " "))
	o = c(trimmed="TRUE",
		numReads = reads,
		numMapped = reads - unaligned,
		numUnmapped = unaligned,
		overallMapRate = as.numeric(ss(y[19], "\\%"))/100,
		concordMapRate = (as.numeric(ss(ss(y[4], "\\(",2), "%"))+as.numeric(ss(ss(y[5], "\\(",2), "%")))/100)	
	}
}
} else {
## all reads unpaired
hisatStats = function(logFile) {
	y = scan(logFile, what = "character", sep= "\n", 
		quiet = TRUE, strip=TRUE)
	o = c(numReads = as.numeric(ss(y[1], " ")),
		numMapped = as.numeric(ss(y[1], " ")) - as.numeric(ss(y[3], " ")),
		numUnmapped = as.numeric(ss(y[3], " ")),
		overallMapRate = as.numeric(ss(y[6], "\\%"))/100)
}
}

logFiles = paste0(MAINDIR, "/HISAT2_out/align_summaries/", pd$SAMPLE_ID, "_summary.txt")
names(logFiles)  = pd$SAMPLE_ID
hiStats = t(sapply(logFiles, hisatStats))

pd = cbind(pd,hiStats)	

### confirm total mapping
pd$totalMapped <- unlist(bplapply(pd$bamFile, getTotalMapped,
    chrs = paste0("chr", c(1:22, 'X', 'Y')), 
    BPPARAM = MulticoreParam(8)))
pd$mitoMapped <- unlist(bplapply(pd$bamFile, getTotalMapped, chrs = 'chrM', 
    BPPARAM = MulticoreParam(8)))
pd$mitoRate <- pd$mitoMapped / (pd$mitoMapped +  pd$totalMapped)


###################################################################

if (hgXX == "hg19") {
	filename = "_Gencode.v25lift37.hg19"
} else if (hgXX == "hg38") {
	filename = "_Gencode.v25.hg38"
}

###############
### gene counts
geneFn = paste0(MAINDIR,"/Counts/gene/",pd$SAMPLE_ID,filename,"_Genes.counts")
names(geneFn) = pd$SAMPLE_ID
all(file.exists(geneFn))

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]

## organize gene map
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
rownames(geneMap) = geneMap$Geneid
geneMap$gencodeID = geneMap$Geneid
geneMap$ensemblID = ss(geneMap$Geneid, "\\.")
geneMap$Geneid = NULL

######### biomart 
if (hgXX=="hg19") {
	# VERSION 75, GRCh37.p13
	ensembl = useMart("ENSEMBL_MART_ENSEMBL", 
		dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
	sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
		values=geneMap$ensemblID, mart=ensembl)
} else if (hgXX=="hg38") {
	# VERSION 85, GRCh38.p7
	ensembl = useMart("ENSEMBL_MART_ENSEMBL",  
		dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
	sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
			values=geneMap$ensemblID, mart=ensembl)
}
#########

geneMap$Symbol = sym$hgnc_symbol[match(geneMap$ensemblID, sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(geneMap$ensemblID, sym$ensembl_gene_id)]

## counts
geneCountList = mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=8)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,pd$SAMPLE_ID] # put in order

# number of reads assigned
geneStatList = lapply(paste0(geneFn, ".summary"), 
	read.delim,row.names=1)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = pd$SAMPLE_ID
pd$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))

# make RPKM
bg = matrix(rep(colSums(geneStats)), nc = nrow(pd), 
	nr = nrow(geneCounts),	byrow=TRUE)
widG = matrix(rep(geneMap$Length), nr = nrow(geneCounts), 
	nc = nrow(pd),	byrow=FALSE)
geneRpkm = geneCounts/(widG/1000)/(bg/1e6)

## save pd
write.csv(pd, file=paste0(MAINDIR,"/annotated_pd.csv"))


###############
### exon counts
exonFn = paste0(MAINDIR,"/Counts/exon/",pd$SAMPLE_ID,filename,"_Exons.counts")
names(exonFn) = pd$SAMPLE_ID
all(file.exists(exonFn))

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$gencodeID = exonMap$Geneid
exonMap$ensemblID = ss(exonMap$Geneid, "\\.")
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$Geneid = NULL

exonMap$Symbol = sym$hgnc_symbol[match(exonMap$ensemblID, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$ensemblID, sym$ensembl_gene_id)]

## counts
exonCountList = mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=8)
exonCounts = do.call("cbind", exonCountList)
rownames(exonCounts) = rownames(exonMap)
exonCounts = exonCounts[,pd$SAMPLE_ID] # put in order

## remove duplicated
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
keepIndex= which(!duplicated(eMap))
exonCounts = exonCounts[keepIndex,]
exonMap = exonMap[keepIndex,]

# number of reads assigned
exonStatList = lapply(paste0(exonFn, ".summary"), 
                      read.delim,row.names=1)
exonStats = do.call("cbind", exonStatList)
colnames(exonStats) = pd$SAMPLE_ID

## make RPKM
bgE = matrix(rep(colSums(exonStats)), nc = nrow(pd), 
	nr = nrow(exonCounts),	byrow=TRUE)
widE = matrix(rep(exonMap$Length), nr = nrow(exonCounts), 
	nc = nrow(pd),	byrow=FALSE)
exonRpkm = exonCounts/(widE/1000)/(bgE/1e6)


#############
##### junctions

## import theJunctions annotation
if (hgXX == "hg19") { 
	#load(paste0(RDIR,"/junction_annotation_hg19_ensembl_v75.rda"))
	load(paste0(RDIR,"/junction_annotation_hg19_gencode_v25lift37.rda"))
} else if (hgXX == "hg38") { 
	#load(paste0(RDIR,"/junction_annotation_hg38_ensembl_v85.rda"))
	load(paste0(RDIR,"/junction_annotation_hg38_gencode_v25.rda"))
}

## via primary alignments only
junctionFiles = paste0(MAINDIR,"/Counts/junction/",pd$SAMPLE_ID,"_junctions_primaryOnly_regtools.count")
all(file.exists(junctionFiles)) #  TRUE

if (PE == TRUE) { 
	juncCounts = junctionCount(junctionFiles, pd$SAMPLE_ID,
		output = "Count", maxCores=8,strandSpecific=TRUE)
} else { 
	source("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq_functions.R")
	juncCounts = junctionCount(junctionFiles, pd$SAMPLE_ID,
		output = "Count", maxCores=8,strandSpecific=FALSE)
}
anno = juncCounts$anno
seqlevels(anno, force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))
#anno = anno[seqnames(anno) %in% paste0("chr", c(1:22,"X","Y","M"))]
#seqlevels(anno) = paste0("chr", c(1:22,"X","Y","M"))

## add additional annotation
anno$inGencode = countOverlaps(anno, theJunctions, type="equal") > 0
anno$inGencodeStart = countOverlaps(anno, theJunctions, type="start") > 0
anno$inGencodeEnd = countOverlaps(anno, theJunctions, type="end") > 0

oo = findOverlaps(anno, theJunctions, type="equal")
anno$gencodeGeneID = NA
anno$gencodeGeneID[queryHits(oo)] = as.character(theJunctions$gencodeID[subjectHits(oo)])
anno$ensemblID = ss(anno$gencodeGeneID, "\\.")
anno$Symbol = NA
anno$Symbol[queryHits(oo)] = theJunctions$symbol[subjectHits(oo)]
anno$gencodeStrand = NA
anno$gencodeStrand[queryHits(oo)] = as.character(strand(theJunctions)[subjectHits(oo)])
anno$gencodeTx = CharacterList(vector("list", length(anno)))
anno$gencodeTx[queryHits(oo)] = theJunctions$tx[subjectHits(oo)]
anno$numTx = elementNROWS(anno$gencodeTx)

# clean up
#anno$Symbol = geneMap$Symbol[match(anno$gencodeGeneID, rownames(geneMap))]

## junction code
anno$code = ifelse(anno$inGencode, "InGen", 
	ifelse(anno$inGencodeStart & anno$inGencodeEnd, "ExonSkip",
	ifelse(anno$inGencodeStart | anno$inGencodeEnd, "AltStartEnd", "Novel")))

## b/w exons and junctions
exonGR = GRanges( exonMap$Chr,	IRanges(exonMap$Start, exonMap$End))
anno$startExon = match(paste0(seqnames(anno),":",start(anno)-1), 
	paste0(seqnames(exonGR), ":", end(exonGR)))
anno$endExon = match(paste0(seqnames(anno),":",end(anno)+1),
	paste0(seqnames(exonGR), ":", start(exonGR)))
g = data.frame(leftGene = exonMap$gencodeID[anno$startExon],
	rightGene = exonMap$gencodeID[anno$endExon],
	leftGeneSym = exonMap$Symbol[anno$startExon],
	rightGeneSym = exonMap$Symbol[anno$endExon],
	stringsAsFactors=FALSE)
g$newGene = NA
g$newGeneSym = NA
g$newGene[which(g$leftGene==g$rightGene)] = 
	g$leftGene[which(g$leftGene==g$rightGene)] 
g$newGeneSym[which(g$leftGene==g$rightGene)] = 
	g$leftGeneSym[which(g$leftGene==g$rightGene)] 
g$newGene[which(g$leftGene!=g$rightGene)] = 
	paste0(g$leftGene,"-",g$rightGene)[which(g$leftGene!=g$rightGene)] 
g$newGeneSym[which(g$leftGene!=g$rightGene)] = 
	paste0(g$leftGeneSym,"-",g$rightGeneSym)[which(g$leftGene!=g$rightGene)] 
g$newGene[which(is.na(g$newGene) & is.na(g$leftGene))] = 
	g$rightGene[which(is.na(g$newGene) & is.na(g$leftGene))] 
g$newGene[which(is.na(g$newGene) & is.na(g$rightGene))] = 
	g$leftGene[which(is.na(g$newGene) & is.na(g$rightGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] = 
	g$rightGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] = 
	g$leftGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] 
g$newGeneSym[g$newGeneSym==""] = NA
g$newGeneSym[g$newGeneSym=="-"] = NA
anno$newGeneID = g$newGene
anno$newGeneSymbol = g$newGeneSym
anno$isFusion = grepl("-", anno$newGeneID)

anno$newGeneSymbol[anno$code =="InGen"] = anno$Symbol[anno$code =="InGen"]
anno$newGeneID[anno$code =="InGen"] = anno$gencodeGeneID[anno$code =="InGen"]
## extract out
jMap = anno
jCounts = juncCounts$countDF
jCounts = jCounts[names(jMap),gsub("-",".",pd$SAMPLE_ID)]

mappedPer10M = sapply(jCounts, sum)/10e6
countsM = DataFrame(mapply(function(x,d) x/d, jCounts , mappedPer10M))
rownames(jCounts) = rownames(countsM) = names(jMap)
jRpkm = as.data.frame(countsM)
rownames(jRpkm) = names(jMap)
colnames(jRpkm)  = colnames(geneRpkm)

## sequence of acceptor/donor sites
left = right = jMap
end(left) = start(left) +1
start(right) = end(right) -1

jMap$leftSeq  = getSeq(Hsapiens, left)
jMap$rightSeq = getSeq(Hsapiens, right)

############################
### add transcript maps ####
if (hgXX == "hg19") { 
	load(paste0(RDIR,"/feature_to_Tx_hg19_gencode_v25lift37.rda"))
} else if (hgXX == "hg38") { 
	load(paste0(RDIR,"/feature_to_Tx_hg38_gencode_v25.rda")) 
}

## gene annotation
geneMap$Class = "InGen"
geneMap$meanExprs = rowMeans(geneRpkm)
mmTx = match(geneMap$gencodeID, names(allTx))
tx = CharacterList(vector("list", nrow(geneMap)))
tx[!is.na(mmTx)] = allTx[mmTx[!is.na(mmTx)]]
geneMap$NumTx = elementNROWS(tx)
geneMap$gencodeTx = sapply(tx,paste0,collapse=";")

## exon annotation
exonMap$Class = "InGen"
exonMap$meanExprs = rowMeans(exonRpkm)
mmTx = match(rownames(exonMap), names(allTx))
tx = CharacterList(vector("list", nrow(exonMap)))
tx[!is.na(mmTx)] = allTx[mmTx[!is.na(mmTx)]]
exonMap$NumTx = elementNROWS(tx)
exonMap$gencodeTx = sapply(tx,paste0,collapse=";")

## junctions
jMap$meanExprs= rowMeans(jRpkm)
#jMap$gencodeGeneID=NULL
#jMap$Symbol=NULL
#colnames(mcols(jMap))[c(8,11,12)] = c("Class", "gencodeID", "Symbol")
colnames(mcols(jMap))[10] = "Class"

### save counts

save(pd, jMap, jCounts, geneCounts, geneMap, exonCounts, exonMap, compress=TRUE,
	file=paste0(MAINDIR,"/rawCounts_",EXPNAME,"_n",N,".rda"))
save(pd, jMap, jRpkm, geneRpkm,	geneMap, exonRpkm, exonMap, compress=TRUE,
	file=paste0(MAINDIR,"/rpkmCounts_",EXPNAME,"_n",N,".rda"))

## write out for coverage
write.table(pd[,c("SAMPLE_ID", "bamFile")], 
	paste0(MAINDIR,"/samples_with_bams.txt"),
	row.names=FALSE, quote=FALSE, sep="\t")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
gotDevtools <- requireNamespace('devtools', quietly = TRUE)
if(gotDevtools) {
    devtools::session_info()
} else {
    sessionInfo()
}

