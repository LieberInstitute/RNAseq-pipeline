##
args = commandArgs(TRUE)
hgXX = args[1]
MAINDIR = args[2]
EXPERIMENT = args[3]
PREFIX = args[4]
PE = args[5]
ERCC = args[6]

source("/users/ajaffe/Lieber/lieber_functions_aj.R")
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")
RDIR="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/junction_txdb"
EXPNAME = paste0(EXPERIMENT,"_",PREFIX)

## change ensembl chr names from default function
getTotalMapped = function(bamFile, mc.cores=1, returnM = TRUE) {
	thecall = paste("samtools idxstats",bamFile)
	tmp = parallel::mclapply(thecall, function(x) {
		cat(".")
		xx = system(x,intern=TRUE)
		xx = do.call("rbind", strsplit(xx, "\t"))
		d = data.frame(chr=xx[,1], L=xx[,2], mapped = xx[,3],
			stringsAsFactors=FALSE)
		d
	},mc.cores=mc.cores)
	
	out = list(totalMapped = sapply(tmp, function(x) sum(as.numeric(x$mapped[x$chr %in% c(1:22,"X","Y")]))),
		mitoMapped = sapply(tmp, function(x) as.numeric(x$mapped[x$chr=="MT"])))
	return(out)
}

## read in pheno	
pd = data.frame(read.table(paste0(MAINDIR,"/SAMPLE_IDs.txt"), as.is=TRUE, header=F))
names(pd)[1] = "SAMPLE_ID"
N = length(pd$SAMPLE_ID)

############################################################ 
###### ercc plots
if (ERCC == TRUE ){
	library(Biostrings)
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
	require(stringr)
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
	require(stringr)
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
libSize = getTotalMapped(pd$bamFile,mc.cores=12)
pd$totalMapped = libSize$totalMapped
pd$mitoMapped = libSize$mitoMapped
pd$mitoRate = pd$mitoMapped / (pd$mitoMapped +  libSize$totalMapped)

###################################################################
library(GenomicRanges)
library(GenomicFeatures)
#library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(biomaRt)


###############
### gene counts
geneFn = paste0(MAINDIR,"/Counts/gene/",pd$SAMPLE_ID,"_Ensembl.rnor6.0.rn6_Genes.counts")
names(geneFn) = pd$SAMPLE_ID
all(file.exists(geneFn))

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]

######### biomart 
# VERSION Rnor_6.0
ensembl = useMart("ensembl")
ensembl = useDataset("rnorvegicus_gene_ensembl",mart=ensembl)
sym = getBM(attributes = c("ensembl_gene_id","rgd_symbol","entrezgene"), 
	values=rownames(geneMap), mart=ensembl)
#########

## organize gene map
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Chr = ifelse(geneMap$Chr=='MT','chrM',paste0('chr',geneMap$Chr))
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
rownames(geneMap) = geneMap$Geneid

geneMap$Symbol = sym$rgd_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]
	
## counts
geneCountList = mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
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
exonFn = paste0(MAINDIR,"/Counts/exon/",pd$SAMPLE_ID,"_Ensembl.rnor6.0.rn6_Exons.counts")
names(exonFn) = pd$SAMPLE_ID
all(file.exists(exonFn))

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$Chr = ifelse(exonMap$Chr=='MT','chrM',paste0('chr',exonMap$Chr))
rownames(exonMap) = paste0("e", rownames(exonMap))

exonMap$Symbol = sym$rgd_symbol[match(exonMap$Geneid, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$Geneid, sym$ensembl_gene_id)]

## counts
exonCountList = mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
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

## via primary alignments only
junctionFiles = paste0(MAINDIR,"/Counts/junction/",pd$SAMPLE_ID,"_junctions_primaryOnly_regtools.count")
all(file.exists(junctionFiles)) #  TRUE

juncCounts = junctionCount(junctionFiles, pd$SAMPLE_ID,
 	output = "Count", maxCores=12,strandSpecific=TRUE)
	
## annotate junctions
load(paste0(RDIR,"/junction_annotation_rn6_ensembl_v86.rda"))

anno = juncCounts$anno
seqlevels(anno,force=TRUE) = c(1:20,"X","Y","MT")
seqlevels(anno) = paste0("chr", c(1:20,"X","Y","M"))

## add additional annotation
anno$inEnsembl = countOverlaps(anno, theJunctions, type="equal") > 0
anno$inEnsemblStart = countOverlaps(anno, theJunctions, type="start") > 0
anno$inEnsemblEnd = countOverlaps(anno, theJunctions, type="end") > 0

oo = findOverlaps(anno, theJunctions, type="equal")
anno$ensemblGeneID = NA
anno$ensemblGeneID[queryHits(oo)] = as.character(theJunctions$ensemblID[subjectHits(oo)])
anno$ensemblSymbol = NA
anno$ensemblSymbol[queryHits(oo)] = theJunctions$symbol[subjectHits(oo)]
anno$ensemblStrand = NA
anno$ensemblStrand[queryHits(oo)] = as.character(strand(theJunctions)[subjectHits(oo)])
anno$ensemblTx = CharacterList(vector("list", length(anno)))
anno$ensemblTx[queryHits(oo)] = theJunctions$tx[subjectHits(oo)]
anno$numTx = elementNROWS(anno$ensemblTx)

# clean up
anno$ensemblSymbol = geneMap$Symbol[match(anno$ensemblGeneID, rownames(geneMap))]

## junction code
anno$code = ifelse(anno$inEnsembl, "InEns", 
	ifelse(anno$inEnsemblStart & anno$inEnsemblEnd, "ExonSkip",
	ifelse(anno$inEnsemblStart | anno$inEnsemblEnd, "AltStartEnd", "Novel")))

## b/w exons and junctions
exonGR = GRanges( exonMap$Chr,	IRanges(exonMap$Start, exonMap$End))
anno$startExon = match(paste0(seqnames(anno),":",start(anno)-1), 
	paste0(seqnames(exonGR), ":", end(exonGR)))
anno$endExon = match(paste0(seqnames(anno),":",end(anno)+1),
	paste0(seqnames(exonGR), ":", start(exonGR)))
g = data.frame(leftGene = exonMap$Geneid[anno$startExon],
	rightGene = exonMap$Geneid[anno$endExon],
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

## extract out
jMap = anno
jCounts = juncCounts$countDF
jCounts = jCounts[names(jMap),paste0('X',gsub("-",".",pd$SAMPLE_ID))]

## rename chr to gencode
names(jMap) = paste0(as.character(seqnames(jMap)),":",
                     start(jMap),"-",end(jMap),"(",as.character(strand(jMap)),")")
rownames(jCounts) = names(jMap)

# MAPPED PER 10 MILLION
mappedPer10M = sapply(jCounts, sum)/10e6
countsM = DataFrame(mapply(function(x,d) x/d, jCounts , mappedPer10M))
rownames(jCounts) = rownames(countsM) = names(jMap)
jRpkm = as.data.frame(countsM)

## sequence of acceptor/donor sites
library(BSgenome.Rnorvegicus.UCSC.rn6)
left = right = anno
end(left) = start(left) +1
start(right) = end(right) -1

jMap$leftSeq  = getSeq(Rnorvegicus, left)
jMap$rightSeq = getSeq(Rnorvegicus, right)

### save counts

save(pd, jMap, jCounts, geneCounts, geneMap, exonCounts, exonMap, compress=TRUE,
	file=paste0(MAINDIR,"/rawCounts_",EXPNAME,"_n",N,".rda"))
save(pd, jMap, jRpkm, geneRpkm,	geneMap, exonRpkm, exonMap, compress=TRUE,
	file=paste0(MAINDIR,"/rpkmCounts_",EXPNAME,"_n",N,".rda"))

## write out for coverage
write.table(pd[,c("SAMPLE_ID", "bamFile")], 
	paste0(MAINDIR,"/samples_with_bams.txt"),
	row.names=FALSE, quote=FALSE, sep="\t")

