## Required libraries
library('derfinder')
library('BiocParallel')
library('GenomicRanges')
library('GenomicFeatures')
library('org.Mm.eg.db')
library('biomaRt')
library('BSgenome.Mmusculus.UCSC.mm10')
library('jaffelab')
library('getopt')

## Specify parameters
spec <- matrix(c(
	'organism', 'o', 2, 'character', 'mm10',
	'maindir', 'm', 1, 'character', 'Main directory',
	'experiment', 'e', 1, 'character', 'Experiment',
	'prefix', 'p', 1, 'character', 'Prefix',
    'paired', 'l', 1, 'logical', 'Whether the reads are paired-end or not',
    'ercc', 'c', 1, 'logical', 'Whether the reads include ERCC or not',
    'cores', 't', 1, 'integer', 'Number of cores to use',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

RDIR="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/junction_txdb"
EXPNAME = paste0(opt$experiment,"_",opt$prefix)


## read in pheno	
metrics <- data.frame(read.table(file.path(opt$maindir, 'SAMPLE_IDs.txt'), as.is=TRUE,
    header = FALSE))
names(metrics)[1] <- "SAMPLE_ID"
metrics$SAMPLE_ID <- basename(metrics$SAMPLE_ID)
N <- length(metrics$SAMPLE_ID)

### add bam file
metrics$bamFile <- file.path(opt$maindir, 'HISAT2_out', paste0(metrics$SAMPLE_ID, '_accepted_hits.sorted.bam'))

### get alignment metrics
if (opt$paired == TRUE) {
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

logFiles = file.path(opt$maindir, 'HISAT2_out', 'align_summaries', paste0(metrics$SAMPLE_ID, '_summary.txt'))
names(logFiles)  = metrics$SAMPLE_ID
hiStats = t(sapply(logFiles, hisatStats))

metrics = cbind(metrics,hiStats)	

### confirm total mapping
metrics$totalMapped <- unlist(bplapply(metrics$bamFile, getTotalMapped,
    chrs = paste0("chr", c(1:19, 'X', 'Y')), 
    BPPARAM = MulticoreParam(opt$cores)))
metrics$mitoMapped <- unlist(bplapply(metrics$bamFile, getTotalMapped, chrs = 'chrM', 
    BPPARAM = MulticoreParam(opt$cores)))
metrics$mitoRate <- metrics$mitoMapped / (metrics$mitoMapped +  metrics$totalMapped)

###################################################################

gencodeGTF = import(con="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCm38_mm10/gencode.vM11.annotation.gtf", format="gtf")
gencodeGENES = mcols(gencodeGTF)[which(gencodeGTF$type=="gene"),c("gene_id","type","gene_type")]
rownames(gencodeGENES) = gencodeGENES$gene_id
rm(gencodeGTF)

###############
### gene counts
geneFn <- file.path(opt$maindir, 'Counts', 'gene', paste0(metrics$SAMPLE_ID, '_Gencode.M11.mm10_Genes.counts'))
names(geneFn) = metrics$SAMPLE_ID
stopifnot(all(file.exists(geneFn)))

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]

######### biomart 
# VERSION M11, GRCm38.p4
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
sym = getBM(attributes = c("ensembl_gene_id","mgi_symbol","entrezgene"), 
	values=rownames(geneMap), mart=ensembl)
#########

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
geneMap$gene_type = gencodeGENES[geneMap$gencodeID,"gene_type"]

geneMap$Symbol = sym$mgi_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]
	
## counts
geneCountList = mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,metrics$SAMPLE_ID] # put in order

# number of reads assigned
geneStatList = lapply(paste0(geneFn, ".summary"), 
                      read.delim,row.names=1)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = metrics$SAMPLE_ID
metrics$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))
# rna Rate
metrics$rRNA_rate = colSums(geneCounts[which(geneMap$gene_type == "rRNA"),])/colSums(geneCounts)


# make RPKM
bg = matrix(rep(colSums(geneStats)), nc = nrow(metrics), 
	nr = nrow(geneCounts),	byrow=TRUE)
widG = matrix(rep(geneMap$Length), nr = nrow(geneCounts), 
	nc = nrow(metrics),	byrow=FALSE)
geneRpkm = geneCounts/(widG/1000)/(bg/1e6)

## save metrics
write.csv(metrics, file = file.path(opt$maindir,
    paste0('read_and_alignment_metrics_', opt$experiment, '_', opt$prefix,
    '.csv')))


###############
### exon counts
exonFn <- file.path(opt$maindir, 'Counts', 'exon', paste0(metrics$SAMPLE_ID, '_Gencode.M11.mm10_Exons.counts'))
names(exonFn) = metrics$SAMPLE_ID
stopifnot(all(file.exists(exonFn)))

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$gencodeID = exonMap$Geneid
exonMap$ensemblID = ss(exonMap$Geneid, "\\.")
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$Geneid = NULL
exonMap$gene_type = gencodeGENES[exonMap$gencodeID,"gene_type"]

exonMap$Symbol = sym$hgnc_symbol[match(exonMap$ensemblID, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$ensemblID, sym$ensembl_gene_id)]

## counts
exonCountList = mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
exonCounts = do.call("cbind", exonCountList)
rownames(exonCounts) = rownames(exonMap)
exonCounts = exonCounts[,metrics$SAMPLE_ID] # put in order

## remove duplicated
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
keepIndex= which(!duplicated(eMap))
exonCounts = exonCounts[keepIndex,]
exonMap = exonMap[keepIndex,]

# number of reads assigned
exonStatList = lapply(paste0(exonFn, ".summary"), 
                      read.delim,row.names=1)
exonStats = do.call("cbind", exonStatList)
colnames(exonStats) = metrics$SAMPLE_ID

## make RPKM
bgE = matrix(rep(colSums(exonStats)), nc = nrow(metrics), 
	nr = nrow(exonCounts),	byrow=TRUE)
widE = matrix(rep(exonMap$Length), nr = nrow(exonCounts), 
	nc = nrow(metrics),	byrow=FALSE)
exonRpkm = exonCounts/(widE/1000)/(bgE/1e6)



#############
##### junctions

## via primary alignments only
junctionFiles <- file.path(opt$maindir, 'Counts', 'junction', paste0(metrics$SAMPLE_ID, '_junctions_primaryOnly_regtools.count'))
stopifnot(all(file.exists(junctionFiles))) #  TRUE

juncCounts = junctionCount(junctionFiles, metrics$SAMPLE_ID,
 	output = "Count", maxCores=12,strandSpecific=TRUE)
	
## annotate junctions
load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_mm10_v79_junction_annotation.rda")

anno = juncCounts$anno
anno = anno[seqnames(anno) %in% paste0("chr", c(1:19,"X","Y","M"))]
seqlevels(anno) = paste0("chr", c(1:19,"X","Y","M"))

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
jCounts = jCounts[names(jMap),gsub("-",".",metrics$SAMPLE_ID)]

# MAPPED PER 10 MILLION
mappedPer10M = sapply(jCounts, sum)/10e6
countsM = DataFrame(mapply(function(x,d) x/d, jCounts , mappedPer10M))
rownames(jCounts) = rownames(countsM) = names(jMap)
jRpkm = as.data.frame(countsM)

## sequence of acceptor/donor sites
left = right = anno
end(left) = start(left) +1
start(right) = end(right) -1

jMap$leftSeq  = getSeq(Mmusculus, left)
jMap$rightSeq = getSeq(Mmusculus, right)

### save counts

save(metrics, jMap, jCounts, geneCounts, geneMap, exonCounts, exonMap, compress=TRUE,
	file= file.path(opt$maindir, paste0('rawCounts_', EXPNAME, '_n', N, '.rda')))
save(metrics, jMap, jRpkm, geneRpkm,	geneMap, exonRpkm, exonMap, compress=TRUE,
	file= file.path(opt$maindir, paste0('rpkmCounts_', EXPNAME, '_n', N, '.rda')))

## write out for coverage
write.table(metrics[,c("SAMPLE_ID", "bamFile")], 
	file.path(opt$maindir, 'samples_with_bams.txt'),
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
