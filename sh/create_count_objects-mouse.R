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
library('devtools')
library('SummarizedExperiment')
library('plyr')

## Specify parameters
spec <- matrix(c(
	'organism', 'o', 2, 'character', 'mm10',
	'maindir', 'm', 1, 'character', 'Main directory',
	'experiment', 'e', 1, 'character', 'Experiment',
	'prefix', 'p', 1, 'character', 'Prefix',
    'paired', 'l', 1, 'logical', 'Whether the reads are paired-end or not',
    'ercc', 'c', 1, 'logical', 'Whether the reads include ERCC or not',
    'cores', 't', 1, 'integer', 'Number of cores to use',
    'stranded', 's', 1, 'character', "Strandedness of the data: Either 'FALSE', 'forward' or 'reverse'",
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}


## For testing
if(FALSE){
    opt <- list('organism' = 'mm10',
        'maindir' = '/dcl01/lieber/ajaffe/Keri/SoLo_032217',
        'experiment' = 'SoLo_032217',
        'prefix' = 'MiSeq',
        'paired' = FALSE,
		'stranded' = FALSE,
        'ercc' = FALSE,
		'cores' = 1
    )
}

stopifnot(opt$stranded %in% c('FALSE', 'forward', 'reverse'))

RDIR="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/junction_txdb"
EXPNAME = paste0(opt$experiment,"_",opt$prefix)


## read in pheno	
manifest <- read.table(file.path(opt$maindir, 'samples.manifest'), sep = '\t',
    header = FALSE, stringsAsFactors = FALSE)
metrics <- data.frame('SAMPLE_ID' = manifest[, ncol(manifest)],
    stringsAsFactors = FALSE)
N <- length(metrics$SAMPLE_ID)


############################################################ 
###### salmon quantification

sampIDs = as.vector(metrics$SAMPLE_ID)

##observed tpm and number of reads
txTpm = bplapply(sampIDs, function(x) {
	read.table(file.path(opt$maindir, "Salmon_tx", x, "quant.sf"),header = TRUE)$TPM }, 
	BPPARAM = MulticoreParam(opt$cores))
txTpm = do.call(cbind,txTpm)

txNumReads = bplapply(sampIDs, function(x) {
	read.table(file.path(opt$maindir, "Salmon_tx", x, "quant.sf"),header = TRUE)$NumReads }, 
	BPPARAM = MulticoreParam(opt$cores))
txNumReads = do.call(cbind,txNumReads)

colnames(txTpm) = colnames(txNumReads) = sampIDs

##get names of transcripts
txNames = read.table(file.path(opt$maindir, "Salmon_tx", sampIDs[1], "quant.sf"),
						header = TRUE)$Name
txNames = as.character(txNames)
txMap = t(ss(txNames, "\\|",c(1,7,2,6,8)))
txMap = as.data.frame(txMap)
rm(txNames)
colnames(txMap) = c("gencodeTx","txLength","gencodeID","Symbol","gene_type")

rownames(txMap) = rownames(txTpm) = rownames(txNumReads) = txMap$gencodeTx


############################################################ 
###### ercc plots
if (opt$ercc == TRUE ){

	##observed kallisto tpm
	erccTPM = sapply(sampIDs, function(x) {
	  read.table(file.path(opt$maindir, "Ercc", x, "abundance.tsv"),header = TRUE)$tpm
	})
	rownames(erccTPM) = read.table(file.path(opt$maindir, "Ercc", sampIDs[1], "abundance.tsv"),
							header = TRUE)$target_id
	#check finiteness / change NaNs to 0s
	erccTPM[which(is.na(erccTPM),arr.ind=T)] = 0
	
	#expected concentration
	spikeIns = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/ercc_actual_conc.txt",
								as.is=TRUE,row.names=2)
	##match row order
	spikeIns = spikeIns[match(rownames(erccTPM),rownames(spikeIns)),]

	pdf(file.path(opt$maindir, 'Ercc', 'ercc_spikein_check_mix1.pdf'),h=12,w=18)
	mypar(4,6)
	for(i in 1:ncol(erccTPM)) {
		plot(log2(10*spikeIns[,"concentration.in.Mix.1..attomoles.ul."]+1) ~ log2(erccTPM[,i]+1),
			xlab="Kallisto log2(TPM+1)", ylab="Mix 1: log2(10*Concentration+1)",
			main = colnames(erccTPM)[i],
			xlim = c(min(log2(erccTPM+1)),max(log2(erccTPM+1))))
		abline(0, 1, lty=2)
	}
	dev.off()

	mix1conc = matrix(rep(spikeIns[,"concentration.in.Mix.1..attomoles.ul."]), 
						nc = ncol(erccTPM), nr = nrow(erccTPM), byrow=FALSE)
	logErr = (log2(erccTPM+1) - log2(10*mix1conc+1))
	metrics$ERCCsumLogErr = colSums(logErr)

	}
############################################################


### add bam file
metrics$bamFile <- file.path(opt$maindir, 'HISAT2_out', paste0(metrics$SAMPLE_ID, '_accepted_hits.sorted.bam'))

### get alignment metrics
if (opt$paired == TRUE) {
    hisatStats = function(logFile) {
    	y = scan(logFile, what = "character", sep= "\n", 
    		quiet = TRUE, strip=TRUE)
		
    	if (as.numeric(ss(ss(y[2], "\\(",2), "%")) == 100) {
        	## 100% of reads paired
        	reads = as.numeric(ss(y[1], " ")) * 2
        	unaligned = as.numeric(ss(y[12], " "))
        	o = data.frame(trimmed = FALSE,
        		numReads = reads,
        		numMapped = reads - unaligned,
        		numUnmapped = unaligned,
        		overallMapRate = as.numeric(ss(y[15], "\\%"))/100,
		        concordMapRate = (as.numeric(ss(ss(y[4], "\\(",2), "%"))+as.numeric(ss(ss(y[5], "\\(",2), "%")))/100,
                stringsAsFactors = FALSE)
    	} else {
        	## Combo of paired and unpaired (from trimming)
        	reads = as.numeric(ss(y[2], " "))*2 + as.numeric(ss(y[15], " "))
        	unaligned = as.numeric(ss(y[12], " ")) + as.numeric(ss(y[16], " "))
        	o = data.frame(trimmed = TRUE,
        		numReads = reads,
        		numMapped = reads - unaligned,
        		numUnmapped = unaligned,
        		overallMapRate = as.numeric(ss(y[19], "\\%"))/100,
        		concordMapRate = (as.numeric(ss(ss(y[4], "\\(",2), "%"))+as.numeric(ss(ss(y[5], "\\(",2), "%")))/100,
                stringsAsFactors = FALSE)
    	}
    }
} else {
    ## all reads unpaired
    hisatStats = function(logFile) {
    	y = scan(logFile, what = "character", sep= "\n", 
    		quiet = TRUE, strip=TRUE)
    	o = data.frame(numReads = as.numeric(ss(y[1], " ")),
    		numMapped = as.numeric(ss(y[1], " ")) - as.numeric(ss(y[3], " ")),
    		numUnmapped = as.numeric(ss(y[3], " ")),
    		overallMapRate = as.numeric(ss(y[6], "\\%"))/100)
    }
}

logFiles = file.path(opt$maindir, 'HISAT2_out', 'align_summaries', paste0(metrics$SAMPLE_ID, '_summary.txt'))
names(logFiles)  = metrics$SAMPLE_ID
hiStats <- do.call(rbind, lapply(logFiles, hisatStats))

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

gencodeEXONS = as.data.frame(gencodeGTF)[which(gencodeGTF$type=="exon"),c("seqnames","start","end","exon_id")]
names(gencodeEXONS) = c("Chr","Start","End","exon_gencodeID")


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
sym = getBM(attributes = c("ensembl_gene_id","mgi_symbol","entrezgene_id"), 
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

geneMap$Symbol = sym$mgi_symbol[match(geneMap$ensemblID, sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene_id[match(geneMap$ensemblID, sym$ensembl_gene_id)]
	
## counts
geneCountList = mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=opt$cores)
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
bg = matrix(rep(as.numeric(geneStats["Assigned",])), nc = nrow(metrics), 
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

exonMap$Symbol = sym$mgi_symbol[match(exonMap$ensemblID, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$ensemblID, sym$ensembl_gene_id)]

## add gencode exon id
exonMap = join(exonMap, gencodeEXONS, type="left", match="first")

## counts
exonCountList = mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=opt$cores)
exonCounts = do.call("cbind", exonCountList)
rownames(exonCounts) = rownames(exonMap)
exonCounts = exonCounts[,metrics$SAMPLE_ID] # put in order

## remove duplicated
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))

## drop runthrough exons with duplicated exons
i = grepl('-', exonMap$Symbol)
j = countOverlaps(eMap[i], eMap[!i], type = 'equal') > 0
dropIndex = which(i)[j]
exonCounts = exonCounts[-dropIndex,]
exonMap = exonMap[-dropIndex,]
eMap = eMap[-dropIndex,]

## drop duplicated exons
keepIndex = which(!duplicated(eMap))
exonCounts = exonCounts[keepIndex,]
exonMap = exonMap[keepIndex,]

## change rownames
exonMap$exon_libdID = rownames(exonMap)
# rownames(exonMap) = rownames(exonCounts) = exonMap$exon_gencodeID

# number of reads assigned
exonStatList = lapply(paste0(exonFn, ".summary"), 
                      read.delim,row.names=1)
exonStats = do.call("cbind", exonStatList)
colnames(exonStats) = metrics$SAMPLE_ID

## make RPKM
bgE = matrix(rep(colSums(exonCounts)), nc = nrow(metrics), 
	nr = nrow(exonCounts),	byrow=TRUE)
widE = matrix(rep(exonMap$Length), nr = nrow(exonCounts), 
	nc = nrow(metrics),	byrow=FALSE)
exonRpkm = exonCounts/(widE/1000)/(bgE/1e6)

## add mean expression
geneMap$meanExprs = rowMeans(geneRpkm)
exonMap$meanExprs = rowMeans(exonRpkm)

## Create gene,exon RangedSummarizedExperiment objects

# ## recount getRPKM version
# getRPKM <- function(rse, length_var = 'Length', mapped_var = NULL) {
    # mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    # bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    # len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    # wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    # assays(rse)$counts / (wid/1000) / (bg/1e6)
# }

gr_genes <- GRanges(seqnames = geneMap$Chr,
    IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand)
names(gr_genes) <- rownames(geneMap)
mcols(gr_genes) <- DataFrame(geneMap[, - which(colnames(geneMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])

rse_gene <- SummarizedExperiment(assays = list('counts' = geneCounts),
    rowRanges = gr_genes, colData = metrics)
save(rse_gene, file = paste0('rse_gene_', EXPNAME, '_n', N, '.Rdata'))

gr_exons <- GRanges(seqnames = exonMap$Chr,
    IRanges(exonMap$Start, exonMap$End), strand = exonMap$Strand)
names(gr_exons) <- rownames(exonMap)
mcols(gr_exons) <- DataFrame(exonMap[, - which(colnames(exonMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])

rse_exon <- SummarizedExperiment(assays = list('counts' = exonCounts),
    rowRanges = gr_exons, colData = metrics)
save(rse_exon, file = paste0('rse_exon_', EXPNAME, '_n', N, '.Rdata'))


###################
##### junctions

## via primary alignments only
junctionFiles <- file.path(opt$maindir, 'Counts', 'junction', paste0(metrics$SAMPLE_ID, '_junctions_primaryOnly_regtools.count'))
stopifnot(all(file.exists(junctionFiles))) #  TRUE

## annotate junctions
load(file.path(RDIR, "junction_annotation_mm10_gencode_vM11.rda"))

if (opt$stranded %in% c('forward', 'reverse')) {
	juncCounts = junctionCount(junctionFiles, metrics$SAMPLE_ID,
		output = "Count", maxCores=opt$cores,strandSpecific=TRUE)
} else {
	juncCounts = junctionCount(junctionFiles, metrics$SAMPLE_ID,
		output = "Count", maxCores=opt$cores,strandSpecific=FALSE)
}
## filter junction counts - drop jxns in <1% of samples
n = max(1, floor(N/100))
jCountsLogical = DataFrame(sapply(juncCounts$countDF, function(x) x > 0))
jIndex = which(rowSums(as.data.frame(jCountsLogical)) >= n) 
juncCounts = lapply(juncCounts, function(x) x[jIndex,])


############ anno/jMap	
anno = juncCounts$anno
seqlevels(anno, pruning.mode="coarse") = paste0("chr", c(1:19,"X","Y","M"))

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

## extract out jMap
jMap = anno
colnames(mcols(jMap))[which(colnames(mcols(jMap))=="code")] = "Class"
rm(anno)

############ jCounts
jCounts = as.matrix(as.data.frame(juncCounts$countDF))
jCounts = jCounts[names(jMap),gsub("-",".",metrics$SAMPLE_ID)] # ensure lines up
colnames(jCounts) = metrics$SAMPLE_ID  # change from '.' to hyphens if needed

############ jRpkm
bgJ = matrix(rep(colSums(jCounts)), nc = nrow(metrics), 
	nr = nrow(jCounts),	byrow=TRUE)
jRpkm = jCounts/(bgJ/10e6)

rownames(jCounts) = rownames(jRpkm) = names(jMap)
colnames(jRpkm)  = metrics$SAMPLE_ID 
jMap$meanExprs = rowMeans(jRpkm)


# ## sequence of acceptor/donor sites
# left = right = jMap
# end(left) = start(left) +1
# start(right) = end(right) -1
# jMap$leftSeq  = getSeq(Hsapiens, left)
# jMap$rightSeq = getSeq(Hsapiens, right)



### save counts

# tosaveCounts = c("metrics", "geneCounts", "geneMap", "exonCounts", "exonMap", "jCounts", "jMap", 
					# "txNumReads", "txMap" )
# tosaveRpkm = c("metrics", "geneRpkm", "geneMap", "exonRpkm", "exonMap", "jRpkm", "jMap", 
					# "txTpm", "txNumReads", "txMap" )

# if (exists("erccTPM")) {
	# tosaveCounts = c("erccTPM", tosaveCounts)
	# tosaveRpkm = c("erccTPM", tosaveRpkm)
# }

# save(list=ls()[ls() %in% tosaveCounts], compress=TRUE,
	# file= file.path(opt$maindir, paste0('rawCounts_', EXPNAME, '_n', N, '.rda')))
# save(list=ls()[ls() %in% tosaveRpkm], compress=TRUE,
	# file= file.path(opt$maindir, paste0('rpkmCounts_', EXPNAME, '_n', N, '.rda')))


## Create RangedSummarizedExperiment objects	
rse_jx <- SummarizedExperiment(assays = list('counts' = jCounts),
    rowRanges = jMap, colData = metrics)
save(rse_jx, file = paste0('rse_jx_', EXPNAME, '_n', N, '.Rdata'))

## transcript
tx = gencodeGTF[which(gencodeGTF$type == "transcript")]
names(tx) = tx$transcript_id
txMap = tx[rownames(txTpm)]
rse_tx = SummarizedExperiment(
	assays = list('tpm' = txTpm, 'counts' = txNumReads),
    colData = metrics, rowRanges = txMap)
save(rse_tx, file = paste0('rse_tx_', EXPNAME, '_n', N, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
