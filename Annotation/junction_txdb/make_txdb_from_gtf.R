#######

library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(rtracklayer)
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")

################################################################
####  hg19  ####################################################
################################################################

## chromosome info
chrInfo = getChromInfoFromUCSC("hg19")
chrInfo$chrom = as.character(chrInfo$chrom)
chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg19"))

##########################
##### ENSEMBL ############

## read in GTF as GRanges
ensembl_v75 = import(con = "../ensembl/Homo_sapiens.GRCh37.75.gtf", format = "gtf")
seqlevels(ensembl_v75,force=TRUE) = c(1:22,"X","Y","MT")
seqlevels(ensembl_v75) = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(ensembl_v75) = si

# get map
ensMap = mcols(ensembl_v75)
ensMap = ensMap[!duplicated(ensMap$transcript_id),
	c("gene_id","gene_name","transcript_id","gene_biotype")]

##  convert to txdb
ensembl_v75_txdb = makeTxDbFromGRanges(ensembl_v75)
#saveDb(ensembl_v75_txdb, file="ensembl_v75.sqlite")

# get introns
introns = intronsByTranscript(ensembl_v75_txdb,use.names=TRUE)
introns = unlist(introns)
introns$TranscriptID = names(introns)
introns$ensemblID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

### make unique junctions
ensembl_v75_junctions = unique(introns)
oo = findOverlaps(ensembl_v75_junctions, introns,type="equal")
xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
ensembl_v75_junctions$tx = CharacterList(xx)
ensembl_v75_junctions$TranscriptID = NULL
ensembl_v75_junctions$numTx = elementLengths(ensembl_v75_junctions$tx)
names(ensembl_v75_junctions) = NULL
theJunctions = ensembl_v75_junctions
save(theJunctions, file="junction_annotation_hg19_ensembl_v75.rda")



##########################
##### GENCODE ############

## read in GTF as GRanges
gencode_v25lift37 = import(con = "../GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf", format = "gtf")
seqlevels(gencode_v25lift37,force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(gencode_v25lift37) = si

# get map
ensMap = mcols(gencode_v25lift37)
ensMap = ensMap[!duplicated(ensMap$transcript_id),
	c("gene_id","gene_name","transcript_id","gene_type")]

##  convert to txdb
gencode_v25lift37_txdb = makeTxDbFromGRanges(gencode_v25lift37)

# get introns
introns = intronsByTranscript(gencode_v25lift37_txdb,use.names=TRUE)
introns = unlist(introns)
introns$TranscriptID = names(introns)
introns$gencodeID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
introns$ensemblID = ss(introns$gencodeID, "\\.")
introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

### make unique junctions
gencode_v25lift37_junctions = unique(introns)
oo = findOverlaps(gencode_v25lift37_junctions, introns,type="equal")
xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
gencode_v25lift37_junctions$tx = CharacterList(xx)
gencode_v25lift37_junctions$TranscriptID = NULL
gencode_v25lift37_junctions$numTx = elementNROWS(gencode_v25lift37_junctions$tx)
names(gencode_v25lift37_junctions) = NULL
theJunctions = gencode_v25lift37_junctions
save(theJunctions, file="junction_annotation_hg19_gencode_v25lift37.rda")


#################################
##### NCBI/RefSeq ###############

## read in GTF as GRanges
refseq_37 = import(con = "/dcl01/lieber/ajaffe/Annotation/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf", format = "gtf")
seqlevels(refseq_37,force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(refseq_37) = si

# get map
refseqMap = mcols(refseq_37)
refseqMap = refseqMap[!duplicated(refseqMap$transcript_id), 
	c("gene_id", "gene_name", "transcript_id", "type")]

##  convert to txdb
refseq_37_txdb = makeTxDbFromGRanges(refseq_37)
#save(refseq_37_txdb, file="refseq_37_txdb.rda")

# get introns
introns = intronsByTranscript(refseq_37_txdb,use.names=TRUE)
introns = unlist(introns)
introns$TranscriptID = names(introns)
introns$refseqID = refseqMap$gene_id[match(introns$TranscriptID, refseqMap$transcript_id)]
introns$symbol = refseqMap$gene_name[match(introns$TranscriptID, 
	refseqMap$transcript_id)]

### make unique junctions
refseq_37_junctions = unique(introns)
oo = findOverlaps(refseq_37_junctions, introns,type="equal")
xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
refseq_37_junctions$tx = CharacterList(xx)
refseq_37_junctions$TranscriptID = NULL
refseq_37_junctions$numTx = elementLengths(refseq_37_junctions$tx)
names(refseq_37_junctions) = NULL
theJunctions = refseq_37_junctions
save(theJunctions, file="junction_annotation_hg19_refseq_grch37.rda")







###################################################################
####  hg38  #######################################################
###################################################################

## chromosome info
chrInfo = getChromInfoFromUCSC("hg38")
chrInfo$chrom = as.character(chrInfo$chrom)
chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg38"))

##########################
##### ENSEMBL ############

## read in GTF as GRanges
ensembl_v85 = import(con = "../ensembl/Homo_sapiens.GRCh38.85.gtf", format = "gtf")
seqlevels(ensembl_v85,force=TRUE) = c(1:22,"X","Y","MT")
seqlevels(ensembl_v85) = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(ensembl_v85) = si

# get map
ensMap = mcols(ensembl_v85)
ensMap = ensMap[!duplicated(ensMap$transcript_id),
	c("gene_id","gene_name","transcript_id","gene_biotype")]
ensMap = ensMap[-is.na(ensMap$transcript_id),]

##  convert to txdb
ensembl_v85_txdb = makeTxDbFromGRanges(ensembl_v85)
#save(ensembl_v85_txdb, file="ensembl_v85_txdb.rda")

# get introns
introns = intronsByTranscript(ensembl_v85_txdb,use.names=TRUE)
introns = unlist(introns)
introns$TranscriptID = names(introns)
introns$ensemblID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

### make unique junctions
ensembl_v85_junctions = unique(introns)
oo = findOverlaps(ensembl_v85_junctions, introns,type="equal")
xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
ensembl_v85_junctions$tx = CharacterList(xx)
ensembl_v85_junctions$TranscriptID = NULL
ensembl_v85_junctions$numTx = elementLengths(ensembl_v85_junctions$tx)
names(ensembl_v85_junctions) = NULL
theJunctions = ensembl_v85_junctions
save(theJunctions, file="junction_annotation_hg38_ensembl_v85.rda")


##########################
##### GENCODE ############

## read in GTF as GRanges
gencode_v25 = import(con = "../GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format = "gtf")
seqlevels(gencode_v25,force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(gencode_v25) = si

# get map
ensMap = mcols(gencode_v25)
ensMap = ensMap[!duplicated(ensMap$transcript_id),
	c("gene_id","gene_name","transcript_id","gene_type")]
ensMap = ensMap[-is.na(ensMap$transcript_id),]

##  convert to txdb
gencode_v25_txdb = makeTxDbFromGRanges(gencode_v25)

# get introns
introns = intronsByTranscript(gencode_v25_txdb,use.names=TRUE)
introns = unlist(introns)
introns$TranscriptID = names(introns)
introns$gencodeID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
introns$ensemblID = ss(introns$gencodeID, "\\.")
introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

### make unique junctions
gencode_v25_junctions = unique(introns)
oo = findOverlaps(gencode_v25_junctions, introns,type="equal")
xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
gencode_v25_junctions$tx = CharacterList(xx)
gencode_v25_junctions$TranscriptID = NULL
gencode_v25_junctions$numTx = elementNROWS(gencode_v25_junctions$tx)
names(gencode_v25_junctions) = NULL
theJunctions = gencode_v25_junctions
save(theJunctions, file="junction_annotation_hg38_gencode_v25.rda")


#################################
##### NCBI/RefSeq ###############

## read in GTF as GRanges
refseq_38 = import(con = "/dcl01/lieber/ajaffe/Annotation/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf", format = "gtf")
seqlevels(refseq_38,force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(refseq_38) = si

# get map
refseqMap = mcols(refseq_38)
refseqMap = refseqMap[!duplicated(refseqMap$transcript_id), 
	c("gene_id", "gene_name", "transcript_id", "type")]

##  convert to txdb
refseq_38_txdb = makeTxDbFromGRanges(refseq_38)
#save(refseq_38_txdb, file="refseq_38_txdb.rda")

# get introns
introns = intronsByTranscript(refseq_38_txdb,use.names=TRUE)
introns = unlist(introns)
introns$TranscriptID = names(introns)
introns$refseqID = refseqMap$gene_id[match(introns$TranscriptID, refseqMap$transcript_id)]
introns$symbol = refseqMap$gene_name[match(introns$TranscriptID, 
	refseqMap$transcript_id)]

### make unique junctions
refseq_38_junctions = unique(introns)
oo = findOverlaps(refseq_38_junctions, introns,type="equal")
xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
refseq_38_junctions$tx = CharacterList(xx)
refseq_38_junctions$TranscriptID = NULL
refseq_38_junctions$numTx = elementLengths(refseq_38_junctions$tx)
names(refseq_38_junctions) = NULL
theJunctions = refseq_38_junctions
save(theJunctions, file="junction_annotation_hg38_refseq_grch38.rda")


