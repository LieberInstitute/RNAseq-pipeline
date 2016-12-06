##

# wrapper for string split and sapply
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

########### for single-end. leo: make work for both
junctionCount = function(junctionFiles, sampleNames=names(junctionFiles), 
	output = c("Count", "Rail"), minOverhang=0, 
	strandSpecific = FALSE, illuminaStranded=FALSE,
	minCount = 1, maxCores=NULL) {
	
	require(GenomicRanges,quietly=TRUE)
	require(parallel,quietly=TRUE)

	if(is.null(maxCores)) {
		maxCores=1
	}
		
	names(junctionFiles) = sampleNames
	cat("Reading in data.\n")
	if(all(is.character(junctionFiles))) {
		theData = mclapply(junctionFiles, function(x) {
			if(output == "Rail") {
				y = read.delim(x, skip = 1, header=FALSE, 
					colClasses = c("character", "integer", 
					"integer", "integer", "integer", "integer"))
				colnames(y) = c("chr", "start", "end", "leftHang", "rightHang", "count")
				y = y[y$count >= minCount,] # filter based on min number
				y = y[y$leftHang > minOverhang & y$rightHang > minOverhang,]
			} else if(output == "Count") {
				y = read.delim(x, skip = 1, header=FALSE, 
				col.names = c("chr", "start","end", "strand", "count"), 
				colClasses = c("character", "integer", "integer", "character","integer"))
				y = y[y$count >= minCount,] # filter based on min number
				y = y[-which(y$strand=="?"),]
			} else stop("Junction formats can only be from Tophat and Rail.\n")
			
			gr = GRanges(y$chr, IRanges(y$start, y$end), 
				strand=y$strand,count = y$count)
			return(gr)
		}, mc.cores=maxCores)
	} else {
		theData = junctionFiles
		stopifnot(all(sapply(theData, class)=="GRanges"))
	}
	cat("Creating master table of junctions.\n")

	## turn into GRangesList
	### THIS STEP IS SLOW...
	grList = GRangesList(theData)

	# each dataset should be checked
	if(illuminaStranded & strandSpecific) {
		grList = GRangesList(mclapply(grList, function(x) {
			strand(x) = ifelse(strand(x)=="+", "-","+")
			return(x)
		},mc.cores=maxCores))
	}
	
	## get into GRanges object of unique junctions
	fullGR = unlist(grList)
	if(!strandSpecific) strand(fullGR) = "*"
	
	fullGR = fullGR[!duplicated(fullGR)] # or unique(fullGR)
	fullGR = sort(fullGR)
	fullGR$count = NULL

	cat(paste("There are", length(fullGR), "total junctions.\n"))
	
	cat("Populating count matrix.\n")

	jNames = paste0(as.character(seqnames(fullGR)),":",
			start(fullGR),"-",end(fullGR),"(",as.character(strand(fullGR)),")")

	## match GRanges
	options(warn=-1)
	mList = mclapply(grList, match, fullGR, 
		ignore.strand = !strandSpecific, mc.cores=maxCores)
	options(warn=0)
	
	countList = mList # initiate 
	M = length(jNames)

	## fill in matrix
	for(i in seq(along=grList)) {
		if(i %% 25 == 0) cat(".")
		cc = rep(0,M)
		cc[mList[[i]]] = theData[[i]]$count
		countList[[i]] = Rle(cc)
	}
	countDF = DataFrame(countList, row.names=jNames)
	
	names(fullGR) = jNames
	## return matrix and GRanges object
	out = list(countDF = countDF, anno = fullGR)
	return(out)
}


### annotate junctions
annotateJunctions = function(juncCounts, build="hg19") {
    library('bumphunter')
	if(build == "hg19") {
		load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_v75_junction_annotation.rda")
	} else stop("Only supports hg19 for now.\n")
	anno = juncCounts$anno

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

	an = annotateNearest(anno, theJunctions)
	anno$nearestSymbol = as.character(theJunctions$symbol[an$matchIndex])
	anno$nearestDist = an$dist

	# put back in
	juncCounts$anno = anno
	return(juncCounts)
}


## junction stats
junctionStats = function(jRpkm, jMap, 
	cuts = c(0,0.5,1,5,10,20,50,100), output="percent") {

	require(GenomicRanges,quietly = TRUE)
	## junction means
	junctionMeans = rowMeans(jRpkm)
	numAboveCut = sapply(cuts, function(x) sum(junctionMeans >= x))

	## junction stats
	theSeq = paste0(jMap$leftSeq, ":", jMap$rightSeq)
	canJ = ifelse(theSeq %in% c("GT:AG", "CT:AC"), "Canonical", "Not")

	juncList = lapply(cuts, function(x) {
		cat(".")
		expIndex=rep(FALSE, length(junctionMeans))
		expIndex[which(junctionMeans > x)] = TRUE
		
		tabIn = table(jMap$inEnsembl, expIndex, dnn = c("Ensembl", "Exprs"))
		tabNear = table((jMap$inEnsemblEnd | jMap$inEnsemblStart) & 
			!(jMap$inEnsemblEnd & jMap$inEnsemblStart) & !jMap$inEnsembl, 
			expIndex, dnn = c("Ensembl", "Exprs"))
		tabNovel = table(jMap$inEnsemblEnd & jMap$inEnsemblStart & !jMap$inEnsembl, 
			expIndex, dnn = c("Ensembl", "Exprs"))
		tabNew = table(!jMap$inEnsemblEnd & !jMap$inEnsemblStart, 
			expIndex, dnn = c("Ensembl", "Exprs"))
		
		canon = mean(canJ[expIndex] == "Canonical")
		list(tabIn = tabIn, tabNear = tabNear,
			tabNovel = tabNovel,tabNew=tabNew,canon = canon)
	})
	if(output == "percent") {
		ensOut = data.frame(numExprs = numAboveCut, knownTrans = sapply(juncList, 
				function(x) x$tabIn["TRUE","TRUE"]/sum(x$tabIn[,"TRUE"])),
			novelTrans = sapply(juncList, 
				function(x) x$tabNovel["TRUE","TRUE"]/sum(x$tabNovel[,"TRUE"])),
			novelJxn = sapply(juncList, 
				function(x) x$tabNear["TRUE","TRUE"]/sum(x$tabNear[,"TRUE"])),
			newJxn = sapply(juncList, 
				function(x) x$tabNew["TRUE","TRUE"]/sum(x$tabNew[,"TRUE"])),
			canon = sapply(juncList, function(x) x$canon))
		ensOut[,-1] = round(ensOut[,-1]*100,2)
	} else if(output == "number") {
		ensOut = data.frame(numExprs = numAboveCut, 
			knownTrans = sapply(juncList, function(x) x$tabIn["TRUE","TRUE"]),
			novelTrans = sapply(juncList, function(x) x$tabNovel["TRUE","TRUE"]),
			novelJxn = sapply(juncList, function(x) x$tabNear["TRUE","TRUE"]),
			newJxn = sapply(juncList, 	function(x) x$tabNew["TRUE","TRUE"]),
			canon = sapply(juncList, function(x) x$canon))
	} else stop("'output' must be 'percent' or 'number'.\n")		
	rownames(ensOut) = cuts
	return(ensOut)
}
