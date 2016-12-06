## Required libraries
library('derfinder')
library('BiocParallel')
library('stringr')
library('jaffelab')

##
args = commandArgs(TRUE)
hgXX = args[1]
MAINDIR = args[2]
EXPERIMENT = args[3]
PREFIX = args[4]

EXPNAME = paste0(EXPERIMENT,"_",PREFIX)

## read in pheno	
pd = data.frame(read.table(paste0(MAINDIR,"/SAMPLE_IDs.txt"), as.is=TRUE, header=F))
names(pd)[1] = "SAMPLE_ID"
N = length(pd$SAMPLE_ID)

### add bam file
pd$bamFile = paste0(MAINDIR, "/HISAT2_out/", pd$SAMPLE_ID, "_accepted_hits.sorted.bam")
pd$bwFile = paste0(MAINDIR, "/Coverage/", pd$SAMPLE_ID, ".bw")

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

## Chrs to use, mitocondrial chromosome has to be the last one for the code
## to work later on
if (hgXX == "rn6") { CHR = c(1:20,"X","Y","MT")
} else if (hgXX == "mm10") { CHR = paste0("chr",c(1:19,"X","Y","M"))
} else { CHR = paste0("chr",c(1:22,"X","Y","M")) }
stopifnot(grepl('M', CHR[length(CHR)]))

### confirm total mapping
pd$totalMapped <- unlist(bplapply(pd$bamFile, getTotalMapped,
    chrs = CHR[-length(CHR)], BPPARAM = MulticoreParam(8)))
pd$mitoMapped <- unlist(bplapply(pd$bamFile, getTotalMapped,
    chrs = CHR[length(CHR)], BPPARAM = MulticoreParam(8)))
pd$mitoRate <- pd$mitoMapped / (pd$mitoMapped +  pd$totalMapped)


###################################################################arg
fullCov <- fullCoverage(files=pd$bamFile, chrs = CHR, mc.cores=12)

save(pd, fullCov, compress=TRUE, file=paste0(MAINDIR,"/fullCoverage_",EXPNAME,"_n",N,".rda"))

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

