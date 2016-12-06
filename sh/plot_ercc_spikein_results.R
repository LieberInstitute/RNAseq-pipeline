## Load required libraries
library('Biostrings')
library('getopt')
library('rafalib')

## Specify parameters
spec <- matrix(c(
	'maindir', 'm', 1, 'character', 'Main directory',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

MAINDIR <- opt$maindir

###
DF = read.table(paste0(MAINDIR,"/SAMPLE_IDs.txt"))
sampIDs = as.vector(DF[,1])


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
	plot(log2(spikeIns[,"concentration.in.Mix.1..attomoles.ul."]) ~ log2(erccTPM[,i]+1),
		xlab="Kallisto log2(TPM+1)", ylab="Mix 1: log2(Concentration)",
		main = colnames(erccTPM)[i],
		xlim = c(min(log2(erccTPM+1)),max(log2(erccTPM+1))))
	abline(0, 1, lty=2)
}
dev.off()

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
