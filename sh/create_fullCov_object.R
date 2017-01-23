## Required libraries
library('derfinder')
library('BiocParallel')
library('jaffelab')
library('getopt')

## Specify parameters
spec <- matrix(c(
	'organism', 'o', 1, 'character', 'Either rn6, mm10 or human',
	'maindir', 'm', 1, 'character', 'Main directory',
	'experiment', 'e', 1, 'character', 'Experiment',
	'prefix', 'p', 1, 'character', 'Prefix',
    'paired', 'l', 1, 'logical', 'Whether the reads are paired-end or not',
    'fullcov', 'f', 1, 'logical', 'Whether to create the full coverage object or not',
    'cores', 'c', 1, 'integer', 'Number of cores to use',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}


hgXX <- opt$organism
MAINDIR <- opt$maindir
EXPERIMENT <- opt$experiment
PREFIX <- opt$prefix
PE <- opt$paired

EXPNAME <- paste0(EXPERIMENT,"_",PREFIX)

if(opt$fullcov) {
    ## read in pheno	
    pd <- data.frame(read.table(file.path(MAINDIR, 'SAMPLE_IDs.txt'), as.is=TRUE,
        header = FALSE))
    names(pd)[1] <- "SAMPLE_ID"
    pd$SAMPLE_ID <- basename(pd$SAMPLE_ID)
    N <- length(pd$SAMPLE_ID)

    ### add bigwig and bam files
    pd$bamFile <- file.path(MAINDIR, 'HISAT2_out', paste0(pd$SAMPLE_ID,
        '_accepted_hits.sorted.bam'))
    pd$bwFile <- file.path(MAINDIR, 'Coverage', paste0(pd$SAMPLE_ID, '.bw'))

    ## Chrs to use, mitocondrial chromosome has to be the last one for the code
    ## to work later on
    if (hgXX == "rn6") { CHR = c(1:20,"X","Y","MT")
    } else if (hgXX == "mm10") { CHR = paste0("chr",c(1:19,"X","Y","M"))
    } else { CHR = paste0("chr",c(1:22,"X","Y","M")) }
    stopifnot(grepl('M', CHR[length(CHR)]))

    ### confirm total mapping
    pd$totalMapped <- unlist(bplapply(pd$bamFile, getTotalMapped,
        chrs = CHR[-length(CHR)], BPPARAM = MulticoreParam(opt$cores)))
    pd$mitoMapped <- unlist(bplapply(pd$bamFile, getTotalMapped,
        chrs = CHR[length(CHR)], BPPARAM = MulticoreParam(opt$cores)))
    pd$mitoRate <- pd$mitoMapped / (pd$mitoMapped +  pd$totalMapped)

    ###################################################################

    fullCov <- fullCoverage(files = pd$bwFile, chrs = CHR, mc.cores = opt$cores)
    save(fullCov, file = file.path(MAINDIR, paste0('fullCoverage_', EXPNAME,
        '_n', N, '.rda')))
}


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
