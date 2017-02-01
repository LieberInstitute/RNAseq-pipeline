## Required libraries
library('getopt')

## Specify parameters
spec <- matrix(c(
    'fastq', 'f', 1, 'character', 'Fastq folder',
    'sampleids', 's', 2, 'character', 'path to SAMPLE_IDs.txt file',
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
if(FALSE) {
    opt <- list(
        fastq = '/dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX',
        sampleids = 'https://raw.githubusercontent.com/nellore/rail/master/ex/dm3_example.manifest'
    )
}

manifest <- read.table(opt$sampleids, sep = '\t', header = FALSE)

## Is the data paired end?
paired <- ncol(manifest) > 3
if(paired) system('touch .paired_end')

## Is merging required?
merged <- length(unique(manifest[, ncol(manifest)])) == nrow(manifest)
if(!merged) system('touch .requires_merging')
    
## Add the fastq folder path if necessary
if(opt$fastq != '') {
    write.table(manifest, file = '.SAMPLE_IDs_original.txt', sep = '\t',
        col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    manifest[, 1] <- file.path(opt$fastq, manifest[, 1])
    if(paired) {
        manifest[, 3] <- file.path(opt$fastq, manifest[, 3])
    }
    write.table(manifest, file = 'SAMPLE_IDs.txt', sep = '\t',
        col.names = FALSE, row.names = FALSE, quote = FALSE)
}

## Find the extension

files <- system(paste0('ls ', manifest[1, 1], '*'), intern = TRUE)

extensions <- c('fastq.gz', 'fq.gz', 'fastq', 'fq')
patterns <- paste0(extensions, '$')
present <- sapply(lapply(patterns, grepl, files), any)
result <- extensions[present]

if(length(result) == 0) {
    error("Unrecognized fastq filename extension.")
}

message(paste(Sys.time(), "the following extensions are available:",
    paste(result, collapse = ', ')))

if(length(result) > 1) result <- result[1]
message(paste(Sys.time(), 'using the following extension', result))

write.table(result, file = '.FILE_extension.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

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
