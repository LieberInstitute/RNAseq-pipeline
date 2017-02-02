## Required libraries
library('getopt')
library('devtools')

## Specify parameters
spec <- matrix(c(
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
        sampleids = 'https://raw.githubusercontent.com/nellore/rail/master/ex/dm3_example.manifest'
    )
}

manifest <- read.table(opt$sampleids, sep = '\t', header = FALSE,
    stringsAsFactors = FALSE)

## Is the data paired end?
paired <- ncol(manifest) > 3
if(paired) system('touch .paired_end')

## Is merging required?
merged <- length(unique(manifest[, ncol(manifest)])) == nrow(manifest)
if(!merged) system('touch .requires_merging')

## Find the file extensions
files <- manifest[, 1]
extensions <- c('fastq.gz', 'fq.gz', 'fastq', 'fq')
patterns <- paste0(extensions, '$')
ext_found <- sapply(files, function(file) {
    extensions[unlist(sapply(patterns, grep, file))[1]]
})
if(any(is.na(ext_found))) {
    stop("Unrecognized fastq filename extension. Should be fastq.gz, fq.gz, fastq or fq")
}

write.table(ext_found, file = '.file_extensions.txt', quote = FALSE,
    row.names = FALSE, col.names = FALSE)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
