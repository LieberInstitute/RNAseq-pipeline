## Required libraries
library('getopt')

## Specify parameters
spec <- matrix(c(
    'file', 'f', 1, 'character', 'Path file without extensions',
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
        file = '/dcl01/lieber/ajaffe/Nina/GSK_PhaseII/data/Sample_R10126_C1BP4ACXX/R10126_C1BP4ACXX_GAGATTCC_L005'
    )
}


files <- system(paste0('ls ', opt$file, '*'), intern = TRUE)

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
