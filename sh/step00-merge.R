## Required libraries
library('getopt')
library('BiocParallel')
library('qualV')

## Specify parameters
spec <- matrix(c(
    'sampleids', 's', 1, 'character', 'Path to the SAMPLE_IDs.txt file',
	'maindir', 'm', 1, 'character', 'Main directory',
    'paired', 'p', 1, 'logical', 'Whether the reads are paired-end or not',
    'extenstion', 'e', 1, 'character', 'The file extension',
    'cores', 'c', 1, 'integer', 'Number of cores to use'
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
        sampleids = '/users/lcollado/SAMPLE_IDs.txt',
        maindir = '/users/lcollado',
        paired = TRUE,
        extension = 'fastq.gz'
    )
}

samples <- read.table(opt$sampleids, stringsAsFactors = FALSE)
outdir <- file.path(opt$maindir, 'merged_fastq')
dir.create(outdir, showWarnings = FALSE)

file_list <- split(samples$V1, samples$V2)

## Get unique names for the samples
new_names <- sapply(file_list, function(common) {
    n <- length(common)
    if(n == 1) return(common)
    chars <- strsplit(common, '')
    result <- chars[[1]]
    j <- 2
    while(j <= n) {
        result <- LCS(result, chars[[j]])$LCS
        j <- j + 1
    }
    result <- paste(result, collapse = '')
    return(result)
})

merge_files <- function(file_names, new_file) {
    message(paste(Sys.time(), 'creating', new_file))
    call <- paste('cat', paste(file_names, collapse = ' '), '>', new_file)
    system(call)
}

res <- bpmapply(function(common, new_name) {
    if(opt$paired) {
        for(read in paste0(c('_R1_001.', '_R2_001.'), opt$extension)) {
            merge_files(paste0(common, read),
                file.path(outdir, paste0(new_name, read)))
        }
    } else {
        read <- paste0('.', opt$extension)
        merge_files(paste0(common, read),
            file.path(outdir, paste0(new_name, read)))
    }
}, file_list, new_names, BPPRAM = MulticoreParam(opt$cores))

message(paste(Sys.time(), 'creating .SAMPLE_IDs_backup.txt'))
system(paste('mv', opt$sampleids,
    paste0('.', gsub('.txt', '_backup.txt', opt$sampleids))))

message(paste(Sys.time(), 'creating the new SAMPLE_IDs.txt file with the merged samples'))
write.table(new_names, file = opt$sampleids, row.names = FALSE, col.names = FALSE, quote = FALSE)

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
