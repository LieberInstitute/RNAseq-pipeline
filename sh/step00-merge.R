## Required libraries
library('getopt')
library('BiocParallel')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'sampleids', 's', 1, 'character', 'Path to the samples.manifest file',
	'outdir', 'o', 1, 'character', 'Full path to directory where the merged fastq files will be saved to',
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

## For testing
testing <- FALSE
if(testing) {
    opt <- list(
        sampleids = 'https://raw.githubusercontent.com/nellore/rail/master/ex/dm3_example.manifest',
        outdir = 'merged_fastq',
        cores = 1
    )
} 

manifest <- read.table(opt$sampleids, sep = '\t', header = FALSE,
    stringsAsFactors = FALSE)
    
if(testing) {
    ## For testing
    manifest[, ncol(manifest)] <- rep('sample1', nrow(manifest))
}

## Is the data paired end?
paired <- ncol(manifest) > 3

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
extensions <- split(ext_found, manifest[, ncol(manifest)])

## Check that extensions are the same per group
if(any(sapply(extensions, function(x) { length(unique(x)) }) != 1)) {
    stop("For each sample name, the extensions of the fastq files to be merged have to be the same")
}

## Create the output directory
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## Split according to the sample names
file_groups <- split(manifest, manifest[, ncol(manifest)])
extensions <- sapply(extensions, '[', 1)

merge_files <- function(file_names, new_file) {
    message(paste(Sys.time(), 'creating', new_file))
    call <- paste('cat', paste(file_names, collapse = ' '), '>', new_file)
    print(call)
    if(!testing) system(call)
}

res <- bpmapply(function(common, new_name, extension) {
    merge_files(common[, 1],
        file.path(opt$outdir, paste0(new_name, '.', extension)))
    if(paired) {
        merge_files(common[, 3],
            file.path(opt$outdir, paste0(new_name, '_read2.', extension)))
    }
}, file_groups, names(file_groups), extensions, 
    BPPARAM = MulticoreParam(opt$cores))

message(paste0(Sys.time(), ' creating .SAMPLE_IDs_backup_', Sys.Date(), '.txt'))
system(paste('mv', opt$sampleids, file.path(dirname(opt$sampleids),
    paste0('.', gsub('.txt', paste0('_backup_', Sys.Date(), '.txt'),
    basename(opt$sampleids))))
))

message(paste(Sys.time(), 'creating the new samples.manifest file with the merged samples'))

if(paired) {
    new_manifest <- data.frame(
        file.path(opt$outdir, paste0(names(file_groups), '.', extensions)),
        rep(0, length(file_groups)),
        file.path(opt$outdir, paste0(names(file_groups), '_read2.',
            extensions)),
        rep(0, length(file_groups)),
        names(file_groups), stringsAsFactors = FALSE
    )
} else {
    new_manifest <- data.frame(
        file.path(opt$outdir, paste0(names(file_groups), '.', extensions)),
        rep(0, length(file_groups)),
        names(file_groups), stringsAsFactors = FALSE
    )
}
## Make names short, in case you want to interactively check the new manifest
colnames(new_manifest) <- paste0('V', seq_len(ncol(new_manifest)))

write.table(new_manifest, file = opt$sampleids, row.names = FALSE,
    col.names = FALSE, quote = FALSE)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
