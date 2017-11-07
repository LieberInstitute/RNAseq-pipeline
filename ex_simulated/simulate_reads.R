# Rscript simulate_reads.R  > simulate_reads_log.txt 2>&1
library('polyester')
library('Biostrings')
library('ShortRead')
library('devtools')

fastapath = system.file("extdata", "chr22.fa", package="polyester")
numtx = count_transcripts(fastapath)
readmat = matrix(20, ncol=10, nrow=numtx)
readmat[1:30, 1:5] = 40

simulate_experiment_countmat(fasta=fastapath, 
    readmat=readmat, outdir='single_end_unstranded', seed=5, paired = FALSE, strand_specific = FALSE)
    
simulate_experiment_countmat(fasta=fastapath, 
    readmat=readmat, outdir='paired_end_unstranded', seed=5, paired = TRUE, strand_specific = FALSE)

simulate_experiment_countmat(fasta=fastapath, 
    readmat=readmat, outdir='single_end_stranded', seed=5, paired = FALSE, strand_specific = TRUE)
    
simulate_experiment_countmat(fasta=fastapath, 
    readmat=readmat, outdir='paired_end_stranded', seed=5, paired = TRUE, strand_specific = TRUE)

# reads <- readFasta('paired_end_stranded/sample_01_1.fasta')
# reads2 <- readFasta('paired_end_stranded/sample_01_2.fasta')
# identical(width(reads), width(reads2))

dirs <- c('single_end_unstranded', 'single_end_stranded')

fasta_files <- dir(dirs, pattern = '.fasta', full.names = TRUE)

write_myfastq <- function(reads, ff) {
    #id <- BStringSet(paste0('read', seq_len(length(reads))))
    fastq <- ShortReadQ(sread(reads), quality = BStringSet(rep(paste(rep('I', max(width(reads))), collapse = ''), length(reads))), id = id(reads))
    writeFastq(fastq, gsub('fasta', 'fastq.gz', ff))
}

## Create fastq files (all with quality I which should be = 40)
xx <- sapply(fasta_files, function(ff) {
    reads <- readFasta(ff)
    reads <- reads[width(reads) == max(width(reads))]
    write_myfastq(reads, ff)
})
unlink(fasta_files)

## now process paired-end reads
dirs <- c('paired_end_unstranded', 'paired_end_stranded')
fasta_files_1 <- dir(dirs, pattern = '_1.fasta', full.names = TRUE)
fasta_files_2 <- dir(dirs, pattern = '_2.fasta', full.names = TRUE)
stopifnot(identical(gsub('_1.fasta', '', fasta_files_1), gsub('_2.fasta', '', fasta_files_2)))

xx <- mapply(function(ff1, ff2) {
    reads1 <- readFasta(ff1)
    reads2 <- readFasta(ff2)
    stopifnot(length(reads1) == length(reads2))
    w1 <- which(width(reads1) == max(width(reads1)))
    w2 <- which(width(reads2) == max(width(reads1)))
    keep <- intersect(w1, w2)
    
    stopifnot(identical(width(reads1[keep]), width(reads2[keep])))
    stopifnot(identical(id(reads1[keep]), id(reads2[keep])))
    
    write_myfastq(reads1[keep], ff1)
    write_myfastq(reads2[keep], ff2)
}, fasta_files_1, fasta_files_2)
unlink(c(fasta_files_1, fasta_files_2))


## Create manifest files
dirs <- c('single_end_unstranded', 'paired_end_unstranded', 'single_end_stranded', 'paired_end_stranded')
xx <- sapply(dirs, function(d) {
    if(grepl('paired', d)) {
        df <- data.frame(path = dir(d, pattern = '_1.fastq'), md5sum = 0, path = dir(d, pattern = '_2.fastq'), md5sum2 = 0, name = gsub('_2.fastq.gz|_', '', dir(d, pattern = '_2.fastq')))
    } else {
        df <- data.frame(path = dir(d, pattern = 'fastq'), md5sum = 0, name = gsub('.fastq.gz|_', '', dir(d, pattern = 'fastq')))
    }
    write.table(df, file = file.path(d, 'samples.manifest'), col.names = FALSE, quote = FALSE, row.names = FALSE, sep = '\t')
    
})

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
