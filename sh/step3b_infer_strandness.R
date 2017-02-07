## In case this step crashes, assume that the data is not stranded
write.table('none', file = 'inferred_strandness_pattern.txt',
    row.names = FALSE, col.names = FALSE, quote = FALSE)

## Infer strandness
strandfiles <- dir('HISAT2_out/infer_strandness', pattern = 'txt',
    full.names = TRUE)
names(strandfiles) <- gsub('.txt', '', dir('HISAT2_out/infer_strandness',
    pattern = 'txt'))

inferred_strandness <- do.call(rbind, lapply(strandfiles, function(sf) {
    message(paste(Sys.time(), 'processing', sf))
    info <- readLines(sf)
    
    if(any(grepl('Unknown Data type', info))) {
        warning(paste('Unknown data type for', sf))
        return(data.frame(
            infer_library = 'Unknown',
            infer_frac_undetermined = NA,
            infer_pattern1 = NA,
            infer_frac_pattern1 = NA,
            infer_pattern2 = NA,
            infer_frac_pattern2 = NA,
            stringsAsFactors = FALSE
        ))
    }
    
    explained <- grep('explained', info)
    data.frame(
        infer_library = gsub('This is | Data', '', info[grep('This is', info)]),
        infer_frac_undetermined = as.numeric(gsub('.*: ', '',
            info[grep('determine', info)])),
        infer_pattern1 = gsub('".*', '', gsub('.*by "', '',
            info[explained[1]])),
        infer_frac_pattern1 = as.numeric(gsub('.*: ', '', info[explained[1]])),
        infer_pattern2 = gsub('".*', '', gsub('.*by "', '',
            info[explained[2]])),
        infer_frac_pattern2 = as.numeric(gsub('.*: ', '', info[explained[2]])),
        stringsAsFactors = FALSE
    ) 
}))

## Print some info
lapply(inferred_strandness[, -grep('frac', colnames(inferred_strandness))], table, useNA = 'ifany')
summary(inferred_strandness[, grep('frac', colnames(inferred_strandness))])

save(inferred_strandness ,
    file = 'HISAT2_out/infer_strandness/inferred_strandness.Rdata')

## Explore visually the results
pdf('HISAT2_out/infer_strandness/inferred_strandness.pdf')
lims <- round(with(inferred_strandness, range(c(infer_frac_pattern1,
    infer_frac_pattern2), na.rm = TRUE)) + c(-0.05, 0.05), 1)
lims[1] <- max(0, lims[1])
lims[2] <- min(1, lims[2])
with(inferred_strandness, plot(infer_frac_pattern1, infer_frac_pattern2,
    xlab = paste('Fraction of reads with pattern 1:',
    names(sort(table(inferred_strandness$infer_pattern1)))[1]),
    ylab = paste('Fraction of reads with pattern 2:',
    names(sort(table(inferred_strandness$infer_pattern2)))[1]),
    ylim = lims, xlim = lims
))
abline(a = 1, b = -1, col = 'red', lwd = 2)
abline(h = 0.5, col = 'grey', lty = 2)
abline(v = 0.5, col = 'grey', lty = 2)
boxplot(inferred_strandness[, grep('frac', colnames(inferred_strandness))],
    ylim = c(0, 1), col = 'light blue', ylab = 'Fraction of reads')
dev.off()

## Determine which pattern to use
observed <- mean(inferred_strandness$infer_frac_pattern1 >
    inferred_strandness$infer_frac_pattern2, na.rm = TRUE)
cutoff <- 0.2
pattern <- ifelse(observed < 0 + cutoff, 
    names(sort(table(inferred_strandness$infer_pattern2)))[1],
    ifelse(observed > 1 - cutoff,
        names(sort(table(inferred_strandness$infer_pattern1)))[1], 'none'))
message(paste(Sys.time(), 'will use the following pattern for bam2wig:',
    pattern))
write.table(pattern, file = 'inferred_strandness_pattern.txt',
    row.names = FALSE, col.names = FALSE, quote = FALSE)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
