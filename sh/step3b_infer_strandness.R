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

save(inferred_strandness ,
    file = 'HISAT2_out/infer_strandness/inferred_strandness.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
