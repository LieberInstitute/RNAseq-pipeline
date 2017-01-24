## 
strandfiles <- dir('HISAT2_out/infer_strandness', pattern = 'txt',
    full.names = TRUE)
names(strandfiles) <- gsub('.txt', '', dir('HISAT2_out/infer_strandness',
    pattern = 'txt'))

inferred_strandness <- do.call(rbind, lapply(strandfiles, function(sf) {
    info <- readLines(sf)
    explained <- grep('explained', info)
    data.frame(
        infer_library = gsub('This is | Data', '', info[grep('This is', info)]),
        infer_frac_undertermined = as.numeric(gsub('.*: ', '',
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
