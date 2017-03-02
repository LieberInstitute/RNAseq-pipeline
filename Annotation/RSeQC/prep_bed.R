## Prepare BED files for http://rseqc.sourceforge.net/#infer-experiment-py
# module load R/3.3.x
# Rscript prep_bed.R -> prep_bed_log.txt 2>&1

library('rtracklayer')
library('devtools')
gtfs <- file.path('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation',
    c('ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf',
    file.path('GENCODE', c('GRCh37_hg19/gencode.v25lift37.annotation.gtf',
    'GRCh38_hg38/gencode.v25.annotationGRCh38.gtf',
    'GRCm38_mm10/gencode.vM11.annotation.gtf'))))
names(gtfs) <- c('rn6', 'hg19', 'hg38', 'mm10')

for(i in seq_along(gtfs)) {
    message(paste(Sys.time(), 'processing', gtfs[i]))
    gr <- import(gtfs[i])
    if('score' %in% colnames(mcols(gr))) {
        mcols(gr) <- mcols(gr)[, -which(colnames(mcols(gr)) == 'score')]
    }
    gr <- gr[gr$type == 'gene']
    message(paste(Sys.time(), 'exporting', paste0(names(gtfs)[i], '.bed')))
    export(gr, paste0(names(gtfs)[i], '.bed'), format = 'bed')
}

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
