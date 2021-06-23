# snakemake variables
transcript_counts = snakemake@input$transcript_counts
annotation_file = snakemake@input$annotation_file
samples_file = snakemake@input$samples_file
converted_gtf_file = snakemake@params$converted_gtf_file
gene_counts = snakemake@output$gene_counts
log_file = snakemake@log[[1]]

# log all console output
log = file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('### Variables used for this analysis:\n')
cat('transcript_counts: ', transcript_counts, '\n', sep = " , ")
cat('annotation_file: "', annotation_file, '"\n', sep = "")
cat('samples_file: "', samples_file, '"\n', sep = "")
cat('converted_gtf_file: "', converted_gtf_file, '"\n', sep = "")
cat('gene_counts: ', gene_counts, '\n', sep = " , ")
cat('log_file: "', log_file, '"\n', sep = "")
cat('\n\n')

cat('### Sessioninfo:\n')
sessionInfo()
cat('\n\n')

suppressMessages(library("tximport", quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("GenomicFeatures", quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("readr", quietly = TRUE, warn.conflicts = FALSE))

if (grepl("\\.gff$", annotation_file) || grepl("\\.gff3$", annotation_file)) {
  if (file.exists(converted_gtf_file)) {
    cat('\n### Producing TxDb object...\n\n')
    txdb = makeTxDbFromGFF(converted_gtf_file)
    cat('\n### TxDb object produced.\n\n')
  }else {
    cat('\n### Annotation file is in GFF format, converting to GTF using GFFRead...\n\n')
    system(paste0("gffread ", annotation_file, " -T -o ", converted_gtf_file, " >> ", log_file, " 2>&1"))
    cat('\n### Annotation GTF convertion finished.\n\n')
    cat('\n### Producing TxDb object...\n\n')
    txdb = makeTxDbFromGFF(converted_gtf_file)
    cat('\n### TxDb object produced.\n\n')
  }
}else {
  cat('\n### Producing TxDb object...\n\n')
  txdb = makeTxDbFromGFF(annotation_file)
  cat('\n### TxDb object produced.\n\n')
}
k = keys(txdb, keytype = "GENEID")
cat('\n### Producing Tx2Gene table...\n\n')
tx2gene = select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene = tx2gene[, 2:1]
cat('\n### Tx2Gene table produced.\n\n')
cat('\n### Tx2Gene table preview:\n')
head(tx2gene)
cat('\n\n')

samples = read.table(samples_file, header = TRUE)
files = transcript_counts
names(files) = paste0(samples$sample_id)
cat('\n### Producing TXi Counts table...\n\n')
txi.salmon = tximport(files, type = "salmon", tx2gene = tx2gene)
cat('\n### TXi Counts table produced.\n\n')

cat('\n### TXi Counts table preview: \n')
head(txi.salmon$counts)
cat('\n\n')

for(i in 1:length(samples$sample_id)) {
  cat("### Writting Gene counts output for: ", samples$sample_id[i], "\n", sep = "")
  write.table("GeneID\tRaw_Counts", file = gene_counts[i], quote = FALSE, sep = '\t', row.names = F, col.names = F )
  write.table(txi.salmon$counts[,c(samples$sample_id[i])], file = gene_counts[i], quote = FALSE, sep = '\t', append = TRUE, col.names = F)
}