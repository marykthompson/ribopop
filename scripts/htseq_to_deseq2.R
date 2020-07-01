#Run gene-level diffexp analysis on the htseq quantification
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#htseq-count-input
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

units_file <- snakemake@params[["units_file"]]
samples_file <- snakemake@params[["samples_file"]]

units <- read.table(units_file, sep = '\t', header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
samples <- read.table(samples_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#add the condition data to the units df
units <- merge(units, samples, by = "sample")
#add the name of the run, this will match the name in the htseq dir
units$run <- paste(units$sample, units$unit, sep = "-")

this_contrast = snakemake@params[["contrast"]]

files <- snakemake@input[["infiles"]]
names(files) <- sapply(files, function(x) unlist(strsplit(x, "/"))[2])
fdf <- data.frame(files)

#add the file names for the htseq results
units <- merge(units, fdf, by.x = "run", by.y = "row.names")

#Subset the data by experment:
exp_df <- units[units$condition %in% this_contrast,]
#only run files that are present
exp_df <- exp_df[exp_df$run %in% names(files)]

condition_df <- exp_df[, c("run", "condition", "files")]
names(condition_df) <- c("sampleName", "condition", "fileName")
condition_df$condition <- factor(condition_df$condition)

#reorder the columns to be as DESeq expects
col_order <- c("sampleName", "fileName", "condition")
condition_df <- condition_df[, col_order]

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = condition_df, design = ~ condition)
dds <- DESeq(ddsHTSeq, sfType = snakemake@params[["sf_method"]])
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file = snakemake@output[["table"]])

# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()
