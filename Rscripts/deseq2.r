# DESeq2 Analysis

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gsubfn))

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

option_list <- list(
  make_option("--counts", type = "character",
              help = "Table that contains the counts"),
  make_option("--design", type = "character",
              help = "File with the design of the experiment."),
  make_option("--out_obj", type = "character",
              help = "DESeq2 object with the result of the analysis."),
  make_option("--organism", type = "character", default = "mouse",
              help = "Organism analyzed. Available = human, mouse. Default = mouse"),
  make_option("--control", type = "character",
              help = "Value from the designs to use as control"),
  make_option("--comparisons", type = "character",
              help = "Values from the designs, comma separated, to compare against control.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load R scripts
source("Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Read the table with the metadata
sampleTableSingle <- read.table(opt$design, fileEncoding = "UTF8")

# Read the table containing the counts
Counts_tab <- read.table(opt$counts, fileEncoding = "UTF8", header = TRUE)

# Move the gene IDs as row names
row.names(Counts_tab) <- Counts_tab$Geneid
Counts_tab$Geneid <- NULL

# Select the columns specified in the provided design and resort the genes
Counts_tab <- Counts_tab[, row.names(sampleTableSingle)]
Counts_tab <- Counts_tab[order(row.names(Counts_tab)), ]

# Add row names as column and subset control samples
sampleTableSingle$rn <- row.names(sampleTableSingle)
control_samples <- sampleTableSingle[sampleTableSingle$Tr1 == opt$control,][['rn']]

# Design model matrix, including batch effect correction
Tr1 <- relevel(factor(sampleTableSingle$Tr1), opt$control)
if (length(levels(factor(sampleTableSingle$Batch))) > 1) {
  design <- model.matrix(~ Tr1 + factor(sampleTableSingle$Batch))
} else {
  design <- model.matrix(~ Tr1)
}

# --------------------------------------------------------------

# Create the experiment from a SummarizedExperiment object
dss <- DESeqDataSetFromMatrix(countData = Counts_tab,
                              colData = sampleTableSingle,
                              design = design)

# Avoid normalization
# sizeFactors(dss) <- 1

# Save the universe (of genes). This is important for downstream analyses
save(dss, file = gsub(".Rda", "_universe.Rda", opt$out_obj, fixed = TRUE))

# filter the counts by number
keep <- rowSums(counts(dss)) >= 15
dss <- dss[keep, ]

# Run the analysis
dds <- DESeq(dss, betaPrior = FALSE)

# --------------------------------------------------------------

# Save the normalized counts
# Get the table from the result
norm_counts <- counts(estimateSizeFactors(dds), normalized = TRUE)

# Get the names of the columns
norm_counts_colnames <- colnames(norm_counts)
# Add the gene names as a new column
norm_counts <- cbind(norm_counts, as.character(mapIds(database, as.character(rownames(norm_counts)), 'SYMBOL', 'ENSEMBL')))
# Rename the columns with the new name
colnames(norm_counts) <- c(norm_counts_colnames, "Genename")

# Save the object both as table and as R object
write.table(norm_counts, file=gsub(".Rda","_norm_counts.tsv", opt$out_obj, fixed = TRUE), sep="\t")
df_norm <- as.data.frame(norm_counts)
save(df_norm, file = gsub(".Rda", "_norm_counts.Rda", opt$out_obj, fixed = TRUE))

# --------------------------------------------------------------

# Get the EnsemblIDs from the counts table
df_norm$EnsGenes <- rownames(df_norm)

# Split the comparisons and run the loop to get each table
for (sample in strsplit(opt$comparisons, ",")[[1]]){
  # Using the names provided in the input as the samples, run this as a loop
  res <- results(dds, name = paste("Tr1", sample, sep = ""), cooksCutoff = FALSE)

  # Save the full result object
  # Contrast name will be replaced by the sample and controls
  res_name <- paste(paste("", sample, opt$control, sep = "_"), "Rda", sep = ".")
  save(res, file = gsub(".Rda", res_name, opt$out_obj, fixed = TRUE))

  # A simple helper function that makes a so-called "MA-plot", i.e. a scatter plot of
  # log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis).
  png(file = gsub(".Rda", sprintf("_%s_%s_MA.png", sample, opt$control), opt$out_obj, fixed = TRUE))
  plotMA(res)
  dev.off()

  # Transform result into data frame
  resdf <- data.frame(res)[complete.cases(data.frame(res)), ]

  # Transform row ensembl IDs into column
  resdf$EnsGenes <- rownames(resdf)

  # Add also gene symbols
  resdf$Genes <- as.character(mapIds(database, as.character(rownames(resdf)),
                                     "SYMBOL", "ENSEMBL"))

  # Save table with all the new names. Replace contrast
  res_tab_name <- paste(paste("", sample, opt$control, sep = "_"), "tsv", sep=".")
  write.table(resdf, file = gsub(".Rda",res_tab_name, opt$out_obj, fixed = TRUE),
              sep = "\t", row.names = FALSE)

  # Subset design table with sample to get vector
  sample_samples <- sampleTableSingle[sampleTableSingle$Tr1 == sample, ][["rn"]]

  # Merge res table with the counts of its samples and controls
  resdf_wcounts <- merge(resdf, df_norm[,c("EnsGenes", control_samples, sample_samples)], by="EnsGenes", all.x=TRUE)

  # filter by normalized counts in order to remove false positives
  # Oder of stuff: Select samples or control columns, transform to numeric with the function up,
  # transform to a data matrix, get the medians, Boolean on who's under 25, select rows
  resdf_wcounts$FLAG <- ifelse((rowMedians(data.matrix(sapply(resdf_wcounts[control_samples], as.numeric))) > 25) | (rowMedians(data.matrix(sapply(resdf_wcounts[sample_samples], as.numeric))) > 25), 'OK', 'WARN: Inconsinstent Counts')

  #Save new table
  res_exp_tab_name = paste(paste("", sample, opt$control, sep='_'), "expanded.tsv", sep="_")
  write.table(resdf_wcounts, file=gsub(".Rda", res_exp_tab_name, opt$out_obj, fixed = TRUE),
              sep = "\t", row.names = FALSE)
}

# --------------------------------------------------------------

# Save the transformed normalized counts, with and without batch effect correction
# Transform the object using Variance Stabilizing Transformation
vsd <- vst(dds)

# Get the counts from the transformed object
tr_counts <- assay(vsd)

# Get the names of the columns
tr_counts_colnames <- colnames(tr_counts)
# Add the gene names as a new column
tr_counts <- cbind(tr_counts, as.character(mapIds(database, as.character(rownames(tr_counts)), 'SYMBOL', 'ENSEMBL')))
# Rename the columns with the new name
colnames(tr_counts) <- c(tr_counts_colnames, "Genename")

# Save the object both as table and as R object
write.table(tr_counts, file=gsub(".Rda","_tr_counts.tsv", opt$out_obj, fixed = TRUE), sep="\t")
df_norm <- as.data.frame(tr_counts)
save(df_norm, file = gsub(".Rda", "_tr_counts.Rda", opt$out_obj, fixed = TRUE))

if (length(levels(factor(sampleTableSingle$Batch))) > 1) {
  # Run limma's Batch effect correction and save again
  assay(vsd) <- removeBatchEffect(assay(vsd), vsd$Batch)

  # Get the counts from the transformed object
  tr_B_counts <- assay(vsd)

  # Get the names of the columns
  tr_B_counts_colnames <- colnames(tr_B_counts)
  # Add the gene names as a new column
  tr_B_counts <- cbind(tr_B_counts, as.character(mapIds(database, as.character(rownames(tr_counts)), 'SYMBOL', 'ENSEMBL')))
  # Rename the columns with the new name
  colnames(tr_B_counts) <- c(tr_B_counts_colnames, "Genename")

  # Save the object both as table and as R object
  write.table(tr_B_counts, file=gsub(".Rda","_tr_B_counts.tsv", opt$out_obj, fixed = TRUE), sep="\t")
  df_norm <- as.data.frame(tr_B_counts)
  save(df_norm, file = gsub(".Rda", "_tr_B_counts.Rda", opt$out_obj, fixed = TRUE))
} else {
  # re-Save the object both as table and as R object
write.table(tr_counts, file=gsub(".Rda","_tr_B_counts.tsv", opt$out_obj, fixed = TRUE), sep="\t")
df_norm <- as.data.frame(tr_counts)
save(df_norm, file = gsub(".Rda", "_tr_B_counts.Rda", opt$out_obj, fixed = TRUE))
}

# Save environment
save.image(file = gsub(".Rda", ".RData", opt$out_obj, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda", "_versions.tsv", opt$out_obj, fixed = TRUE))
