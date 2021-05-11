# DESeq2 Analysis

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gsubfn))

option_list = list(
  make_option("--counts", type="character",
              help="Table that contains the counts"),
  make_option("--design", type="character",
              help="File with the design of the experiment."),
  make_option("--out_obj", type="character",
              help="DESeq2 object with the result of the analysis."),
  make_option("--organism", type="character", default= "mouse",
              help="Organism analyzed. Available = human, mouse. Default = mouse"),
  make_option("--control", type="character",
              help="Value from the designs to use as control"),
  make_option("--comparisons", type="character",
              help="Values from the designs, comma separated, to compare against control.")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Read the table with the metadata
sampleTableSingle = read.table(opt$design, fileEncoding = "UTF8")

# Read the table containing the counts
Counts_tab = read.table(opt$counts, fileEncoding = "UTF8", header=TRUE)
row.names(Counts_tab) <- Counts_tab$Geneid
Counts_tab$Geneid = NULL
Counts_tab <- Counts_tab[,row.names(sampleTableSingle)]
Counts_tab <- Counts_tab[order(row.names(Counts_tab)),]

# Add row names as column and subset control samples
sampleTableSingle$rn <- row.names(sampleTableSingle)
control_samples <- sampleTableSingle[sampleTableSingle$Tr1 == opt$control,][['rn']]

# Design model matrix
#design <- model.matrix( ~ as.character(sampleTableSingle[,1]) + as.character(sampleTableSingle[,2]))
Tr1 = relevel(sampleTableSingle[,1], opt$control)
design <- model.matrix( ~ Tr1)

# Create the experiment from a SummarizedExperiment object
dss <- DESeqDataSetFromMatrix(countData = Counts_tab,
                              colData = sampleTableSingle,
                              design = design)

# Save the universe (of genes)
save(dss,file=gsub(".Rda","_universe.Rda",opt$out_obj, fixed = TRUE))

# filter the counts
keep <- rowSums(counts(dss)) >= 25
dss <- dss[keep,]

# Run the analysis
dds <- DESeq(dss, betaPrior=FALSE)

# Check the names of the main contrasts
# print(resultsNames(dds))

# Plot the dispersion of the set
# plotDispEsts(dds)

# Save the normalized counts
# <TO_DO>: The header is odd in the file, so check the samples or load the whole object when working with the normalized counts.
norm_counts <- counts(estimateSizeFactors(dds), normalized = T)
norm_counts_colnames <- colnames(norm_counts)
norm_counts <- cbind(norm_counts, as.character(mapIds(database, as.character(rownames(norm_counts)), 'SYMBOL', 'ENSEMBL')))
colnames(norm_counts) <- c(norm_counts_colnames, "Genename")

write.table(norm_counts, file=gsub(".Rda","_norm_counts.tsv",opt$out_obj, fixed = TRUE),sep="\t")
df_norm <- as.data.frame(norm_counts)
save(df_norm, file=gsub(".Rda","_norm_counts.Rda",opt$out_obj, fixed = TRUE))

df_norm$EnsGenes <- rownames(df_norm)

#<TO_DO>: Now, this is doubtful, because now would be time to do the contrasts.
# To do these accuratelz, we have to consider, number of factors, levels of each factor,
#  interaction of the factors, and which factors interact.
# Maybe define the design via input? (Galaxy does it like this)

# Anyway, the output should be a filtered DESeq2 result object and a table with
# the results (named after the contrast used)

# Contrast may vary. Probably need to use loop for all contrasts,
# even more if there is interaction. <TO_DO>: Also, use given threshold
# res <- results(dds, alpha = 0.001, contrast=c("Treatment","REC_FUL","HEALTHY"))

# Split the comparisons and run the loop of tables
for (sample in strsplit(opt$comparisons, ",")[[1]]){
  # Using the names provided in the input as the samples, run this as a loop
  res <- results(dds, name = paste("Tr1",sample,sep=""), cooksCutoff=FALSE)
  # Save the full result object
  # Contrast name will be replaced by the sample and controls
  res_name = paste(paste("", sample, opt$control, sep='_'),"Rda", sep=".")
  save(res,file=gsub(".Rda",res_name,opt$out_obj, fixed = TRUE))
  
  # A simple helper function that makes a so-called "MA-plot", i.e. a scatter plot of
  # log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis).
  png(file=gsub(".Rda","_MA.png",opt$out_obj, fixed = TRUE))
  plotMA(res)
  dev.off()
  
  # Transform result into data frame
  resdf <- data.frame(res)[complete.cases(data.frame(res)),]
  
  # Transform row ensembl IDs into column
  resdf$EnsGenes <- rownames(resdf)
  
  # Add also gene symbols
  resdf$Genes <- as.character(mapIds(database, as.character(rownames(resdf)),
                                     'SYMBOL', 'ENSEMBL'))
  
  # Save table with all the new names. Replace contrast
  res_tab_name = paste(paste("", sample, opt$control, sep='_'),"tsv", sep=".")
  write.table(resdf, file=gsub(".Rda",res_tab_name,opt$out_obj, fixed = TRUE),
              sep="\t", row.names = FALSE)
  
  # Subset design table with sample to get vector
  sample_samples <- sampleTableSingle[sampleTableSingle$Tr1 == sample,][['rn']]
  
  # Merge res table with the counts of its samples
  resdf_wcounts <- merge(resdf, df_norm[,c("EnsGenes", control_samples, sample_samples)], by='EnsGenes', all.x=TRUE)
  
  #Save new table
  res_exp_tab_name = paste(paste("", sample, opt$control, sep='_'),"expanded.tsv", sep="_")
  write.table(resdf_wcounts, file=gsub(".Rda",res_exp_tab_name,opt$out_obj, fixed = TRUE),
              sep="\t", row.names = FALSE)
}
# Save environment
save.image(file=gsub(".Rda",".RData",opt$out_obj, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda","_versions.tsv",opt$out_obj, fixed = TRUE))
