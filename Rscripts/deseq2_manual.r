# DESeq2 Analysis
library(DESeq2)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism("mouse")

# Read the table with the metadata
sampleTableSingle = read.table("/VAULT/Human_data/Macrophage_stim/design.txt", fileEncoding = "UTF8")

# Read the table containing the counts
Counts_tab = read.table("/VAULT/Human_data/Macrophage_stim/counts.tsv", fileEncoding = "UTF8", header=TRUE)
row.names(Counts_tab) <- Counts_tab$Geneid
Counts_tab$Geneid = NULL
Counts_tab <- Counts_tab[,row.names(sampleTableSingle)]
Counts_tab <- Counts_tab[order(row.names(Counts_tab)),]

# Design model matrix (STILL MANUAL)
Tr1 = relevel(sampleTableSingle[,1],"BMDM-Control")
Tr2 = relevel(sampleTableSingle[,2], "mock")
design <- model.matrix( ~ Tr1)

# Create the experiment from a SummarizedExperiment object
dss <- DESeqDataSetFromMatrix(countData = Counts_tab,
                              colData = sampleTableSingle,
                              design = design)

# Keep the universe
save(dss,file="/VAULT/Human_data/Macrophage_stim/universe.rda")

# Get the genomic ranges
Grang_tab = read.table("/DATA/DSS_rec_evolution/DSS_rec_evol.counts.ranges.tsv", fileEncoding = "UTF8", header=TRUE)
Grlist <- makeGRangesListFromDataFrame(Grang_tab, split.field = "Geneid")
Grlist_filt <- Grlist[names(Grlist) %in% rownames(assay(dss))]
rowRanges(dss) <- Grlist_filt

# filter the counts
keep <- rowSums(counts(dss)) >= 25
dss <- dss[keep,]

# Get and save the fpkm
fpkm_df <- as.data.frame(fpkm(dss))
save(fpkm_df, file="/DATA/DSS_rec_evolution/DSS_rec_evol.fpkm.Rda")

# Run the analysis
dds <- DESeq(dss)

# Check the names of the main contrasts
resultsNames(dds)

# Plot the dispersion of the set
plotDispEsts(dds)

# Save the normalized counts
# <TO_DO>: The header is odd in the file, so check the samples or load the whole object when working with the normalized counts.
norm_counts <- counts(estimateSizeFactors(dds), normalized = T)
write.table(norm_counts, file="/VAULT/Human_data/Macrophage_stim/norm_counts.tsv",sep="\t")
df_norm <- as.data.frame(norm_counts)
save(df_norm, file="/VAULT/Human_data/Macrophage_stim/norm_counts.Rda")

#<TO_DO>: Now, this is doubtful, because now would be time to do the contrasts.
# To do these accuratelz, we have to consider, number of factors, levels of each factor,
#  interaction of the factors, and which factors interact.
# Maybe define the design via input? (Galaxy does it like this)

# Anyway, the output should be a filtered DESeq2 result object and a table with
# the results (named after the contrast used)

# Contrast may vary. Probably need to use loop for all contrasts,
# even more if there is interaction. <TO_DO>: Also, use given threshold
res <- results(dds, alpha = 0.001, name = 'Tr1BMDM.IL4.12h')

# Save the full result object
# <TO_DO>: Change 'contrast' for actual contrast name
save(res,file="/VAULT/Human_data/Macrophage_stim/BMDM.IL4.12h.rda")

# A simple helper function that makes a so-called "MA-plot", i.e. a scatter plot of
# log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis).
png(file=gsub(".Rda","_MA.png",opt$obj_out, fixed = TRUE))
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
write.table(resdf, file="/VAULT/Human_data/Macrophage_stim/BMDM.IL4.12h.tsv",
            sep="\t", row.names = FALSE)

# Save environment
save.image(file=gsub(".Rda",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda","_versions.tsv",opt$obj_out, fixed = TRUE))
