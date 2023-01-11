# DESeq2 Analysis
library(DESeq2)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism("human")

# Read the table with the metadata
sampleTableSingle = read.table("/VAULT/20220126_GSE90830_CancerCellLines/design.txt", fileEncoding = "UTF8")

# Read the table containing the counts
Counts_tab = read.table("/VAULT/20220126_GSE90830_CancerCellLines/counts.tsv", fileEncoding = "UTF8", header=TRUE)
row.names(Counts_tab) <- Counts_tab$Geneid
Counts_tab$Geneid = NULL
Counts_tab <- Counts_tab[,row.names(sampleTableSingle)]
Counts_tab <- Counts_tab[order(row.names(Counts_tab)),]

# Add row names as column and subset control samples
sampleTableSingle$rn <- row.names(sampleTableSingle)
control_samples <- sampleTableSingle[sampleTableSingle$Tr1 == 'CACO2',][['rn']]

# Design model matrix
#design <- model.matrix( ~ as.character(sampleTableSingle[,1]) + as.character(sampleTableSingle[,2]))
Tr1 = relevel(sampleTableSingle[,1], 'CACO2')
design <- model.matrix( ~ Tr1)

# Create the experiment from a SummarizedExperiment object
dss <- DESeqDataSetFromMatrix(countData = Counts_tab,
                              colData = sampleTableSingle,
                              design = design)

# Keep the universe
save(dss,file="/VAULT/20191216_Marta_Request/universe.rda")

# Get the genomic ranges
Grang_tab = read.table("/VAULT/20220126_GSE90830_CancerCellLines/counts.ranges.fix.tsv", fileEncoding = "UTF8", header=TRUE)
Grlist <- makeGRangesListFromDataFrame(Grang_tab, split.field = "Geneid")
Grlist_filt <- Grlist[names(Grlist) %in% rownames(assay(dss))]

keep <- rownames(assay(dss)) %in% names(Grlist)
dss <- dss[keep,]

rowRanges(dss) <- Grlist_filt

# filter the counts
keep <- rowSums(counts(dss)) >= 25
dss <- dss[keep,]

# Get and save the fpkm
fpkm_df <- as.data.frame(fpkm(dss))

fpkm_df_colnames <- colnames(fpkm_df)
fpkm_df <- cbind(fpkm_df, as.character(mapIds(database, as.character(rownames(fpkm_df)), 'SYMBOL', 'ENSEMBL')))
colnames(fpkm_df) <- c(fpkm_df_colnames, "Genename")

save(fpkm_df, file="/VAULT/20220126_GSE90830_CancerCellLines/detables/norm_counts.fpkm.Rda")
write.table(fpkm_df, file='/VAULT/20220126_GSE90830_CancerCellLines/detables/norm_counts.fpkm.tsv',sep="\t")

# Run the analysis
dds <- DESeq(dss)

# Only for cases with high variation
dds <- DESeq(dss, betaPrior=FALSE)

# Check the names of the main contrasts
resultsNames(dds)

# Plot the dispersion of the set
plotDispEsts(dds)

# Save the normalized counts
# <TO_DO>: The header is odd in the file, so check the samples or load the whole object when working with the normalized counts.
norm_counts <- counts(estimateSizeFactors(dds), normalized = T)
write.table(norm_counts, file="/VAULT/20191216_Marta_Request/norm_counts.tsv",sep="\t")
df_norm <- as.data.frame(norm_counts)
save(df_norm, file="/VAULT/20191216_Marta_Request/norm_counts.Rda")

#<TO_DO>: Now, this is doubtful, because now would be time to do the contrasts.
# To do these accuratelz, we have to consider, number of factors, levels of each factor,
#  interaction of the factors, and which factors interact.
# Maybe define the design via input? (Galaxy does it like this)

# Anyway, the output should be a filtered DESeq2 result object and a table with
# the results (named after the contrast used)

# Contrast may vary. Probably need to use loop for all contrasts,
# even more if there is interaction. <TO_DO>: Also, use given threshold
res <- results(dds, name = 'Tr1B6mTEC')

# High variation cases
res <- results(dds, name = 'Tr1O12dc', cooksCutoff=FALSE)

# Save the full result object
# <TO_DO>: Change 'contrast' for actual contrast name
save(res,file="/VAULT/20191216_Marta_Request/B6mTEC_v_B6cTEC.rda")

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
write.table(resdf, file="/VAULT/20191216_Marta_Request/B6mTEC_v_B6cTEC.tsv",
            sep="\t", row.names = FALSE)

# Save environment
save.image(file=gsub(".Rda",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda","_versions.tsv",opt$obj_out, fixed = TRUE))
