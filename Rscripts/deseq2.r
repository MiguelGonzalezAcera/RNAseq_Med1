# DESeq2 Analysis

library(DESeq2)
library(optparse)
library(data.table)

option_list = list(
  make_option("--counts", type="character",
              help="Table that contains the counts"),
  make_option("--ranges", type="character",
              help="Table that contains the ranges"),
  make_option("--bamfiles", type="character",
              help="list of bam files of the assay, with metadata."),
  make_option("--out_obj", type="character",
              help="DESeq2 object with the result of the analysis."),
  make_option("--organism", type="character", default= "human",
              help="Organism analyzed. Available = human, mouse. Default = Human")
)

opt_parser = OprtionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Read the table with the metadata
sampleTableSingle = read.table(opt$bamfiles, fileEncoding = "UTF8")

# Read the table containing the counts
Counts_tab = read.table(opt$counts, fileEncoding = "UTF8", header=TRUE)
row.names(Counts_tab) <- Counts_tab$Geneid
Counts_tab$Geneid = NULL
Counts_tab <- Counts_tab[,row.names(sampleTableSingle)]
Counts_tab <- Counts_tab[order(row.names(Counts_tab)),]

# Design model matrix (STILL MANUAL)
design <- model.matrix( ~ as.character(sampleTableSingle[,1]) + as.character(sampleTableSingle[,2]))

# Create the experiment from a SummarizedExperiment object
dss <- DESeqDataSetFromMatrix(countData = Counts_tab,
                              colData = sampleTableSingle,
                              design = design)

# Get the genomic ranges
Grang_tab = read.table(opt$ranges, fileEncoding = "UTF8", header=TRUE)
Grlist <- makeGRangesListFromDataFrame(Grang_tab, split.field = "Geneid")
Grlist_filt <- Grlist[names(Grlist) %in% rownames(assay(dss))]
rowRanges(dss) <- Grlist_filt

# filter the counts
keep <- rowSums(counts(dss)) >= 50
dss <- dss[keep,]

# Get and save the fpkm
fpkm_df <- as.data.frame(fpkm(dss))
save(fpkm_df, file=gsub(".Rda","_fpkm.Rda",opt$obj_out, fixed = TRUE))

# Run the analysis
dds <- DESeq(dss)

# Check the names of the main contrasts
resultsNames(dds)

# Plot the dispersion of the set
plotDispEsts(dds)

# Save the normalized counts
# <TO_DO>: The header is odd in the file, so check the samples or load the whole object when working with the normalized counts.
norm_counts <- counts(estimateSizeFactors(dds), normalized = T)
write.table(norm_counts, file=gsub(".Rda","_norm_counts.tsv",opt$obj_out, fixed = TRUE),sep="\t")
df_norm <- as.data.frame(norm_counts)
save(df_norm, file=gsub(".Rda","_norm_counts.Rda",opt$obj_out, fixed = TRUE))

#<TO_DO>: Now, this is doubtful, because now would be time to do the contrasts.
# To do these accuratelz, we have to consider, number of factors, levels of each factor,
#  interaction of the factors, and which factors interact.
# Maybe define the design via input? (Galaxy does it like this)

# Anyway, the output should be a filtered DESeq2 result object and a table with
# the results (named after the contrast used)

# Contrast may vary. Probably need to use loop for all contrasts,
# even more if there is interaction. <TO_DO>: Also, use given threshold
res <- results(dds, alpha = 0.001, contrast=c("Treatment","REC_FUL","HEALTHY"))

# Save the full result object
# <TO_DO>: Change 'contrast' for actual contrast name
save(res,file=gsub(".Rda","_contrast.Rda",opt$obj_out, fixed = TRUE))

#<TO_DO>: Use given threshold
sel = res$padj < 0.1

# Filter by most DE according to padj
sel = replace(sel, is.na(sel), FALSE)
res <- res[sel,]

# Save the list of genes
entrezgeneids <- as.character(mapIds(database, as.character(rownames(resdf)),
                                     'ENTREZID', 'ENSEMBL'))
fwrite(list(entrezgeneids), file = gsub(".Rda","_entrezgeneids.txt",opt$obj_out,
                                        fixed = TRUE))

# Save the complete list of genes.
universeids <- unique(as.character(mapIds(database,
                                          as.character(rownames(assay(Test_experiment))),
                                          'ENTREZID', 'ENSEMBL')))
fwrite(list(universeids), file = gsub(".Rda","_universeids.txt",opt$obj_out,
                                        fixed = TRUE))

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
write.table(resdf, file=gsub(".Rda","_contrast.tsv",opt$obj_out, fixed = TRUE),
            sep="\t", row.names = FALSE)

# Save environment
save.image(file=gsub(".Rda",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda","_versions.tsv",opt$obj_out, fixed = TRUE))
