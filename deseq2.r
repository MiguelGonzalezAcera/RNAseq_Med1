# DESeq2 Analysis

library(DESeq2)
library(optparse)
library(data.table)

option_list = list(
  make_option("--robj", type="character",
              help="R object that contains the counts"),
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
source("D:/Documentos/LATESIS/Scripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Read the table with the metadata
sampleTableSingle = read.table(opt$bamfiles, fileEncoding = "UTF8")

# Load the object containing the counts
load(opt$robj)

# Design model matrix (STILL UNKNOWN HOW TO DO)
design <- model.matrix(~UNKNOWN)

# Create the experiment from a SummarizedExperiment object
dss <- DESeqDataSet(Test_experiment, design = design)

# filter the counts
keep <- rowSums(counts(dss)) >= 50
dss <- dss[keep,]

# Run the analysis
dds <- DESeq(dss)

# Check the names of the main contrasts
resultsNames(dds)

# Plot the dispersion of the set
plotDispEsts(dds)

# Save the normalized counts
norm_counts <- counts(estimateSizeFactors(dds), normalized = T)
write.table(norm_counts, file=gsub(".Rda","_norm_counts.tsv",opt$obj_out, fixed = TRUE),
            sep="\t", row.names = FALSE)

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