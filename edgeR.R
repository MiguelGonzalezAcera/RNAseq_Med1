# EdgeR Analysis

library(edgeR)
library(optparse)
library(data.table)

option_list = list(
  make_option("--robj", type="character",
              help="R object that contains the counts"),
  make_option("--bamfiles", type="character",
              help="list of bam files of the assay, with metadata."),
  make_option("--out_obj", type="character",
              help="EdgeR object with the result of the analysis."),
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

# Filter genes that have not at least 1cpm in at leats two samples
#<TO_DO>: Figure out how to insert the grouping

# Create the edgeR object
dge.n = DGEList(counts=assay(Test_experiment),
                group=colData(Test_experiment)[,"Treatment"])

# Get the counts per million
countsPerMillion <- cpm(dge.n)
summary(countsPerMillion)

# Obtain the rows that have at least 1 cpm in two samples
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dge.n <- dge.n[keep,]
summary(cpm(dge.n))

# Normalize the counts
# Calculate normalization factors
# Possible methods: "TMM","TMMwzp","RLE","upperquartile","none". Read details for references.
dge.n = calcNormFactors(dge.n, method = "TMM")

# Calculate the normalized counts once the normalization factors have been set.
dge.n = estimateCommonDisp(dge.n)

# Extract the normalized counts.
norm_counts = dge.n$pseudo.counts

# Re-create the object with the normalized counts
#<TO_DO>: Figure out how to insert the grouping
dge = DGEList(counts = norm_counts,
              group = colData(Test_experiment)[,"Treatment"],
              genes=rownames(norm_counts))

# Create the design matrix
#<TO_DO>: Set the treatments inputs correctly
# Extracted from https://support.bioconductor.org/p/69802/
design <- model.matrix(~Treatment)

# Estimate a common negative binomial dispersion parameter
dge <- estimateGLMCommonDisp(dge, design, verbose=TRUE)

# Estimate the abundance-dispersion trend by Cox-Reid approximate profile likelihood
dge <- estimateGLMTrendedDisp(dge, design)

# Compute an empirical Bayes estimate of the negative binomial dispersion parameter
# for each tag, with expression levels specified by a log-linear model.
dge <- estimateGLMTagwiseDisp(dge, design)

# Plot the biological coeficient of variation.
# x: log counts per million
# y: BCV
# If trend is too high, there is too much variation between samples.
# Less expressed genes always have a lot of variation.
plotBCV(dge)

# Once studied the dispersions, adjust to model and do DE anaysis
#<TO_DO>: the design thingy
# Fit a negative binomial generalized log-linear model to the read counts for each gene
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, contrast = c(0,0,1,0,0))

# Plot the genes over and under the threshold
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags = deGenes)
abline(h=c(-1, 1), col=2)

# Transform first 500 results into formatted table
resdf <- as.data.frame(topTags(lrt, n=500, adjust.method = "BH", p.value = 0.1))
colnames(resdf) <- c("EnsGenes","log2FoldChange","logCPM","LR","pvalue","padj")
resdf$EnsGenes <- rownames(resdf)
resdf$Genes <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(resdf)),
                                   'SYMBOL', 'ENSEMBL'))

# Save table with all the new names. Replace contrast
write.table(resdf, file=gsub(".Rda","_contrast.tsv",opt$obj_out, fixed = TRUE),
            sep="\t", row.names = FALSE)

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

# Save environment
save.image(file=gsub(".Rda",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda","_versions.tsv",opt$obj_out, fixed = TRUE))