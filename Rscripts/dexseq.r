# DEXeq Analysis. Tutorial from https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-transcript-usage/dtu-using-dexseq#stager-procedure

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(stageR))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(DRIMSeq))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gsubfn))
suppressPackageStartupMessages(library(IHW))

option_list <- list(
  make_option("--salmon_counts", type = "character",
              help = "Folder with the results of the salmon mapping."),
  make_option("--annotation", type = "character",
              help = "Location of the table with the annotation per transcript."),
  make_option("--design", type = "character",
              help = "File with the design of the experiment."),
  make_option("--out_obj", type = "character",
              help = "DEXseq object with the result of the analysis."),
  make_option("--organism", type = "character", default = "mouse",
              help = "Organism analyzed. Available = human, mouse. Default = mouse"),
  make_option("--control", type = "character",
              help = "Value from the designs to use as control"),
  make_option("--comparison", type = "character",
              help = "Values from the designs to compare against control.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load R scripts
source("Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Read the table with the metadata
sampleTableSingle <- read.table(opt$design, fileEncoding = "UTF8")

# Add row names as column and subset control and experimental samples
sampleTableSingle$rn <- row.names(sampleTableSingle)
control_samples <- sampleTableSingle[sampleTableSingle$Tr1 == opt$control,][['rn']]
experim_samples <- sampleTableSingle[sampleTableSingle$Tr1 == opt$comparison,][['rn']]

# Relevel the design table to the control
sampleTableSingle$Tr1 <- relevel(factor(sampleTableSingle$Tr1), opt$control)

# Select only the samples requested
sampleTableSingle <- sampleTableSingle[c(control_samples, experim_samples),]

# Load the annotation table
annotTab <- read.table(opt$annotation, sep = '\t', header = TRUE,fileEncoding = "UTF8")

# Generate columns of the ensembl id with the version
annotTab$GeneID <- paste(annotTab$gene_id, annotTab$gene_version, sep=".")
annotTab$TranscriptID <- paste(annotTab$transcript_id, annotTab$transcript_version, sep=".")

# Get a table that relates transcript ID and gene ID, with versions
tx2gene <- annotTab[annotTab$item == 'transcript',][c('TranscriptID', 'GeneID')]

colnames(tx2gene) <- c('TXNAME','GENEID')

# add the chromosome information
tx2gene_chr <- annotTab[annotTab$item == 'chromosome',][c('chr','chr')]

colnames(tx2gene_chr) <- c('TXNAME','GENEID')

# concatenate
tx2gene <- rbind(tx2gene, tx2gene_chr)

# -----------------------------------------------------

#Get the sf files from the provided folder
files <- list.files(opt$salmon_counts, pattern="*.sf", full.names=TRUE)
names(files) <- gsub(".sf", "", list.files(opt$salmon_counts, pattern="*.sf"))

# Select just the files in the design
files <- files[row.names(sampleTableSingle)]

# -----------------------------------------------------

# Read the salmon files with the DTU moniker
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txIn = TRUE, txOut = TRUE, countsFromAbundance = "dtuScaledTPM")

# Save the counts tables for registries and possible reanalyses. Notice the file per each comparison.
write.table(txi.salmon$counts, file=paste(opt$salmon_counts, sprintf("%s_%s_salmon_dtuScaledTPM.tsv", opt$comparison, opt$control), sep = "/"), sep="\t")
txi.salmon.counts <- as.data.frame(txi.salmon$counts)
save(txi.salmon.counts, file = paste(opt$salmon_counts,sprintf("%s_%s_salmon_dtuScaledTPM.Rda", opt$comparison, opt$control), sep = "/"))

# Build a dataframe with just the counts and filter the empty genes
cts = txi.salmon$counts
cts = cts[rowSums(cts) > 0,]

# ------------------------------------------------------

# Pre-filter with DRIMSeq dmFilter

# Create a coutns dataframe
txdf.sub <- tx2gene[match(rownames(cts),tx2gene$TXNAME),]
counts = data.frame(gene_id = txdf.sub$GENEID, feature_id = txdf.sub$TXNAME, cts)

# Make a version of the design table with a different layout (DRIMSeq is picky)
samps = data.frame(sample_id = rownames(sampleTableSingle), group = factor(sampleTableSingle$Tr1), batch = factor(sampleTableSingle$Batch))

# Create a dmDSdata object
d = DRIMSeq::dmDSdata(counts = counts, samples = samps)

# filter (parameters from the tutorial)
n = nrow(samps)
n.small = min(table(samps$group))

d = DRIMSeq::dmFilter(d,
                      min_samps_feature_expr = n.small, min_feature_expr = 10,
                      min_samps_feature_prop = n.small, min_feature_prop = 0.1,
                      min_samps_gene_expr = n, min_gene_expr = 10)

# ------------------------------------------------------

# DEXSeq

# get the counts matrix filtered and build the DEXSeq object
countData = round(as.matrix(counts(d)[,-c(1:2)]))

# Create the basic dexseq object
dxd = DEXSeqDataSet(countData = countData, sampleData = samps, 
                    design = ~sample + exon + group:exon, 
                    featureID = counts(d)$feature_id, groupID = counts(d)$gene_id)

# Perform DEXSeq. Consider potential batch effects (this one comes from the DEXSeq vignette)
if (length(levels(factor(sampleTableSingle$Batch))) > 1) {
    # Generate the model to employ
    formulaFullModel = ~ sample + exon + batch:exon + group:exon
    formulaReducedModel = ~ sample + exon + batch:exon

    # Run DEXSeq with the considerations
    dxd = estimateSizeFactors(dxd)
    dxd = estimateDispersions(dxd, formula = formulaFullModel)
    dxd = testForDEU(dxd, reducedModel = formulaReducedModel, fullModel = formulaFullModel)
} else {
    # Run DEXSeq raw
    dxd = estimateSizeFactors(dxd)
    dxd = estimateDispersions(dxd)
    dxd = testForDEU(dxd)
}

# Get the results
dxr = DEXSeqResults(dxd, independentFiltering = FALSE)

# Get the filtered ones for flagging later
dxr_filt = DEXSeqResults(dxd, independentFiltering = TRUE)

dxr_filt_lst = setdiff(dxr_filt$featureID, dxr$featureID)

# Get a per-gene adjusted p-val, aggregating evidence for multiple tests per gene
# the important table is dxr.g, which can be merged in the result table.
qval = perGeneQValue(dxr)
dxr.g = data.frame(gene = names(qval), qval)
dxr.t = as.data.frame(dxr[, c("featureID","groupID","pvalue")])

# ------------------------------------------------------

# StageR analysis
# (I don't quite know what this is for, but I think is for p value adjustment in two-stages tests??)

# Construct an object with the per-gene p-values for the screening
pScreen = qval
names(pScreen) <- gsub("[.].*$", "", as.character(names(pScreen)), perl = TRUE)

# Get the per-transcript confirmation p-values
pConfirmation = matrix(dxr.t$pvalue, ncol=1)
dimnames(pConfirmation) = list(gsub("[.].*$", "", as.character(dxr.t$featureID), perl = TRUE),"transcript")

# Make a conversor table with and without the ensembl version
tx2gene_version = data.frame(dxr.t[,c("featureID", "groupID")], 
                     dxr.t[,c("featureID", "groupID")])
for (i in 1:2) tx2gene_version[,i] = gsub("[.].*$", "", as.character(tx2gene_version[,i]), perl = TRUE)

# Construct the object and perform the stageR analysis
# qval used in pScreen, hence pScreenAdjusted=TRUE
stageRObj = stageRTx(pScreen = pScreen, 
                     pConfirmation = pConfirmation, 
                     pScreenAdjusted = TRUE, 
                     tx2gene = tx2gene_version[1:2]
                     )

# Note: There must be a way to keep all the genes, instead of just the significant ones
stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05, allowNA=TRUE)

dex.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)

dex.padj = merge(tx2gene_version, dex.padj, by.x = c("groupID","featureID"), by.y = c("geneID","txID"))

# Change the colnames for easier merging later
colnames(dex.padj) <- c("groupID", "featureID", "featureID.1", "groupID.1", "gene_pval", "transcript_pval")

# ------------------------------------------------------

# Export the data

# Get the normalized counts per transcript
# The /2 in the rows is because DEXSeq produces double rows, one set with the counts in the exons (per its design)
# and another with the /rest/ of the counts. Since I'm only interested in the counts per transcript, i have to chop it.
dex.norm = cbind(as.data.frame(stringr::str_split_fixed(rownames(counts(dxd)), ":", 2)), as.data.frame(counts(dxd, normalized = TRUE))[,1:(ncol(counts(dxd))/2)])
colnames(dex.norm) = c("groupID", "featureID", as.character(colData(dxd)$sample_id)[1:(ncol(counts(dxd))/2)])
row.names(dex.norm) = NULL

# Per-group normalized mean
dex.mean = as.data.frame(sapply( levels(samps$group), function(lvl) rowMeans(dex.norm[, 3:ncol(dex.norm)][, samps$group == lvl, drop = FALSE])))

# log2FC of the expression
dex.log2fc = log2(dex.mean[opt$comparison]/dex.mean[opt$control])
colnames(dex.log2fc) = "log2fc"
rownames(dex.log2fc) = dex.norm$featureID

# Bind all columns in order in a single table
dexData = cbind(dex.norm[,1:2], dex.mean)
dexData = merge(dexData, dex.padj[,c("featureID.1","groupID.1","gene_pval","transcript_pval")], by.x = c("groupID","featureID") , by.y = c("groupID.1", "featureID.1"))
dexData = merge(dexData, dex.log2fc, by.x = "featureID", by.y = "row.names")
dexData = merge(annotTab[annotTab$item == 'transcript',][c("chr","start","end","strand","gene_name","gene_biotype","transcript_name","transcript_biotype","ccds_id","tag","GeneID","TranscriptID")], dexData, by.x = c("GeneID","TranscriptID"), by.y = c("groupID","featureID"))
dexData = merge(dexData, dex.norm, by.x = c("GeneID","TranscriptID"), by.y = c("groupID","featureID"))

# Flag the flags
# filter by normalized counts in order to remove false positives
# Oder of stuff: Select samples or control columns, transform to numeric with the function up,
# transform to a data matrix, get the medians, Boolean on who's under 25, select rows
dexData$FLAG <- ifelse(
  dexData$TranscriptID %in% dxr_filt_lst,
  'FAIL: Filtered by indFilt',
  ifelse(
    (rowMedians(data.matrix(sapply(dexData[control_samples], as.numeric))) > 25) | (rowMedians(data.matrix(sapply(dexData[experim_samples],as.numeric))) > 25),
    ifelse(
      (
        (
          rowSds(data.matrix(sapply(dexData[control_samples], as.numeric)))*2 > rowMeans(data.matrix(sapply(dexData[control_samples], as.numeric)))
        ) & (
          rowSums(data.matrix(sapply(dexData[control_samples], as.numeric))) > 25
        )
      ) | (
        (
          rowSds(data.matrix(sapply(dexData[experim_samples], as.numeric)))*2 > rowMeans(data.matrix(sapply(dexData[experim_samples], as.numeric)))
        ) & (
          rowSums(data.matrix(sapply(dexData[experim_samples], as.numeric))) > 25
        )
      ),
      "INFO: High variation in condition",
      'OK'
    ),
    'WARN: Inconsinstent Counts'
  )
)

# Save it
# Save as R object
res_exp_name = paste(paste("", opt$comparison, opt$control, sep='_'), ".Rda", sep="")
save(dexData, file = gsub(".Rda", res_exp_name, opt$out_obj, fixed = TRUE))

#Save new table
res_exp_tab_name = paste(paste("", opt$comparison, opt$control, sep='_'), ".tsv", sep="")
write.table(dexData, file=gsub(".Rda", res_exp_tab_name, opt$out_obj, fixed = TRUE),
            sep = "\t", row.names = FALSE)

# Save environment
save.image(file = gsub(".Rda", ".RData", opt$out_obj, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda", "_versions.tsv", opt$out_obj, fixed = TRUE))