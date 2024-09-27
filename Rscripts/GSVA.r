suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(cluster))

# Create options
option_list <- list(
  make_option("--heatmap", type = "character",
              help = "Heatmap of the coverage of the genes over the samples"),
  make_option("--counts", type = "character",
              help = "An r object with the normalized counts. Produced in the DE script."),
  make_option("--design", type = "character",
              help = "File with the design of the experiment."),
  make_option("--out_obj", type = "character",
              help = "DESeq2 object with the result of the analysis."),
  make_option("--organism", type = "character", default = "mouse",
              help = "Organism analyzed. Available = human, mouse. Default = mouse"),
  make_option("--control", type = "character",
              help = "Value from the designs to use as control"),
  make_option("--comparisons", type = "character",
              help = "Values from the designs, comma separated, to compare against control."),
  make_option("--genesets", type = "character", default = '',
              help = "List of files with the sets of genes desired, ensemblIDs, comma separated."),
  make_option("--dims", type = "character", default = "3500,3500",
              help = "Dimensions of the plot in pixels. Default = 3500,3500"),
  make_option("--colors", type = "character", default = "purple,white,orange",
              help = "Colors for the heatmap, from lower to higher. Default = purple,white,orange"),
  make_option("--limits", type = "character", default = "-1,0,1",
              help = "Limits and center for the color scale. Default = -1,0,1")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load R scripts
source("Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Load the r object containing the data.
load(opt$counts)

# Remove the genename column
wdf_norm <- df_norm[, !names(df_norm) %in% c("Genename")]

# Create the lists of genes
genes <- list()

# Check if there are genelists to be used
if (opt$genesets == '') {
  # load the respective set of gene markers according to organism
  if (opt$organism == "mouse") {
    load("Static/Mouse_genemarkers_ensembl.Rda")
  } else if (opt$organism == "human") {
    load("Static/Human_genemarkers_ensembl.Rda")
  }
} else {
  for (genefile in strsplit(opt$genesets, ",")[[1]]) {
    genefile_name = gsub("_ensembl.txt","", tail(strsplit(genefile, "/")[[1]], n=1), fixed = TRUE)

    genes[[genefile_name]] <- readLines(genefile)
  }
}

# transform the data to numeric matrix
options(digits = 20)
numwdf_norm <- matrix(as.numeric(as.matrix(wdf_norm)), ncol = ncol(as.matrix(wdf_norm)))
rownames(numwdf_norm) <- rownames(wdf_norm)
colnames(numwdf_norm) <- colnames(wdf_norm)

# Run the gsva ove the set of counts
gsvaPar <- gsvaParam(numwdf_norm, genes)
gsva.es <- gsva(gsvaPar, verbose = FALSE)

# Save the GSVA table as R object and table
save(gsva.es, file = opt$out_obj)
write.table(gsva.es, file = gsub(".Rda",".tsv", opt$out_obj, fixed = TRUE), sep = "\t")

# -------------------------------------------------------------------------

# Make the trnasformed counts heatmap
# Establish colors
color <- colorRamp2(
  c(
    as.integer(strsplit(opt$limits, ",")[[1]][1]),
    as.integer(strsplit(opt$limits, ",")[[1]][2]),
    as.integer(strsplit(opt$limits, ",")[[1]][3])
  ),
  c(
    strsplit(opt$colors, ",")[[1]][1],
    strsplit(opt$colors, ",")[[1]][2],
    strsplit(opt$colors, ",")[[1]][3]
  )
)

# Create a heatmap of the values
png(
  file = opt$heatmap,
  width = as.integer(strsplit(opt$dims, ",")[[1]][1]),
  height = as.integer(strsplit(opt$dims, ",")[[1]][2]),
  res = 600
)
Heatmap(gsva.es, cluster_columns = FALSE,
        col = color, column_dend_height = unit(5, "cm"),
        row_dend_width = unit(2, "cm"),
        heatmap_legend_param = list(
          title = "Relative counts",
          at = c(
            as.integer(strsplit(opt$limits, ",")[[1]][1]) * 2,
            as.integer(strsplit(opt$limits, ",")[[1]][1]),
            as.integer(strsplit(opt$limits, ",")[[1]][2]),
            as.integer(strsplit(opt$limits, ",")[[1]][3]),
            as.integer(strsplit(opt$limits, ",")[[1]][3]) * 2
            )
        )
        )
dev.off()

# ----------------------------------------------------------------

# Run a differential expression analysis
# diff expression with limma
sampleTableSingle <- read.table(opt$design, fileEncoding = "UTF8")

# Design model matrix
Tr1 <- relevel(factor(sampleTableSingle[, 1]), opt$control)
design <- model.matrix(~ Tr1)

# Run limma
fit <- lmFit(gsva.es, design)
fit <- eBayes(fit)
res <- decideTests(fit, p.value = 0.1)

# Save each of the results by sample
for (sample in strsplit(opt$comparisons, ",")[[1]]){
  # Get result of the diff expression
  restab <- topTable(fit, coef = paste("Tr1", sample, sep = ""), number = 100)

  # Save the table
  restab_name <- paste(paste("", sample, opt$control, sep= "_"), "tsv", sep = ".")
  write.table(restab, file = gsub(".Rda", restab_name, opt$out_obj, fixed = TRUE),
              sep = "\t")
}

# Save environment
save.image(file = gsub(".Rda", ".RData", opt$out_obj, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda", "_versions.tsv", opt$out_obj, fixed = TRUE))
