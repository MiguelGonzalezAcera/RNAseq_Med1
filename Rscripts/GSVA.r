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
  make_option("--dims", type = "character", default = "3500,3500",
              help = "Dimensions of the plot in pixels. Default = 3500,3500"),
  make_option("--colors", type = "character", default = "blue,white,red",
              help = "Colors for the heatmap, from lower to higher. Default = blue,white,red"),
  make_option("--limits", type = "character", default = "-1,0,1",
              help = "Limits and center for the color scale. Default = -1,0,1"),
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

# load the respective set of gene markers according to organism
if (opt$organism == "mouse") {
  source("Static/Mouse_genemarkers_ensembl.Rda")
} else if (opt$organism == "human") {
  source("Static/Human_genemarkers_ensembl.Rda")
}

# transform the data to numeric matrix
options(digits = 20)
numwdf_norm <- matrix(as.numeric(as.matrix(wdf_norm)), ncol = ncol(as.matrix(wdf_norm)))
rownames(numwdf_norm) <- rownames(wdf_norm)
colnames(numwdf_norm) <- colnames(wdf_norm)

# Run the gsva ove the set of counts
gsva.es <- gsva(numwdf_norm, genes, verbose = FALSE)

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
        row_dend_width = unit(2, "cm"))
dev.off()

# ----------------------------------------------------------------

# Run a differential expression analysis and make a heatmap with it
# diff expression with limma
sampleTableSingle <- read.table(opt$design, fileEncoding = "UTF8")

# Design model matrix
Tr1 <- relevel(factor(sampleTableSingle[, 1]), opt$control)
design <- model.matrix(~ Tr1)

# Run limma
fit <- lmFit(gsva.es, design)
fit <- eBayes(fit)
res <- decideTests(fit, p.value = 0.1)

# Get a dataframe with the columns of the fold change of all the samples
clust_df <- NULL
pval_df <- NULL
filenames <- c()

# Save each of the results by sample
for (sample in strsplit(opt$comparisons, ",")[[1]]){
  # Get result of the diff expression
  restab <- topTable(fit, coef = paste("Tr1", sample, sep = ""), number = 100)

  # Save the table
  restab_name <- paste(paste("", sample, opt$control, sep= " "), "tsv", sep = ".")
  write.table(restab, file = gsub(".Rda", restab_name, opt$out_obj, fixed = TRUE),
              sep = "\t")

  # Transform into dataframe
  full_df <- as.data.frame(restab)

  # Get the fold change and the pvalue column
  FC_df <- full_df["logFC"]
  pv_df <- full_df["P.Value"]

  # Make the filename and rename the columns
  filename <- paste(sample, opt$control, sep = "_")
  # Name the column as the file and save the name
  colnames(FC_df) <- c(filename)
  colnames(pv_df) <- c(filename)
  filenames <- c(filenames, filename)

  if (is.null(clust_df) == TRUE) {
    # If the final df is empty, fill it with one column
    clust_df <- FC_df
    pval_df <- pv_df
  } else {
    # If not, add the column to the df
    clust_df <- merge(clust_df, FC_df, by = 0, all = TRUE)
    rownames(clust_df) <- clust_df$Row.names
    clust_df$Row.names <- NULL

    pval_df <- merge(pval_df, pv_df, by = 0, all = TRUE)
    rownames(pval_df) <- pval_df$Row.names
    pval_df$Row.names <- NULL
  }
}

# Drop Na values
clust_df <- clust_df[complete.cases(clust_df), ]
pval_df <- pval_df[complete.cases(pval_df), ]

# Keep only genes with valid pvalues
clust_df <- clust_df[rownames(clust_df) %in% rownames(pval_df), drop = FALSE]

# Change column names
colnames(clust_df) <- filenames
colnames(pval_df) <- filenames

# Failsafe for clusterings with low instances
if (length(rownames(clust_df)) < 2) {
  print("GSVA doesn\'t have the necessary length to do the clustering in this group of samples")
} else {
  # Perform the clustering analysis over the table
  # Tree construction (rows)
  hr <- hclust(as.dist(1 - cor(t(data.matrix(clust_df)),
                            method = "pearson")), method = "average")

  # Tree cutting
  mycl <- cutree(hr, h = max(hr$height) / 1.3)

  # Clustering boxes
  mycolhc <- rainbow(length(unique(mycl)), start = 0.1, end = 0.9)
  mycolhc <- mycolhc[as.vector(mycl)]

  # change heatmap filename for the Fc one
  FC_hmap_path <- gsub(".png", "_FC.png", opt$heatmap, fixed = TRUE)

  png(
    file = FC_hmap_path,
    width = int(strsplit(opt$dims, ",")[0]),
    height = int(strsplit(opt$dims, ",")[1]),
    res = 600
  )

  # Mount the heatmap
  row_den <- color_branches(hr, h = max(hr$height) / 1.5)
  Heatmap(
    data.matrix(clust_df), cluster_rows = as.dendrogram(hr),
    cluster_columns = FALSE, col = color, row_dend_width = unit(3, "cm"),
    row_names_gp = gpar(fontsize = (90 / length(rownames(clust_df)) + 5)),
    column_names_gp = gpar(fontsize = (90 / length(rownames(clust_df)) + 5) + 2),
    column_names_max_height = unit(8, "cm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        sprintf("%.2f", data.matrix(pval_df)[i, j]), x, y,
        gp = gpar(fontsize = (80 / length(rownames(clust_df)) + 3))
      )
    }
  )
  dev.off()
}

# Save environment
save.image(file = gsub(".Rda", ".RData", opt$out_obj, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda", "_versions.tsv", opt$out_obj, fixed = TRUE))
