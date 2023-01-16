suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(optparse))


option_list = list(
  make_option("--Rdata", type="character", default = "",
              help="Files with the data. Useful when comparing projects."),
  make_option("--colnames", type="character", default = "",
              help="Names of the columns, in the same order as Rdata. Useful when comparing projects."),
  make_option("--heatmap", type="character",
              help="Heatmap of the coverage of the genes over the samples"),
  make_option("--organism", type="character", default= "mouse",
              help="Organism analyzed. Available = human, mouse. Default = mouse")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Read this from the metadata file.
# These are the full DE analysis results of all the samples
treats = strsplit(opt$Rdata, ",")[[1]]
colnames = strsplit(opt$colnames, ",")[[1]]

# 1.- FOLD CHANGE
# Get a dataframe with the columns of the fold change of all the samples
clust_df <- NULL
pval_df <- NULL
for (filename in treats){
  # Read each file
  load(filename)
  full_df <- as.data.frame(restab)
  
  # Get the fold change and the pvalue column
  FC_df <- full_df["logFC"]
  pv_df <- full_df["P.Value"]
  
  filename_col <- gsub(".Rda","", tail(strsplit(filename, "/"), n=1), fixed = TRUE)
  # Name the column as the file
  colnames(FC_df) <- c(filename)
  colnames(pv_df) <- c(filename)
  
  if (is.null(clust_df) == T){
    # If the final df is empty, fill it with one column
    clust_df <- FC_df
    pval_df <- pv_df
  } else {
    # If not, add the column to the df
    clust_df <- merge(clust_df, FC_df, by=0,all=T)
    rownames(clust_df) <- clust_df$Row.names
    clust_df$Row.names <- NULL
    
    pval_df <- merge(pval_df, pv_df, by=0,all=T)
    rownames(pval_df) <- pval_df$Row.names
    pval_df$Row.names <- NULL
  }
}
# Drop Na values
clust_df <- clust_df[complete.cases(clust_df), ]
pval_df <- pval_df[complete.cases(pval_df), ]

# Keep only genes with valid pvalues
clust_df <- clust_df[rownames(clust_df) %in% rownames(pval_df), ,drop=FALSE]

# Change column names
colnames(clust_df) <- colnames
colnames(pval_df) <- colnames

# Failsafe for clusterings with low instances
if (length(rownames(clust_df)) < 2) {
  print('GSVA doesn\'t have the necessary length to do the clustering in this group of samples')
  quit()
}

# Perform the clustering analysis over the table
# Tree construction (rows and columns)
hr <- hclust(as.dist(1-cor(t(data.matrix(clust_df)),
                           method="pearson")), method="average")
hc <- hclust(as.dist(1-cor(data.matrix(clust_df),
                           method="pearson")), method="average")

# Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.3)

# Clustering boxes
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
mycolhc <- mycolhc[as.vector(mycl)]

# Establish colors
color <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

png(file=opt$heatmap, width = 2000, height = 2000, res = 300)
# Mount the heatmap
#<TO_DO>: Add the title of the plot, according to whatever
row_den = color_branches(hr, h = max(hr$height)/1.5)
Heatmap(data.matrix(clust_df), cluster_rows = as.dendrogram(hr),
        #as.dendrogram(row_den),
        cluster_columns = FALSE,
        col=color, column_dend_height = unit(5, "cm"),
        #row_dend_width = unit(3, "cm"),
        row_names_gp = gpar(fontsize = (90/length(rownames(clust_df))+5)),
        column_names_gp = gpar(fontsize = (90/length(rownames(clust_df))+5) + 2),
        column_names_max_height = unit(8, "cm"),
        #split = max(mycl), gap = unit(2, "mm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", data.matrix(pval_df)[i, j]), x, y, gp = gpar(fontsize = (80/length(rownames(clust_df))+3)))
        }
)
dev.off()
