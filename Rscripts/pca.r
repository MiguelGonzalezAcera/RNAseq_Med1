# Run the counts form the bamfile

library(scatterplot3d)
library(optparse)

option_list = list(
  make_option("--counts", type="character",
              help="Table that contains the counts"),
  make_option("--bamfiles", type="character",
              help="list of bam files of the assay, with metadata."),
  make_option("--out_plot", type="character",
              help="file that contains the plot of the PCA.")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load functions
source("/DATA/RNAseq_test/Scripts/Rfunctions.R")

# Read the table with the metadata
sampleTableSingle = read.table(opt$bamfiles, fileEncoding = "UTF8")

# Read the table containing the counts
Counts_tab = read.table(opt$counts, fileEncoding = "UTF8", header=TRUE)
row.names(Counts_tab) <- Counts_tab$Geneid
Counts_tab$Geneid = NULL
Counts_tab <- Counts_tab[,row.names(sampleTableSingle)]

# do principal components analyses
# as written in monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/
# RNAseq_DE_analysis_with_R.html

# Select the 1000 most highly expressed genes.
select = order(rowMeans(Counts_tab), decreasing = TRUE)[1:1000]
highexprgenes_counts <- Counts_tab[select,]

# Get all treatment columns into one:
#<TO_DO> Adjust the format of the table (first names, then treatments)
Treatment <- do.call(paste, c(sampleTableSingle[colnames(
  sampleTableSingle)[(1:length(colnames(sampleTableSingle)))]], sep = "_"))

# Annotate the columns with the treatment.
colnames(highexprgenes_counts) <- factor(Treatment)

# Get the genes as columns
data_for_PCA <- t(highexprgenes_counts)

# Matrix of dissimilarities. Return also the eigenvalues.
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)

# Transform the eigenvalues into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)

# Plot the eigenvalues
# NOTE: This plot should not be activated in the final script, since is only for checking.
#<TO_DO>: Subsitute for an if 2d/3d using a threshold (ex. 10%)
# barplot(eig_pc, las=1, xlab="Dimensions", ylab="prop of explained value", y.axis=NULL)

# Translate this into parameter
dimthr <- 10
if (eig_pc[3] < dimthr) {
  # Do 2d plot only
  mds <- cmdscale(dist(data_for_PCA), k=2)

  # Transform the data to a data frame and keep the row names
  mdsdf = as.data.frame(mds)
  mdsdf$Treatment <- rownames(mds)

  # Do the 2d graph
  pca2d(mdsdf,1,2,eig_pc,opt$out_plot)

} else if (eig_pc[3] >= dimthr) {
  # Do 2d plots for each dimension
  mds <- cmdscale(dist(data_for_PCA), k=3)

  # Transform the data to a data frame and keep the row names
  mdsdf = as.data.frame(mds)
  mdsdf$Treatment <- rownames(mds)

  # Do the 2d graphs
  pca2d(mdsdf,1,2,eig_pc,opt$out_plot)
  pca2d(mdsdf,2,3,eig_pc,opt$out_plot)
  pca2d(mdsdf,1,3,eig_pc,opt$out_plot)

  # Do the 3d graph
  # Shapes of the points by treatment
  shapes = c(1:length(levels(factor(mdsdf$Treatment)))) + 15
  shape <- shapes[as.numeric(factor(mdsdf$Treatment))]

  # Colors of the points by treatment
  color = rainbow(length(levels(factor(mdsdf$Treatment))))
  colore <- color[as.numeric(factor(mdsdf$Treatment))]

  # Produce many images in order to do a gif for web visualization
  for (ang in c(seq(1, 180))){
    # Create the plot
    png(file=gsub(".png",sprintf("_3d_%s.png",ang),opt$out_plot, fixed = TRUE))
    s3d <- scatterplot3d(mdsdf[,1], mdsdf[,2], mdsdf[,3], pch=shape, color=colore,
                         cex.symbols=1.5, angle=ang, label.tick.marks = FALSE,
                         xlab = sprintf("Dimension 1 (%s %%)",eig_pc[1]),
                         ylab = sprintf("Dimension 2 (%s %%)",eig_pc[2]),
                         zlab = sprintf("Dimension 3 (%s %%)",eig_pc[3]))
    # Add the legend
    legend("bottomright", legend=levels(factor(mdsdf$Treatment)),
           col = color,pch=shapes, inset = -0.1, xpd=TRUE)
    dev.off()
  }
}

# Save environment
save.image(file=gsub(".png",".RData",opt$out_plot, fixed = TRUE))

# Save versions
get_versions(gsub(".png","_versions.tsv",opt$out_plot, fixed = TRUE))
