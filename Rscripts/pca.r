# Run the counts form the bamfile

library(scatterplot3d)
library(optparse)
library(ggplot2)

option_list <- list(
  make_option("--counts", type = "character",
              help = "Table that contains the counts"),
  make_option("--design", type = "character",
              help = "list of bam files of the assay, with metadata."),
  make_option("--out_dir", type = "character",
              help = "Folder to contain the plots of the PCA.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load R scripts
source("Rscripts/Rfunctions.R")

# Read the table with the metadata
sampleTableSingle <- read.table(opt$design, fileEncoding = "UTF8")

# Read the table containing the counts
Counts_tab <- read.table(opt$counts, fileEncoding = "UTF8", header = TRUE)
row.names(Counts_tab) <- Counts_tab$Geneid
Counts_tab$Geneid <- NULL
Counts_tab <- Counts_tab[, row.names(sampleTableSingle)]

# do principal components analyses
# as written in monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/
# RNAseq_DE_analysis_with_R.html

# Select the 1000 most highly expressed genes.
select <- order(rowMeans(Counts_tab), decreasing = TRUE)[1:1000]
select <- order(rowMeans(Counts_tab), decreasing = TRUE)
highexprgenes_counts <- Counts_tab[select, ]

# Get all treatment columns into one:
#<TO_DO> Adjust the format of the table (first names, then treatments)
lenCol <- length(colnames(sampleTableSingle))

Treatment <- do.call(paste, c(sampleTableSingle[colnames(
  sampleTableSingle)[(1:lenCol)]], sep = "_"))

# Annotate the columns with the treatment.
sampleids <- paste(colnames(highexprgenes_counts), factor(Treatment), sep = " - ")
colnames(highexprgenes_counts) <- factor(Treatment)

# Get the genes as columns
data_for_PCA <- t(highexprgenes_counts)

# Matrix of dissimilarities. Return also the eigenvalues.
mds <- cmdscale(dist(data_for_PCA), k = 3, eig = TRUE)

# Transform the eigenvalues into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)

# Plot the eigenvalues
# NOTE: This plot should not be activated in the final script, since is only for checking.
#<TO_DO>: Subsitute for an if 2d/3d using a threshold (ex. 10%)
# barplot(eig_pc, las=1, xlab="Dimensions", ylab="prop of explained value", y.axis=NULL)

# Create 2d plotting function
pca2d <- function(tab, d1, d2, perc, out) {
  # Takes a table and its columns and does a plot
  # Opt 1: Regular graph engine
  # png(file=paste(out,sprintf("pca_2d_d%s_d%s.png",d1,d2), sep = "/"))
  # plot(tab[,d1], tab[,d2], type="n",
  #      xlab=sprintf("Dimension %s (%s %%)",d1,perc[d1]),
  #      ylab=sprintf("Dimension %s (%s %%)",d2,perc[d2]), main="")
  # text(tab[,d1], tab[,d2], rownames(tab), cex=0.8)
  # dev.off()
  
  # Opt 2: ggplot
  tab_plot <- ggplot(tab, aes(x = get(sprintf("V%s", d1)), y = get(sprintf("V%s", d2)), color = Treatment, size = 20, label = SampleID)) + 
    geom_point(size = 2) +
    xlab(sprintf("Dimension %s (%s %%)",d1,perc[d1])) + ylab(sprintf("Dimension %s (%s %%)",d2,perc[d2])) + 
    #geom_text(hjust = 0, vjust = 0, size = 2) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    )
  ggsave(file = paste(out,sprintf("pca_2d_d%s_d%s.png", d1, d2), sep = "/"), tab_plot, device = "png", bg = "white", width = 14, height = 10, units = 'cm', dpi=600)
  ggsave(file = paste(out,sprintf("pca_2d_d%s_d%s_transparent.png", d1, d2), sep = "/"), tab_plot, device = "png", bg = "transparent", width = 14, height = 10, units = 'cm', dpi=600)
  ggsave(file = paste(out,sprintf("pca_2d_d%s_d%s.svg", d1, d2), sep = "/"), tab_plot, device = "svg", bg = "white", width = 14, height = 10, units = 'cm', dpi=600)
  ggsave(file = paste(out,sprintf("pca_2d_d%s_d%s_transparent.svg", d1, d2), sep = "/"), tab_plot, device = "svg", bg = "transparent", width = 14, height = 10, units = 'cm', dpi=600)
}

# Translate this into parameter
dimthr <- 10

# Do main 2d plot
mds <- cmdscale(dist(data_for_PCA), k = 3)

# Transform the data to a data frame and keep the row names
mdsdf <- as.data.frame(mds)
mdsdf$Treatment <- rownames(mds)
mdsdf$SampleID <- sampleids

# Do the 2d graph
pca2d(mdsdf, 1, 2, eig_pc, opt$out_dir)

# If 3rd eigenvalue is higher than 10, produce the dimension combination and the gif
if (eig_pc[3] >= dimthr) {
  # Do the rest of the 2d graphs
  pca2d(mdsdf, 2, 3, eig_pc, opt$out_dir)
  pca2d(mdsdf, 1, 3, eig_pc, opt$out_dir)

  # Do the 3d graph
  # Shapes of the points by treatment
  shapes <- rep_len(c(1:14), length.out = length(levels(factor(mdsdf$Treatment))))
  shape <- shapes[as.numeric(factor(mdsdf$Treatment))]

  # Colors of the points by treatment
  color <- rainbow(length(levels(factor(mdsdf$Treatment))))
  colore <- color[as.numeric(factor(mdsdf$Treatment))]

  # Produce many images in order to do a gif for web visualization
  for (ang in c(seq(1, 180))){
    # Create the plot
    png(file = paste(opt$out_dir, sprintf("pca_3d_%03d.png", ang), sep = "/"))
    #png(file=sprintf("/DATA/Thesis_proj/SEPIA_The_Paper/FullPCA/pca_FC_full_3d_%03d.png",ang), width = 1200, height = 1200, res = 150)
    s3d <- scatterplot3d(mdsdf[, 1], mdsdf[, 2], mdsdf[, 3], pch = shape, color = colore,
                         cex.symbols = 1.5, angle = ang,
                         #label.tick.marks = TRUE,
                         xlab = sprintf("Dimension 1 (%s %%)", eig_pc[1]),
                         ylab = sprintf("Dimension 2 (%s %%)", eig_pc[2]),
                         zlab = sprintf("Dimension 3 (%s %%)", eig_pc[3]))
    # Add the legend
    legend("bottomright", legend = levels(factor(mdsdf$Treatment)),
           col = color, pch = shapes, inset = -0.1, xpd = TRUE, pt.cex = 1, cex = 0.5)
    dev.off()
  }
}

# Save the table with the eigenvalues for alternative displays
write.table(mdsdf, file = paste(opt$out_dir, "pca_eigenvalues.tsv", sep = "/"), sep = "\t")

# Save environment
save.image(file = paste(opt$out_dir, "pca.RData", sep = "/"))

# Save versions
get_versions(paste(opt$out_dir, "pca_versions.tsv", sep = "/"))
