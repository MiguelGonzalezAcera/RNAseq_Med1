library(DESeq2)
library(optparse)
library(calibrate)

# Create options
option_list = list(
  make_option("--res", type="character",
              help="R object with the result of the analysis"),
  make_option("--out_plot", type="character",
              help="file that contains the volcano."),
  make_option("--organism", type="character", default="mouse",
              help="Organism for the genenames")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load Rfunctions
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select database
database <- select.organism(opt$organism)

# Read data table
load(opt$res)
resdf <- data.frame(res)[complete.cases(data.frame(res)),]

# add genenames to the table
resdf$Genes <- as.character(mapIds(database, as.character(rownames(resdf)),
                                   'SYMBOL', 'ENSEMBL'))

png(file=opt$out_plot, width = 3000, height = 3000, res = 600)
# Create scatterplot for the volcano
with(resdf, plot(log2FoldChange, -log10(pvalue), pch=20))

# Color the points using thresholds
#<TO_DO>: replace the thresholds by parameters
with(subset(resdf, padj<.001), points(log2FoldChange, -log10(pvalue),
                                      pch=20, col="red"))
with(subset(resdf, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue),
                                                  pch=20, col="orange"))
with(subset(resdf, padj<.001 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue),
                                                              pch=20, col="green"))

resdf$Genes_filt <- rapply(as.list(resdf$Genes),function(x) ifelse(startsWith(x,"Gsdmc"),x,""), how = "replace")
resdf$Genes_filt[is.na(resdf$Genes_filt)] <- ""

# Annotate significant points with genenames
with(subset(resdf, padj<.001 & abs(log2FoldChange)>1), textxy(log2FoldChange*0.9, -log10(pvalue)*1.02,
                                                              labs=Genes, cex = .5))

dev.off()

# Save environment
save.image(file=gsub(".png",".RData",opt$out_plot, fixed = TRUE))

# Save versions
get_versions(gsub(".png","_versions.tsv",opt$out_plot, fixed = TRUE))