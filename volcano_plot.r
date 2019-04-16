library(DESeq2)
library(optparse)
library(calibrate)

# Create options
option_list = list(
  make_option("--tab", type="character",
              help="Table with the result of the analysis"),
  make_option("--out_plot", type="character",
              help="file that contains the volcano.")
)

opt_parser = OprtionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load Rfunctions
source("D:/Documentos/LATESIS/Scripts/Rfunctions.R")

# Read data table
resdf = read.table(opt$tab, fileEncoding = "UTF8")

png(file=opt$out_plot)
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

# Annotate significant points with genenames
with(subset(resdf, padj<.001 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue),
                                                              labs=Genes))

dev.off()

# Save environment
save.image(file=gsub(".png",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".png","_versions.tsv",opt$obj_out, fixed = TRUE))