library(DESeq2)
library(optparse)
library(calibrate)

# Create options
option_list = list(
  make_option("--res", type="character",
              help="R object with the result of the analysis"),
  make_option("--out_plot", type="character",
              help="file that contains the volcano.")
)

opt_parser = OprtionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load Rfunctions
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Read data table
load("/DATA/IECs_rhoa_Rocio/Tr1KO.rda")
resdf <- data.frame(res)[complete.cases(data.frame(res)),]

# add genenames to the table
resdf$Genes <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(resdf)),
                                   'SYMBOL', 'ENSEMBL'))

png(file="/DATA/IECs_rhoa_Rocio/Tr1KO_volcano.png", width = 8000, height = 6000, res = 600)
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