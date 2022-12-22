suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(calibrate))

# Create options
option_list = list(
  make_option("--res", type="character",
              help="R object with the result of the analysis"),
  make_option("--out_plot", type="character",
              help="file that contains the volcano."),
  make_option("--organism", type="character", default="mouse",
              help="Organism for the genenames"),
  make_option("--genelist", type="character", default="",
              help="List of genes to be included in the plot, in txt format. Ensembl IDs only"),
  make_option("--labels", type="character", default=FALSE,
              help="Add labels of the genes")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load Rfunctions
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select database
database <- select.organism(opt$organism)

# Read data table
load(opt$res)
#load("/VAULT/Thesis_proj/detables/Mouse_models_cDSSdc_Cerldc.Rda")
resdf <- data.frame(res)[complete.cases(data.frame(res)),]

if (opt$genelist != "") {
  #Load list of genes
  genes = readLines(opt$genelist)
  
  resdf$Genename <- rownames(resdf)
  
  # Get rows in the list of genes
  resdf <- resdf[resdf$Genename %in% genes, ,drop=FALSE]
  resdf$Genename = NULL
}

# add genenames to the table
resdf$Genes <- as.character(mapIds(database, as.character(rownames(resdf)),
                                   'SYMBOL', 'ENSEMBL'))

# Make bool column for coloring
resdf$Colour <- (resdf$log2FoldChange < -1 | resdf$log2FoldChange > 1) & resdf$pvalue < 0.05

#Plot with ggplot
ggplot(resdf, aes(log2FoldChange, -log10(pvalue), color = Colour)) + 
  geom_rect(aes(xmin=-Inf, xmax=-1, ymax=Inf, ymin=-log10(0.05)), fill="blue", alpha=0.5) + 
  geom_point()




png(file=opt$out_plot, width = 3000, height = 3000, res = 600)
#png(file="/DATA/Thesis_proj/Requests_n_stuff/20200217_Le_request/cDSS_volcano_plot.png", width = 3000, height = 3000, res = 600)
# Create scatterplot for the volcano
with(resdf, plot(log2FoldChange, -log10(pvalue), pch=20, cex = 0.6))

rect(min(resdf$log2FoldChange)-0.5, -log10(min(resdf$pvalue))*-1.5, -1, -log10(min(resdf$pvalue))*1.5, col=rgb(0,0,1.0,alpha=0.5))

# Color the points using thresholds
#<TO_DO>: replace the thresholds by parameters
with(subset(resdf, padj<.05), points(log2FoldChange, -log10(pvalue),
                                     pch=20, col="red", cex = 0.6))
with(subset(resdf, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue),
                                                  pch=20, col="orange", cex = 0.6))
with(subset(resdf, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue),
                                                             pch=20, col="green", cex = 0.6))

if (opt$genelist != "" & opt$labels == TRUE){
  genelist_ens <- scan(opt$genelist, character(), quote="")
  genelist <- as.character(mapIds(database, as.character(genelist_ens),
                                  'SYMBOL', 'ENSEMBL'))
  
  with(subset(resdf, padj<.05 & abs(log2FoldChange)>1 & Genes %in% genelist), points(log2FoldChange, -log10(pvalue),pch=20, col="cyan", cex = 1.2))
  resdf$Genes_filt <- rapply(as.list(resdf$Genes),function(x) ifelse(x %in% genelist,x,""), how = "replace")
  resdf$Genes_filt[is.na(resdf$Genes_filt)] <- ""
  with(subset(resdf, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange*0.9, -log10(pvalue)*1, labs=Genes_filt, cex = .6)) 
  
}

# Annotate significant points with genenames
# with(subset(resdf, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange*0.9, -log10(pvalue)*1.02,
#                                                               labs=Genes, cex = .5))

dev.off()

# Save environment
save.image(file=gsub(".png",".RData",opt$out_plot, fixed = TRUE))

# Save versions
get_versions(gsub(".png","_versions.tsv",opt$out_plot, fixed = TRUE))