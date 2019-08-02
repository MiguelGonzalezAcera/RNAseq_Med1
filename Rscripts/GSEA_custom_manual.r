library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

# Load the R object with the result
load("/DATA/ILC_cell_lines_Stefan/Wt_ko.rda")

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism("mouse")

# Transform the object into a named genelist (all of the genes)
geneList <- res$log2FoldChange
names(geneList) <- as.character(mapIds(database, as.character(rownames(res)),
                                       'ENTREZID', 'ENSEMBL'))

# Read the custom groups (column 1: name of the group, and column 2, gene id in entrez)
groups = read.table("/DATA/ILC_cell_lines_Stefan/genelist_cytokine_prod/genegroups.txt", fileEncoding = "UTF8")

# Do the gene set enrichment analysis
z = GSEA(sort(geneList,decreasing=T), TERM2GENE = groups, pvalueCutoff = 1, minGSSize = 5)

# Plot the result of the GSA
if (length(rownames(as.data.frame(z))) >= 5){
  len = 5
} else {
  len = length(rownames(as.data.frame(z)))
}

# Save enrichment table
write.table(as.data.frame(z), file="/DATA/ILC_cell_lines_Stefan/genelist_interferon_gamma/GSEA_citokine_groups_Wt_ko.tsv",sep="\t")

png(file="/DATA/ILC_cell_lines_Stefan/genelist_interferon_gamma/GSEA_citokine_groups_Wt_ko.png", width = 8000, height = 8000, res = 600)
gseaplot2(z, geneSetID = 1:len, pvalue_table = TRUE)
dev.off()

# Barplot
png(file="", width = 8000, height = 8000, res = 600)
barplot(x, showCategory=)
dev.off()

# Enrichment map
png(file="", width = 8000, height = 8000, res = 600)
emapplot(x)
dev.off()
