# Resctome enrichment

library(clusterProfiler)
library(pathview)
library(DESeq2)
library(optparse)
library(org.Mm.eg.db)

option_list = list(
  make_option("--out_tab", type="character",
              help="Table with the result of the analysis."),
  make_option("--genelist", type="character", default = c(),
              help="A genelist of entrez genes with the gene group
              to use for the enrichment. Default = None"),
  make_option("--organism", type="character", default= "human",
              help="Organism analyzed. Available = human, mouse. Default = Human")
)

opt_parser = OprtionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R object
load()

# Load R scripts
source("D:/Documentos/LATESIS/Scripts/Rfunctions.R")

# Select organism
database <- select.organism("mouse")

# Generate named list of FC
geneList <- res$log2FoldChange
names(geneList) <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(res)),
                                       'ENTREZID', 'ENSEMBL'))

# Obtain names
entrezgeneids <- names(geneList)
genes <- scan("/DATA/DSS_rec_evolution/genelist_pg.txt", character(), quote="")
entrezgeneids <- as.character(mapIds(org.Mm.eg.db, genes, 'ENTREZID', 'SYMBOL'))

# Do the GSEA
#<TO_DO>: Change parameters cutoff and organism.
hgCutoff <- 0.05
x <- enrichKEGG(entrezgeneids, organism="mmu", pvalueCutoff=hgCutoff, pAdjustMethod="BH",
                qvalueCutoff=0.1)

# Transform the result into a data frame
KEGGtable <- as.data.frame(x)

# Save the data frame
write.table(KEGGtable, file=opt$out_tab, sep="\t", row.names = FALSE)

# Obtain plots
# barplot
png(file=gsub(".tsv","_barplot.png",opt$out_plot, fixed = TRUE))
barplot(x, showCategory=16)
dev.off()

# Enrichment map
png(file=gsub(".tsv","_emap.png",opt$out_plot, fixed = TRUE))
emapplot(x)
dev.off()

# Show the main pathways
top_pathways <- rownames(KEGGtable)[1:10]

for (pway in top_pathways) {
  pathway <- pathview(gene.data=geneList, pathway.id = pway, species = "mmu", kegg.dir = "/DATA/DSS_rec_evolution/", out.suffix = "pathview_02inf_hi")
  wd <- paste(c(getwd(),paste(c(pway,'pathview_02inf_hi.png'), collapse = '.')), collapse = '/')
  
  file.copy(wd, "/DATA/DSS_rec_evolution/")
  file.remove(wd)
}

# Save environment
save.image(file=gsub(".tsv",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".tsv","_versions.tsv",opt$obj_out, fixed = TRUE))