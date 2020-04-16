# Resctome enrichment

library(clusterProfiler)
library(pathview)
library(DESeq2)
library(optparse)
library(enrichplot)

option_list = list(
  make_option("--out_tab", type="character",
              help="Table with the result of the analysis. TSV extension"),
  make_option("--in_obj", type="character",
              help="Robject with the DE analysis. Rda extension"),
  make_option("--id", type="character",
              help="Sample or group name for the output. STR"),
  make_option("--genelist", type="character", default = "",
              help="A genelist of genes with the gene group
              to use for the enrichment. TXT extension. Default = None"),
  make_option("--organism", type="character", default= "mouse",
              help="Organism analyzed. STR. Available = human, mouse. Default = Mouse")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R object
load(opt$in_obj)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)
if (opt$organism == "human") {
  org_db = "hsa"
} else if (opt$organism == "mouse") {
  org_db = "mmu"
}

# Generate named list of FC
geneList <- res$log2FoldChange
names(geneList) <- as.character(mapIds(database, as.character(rownames(res)),
                                       'ENTREZID', 'ENSEMBL'))

# Obtain genelist
if (opt$genelist == ""){
  entrezgeneids <- (as.character(mapIds(database, as.character(rownames(res[which((res$log2FoldChange < -1 | res$log2FoldChange > 1) & (res$padj < 0.05)),])), 'ENTREZID', 'ENSEMBL')))
} else {
  genes <- scan(opt$genelist, character(), quote="")

  # Transform the ensembl names into gene symbol. NOTE that the name of the variable must change.
  entrezgeneids <- as.character(mapIds(database, as.character(genes), 'ENTREZID', 'ENSEMBL'))
}

# Do the GSEA
#<TO_DO>: Change parameters cutoff and organism.
hgCutoff <- 0.1
x <- enrichKEGG(entrezgeneids, organism=org_db, pvalueCutoff=hgCutoff, pAdjustMethod="BH",
                qvalueCutoff=0.1)

# Transform the result into a data frame
KEGGtable <- as.data.frame(x)

# Save the data frame
write.table(KEGGtable, file=opt$out_tab, sep="\t", row.names = FALSE)

# Obtain plots
# barplot
png(file=gsub(".tsv","_barplot.png",opt$out_tab, fixed = TRUE), width = 8000, height = 6000, res = 600)
barplot(x, showCategory=16)
dev.off()

# Enrichment map
png(file=gsub(".tsv","_emap.png",opt$out_tab, fixed = TRUE), width = 8000, height = 6000, res = 600)
emapplot(x)
dev.off()

# Get final path
pathlist <- unlist(strsplit(opt$out_tab,"/"))
path = paste(pathlist[1:length(pathlist)-1], collapse = "/")

# Show the main pathways
if (length(rownames(KEGGtable)) >= 50){
  top_pathways <- rownames(KEGGtable)[1:50]
} else {
  top_pathways <- rownames(KEGGtable)[1:length(rownames(KEGGtable))]
  }

for (pway in top_pathways) {
  pathway <- pathview(gene.data=geneList, pathway.id = pway, species = org_db, kegg.dir = "/DATA/tmp/", out.suffix = opt$id, limit=list(gene=2, cpd=0.25))
  wd <- paste(c(getwd(), paste(c(pway, opt$id, "png"), collapse = '.')), collapse = '/')
  print(wd)

  file.copy(wd, path)
  file.remove(wd)
}

# Save environment
save.image(file=gsub(".tsv",".RData",opt$out_tab, fixed = TRUE))

# Save versions
get_versions(gsub(".tsv","_versions.tsv",opt$out_tab, fixed = TRUE))
