# Resctome enrichment

library(clusterProfiler)
library(pathview)
library(DESeq2)
library(optparse)

option_list = list(
  make_option("--out_tab", type="character",
              help="Table with the result of the analysis. TSV extension"),
  make_option("--in_obj", type="character",
              help="Robject with the DE analysis. Rda extension"),
  make_option("--id", type="character",
              help="Sample or group name for the output. STR"),
  make_option("--genelist", type="character", default = c(),
              help="A genelist of entrez genes with the gene group
              to use for the enrichment. TXT extension. Default = None"),
  make_option("--organism", type="character", default= "mouse",
              help="Organism analyzed. STR. Available = human, mouse. Default = Mouse")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R object
load(opt$in_obj)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Generate named list of FC
geneList <- res$log2FoldChange
names(geneList) <- as.character(mapIds(database, as.character(rownames(res)),
                                       'ENTREZID', 'ENSEMBL'))

# Obtain names
genes <- scan(opt$genelist, character(), quote="")
entrezgeneids <- as.character(mapIds(database, genes, 'ENTREZID', 'SYMBOL'))

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
png(file=gsub(".tsv","_barplot.png",opt$out_tab, fixed = TRUE))
barplot(x, showCategory=16)
dev.off()

# Enrichment map
png(file=gsub(".tsv","_emap.png",opt$out_tab, fixed = TRUE))
emapplot(x)
dev.off()

# Get final path
pathlist <- unlist(strsplit(opt$out_tab,"/"))
path = paste(pathlist[1:length(pathlist)-1], collapse = "/")

# Show the main pathways
top_pathways <- rownames(KEGGtable)[1:10]

for (pway in top_pathways) {
  pathway <- pathview(gene.data=geneList, pathway.id = pway, species = "mmu", kegg.dir = "/DATA/tmp/", out.suffix = opt$id)
  wd <- paste(c(getwd(), paste(c(pway, opt$id, "png"), collapse = '.')), collapse = '/')
  print(wd)

  file.copy(wd, path)
  file.remove(wd)
}

# Save environment
save.image(file=gsub(".tsv",".RData",opt$out_tab, fixed = TRUE))

# Save versions
get_versions(gsub(".tsv","_versions.tsv",opt$out_tab, fixed = TRUE))
