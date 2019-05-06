# Resctome enrichment

library(clusterProfiler)
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
database <- select.organism(opt$organism)

# Read genelist
entrezgeneids <- scan(opt$genelist, character(), quote="")

# Do the GSEA
#<TO_DO>: Change parameters cutoff and organism.
hgCutoff <- 0.05
x <- enrichPathway(gene = entrezgeneids, pvalueCutoff = hgCutoff, readable = T,
              organism = "mouse")

# Transform the result into a data frame
Reactometable <- as.data.frame(x)

# Save the data frame
write.table(Reactometable, file=opt$out_tab, sep="\t", row.names = FALSE)

# Obtain plots
# barplot
png(file=gsub(".tsv","_barplot.png",opt$out_plot, fixed = TRUE))
barplot(x, showCategory=16)
dev.off()

# Enrichment map
png(file=gsub(".tsv","_emap.png",opt$out_plot, fixed = TRUE))
emapplot(x)
dev.off()

# Gene-Concept Network
# plot linkages of genes and enriched concepts (e.g. GO categories, KEGG pathways)
png(file=gsub(".tsv","_cnet.png",opt$out_plot, fixed = TRUE))
cnetplot(x, categorySize="pvalue", foldChange = entrezgeneids)
dev.off()

# Save environment
save.image(file=gsub(".tsv",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".tsv","_versions.tsv",opt$obj_out, fixed = TRUE))