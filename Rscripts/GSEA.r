suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(optparse))


option_list = list(
  make_option("--genegroup", type="character",
              help="Path for the groups of genes, in particular entrez table format"),
  make_option("--in_obj", type="character",
              help="Robject with the DE analysis. Rda extension"),
  make_option("--gseaplot", type="character",
              help="Path for the GSEAplot. png extension"),
  make_option("--organism", type="character", default= "mouse",
              help="Organism analyzed. STR. Available = human, mouse. Default = Mouse")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load the R object with the result
load(opt$in_obj)

# Load R scripts
source("Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Transform the object into a named genelist (all of the genes)
geneList <- res$log2FoldChange
names(geneList) <- as.character(mapIds(database, as.character(rownames(res)),
                                       'ENTREZID', 'ENSEMBL'))

# Read the custom groups (column 1: name of the group, and column 2, gene id in entrez)
groups = read.table(opt$genegroup, fileEncoding = "UTF8")

# Do the gene set enrichment analysis
z = GSEA(sort(geneList,decreasing=T), TERM2GENE = groups, pvalueCutoff = 1, minGSSize = 5)

# Plot the result of the GSA
if (length(rownames(as.data.frame(z))) >= 10){
  len = 10
} else {
  len = length(rownames(as.data.frame(z)))
}

# Save enrichment table
write.table(as.data.frame(z), file=gsub(".png",".tsv",opt$gseaplot, fixed = TRUE),sep="\t",row.names = FALSE)

png(file=opt$gseaplot, width = 2000, height = 2000, res = 200)
gseaplot2(z, geneSetID = 1:len, pvalue_table = FALSE, base_size = 24)
dev.off()

# Barplot
# png(file=gsub(".png","_barplot.png",opt$gseaplot, fixed = TRUE), width = 8000, height = 8000, res = 600)
# barplot(z)
# dev.off()

# Enrichment map
# png(file=gsub(".png","_emapplot.png",opt$gseaplot, fixed = TRUE), width = 8000, height = 8000, res = 600)
# emapplot(z)
# dev.off()
