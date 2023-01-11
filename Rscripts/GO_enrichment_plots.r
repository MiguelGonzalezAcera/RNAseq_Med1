# GO enrichment plots

suppressPackageStartupMessages(library(GOstats))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option("--out_tab", type="character",
              help="Table object with the result of the analysis."),
  make_option("--organism", type="character", default= "human",
              help="Organism analyzed. Available = human, mouse. Default = Human"),
  make_option("--geneids", type="character",
              help="List of genes used in the GO enrichment analysis")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load Rfunctions
source("Rscripts/Rfunctions.R")

# Load R object
load(opt$out_tab)

# Load geneids
load(opt$geneids)

# Select organism
database <- select.organism(opt$organism)

# Obtain plots
# barplot
png(file=sprintf("%s_barplot.png",
                 gsub(".rda","",opt$out_tab, fixed=TRUE)), width = 8000, height = 6000, res = 600)
barplot(x, showCategory=16)
dev.off()

# dotplot
png(file=sprintf("%s_dotplot.png",
                 gsub(".rda","",opt$out_tab, fixed=TRUE)), width = 8000, height = 6000, res = 600)
dotplot(x, showCategory=16)
dev.off()

# Enrichment map
png(file=sprintf("%s_emap.png",
                 gsub(".rda","", opt$out_tab, fixed=TRUE)), width = 8000, height = 6000, res = 600)
emapplot(x)
dev.off()

# Gene-Concept Network
# plot linkages of genes and enriched concepts (e.g. GO categories, KEGG pathways)
png(file=sprintf("%s_cnet.png",
                 gsub(".rda","", opt$out_tab, fixed=TRUE)), width = 8000, height = 6000, res = 600)
cnetplot(x, categorySize="pvalue", foldChange = entrezgeneids)
dev.off()
