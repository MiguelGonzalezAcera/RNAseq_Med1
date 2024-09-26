# GO enrichment

suppressPackageStartupMessages(library(GOstats))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(enrichplot))

option_list <- list(
  make_option("--out_tab", type = "character",
              help = "Table with the result of the analysis."),
  make_option("--obj", type = "character",
              help = "A genelist of entrez genes with the gene group
              to use for the enrichment"),
  make_option("--universe", type = "character",
              help = "The r object containing the complete genelist of entrez genes in the assay. It is saved in the differential expression analysis"),
  make_option("--organism", type = "character", default = "human",
              help = "Organism analyzed. Available = human, mouse. Default = Human"),
  make_option("--genelist", type = "character", default = "",
              help = "List of genes to analyze. Optional.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load Rfunctions
source("Rscripts/Rfunctions.R")

# Load R object
load(opt$obj)

# Filter the object by fold change
res <- res[which((res$log2FoldChange < -1 | res$log2FoldChange > 1) & (res$padj < 0.05)), ]

# Load the universe
load(opt$universe)

# Select organism
database <- select.organism(opt$organism)

# Obtain genelist
if (opt$genelist == "") {
  entrezgeneids <- tryCatch(
    {
      as.character(mapIds(database, as.character(rownames(res)), "ENTREZID", "ENSEMBL"))
    },
    error = function(cond) {
      message(paste("Error with gene set:", rownames(res)))
      message("Error message:")
      message(cond)
      message()
      quit()
    },
    warning = function(cond) {
      message("Database had a warning:")
      message(cond)
      message()
    },
    finally = {
    }
  )
} else {
  genes <- readLines(opt$genelist)

  # Transform the ensembl names into gene symbol. NOTE that the name of the variable must change.
  entrezgeneids <- tryCatch(
    {
      entrezgeneids <- as.character(mapIds(database, as.character(genes), "ENTREZID", "ENSEMBL"))
    },
    error = function(cond) {
      message(paste("Error with gene set:", genes))
      message("Error message:")
      message(cond)
      quit()
    },
    warning = function(cond) {
      message("Database had a warning:")
      message(cond)
    },
    finally = {
    }
  )
}

# Obtain universe ids
universeids <- unique(as.character(mapIds(database, as.character(rownames(counts(dss))), "ENTREZID", "ENSEMBL")))

# Define ontologies
ontologies <- c("BP", "MF", "CC")

# Loop for each ontology
for (ont in ontologies) {
  # Do the GSEA
  hgCutoff <- 0.1
  x <- enrichGO(entrezgeneids, database, ont=ont, pvalueCutoff = hgCutoff, readable = TRUE,
                pAdjustMethod = "BH", universe = universeids, minGSSize = 1, maxGSSize = 1000)

  # Transform the result into a data frame
  GOtable <- as.data.frame(x)

  # Save the data frame
  write.table(GOtable, file = sprintf("%s_%s.tsv",
              gsub(".tsv", "", opt$out_tab, fixed = TRUE), ont),
              sep = "\t", row.names = FALSE)

  # Obtain plots
  # barplot
  png(file = sprintf("%s_%s_barplot.png",
                  gsub(".tsv", "", opt$out_tab, fixed = TRUE), ont),
                  width = 8000, height = 6000, res = 600)
  # when plotting inside a for loop, you have to explicitely print the plot
  bplot <- tryCatch(
    {
      barplot(x, showCategory = 16)
    },
    error = function(cond) {
      message(paste("Error doing barplot of: ", ont))
      message("Error message:")
      message(cond)
      message()
      quit()
    },
    warning = function(cond) {
      message("Database had a warning:")
      message(cond)
      message()
    },
    finally = {
    }
  )
  print(bplot)
  dev.off()

  # dotplot
  png(file = sprintf("%s_%s_dotplot.png",
                  gsub(".tsv", "", opt$out_tab, fixed = TRUE), ont),
                  width = 8000, height = 6000, res = 600)
  dplot <- tryCatch(
    {
      dotplot(x, showCategory = 16)
    },
    error = function(cond) {
      message(paste("Error doing dotplot of: ", ont))
      message("Error message:")
      message(cond)
      message()
      quit()
    },
    warning = function(cond) {
      message("Database had a warning:")
      message(cond)
      message()
    },
    finally = {
    }
  )
  print(dplot)
  dev.off()

  # Enrichment map
  png(file = sprintf("%s_%s_emap.png",
                  gsub(".tsv", "", opt$out_tab, fixed = TRUE), ont),
                  width = 8000, height = 6000, res = 600)
  emplot <- tryCatch(
    {
      x2 <- pairwise_termsim(x)
      emapplot(x2)
    },
    error = function(cond) {
      message(paste("Error doing barplot of: ", ont))
      message("Error message:")
      message(cond)
      message()
      quit()
    },
    warning = function(cond) {
      message("Database had a warning:")
      message(cond)
      message()
    },
    finally = {
    }
  )
  print(emplot)
  dev.off()

  # Gene-Concept Network
  # plot linkages of genes and enriched concepts
  # (e.g. GO categories, KEGG pathways)
  png(file = sprintf("%s_%s_cnet.png",
                  gsub(".tsv", "", opt$out_tab, fixed = TRUE), ont),
                  width = 8000, height = 6000, res = 600)
  # Cnet misses sometimes, so we can put it in a error structure
  cnet <- tryCatch(
    {
      cnetplot(x, categorySize = "pvalue", foldChange = entrezgeneids)
    },
    error = function(cond) {
      message(paste("Error doing cnet of: ", ont))
      message("Error message:")
      message(cond)
      message()
      quit()
    },
    warning = function(cond) {
      message("Database had a warning:")
      message(cond)
      message()
    },
    finally = {
    }
  )
  print(cnet)
  dev.off()

  # Save the objects
  save(x, file = sprintf("%s_%s.Rda",
                       gsub(".tsv", "", opt$out_tab, fixed = TRUE), ont))

}

# Save environment
save.image(file = gsub(".tsv", ".RData", opt$out_tab, fixed = TRUE))

# Save versions
get_versions(gsub(".tsv", "_versions.tsv", opt$out_tab, fixed = TRUE))
