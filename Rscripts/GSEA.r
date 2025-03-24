suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(rjson))

option_list <- list(
  make_option("--pathways", type = "character",
              help = "Path for the selected pathways, in json format"),
  make_option("--in_obj", type = "character",
              help = "Robject with the DE analysis. Rda extension"),
  make_option("--gseaplot", type = "character",
              help = "Path for the GSEAplot. png extension"),
  make_option("--dims", type = "character", default = "2000,2000",
              help = "Dimensions of the plot in pixels. Default = 2000,2000"),
  make_option("--organism", type = "character", default = "mouse",
              help = "Organism analyzed. STR. Available = human, mouse. Default = Mouse")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load the R object with the result
load(opt$in_obj)

# Load R scripts
source("Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Transform the object into a named genelist (all of the genes)
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)

# Sort the genelist by value
geneList <- geneList[order(geneList)]

# read the gene lists from the json
pways_data <- fromJSON(file = opt$pathways)

# Do the gene set ernichment analysis
z <- fgsea(pathways = pways_data, stats = geneList, minSize = 5, maxSize = 500)

# fix the gene cells for the table
z$leadingEdge<- unlist(lapply(z$leadingEdge, function(x) {paste(x, collapse = ",")}))

# Save enrichment table
write.table(as.data.frame(z), file = gsub(".svg", "_GSEA.tsv", opt$gseaplot, fixed = TRUE), sep = "\t", row.names = FALSE)

# Make the plots
if (nrow(as.data.frame(z)) != 0) {
  for (pway in as.data.frame(z)[['pathway']][1:min(10:length(as.data.frame(z)[['pathway']]))]) {
    # Make the original plot
    p <- plotEnrichment(pways_data[[pway]], geneList) + labs(title = pway)
    print(pway)
    # Replace the line color with a chosen one
    if (as.data.frame(z)[as.data.frame(z)["pathway"] == pway, ][["ES"]] <= 0) {
      p$layers[[1]]$aes_params$colour <- "#7f00ff"
    } else {
      p$layers[[1]]$aes_params$colour <- "#ff8000"
    }

    # Save the plot
    png(
      file = gsub(".svg", sprintf("_%s_GSEA.png", pway), opt$gseaplot, fixed = TRUE),
      width = as.integer(strsplit(opt$dims, ",")[[1]][1]),
    height = as.integer(strsplit(opt$dims, ",")[[1]][2]),
      res = 300)
    print(p)
    dev.off()

    svg(
      file = gsub(".svg", sprintf("_%s_GSEA.svg", pway), opt$gseaplot, fixed = TRUE),
      width = as.integer(strsplit(opt$dims, ",")[[1]][1]),
    height = as.integer(strsplit(opt$dims, ",")[[1]][2])
      )
    print(p)
    dev.off()
  }
}
