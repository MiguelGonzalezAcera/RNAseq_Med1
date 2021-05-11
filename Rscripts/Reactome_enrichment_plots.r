library(clusterProfiler)
library(DESeq2)
library(optparse)
library(ReactomePA)

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

# Transform result to table
Reactometable <- as.data.frame(x)

if (length(rownames(Reactometable)) >= 10){
  top_pathways <- rownames(Reactometable)[1:10]
} else {
  top_pathways <- rownames(Reactometable)[1:length(rownames(Reactometable))]
}

# Show routes
for (pway in top_pathways) {
  print(gsub(".rda",sprintf("_%s_pway.png",pway),opt$out_tab, fixed = TRUE))
  png(file=gsub(".rda",sprintf("_%s_pway.png",pway),opt$out_tab, fixed = TRUE), width = 8000, height = 6000, res = 600)
  tryCatch(
    {
      plot_res <- viewPathway(Reactometable[pway,"Description"], organism = opt$organism, readable=TRUE, foldChange=geneList)
    },
    error=function(cond) {
      message(paste("There is an issue with this pathway:", Reactometable[pway,"Description"]))
      message("Here's the original error message:")
      message(cond)
    }
  )
  print(plot_res)
  dev.off()
}
