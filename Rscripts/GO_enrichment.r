# GO enrichment

library(GOstats)
library(clusterProfiler)
library(DESeq2)
library(optparse)
library(org.Mm.eg.db)

option_list = list(
  make_option("--out_tab", type="character",
              help="Table with the result of the analysis."),
  make_option("--obj", type="character",
              help="A genelist of entrez genes with the gene group
              to use for the enrichment"),
  make_option("--universe", type="character",
              help="The r object containing the complete genelist of entrez genes in the assay. It is saved in the differential expression analysis"),
  make_option("--organism", type="character", default= "human",
              help="Organism analyzed. Available = human, mouse. Default = Human"),
  make_option("--genelist", type="character", default='',
              help='List of genes to analyze. Optional.')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load Rfunctions
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Load R object
load(opt$obj)

# Load the universe
load(opt$universe)

# Select organism
database <- select.organism(opt$organism)

# Obtain genelist
if (opt$genelist == ""){
  entrezgeneids <- (as.character(mapIds(database, as.character(rownames(res)), 'ENTREZID', 'ENSEMBL')))
} else {
  genes = readLines(opt$genelist)
  
  # Transform the ensembl names into gene symbol. NOTE that the name of the variable must change.
  entrezgeneids <- as.character(mapIds(database, as.character(genes), 'ENTREZID', 'SYMBOL'))
}

# Obtain universe ids
universeids <- unique(as.character(mapIds(database, as.character(rownames(counts(dss))), 'ENTREZID', 'ENSEMBL')))

# Define ontologies
ontologies <- c("BP","MF","CC")

##OPTION 1
#<TO_DO>: Change parameters cutoff and database.
#<TO_DO>: Run this three times, one for each ontology (BP, MF, CC)
for (ont in ontologies) {
  # Run the hypergeometric test
  hgCutoff <- 0.05
  params <- new("GOHyperGParams", annotation=database$packageName, geneIds=entrezgeneids,
                universeGeneIds=universeids, ontology=ont, pvalueCutoff=hgCutoff,
                conditional=FALSE, testDirection="over")
  
  hg <- hyperGTest(params)
  
  # Get the pvalues
  hg.pv <- pvalues(hg)
  
  # Adjust the p values for comparisons
  hg.pv.fdr <- p.adjust(hg.pv,'fdr')
  
  # Get the significately enriched GO terms
  sigGO.ID <- names(hg.pv.fdr[hg.pv.fdr < hgCutoff])
  ## length(sigGO.ID)
  
  # Select the GO terms from the table
  df <- summary(hg)
  GOannot.table <- df[df[,1] %in% sigGO.ID,]
  
  # Save the enrichment table
  write.table(GOannot.table, file=sprintf("%s_%s.tsv", 
                                          gsub(".tsv","", opt$out_tab, fixed=TRUE),ont),
              sep="\t", row.names = FALSE)
}
  
##OPTION 2
# Also loop for each ontology
for (ont in ontologies) {
  # Do the GSEA
  ont = "CC"
  hgCutoff <- 0.05
  x <- enrichGO(entrezgeneids, database, ont=ont, pvalueCutoff = hgCutoff, readable = T,
                pAdjustMethod = "BH", universe = universeids)
  
  # Transform the result into a data frame
  GOtable <- as.data.frame(x)
  
  # Save the data frame
  write.table(GOtable, file=sprintf("%s_%s.tsv", 
                                    gsub(".tsv","", opt$out_tab, fixed=TRUE),ont),
              sep="\t", row.names = FALSE)
  
  # Obtain plots
  # barplot
  png(file=sprintf("%s_%s_barplot.png", 
                   gsub(".tsv","",opt$out_tab, fixed=TRUE),ont), width = 8000, height = 6000, res = 600)
  barplot(x, showCategory=16)
  dev.off()
  
  # Enrichment map
  png(file=sprintf("%s_%s_emap.png", 
                   gsub(".tsv","", opt$out_tab, fixed=TRUE),ont), width = 8000, height = 6000, res = 600)
  emapplot(x)
  dev.off()
  
  # Gene-Concept Network
  # plot linkages of genes and enriched concepts (e.g. GO categories, KEGG pathways)
  png(file=sprintf("%s_%s_cnet.png", 
                   gsub(".tsv","", opt$out_tab, fixed=TRUE),ont), width = 8000, height = 6000, res = 600)
  cnetplot(x, categorySize="pvalue", foldChange = entrezgeneids)
  dev.off()
}

# Save environment
save.image(file=gsub(".tsv",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".tsv","_versions.tsv",opt$obj_out, fixed = TRUE))