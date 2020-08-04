# Clustering meddling

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(optparse))

# NOTE: this chunk of code shall be replaced by a genelist input.
# Might be a table from clusterProfiler, the raw genelist, or might be other source
option_list = list(
        make_option("--heatmap", type="character",
                    help="Heatmap of the coverage of the genes over the samples"),
        make_option("--counts", type="character",
                    help="An r object with the normalized counts. Produced in the DE script."),
        make_option("--genelist", type="character",
                    help='List of genes to analyze. Ensembl coding'),
        make_option("--organism", type="character", default= "mouse",
                    help="Organism analyzed. Available = human, mouse. Default = mouse")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Load the r object containing the data. 
load(opt$counts)
#df_norm <- subset(df_norm, select=c("Mock_1", "Mock_2", "Mock_3", "IL13_1", "IL13_2", "IL13_3"))

genes = readLines(opt$genelist)

# Transform the ensembl names into gene symbol.
df_norm$Genename <- rownames(df_norm)

# Get rows in the list of genes
clust_df <- df_norm[df_norm$Genename %in% genes, ,drop=FALSE]
clust_df$Genename = NULL

rows_hm <- as.character(mapIds(database, as.character(rownames(clust_df)),
                               'SYMBOL', 'ENSEMBL'))
new <- 1000:2000
rows_hm[is.na(rows_hm)] <- paste("Unk",new[1:sum(is.na(rows_hm))], sep="")
rownames(clust_df) <- rows_hm

# Perform the clustering analysis over the table
# Tree construction (rows and columns)
hr <- hclust(as.dist(1-cor(t(data.matrix(clust_df)),
                           method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(log(data.matrix(clust_df) + 1 ),
                           method="pearson")), method="complete") 

# Establish colors
color <- colorRamp2(c(0, 2), c("white", "red"))

png(file=opt$heatmap, width = 3500, height = 10000, res = 600)
# Mount the heatmap
#<TO_DO>: Add the title of the plot, according to whatever
Heatmap(t(scale(t(data.matrix(clust_df)))), cluster_rows = as.dendrogram(hr),
        cluster_columns = FALSE,
        col=color, column_dend_height = unit(5, "cm"),
        row_dend_width = unit(2, "cm"))
dev.off()

# Save environment
save.image(file=gsub(".png",".RData",opt$heatmap, fixed = TRUE))

# Save versions
get_versions(gsub(".png","_versions.tsv",opt$heatmap, fixed = TRUE))