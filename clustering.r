# Clustering meddling

BiocManager::install(c("ComplexHeatmap","dendextend"))

library(ComplexHeatmap)
library(dendextend)
library(cluster)
library(gplots)
library(org.Mm.eg.db)

# NOTE: this chunk of code shall be replaced by a genelist input.
# Might be a table from clusterProfiler, the raw genelist, or might be other source

# Load the r object containing the data. Can be fpkm or normalized counts.
load("/DATA/DSS_rec_evolution/DSS_rec_evol.fpkm.Rda")
load("/DATA/DSS_rec_evolution/DSS_rec_evol.norm_counts.Rda")

genes = readLines("/DATA/DSS_rec_evolution/genelist_markers.txt")

# Transform the ensembl names into gene symbol. NOTE that the name of the variable must change.
fpkm_df$Genenames <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(fpkm_df)),
                                  'SYMBOL', 'ENSEMBL'))

# Get rows in the list of genes
clust_df <- fpkm_df[fpkm_df$Genenames %in% genes, ,drop=FALSE]

# Genenames as Gene symbol
rows_hm <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(clust_df)),
                               'SYMBOL', 'ENSEMBL'))
rows_hm[is.na(rows_hm)] <- "Unk"
rownames(clust_df) <- rows_hm
clust_df$Genenames = NULL

# Repeat process with the complete data for the column tree. Take into account the origin
df_norm$Genenames = NULL

# Perform the clustering analysis over the table
# Tree construction (rows and columns)
hr <- hclust(as.dist(1-cor(t(data.matrix(clust_df)),
                           method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(data.matrix(df_norm),
                           method="pearson")), method="average") 

# Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5)

# Clustering boxes
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
mycolhc <- mycolhc[as.vector(mycl)] 

# Establish colors
color <- colorpanel(100, "blue", "white", "red")

png(file="/DATA/DSS_rec_evolution/Test_heatmap.png", width = 8000, height = 8000, res = 600)
# Mount the heatmap
#<TO_DO>: Add the title of the plot, according to whatever
row_den = color_branches(hr, h = max(hr$height)/1.5) 
Heatmap(t(scale(t(data.matrix(clust_df)))), cluster_rows = as.dendrogram(row_den),
        cluster_columns = as.dendrogram(hc), 
        col=color, column_dend_height = unit(5, "cm"),
        row_dend_width = unit(10, "cm"), 
        row_names_gp = gpar(fontsize = (150/length(genelist)+5)),
        split = max(mycl), gap = unit(2, "mm"))
dev.off()