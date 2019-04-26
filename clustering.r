# Clustering meddling

BiocManager::install(c("ComplexHeatmap","dendextend"))

library(ComplexHeatmap)
library(dendextend)
library(cluster)
library(gplots)
library(org.Mm.eg.db)

# NOTE: this chunk of code shall be replaced by a genelist input.
# Might be a table from clusterProfiler, the raw genelist, or might be other source

load("/DATA/DSS_rec_evolution/DSS_rec_evol.norm_counts.Rda")

genes = readLines("/DATA/DSS_rec_evolution/genelist_pg.txt")
genes = rownames(df_norm)

df_norm$Genenames <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(df_norm)),
                                  'SYMBOL', 'ENSEMBL'))

genelist = df_norm$Genenames[df_norm$Genenames %in% genes]

# Get rows in the list of genes
clust_df <- df_norm[df_norm$Genenames %in% genelist, ,drop=FALSE]

# Genenames as Gene symbol
rows_hm <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(clust_df)),
                               'SYMBOL', 'ENSEMBL'))
rows_hm[is.na(rows_hm)] <- "Unk"
rownames(clust_df) <- rows_hm
clust_df$Genenames = NULL

# Perform the clustering analysis over the table
# Tree construction (rows and columns)
hr <- hclust(as.dist(1-cor(t(data.matrix(clust_df)),
                           method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(data.matrix(clust_df),
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

# Option 1.- Use Heatmap2
heatmap.2(data.matrix(clust_df),main="A title",
          Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=color,
          scale="row", density.info="none", trace="none", RowSideColors=mycolhc,
          margins = c(10,5))

# Option 2.- Use ComplexHeatmap package
#<TO_DO>: Scale the fontsize according to the number of genes and the number of gaps
row_den = color_branches(hr, h = max(hr$height)/1.5) 
Heatmap(t(scale(t(data.matrix(clust_df)))), cluster_rows = as.dendrogram(row_den),
        cluster_columns = as.dendrogram(hc), 
        col=color, column_dend_height = unit(5, "cm"),
        row_dend_width = unit(10, "cm"), 
        row_names_gp = gpar(fontsize = (150/length(genelist)+2)),
        split = max(mycl), gap = unit(2, "mm"))
dev.off()

# This must be done after checking the heatmap
# Transform the subtrees into factor
subtrees <- factor(mycl)

# Get the ID from a gene in the group and extract all of the genes from the group
#<TO_DO>: Geneneame as input. Probably user input. This maybe is another script.
genename <- "Csf1"
subtree_id <- as.character(subtrees[names(subtrees) == genename])
subtree_genelist <- names(subtrees[subtrees == subtree_id])

# Transform the list into entrezids for enrichment by other methods.
subtree_entrez <- as.character(mapIds(org.Mm.eg.db, subtree_genelist,
                                      'ENTREZID', 'SYMBOL'))