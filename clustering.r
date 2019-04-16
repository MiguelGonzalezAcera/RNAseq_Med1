# Clustering meddling

BiocManager::install(c("ComplexHeatmap","dendextend"))

library(ComplexHeatmap)
library(dendextend)
library(cluster)
library(gplots)

# NOTE: this chunk of code shall be replaced by a genelist input.
# Might be a table from clusterProfiler, the raw genelist, or might be other source
rows <- c(1)
genelist <- c()
for (i in rows) {
  genelist <- c(genelist,as.character(mapIds(org.Mm.eg.db, 
                                  unlist(strsplit(as.data.frame(x)$geneID[i], "/")),
                                  'ENSEMBL', 'SYMBOL')))
  genelist <- unique(genelist)
}
length(genelist)
genelist <- genelist[1:50]

# Read this from the metadata file. 
# These are the full DE analysis results of all the samples
treats <- c("test_cluster_INF_MI.Rda","test_cluster_INF_HI.Rda",
            "test_cluster_REC_MOD.Rda","test_cluster_REC_FUL.Rda")

# 1.- FOLD CHANGE
# Get a dataframe with the columns of the fold change of all the samples
clust_df <- NULL
for (filename in treats){
  # Read each file
  load(filename)
  full_df <- as.data.frame(res)
  
  # Sort the names by rowname (ENSEMBLID)
  full_df <- full_df[order(row.names(full_df)),]
  
  # Get the fold change column
  FC_df <- full_df[genelist[order(genelist)],"log2FoldChange",drop=FALSE]
  
  # Name the column as the file
  colnames(FC_df) <- c(filename)
  
  if (is.null(clust_df) == T){
    # If the final df is empty, fill it with one column
    clust_df <- FC_df
  } else {
    # If not, add the column to the df
    clust_df <- merge(clust_df, FC_df, by=0,all=T)
    rownames(clust_df) <- clust_df$Row.names
    clust_df$Row.names <- NULL
    head(clust_df)
  }
}

# 2.- EXPRESSION
# Get the normalized counts
#<TO_DO>: Load dss object
clust_df <- as.data.frame(counts(estimateSizeFactors(dds), normalized = T))
clust_df <- clust_df[genelist[order(genelist)], ,drop=FALSE]

# Genenames as Gene symbol
rows_hm <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(clust_df)),
                               'SYMBOL', 'ENSEMBL'))
rows_hm[is.na(rows_hm)] <- "Unk"
rownames(clust_df) <- rows_hm

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

png(file="Test_heatmap.png", width = 8000, height = 8000, res = 600)
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