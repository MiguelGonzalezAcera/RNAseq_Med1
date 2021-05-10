suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RMariaDB))


option_list = list(
  make_option("--project", type="character", default = "",
              help="Name of the project used for the heatmap. Must be processed and included in the mysql database"),
  make_option("--Rdata", type="character", default = "",
              help="Files with the data. Useful when comparing projects. Only used in case the project name is empty."),
  make_option("--colnames", type="character", default = "",
              help="Names of the columns, in the same order as Rdata. Useful when comparing projects. Only used in case the project name is empty."),
  make_option("--genelist", type="character",
              help="List of genes to be included in the plot, in txt format. Ensembl IDs only"),
  make_option("--heatmap", type="character",
              help="Heatmap of the coverage of the genes over the samples"),
  make_option("--organism", type="character", default= "mouse",
              help="Organism analyzed. Available = human, mouse. Default = mouse")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Read this from the metadata file. 
# These are the full DE analysis results of all the samples

## UNUSED
# treats <- strsplit(opt$DE, ",", fixed = TRUE)[[1]]
# colnames <- strsplit(opt$names, ",", fixed = TRUE)[[1]]
if (opt$project != "") {
  # For the obtaining of the sample data, we shall query the database, as it contains the location of every file processed.
  # Connect to the database
  storiesDb <- dbConnect(RMariaDB::MariaDB(), user='root', password="Plater1a", dbname='Projects', host='localhost')
  # dbListTables(storiesDb)
  
  # Create the query
  query <- sprintf("select * from %s order by Control ASC, Sample ASC;", opt$project)
  
  # Execute and retriece the query
  rsInsert <- dbSendQuery(storiesDb, query)
  dbRows<-dbFetch(rsInsert)
  
  # Select the needed data
  treats = dbRows["Robj_path"][[1]]
  colnames = dbRows["Comparison"][[1]]

  # Disconnect from the database  
  dbDisconnect(storiesDb)
} else {
  treats = strsplit(opt$Rdata, ",")[[1]]
  colnames = strsplit(opt$colnames, ",")[[1]]
}

genelist = readLines(opt$genelist)

#genelist <- unique(as.character(mapIds(database, genes, 'ENSEMBL', 'SYMBOL')))
#genelist <- genelist[!is.na(genelist)]

# 1.- FOLD CHANGE
# Get a dataframe with the columns of the fold change of all the samples
clust_df <- NULL
pval_df <- NULL
for (filename in treats){
  # Read each file
  load(filename)
  full_df <- as.data.frame(res)
  
  # Sort the names by rowname (ENSEMBLID)
  full_df <- full_df[order(row.names(full_df)),]
  
  # Get the fold change and the pvalue column
  FC_df <- full_df[genelist[order(genelist)],"log2FoldChange",drop=FALSE]
  pv_df <- full_df[genelist[order(genelist)],"pvalue",drop=FALSE]
  
  filename_col <- gsub(".Rda","", tail(strsplit(filename, "/"), n=1), fixed = TRUE)
  # Name the column as the file
  colnames(FC_df) <- c(filename)
  colnames(pv_df) <- c(filename)
  
  if (is.null(clust_df) == T){
    # If the final df is empty, fill it with one column
    clust_df <- FC_df
    pval_df <- pv_df
  } else {
    # If not, add the column to the df
    clust_df <- merge(clust_df, FC_df, by=0,all=T)
    rownames(clust_df) <- clust_df$Row.names
    clust_df$Row.names <- NULL
    
    pval_df <- merge(pval_df, pv_df, by=0,all=T)
    rownames(pval_df) <- pval_df$Row.names
    pval_df$Row.names <- NULL
  }
}
# Drop Na values
clust_df <- clust_df[complete.cases(clust_df), ]
pval_df <- pval_df[complete.cases(pval_df), ]

# Keep only genes with valid pvalues

clust_df <- clust_df[rownames(clust_df) %in% rownames(pval_df), ,drop=FALSE]

# Change column names
colnames(clust_df) <- colnames
colnames(pval_df) <- colnames

# Sort by original genelist
#clust_df <- clust_df[match(genelist, rownames(clust_df)),]

# Genenames as Gene symbol
rows_hm <- as.character(mapIds(database, as.character(rownames(clust_df)),
                               'SYMBOL', 'ENSEMBL'))
# new <- 1000:2000
rows_hm[is.na(rows_hm)|duplicated(rows_hm)] <- rownames(clust_df)[is.na(rows_hm)|duplicated(rows_hm)]
#paste("Unk",new[1:sum(is.na(rows_hm))], sep="")
rownames(clust_df) <- rows_hm

if (length(rownames(clust_df)) < 2) {
  print('Genelist doesn\'t have the necessary length to do the clustering in this group of samples')
  quit()
}

# Perform the clustering analysis over the table
# Tree construction (rows and columns)
hr <- hclust(as.dist(1-cor(t(data.matrix(clust_df)),
                           method="pearson")), method="average")
hc <- hclust(as.dist(1-cor(data.matrix(clust_df),
                           method="pearson")), method="average") 

# Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.3)

# Clustering boxes
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
mycolhc <- mycolhc[as.vector(mycl)] 

# Establish colors
color <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

png(file=opt$heatmap, width = 1000, height = 4000, res = 300)
# Mount the heatmap
#<TO_DO>: Add the title of the plot, according to whatever
row_den = color_branches(hr, h = max(hr$height)/1.5) 
Heatmap(data.matrix(clust_df), cluster_rows = as.dendrogram(hr),
        #as.dendrogram(row_den),
        cluster_columns = FALSE, 
        col=color, column_dend_height = unit(5, "cm"),
        #row_dend_width = unit(3, "cm"), 
        row_names_gp = gpar(fontsize = (90/length(genelist)+5)),
        column_names_gp = gpar(fontsize = (90/length(genelist)+5) + 2),
        column_names_max_height = unit(8, "cm"),
        #split = max(mycl), gap = unit(2, "mm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", data.matrix(pval_df)[i, j]), x, y, gp = gpar(fontsize = (80/length(genelist)+3)))
        }
        )
dev.off()

