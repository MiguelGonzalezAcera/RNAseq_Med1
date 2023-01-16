# Clustering meddling

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(optparse))

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

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
                    help="Organism analyzed. Available = human, mouse. Default = mouse"),
        make_option("--dims", type="character", default= "2000,2000",
                    help="Dimensions of the plot in pixels. Default = 2000,2000"),
        make_option("--design", type="character", default = "",
                    help="Organism analyzed. Available = human, mouse. Default = mouse")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R scripts
source("Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Load the r object containing the data.
load(opt$counts)

# Select the columns for the design if there were any
if (opt$design != "") {
        sampleTableSingle = read.table(opt$design, fileEncoding = "UTF8")

        df_norm <- df_norm[c(rownames(sampleTableSingle) ,'Genename')]
}

# Read the gene file
genes = readLines(opt$genelist)

# Transform the ensembl names into gene symbol.
df_norm$Genename <- rownames(df_norm)

# Get rows in the list of genes
m <- match(df_norm$Genename, genes)
clust_df <- df_norm[!is.na(m),][order(na.omit(m)),]
#clust_df <- df_norm[df_norm$Genename %in% genes, ,drop=FALSE]
clust_df$Genename = NULL

if (dim(clust_df)[1] <= 1) {
        print(sprintf("Could not find enough genes expressed for marker %s.", opt$genelist))
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))
        quit(save = "no")
}

rows_hm <- as.character(mapIds(database, as.character(rownames(clust_df)),
                               'SYMBOL', 'ENSEMBL'))
# new <- 1000:2000
rows_hm[is.na(rows_hm)|duplicated(rows_hm)] <- rownames(clust_df)[is.na(rows_hm)|duplicated(rows_hm)]
        #paste("Unk",new[1:sum(is.na(rows_hm))], sep="")
rownames(clust_df) <- rows_hm

# Transform matrix to numeric
cdf = sapply(clust_df, as.numeric)
rownames(cdf) <- rows_hm

# Quick fix for column clustering in case some values are equal
a <- cor(log(cdf + 1), method='pearson')
a[is.na(a)] = 0

# Quick fix for the row clustering in case some values are equal
b <- cor(log(t(cdf)), method='pearson')
b[is.na(b)] = 0

# Perform the clustering analysis over the table
# Tree construction (rows and columns)
hr <- hclust(as.dist(1-b), method="complete")
hc <- hclust(as.dist(1-a), method="complete")

# Establish colors
color <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

png(file=opt$heatmap, width = int(strsplit(opt$dims, ',')[0]), height = int(strsplit(opt$dims, ',')[1]), res = 300)
# Mount the heatmap
#<TO_DO>: Add the title of the plot, according to whatever
Heatmap(t(scale(t(log(cdf + 1)))), cluster_rows = FALSE,
        cluster_columns = FALSE,
        col=color, column_dend_height = unit(5, "cm"),
        row_names_gp = gpar(fontsize = (90/length(genes)+5)),
        row_dend_width = unit(2, "cm"), show_row_names = TRUE)
# dev.off()

# Save environment
save.image(file=gsub(".png",".RData",opt$heatmap, fixed = TRUE))

# Save versions
get_versions(gsub(".png","_versions.tsv",opt$heatmap, fixed = TRUE))
