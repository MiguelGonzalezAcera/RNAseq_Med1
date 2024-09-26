# Clustering meddling

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(optparse))

# Might be a table from clusterProfiler, the raw genelist, or might be other source
option_list <- list(
        make_option("--heatmap", type = "character",
                    help = "Heatmap of the coverage of the genes over the samples"),
        make_option("--counts", type = "character",
                    help = "An r object with the normalized counts. Produced in the DE script."),
        make_option("--genelist", type = "character",
                    help = "List of genes to analyze. Ensembl coding"),
        make_option("--organism", type = "character", default = "mouse",
                    help = "Organism analyzed. Available = human, mouse. Default = mouse"),
        make_option("--dims", type = "character", default = "2000,2000",
                    help = "Dimensions of the plot in pixels. Default = 2000,2000"),
        make_option("--colors", type = "character", default = "#7F00FF,#faebd7,#FF8000",
                    help = "Colors for the heatmap, from lower to higher. Default = #7F00FF,#faebd7,#FF8000"),
        make_option("--limits", type = "character", default = "-1,0,1",
                    help = "Limits and center for the color scale. Default = -1,0,1"),
        make_option("--cluster_cols", type = "character", default = "FALSE",
                    help = "Enable column clustering. Options = TRUE, FALSE. Default = FALSE"),
        make_option("--cluster_rows", type = "character", default = "FALSE",
                    help = "Enable row clustering. Options = TRUE, FALSE. Default = FALSE"),
        make_option("--design", type = "character", default = "",
                    help = "Organism analyzed. Available = human, mouse. Default = mouse")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load R scripts
source("Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Load the r object containing the data.
load(opt$counts)

# Select the columns for the design if there were any
if (opt$design != "") {
        sampleTableSingle <- read.table(opt$design, fileEncoding = "UTF8")

        df_norm <- df_norm[c(rownames(sampleTableSingle), "Genename")]
}

# Read the gene file
genes <- readLines(opt$genelist)

# Transform the ensembl names into gene symbol.
df_norm$Genename <- rownames(df_norm)

# Get rows in the list of genes
m <- match(df_norm$Genename, genes)
clust_df <- df_norm[!is.na(m), ][order(na.omit(m)), ]
clust_df$Genename <- NULL

# Check if enough the selected genes are in the provided table
if (dim(clust_df)[1] <= 1) {
        print(sprintf("Could not find enough genes expressed for marker %s.", opt$genelist))
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))
        quit(save = "no")
}

# Change names to gene symbol
rows_hm <- as.character(mapIds(database, as.character(rownames(clust_df)),
                               "SYMBOL", "ENSEMBL"))

# Replace those names without symbol in the databse package with the ensemblid
rows_hm[is.na(rows_hm)|duplicated(rows_hm)] <- rownames(clust_df)[is.na(rows_hm)|duplicated(rows_hm)]

# Transform matrix to numeric
# Copy just in case
cdf <- clust_df

# Columns of the dataframe are factors, so just transforming into matrix will fuck up the data.
# Solution from: https://stackoverflow.com/questions/27528907/how-to-convert-data-frame-column-from-factor-to-numeric
# https://stackoverflow.com/questions/55320327/how-to-apply-a-function-in-each-column-of-a-data-frame
cdf <- sapply(cdf, function(x) as.numeric(as.character(x)))

# Establish the symbols as rows
rownames(cdf) <- rows_hm

# #Wee trick for now
# rows_hm[!(rows_hm %in% c('Arg2','Dio1','Rdh9','Rdh16','Aldh1l1','Ppargc1a'))] <- ""

# Cluster columns if requested
if (as.logical(opt$cluster_cols)) {
        # Quick fix for column clustering in case some values are equal
        a <- cor(log(cdf + 1), method = "pearson")
        a[is.na(a)] <- 0

        # Perform the clustering analysis over the table
        # Tree construction
        cclust <- hclust(as.dist(1 - a), method = "complete")
} else {
        # Assign the FALSe value to the variable
        cclust <- as.logical(opt$cluster_cols)
}

# Cluster rows if requested
if (as.logical(opt$cluster_rows)) {
        # Quick fix for the row clustering in case some values are equal
        b <- cor(log(t(cdf)), method = "pearson")
        b[is.na(b)] <- 0

        # Perform the clustering analysis over the table
        # Tree construction
        rclust <- hclust(as.dist(1 - b), method = "complete")
} else {
        rclust <- as.logical(opt$cluster_rows)
}

# Establish colors and limits of the gradient
color <- colorRamp2(
        c(
                as.integer(strsplit(opt$limits, ",")[[1]][1]),
                as.integer(strsplit(opt$limits, ",")[[1]][2]),
                as.integer(strsplit(opt$limits, ",")[[1]][3])
        ),
        c(
                strsplit(opt$colors, ",")[[1]][1],
                strsplit(opt$colors, ",")[[1]][2],
                strsplit(opt$colors, ",")[[1]][3]
        )
)

# Open the image with the given parameters
png(
        file = opt$heatmap,
        width = as.integer(strsplit(opt$dims, ",")[[1]][1]),
        height = as.integer(strsplit(opt$dims, ",")[[1]][2]),
        res = 300
)
# Mount the heatmap with the respective transformations
Heatmap(t(scale(t(log(cdf + 1)))), cluster_rows = rclust,
        cluster_columns = cclust,
        col = color, column_dend_height = unit(5, "cm"),
        row_labels = rows_hm,
        row_names_gp = gpar(fontsize = (90 / length(genes) + 5)),
        row_dend_width = unit(1, "cm"),
        show_row_names = TRUE,
        #show_row_names = FALSE,
        show_column_names = TRUE,
        heatmap_legend_param = list(
          title = "Fold",
          at = c(
            as.integer(strsplit(opt$limits, ",")[[1]][1]) * 2,
            as.integer(strsplit(opt$limits, ",")[[1]][1]),
            as.integer(strsplit(opt$limits, ",")[[1]][2]),
            as.integer(strsplit(opt$limits, ",")[[1]][3]),
            as.integer(strsplit(opt$limits, ",")[[1]][3]) * 2
            )
        )
        )
dev.off()

#Save as svg also
svg(
        file = gsub(".png", ".svg", opt$heatmap, fixed = TRUE),
        width = as.integer(strsplit(opt$dims, ",")[[1]][1]) / 500,
        height = as.integer(strsplit(opt$dims, ",")[[1]][2]) / 500
)
# Mount the heatmap with the respective transformations
Heatmap(t(scale(t(log(cdf + 1)))), cluster_rows = rclust,
        cluster_columns = cclust,
        col = color, column_dend_height = unit(5, "cm"),
        row_names_gp = gpar(fontsize = (90 / length(genes) + 5)),
        row_dend_width = unit(1, "cm"),
        #show_row_names = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(
          title = "Fold",
          at = c(
            as.integer(strsplit(opt$limits, ",")[[1]][1]) * 2,
            as.integer(strsplit(opt$limits, ",")[[1]][1]),
            as.integer(strsplit(opt$limits, ",")[[1]][2]),
            as.integer(strsplit(opt$limits, ",")[[1]][3]),
            as.integer(strsplit(opt$limits, ",")[[1]][3]) * 2
            )
        )
        )
dev.off()

# Save environment
save.image(file = gsub(".png", ".RData", opt$heatmap, fixed = TRUE))

# Save versions
get_versions(gsub(".png", "_versions.tsv", opt$heatmap, fixed = TRUE))
