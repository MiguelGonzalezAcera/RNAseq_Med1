suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(circlize))

# Create options
option_list = list(
  make_option("--heatmap", type="character",
              help="Heatmap of the coverage of the genes over the samples"),
  make_option("--counts", type="character",
              help="An r object with the normalized counts. Produced in the DE script."),
  make_option("--design", type="character",
              help="File with the design of the experiment."),
  make_option("--out_obj", type="character",
              help="DESeq2 object with the result of the analysis."),
  make_option("--organism", type="character", default= "mouse",
              help="Organism analyzed. Available = human, mouse. Default = mouse"),
  make_option("--control", type="character",
              help="Value from the designs to use as control"),
  make_option("--comparisons", type="character",
              help="Values from the designs, comma separated, to compare against control.")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load R scripts
source("/DATA/RNAseq_test/Scripts/Rscripts/Rfunctions.R")

# Select organism
database <- select.organism(opt$organism)

# Load the r object containing the data.
load(opt$counts)

# Remove the genename column
wdf_norm <- df_norm[ , !names(df_norm) %in% c("Genename")]

# Create the lists of genes
genes <- list()

if (opt$organism == 'mouse'){
  genes[['B_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/B_cells_ensembl.txt")
  genes[['Endothelial']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Endothelial_ensembl.txt")
  genes[['EntericGlial']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/EntericGlial_ensembl.txt")
  genes[['EntericNeuron']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/EntericNeuron_ensembl.txt")
  genes[['EnterocyteDist']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/EnterocyteDist_ensembl.txt")
  genes[['EnterocyteProx']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/EnterocyteProx_ensembl.txt")
  genes[['Enteroendocrine']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Enteroendocrine_ensembl.txt")
  genes[['Fibroblasts']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Fibroblasts_ensembl.txt")
  genes[['Goblet']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Goblet_ensembl.txt")
  genes[['Mast_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Mast_cells_ensembl.txt")
  genes[['M_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/M_cells_ensembl.txt")
  genes[['MO_DC']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/MO_DC_ensembl.txt")
  genes[['Neutrophils']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Neutrophils_ensembl.txt")
  genes[['NK_ILC1']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/NK_ILC1_ensembl.txt")
  genes[['Paneth']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Paneth_ensembl.txt")
  genes[['Plasma_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Plasma_cells_ensembl.txt")
  genes[['Smooth_muscle']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Smooth_muscle_ensembl.txt")
  genes[['StemProg']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/StemProg_ensembl.txt")
  genes[['TAProg']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/TAProg_ensembl.txt")
  genes[['T_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/T_cells_ensembl.txt")
  genes[['Tuft']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Tuft_ensembl.txt")
} else if (opt$organism == 'human'){
  genes[['B_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/B_cells_human_ensembl.txt")
  genes[['Endothelial']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Endothelial_human_ensembl.txt")
  genes[['EntericGlial']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/EntericGlial_human_ensembl.txt")
  genes[['EntericNeuron']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/EntericNeuron_human_ensembl.txt")
  genes[['EnterocyteDist']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/EnterocyteDist_human_ensembl.txt")
  genes[['EnterocyteProx']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/EnterocyteProx_human_ensembl.txt")
  genes[['Enteroendocrine']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Enteroendocrine_human_ensembl.txt")
  genes[['Fibroblasts']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Fibroblasts_human_ensembl.txt")
  genes[['Goblet']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Goblet_human_ensembl.txt")
  genes[['Mast_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Mast_cells_human_ensembl.txt")
  genes[['M_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/M_cells_human_ensembl.txt")
  genes[['MO_DC']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/MO_DC_human_ensembl.txt")
  genes[['Neutrophils']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Neutrophils_human_ensembl.txt")
  genes[['NK_ILC1']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/NK_ILC1_human_ensembl.txt")
  genes[['Paneth']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Paneth_human_ensembl.txt")
  genes[['Plasma_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Plasma_cells_human_ensembl.txt")
  genes[['Smooth_muscle']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Smooth_muscle_human_ensembl.txt")
  genes[['StemProg']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/StemProg_human_ensembl.txt")
  genes[['TAProg']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/TAProg_human_ensembl.txt")
  genes[['T_cells']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/T_cells_human_ensembl.txt")
  genes[['Tuft']] <- readLines("/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Tuft_human_ensembl.txt")
}

# transform the data to numeric matrix
options(digits=20)
numwdf_norm <- matrix(as.numeric(as.matrix(wdf_norm)),  ncol = ncol(as.matrix(wdf_norm)))
rownames(numwdf_norm) <- rownames(wdf_norm)
colnames(numwdf_norm) <- colnames(wdf_norm)

# Run the gsva ove the set of counts
gsva.es <- gsva(numwdf_norm, genes, verbose=FALSE)

# Save the GSVA table as R object and table
save(gsva.es, file=opt$out_obj)
write.table(gsva.es, file=gsub(".Rda",".tsv",opt$out_obj, fixed = TRUE),sep="\t")

# Establish colors
color <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create a heatmap of the values
png(file=opt$heatmap, width = 3500, height = 3500, res = 600)
Heatmap(gsva.es,cluster_columns = FALSE,
        col=color, column_dend_height = unit(5, "cm"),
        row_dend_width = unit(2, "cm"))
dev.off()

# diff expression with limma
sampleTableSingle = read.table(opt$design, fileEncoding = "UTF8")

# Design model matrix
#design <- model.matrix( ~ as.character(sampleTableSingle[,1]) + as.character(sampleTableSingle[,2]))
Tr1 = relevel(factor(sampleTableSingle[,1]), opt$control)
design <- model.matrix( ~ Tr1)

# Run limma
fit <- lmFit(gsva.es, design)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.1)

# Save each of the results by sample
for (sample in strsplit(opt$comparisons, ",")[[1]]){
  # Get result of the diff expression
  restab <- topTable(fit, coef=paste("Tr1",sample,sep=""), number = 100)
  
  # Save the full result object
  # Contrast name will be replaced by the sample and controls
  res_name = paste(paste("", sample, opt$control, sep='_'),"Rda", sep=".")
  save(restab,file=gsub(".Rda",res_name,opt$out_obj, fixed = TRUE))
  
  # Save the table
  restab_name = paste(paste("", sample, opt$control, sep='_'),"tsv", sep=".")
  write.table(restab, file=gsub(".Rda",restab_name,opt$out_obj, fixed = TRUE),
              sep="\t")
}

# Save environment
save.image(file=gsub(".Rda",".RData",opt$out_obj, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda","_versions.tsv",opt$out_obj, fixed = TRUE))
