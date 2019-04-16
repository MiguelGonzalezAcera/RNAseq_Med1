# Run the counts form the bamfile

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(optparse)

# Parse all arguments into an object
option_list = list(
  make_option("--annotation", type="character",
              help="Annotation file, gtf format."),
  make_option("--bamfiles", type="character",
              help="list of bam files of the assay, with metadata."),
  make_option("--obj_out", type="character",
              help="SummarizedExperiment object path with the counts")
)

opt_parser = OprtionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load Rfunctions
source("D:/Documentos/LATESIS/Scripts/Rfunctions.R")

# Read the annotation file
gtfFile = opt$annotation 
txdb = makeTxDbFromGFF(gtfFile, format="gtf")
genes = exonsBy(txdb, by="gene")

# REad the bam file list
sampleTableSingle = read.table(opt$bamfiles, fileEncoding = "UTF8")
fls = sampleTableSingle[,1]
bamLst = BamFileList(fls, index=character(),yieldSize=100000,obeyQname=TRUE)

# Generate the counts for the bamfiles
#<TO_DO>: This step is probably going to be done either externally or 
# internally with featurecounts. Install featurecounts.
Test_experiment = summarizeOverlaps(features = genes,read=bamLst,
                                    mode="Union",
                                    singleEnd=TRUE,
                                    ignore.strand=TRUE,
                                    fragments=FALSE)

# Add the data of the experiment and rename the coumns of the assay
colData(Test_experiment) <- sampleTableSingle

colnames(Test_experiment) <- fls

# Save the R object into a file
save(Test_experiment,file=opt$obj_out)

# Save environment
save.image(file=gsub(".Rda",".RData",opt$obj_out, fixed = TRUE))

# Save versions
get_versions(gsub(".Rda","_versions.tsv",opt$obj_out, fixed = TRUE))