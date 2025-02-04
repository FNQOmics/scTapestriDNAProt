library(data.table)
if(!require(rhdf5)){install.packages("rhdf5")}

# Set working directory and file paths
args <- commandArgs(trailingOnly = TRUE)
# Check if arguments are provided
if (length(args) == 0) {
  stop("No arguments provided; requires h5_object and optional output_directory")
}
if (length(args) ==1) {
  outdir <- getwd()
}
if (length(args) ==2) {
  outdir <- args[2]
}


file_h5 <- args[1]

# Check if the file exists
if (!file.exists(file_h5)) {
  stop("HDF5 file not found at the specified path.")
}

# Save the data as a CSV file
if (!dir.exists(outdir)) {
  stop("Output directory does not exist.")
}

setwd(outdir)

h5_file_only <- basename(file_h5)
h5_tsv = paste0(outdir,"/",h5_file_only)
h5_tsv <- gsub(".h5",".tsv", h5_tsv)
h5_tsv <- gsub(".hdf5",".tsv", h5_tsv)
print(h5_tsv)

# Close any open HDF5 connections
h5closeAll()

# List all groups/datasets in the file for inspection
cat("Listing contents of the HDF5 file:\n")
print(h5ls(file_h5))

h5_obj <- h5read(file_h5, name="/assays/")


# Read relevant datasets from the HDF5 file
cat("\nReading datasets...\n")
chrom <- h5read(file_h5, name = "/assays/dna_variants/ca/CHROM")
pos <- h5read(file_h5, name = "/assays/dna_variants/ca/POS")
ref <- h5read(file_h5, name = "/assays/dna_variants/ca/REF")
alt <- h5read(file_h5, name = "/assays/dna_variants/ca/ALT")
id_data <- h5read(file_h5, name = "/assays/dna_variants/ca/id")
af <- h5read(file_h5, name = "/assays/dna_variants/layers/AF")
ngt <- h5read(file_h5, name = "/assays/dna_variants/layers/NGT")
dp <- h5read(file_h5, name = "/assays/dna_variants/layers/DP")
gq <- h5read(file_h5, name = "/assays/dna_variants/layers/GQ")
cell_barcodes <- h5read(file_h5, name = "/assays/dna_variants/ra/barcode")

# Validate dimensions
if (nrow(af) != length(chrom)) {
  stop("Mismatch between the number of variants and AF data rows.")
}
if (ncol(af) != length(cell_barcodes)) {
  stop("Mismatch between the number of cells and AF data columns.")
}

# Create a data frame for variants
variant_info <- data.frame(
  CHROM = chrom,
  POS = pos,
  ID = id_data,
  REF = ref,
  ALT = alt
)

# Create per-cell data with actual cell barcodes as column identifiers
cat("\nCreating per-cell data...\n")
cell_data_list <- list()

for (i in seq_along(cell_barcodes)) {
  cell_data_list[[i]] <- paste0("DP=",dp[,i],';AF=',af[,i],";GQ=",gq[,i],";NGT=",ngt[,i])
}

# Combine all cell data into one data frame
cat("\nCombining data...\n")
all_cells_data <- do.call(cbind, cell_data_list)
colnames(all_cells_data) <- cell_barcodes
# Merge variant and cell data
final_data <- cbind(variant_info, all_cells_data)


cat("\nSaving final data as TSV...\n")
write.table(
  final_data,
  file = h5_tsv,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

cat("\nProcess completed successfully!\n")

