library(Matrix)
library(dplyr)
library(ggplot2)
library(rhdf5)
library(Seurat)
library(ggplot2)
library(egg)
library(ggpubr)
library(RColorBrewer)
library(pdftools)
library(data.table)
library(tidyr)
rm(list = ls())

setwd("/Users/jeromesamir/Documents/Jerome/Kirby/Tapestri_pipeline/")
Sys.setenv('R_MAX_VSIZE'=32000000000)
resolutions_to_test = c(0.2, 0.35, 0.5, 0.8, 1, 1.2)
data_dir = file.path("Data")

# Get list of tapestri files to run data on
files = list.files(data_dir, pattern = ".h5$", recursive = T)
files = files[! files %like% "OLD"]
files = file.path(data_dir, files)
all_samples = gsub(".*/", "", gsub(".dna[+]protein.h5", "", sort(files)))
names(files) = all_samples

# Run clustering for each sample and save results
# Exrtract protein data from each tapestri directory
# scale and mormalise data, calculate PCA and UMAP and create metadata containing this information

if (! file.exists('tapestri_data.Rdata')) {
  assays = list()
  metadatas = list()
} else {
  files_new = files
  all_samples_new = all_samples
  load('tapestri_data.Rdata')
  files = files_new
  all_samples = all_samples_new
  remove(all_samples_new, files_new)
}

redo = T # Change this to TRUE if you want to redo all the already processed samples
for (file in all_samples) {
  if (! file %in% names(assays) | ! file %in% names(metadatas) | redo == T) {
    print(file)
    file_assays = suppressWarnings(h5read(files[[file]],"assays"))
    h5closeAll()
    protein_count_matrix = file_assays$protein_read_counts$layers$read_counts
    rownames(protein_count_matrix) = file_assays$protein_read_counts$ca$id
    colnames(protein_count_matrix) = file_assays$protein_read_counts$ra$barcode
    remove(file_assays) # large variable
    gc()
    
    seurat_object = CreateSeuratObject(counts = protein_count_matrix)
    seurat_object = NormalizeData(seurat_object, normalization.method = "CLR",margin=2)
    norm_count_protein_mat = t(seurat_object@assays$RNA@data)
    # norm_count_protein_mat = exp(norm_count_protein_mat) - 1
    
    # Run PCA using variable genes and select the number of PCs to use for clustering using an elbow plot
    seurat_object@assays$RNA@scale.data = as.matrix(seurat_object@assays$RNA@data)
    
    seurat_object = RunPCA(seurat_object, features = rownames(seurat_object))
    # Take only PCs which have a change in variation is > 0.1% from the next -- OTHERWISE, can just use 7
    elb_plot = ElbowPlot(seurat_object, ndims=30)$data
    pct = elb_plot$stdev
    max_pca_dim = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    elb_plot = ggplot() + 
      geom_point(data = elb_plot %>% filter(dims <= max_pca_dim), mapping = aes(dims, stdev), size = 5, color = 'green') +
      geom_point(data = elb_plot %>% filter(dims > max_pca_dim), mapping = aes(dims, stdev), size = 5, color = 'red') +
      ggtitle(file)
    
    # set.seed(888)
    seurat_object = RunUMAP(seurat_object, reduction = "pca", dims = 1:max_pca_dim, min.dist = 0,random.seed=888)
    # set.seed(888)
    seurat_object = FindNeighbors(seurat_object, reduction = "pca", dims = 1:max_pca_dim,k.param=max_pca_dim,random.seed=888)
    
    for (res in resolutions_to_test) {
      # set.seed(888)
      seurat_object = FindClusters(seurat_object, resolution = res,random.seed=888)
    }
    seurat_object@meta.data$sample = file
    
    seurat_object@meta.data = cbind(seurat_object@meta.data, as.data.frame(seurat_object@reductions$pca@cell.embeddings)[rownames(seurat_object@meta.data),c(1:30)])
    seurat_object@meta.data = cbind(seurat_object@meta.data, as.data.frame(seurat_object@reductions$umap@cell.embeddings)[rownames(seurat_object@meta.data),c(1:2)])
    assays[[file]] = norm_count_protein_mat
    metadatas[[file]] = seurat_object@meta.data %>% select(-seurat_clusters)
  }
}
assays = assays[all_samples]
metadatas = metadatas[all_samples]
files = files[all_samples]
save(all_samples,assays,metadatas,files,file='tapestri_data.Rdata')
