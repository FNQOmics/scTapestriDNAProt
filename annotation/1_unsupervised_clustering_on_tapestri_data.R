try(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
rm(list = ls())
library(Matrix)
library(dplyr)
library(ggplot2)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pdftools)
library(data.table)
library(tidyr)
setwd("/Users/jeromesamir/Documents/Jerome/Kirby/Tapestri_pipeline")
source("~/Documents/Jerome/Kirby/Scripts/jeromefuncs.R")
data_dir = "Data"
# Set themes and colours
theme_set(theme_classic())
colours = set1
saved_colours = colours
colours = c(colours, 'cyan', brewer.pal(8, 'Set3')[c(1,3,6,7)], "orange")
names(colours) = c("CD4-CD8- T","CD4 T","B cell","CD8 T","CD3- T (CD7-)","Myeloid","Doublet","Dead","CD4+CD8+ T", 'Unknown', "NK", "DC/Monocyte", "CD3- T (CD7+)", "Fibroblast", "CD3- T")
proteins_to_show = c("CD19","CD45","CD3","CD4","CD8","CD326","CD45RA","CD56", "CD57","CD7","CD14", "CD16", "CD11b", "CD11c","HLA_DR", "CD38", "CD103", "CD90", "CD64", "CD13")
percentage_nonzero <- function(x) {
  sum(x != 0, na.rm = TRUE) / length(x) * 100
}
addSmallLegend <- function(myPlot, pointSize = 5, textSize = 8, spaceLegend = 0.7) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


# Load tapestri data which was created in v0
load('tapestri_data.Rdata')
samples_to_run_unsupervised_clustering_on = all_samples


# Plot proteins for different clusters at different resolutions
dotplots = list()
for (s in samples_to_run_unsupervised_clustering_on) {
  print(s)
  met = metadatas[[s]]
  all_proteins = gsub("[-/]", "_", colnames(assays[[s]]))
  colnames(assays[[s]]) = all_proteins
  all_proteins = all_proteins[all_proteins %like% "^[A-Za-z0-9_.]*$"]
  assays[[s]] = assays[[s]][,all_proteins]
  met = cbind(met, assays[[s]][rownames(met),])
  dotplots[[s]] = data.frame()
  
  # Just take highest resolution by default except in some few cases
  res = 1
  if (s %like% "AP_21_RA") {
    res = 0.8 # or 0.5
  } else if (s == 'AP23-RA_-_2_lanes') {
    res = 0.8
  } else if (s == 'IBM-GMR-muscle') {
    res = 0.5 # or 0.35
  } else if (s == 'KASH-_synovium_fluid') {
    res = 0.5
  }
  
  formula = as.formula(paste0('cbind(', paste(all_proteins, collapse = ","), ") ~ RNA_snn_res.", res))  
  result = aggregate(formula, data = met, FUN = mean, na.rm = TRUE) %>% 
    pivot_longer(cols = 2:(length(all_proteins) + 1),
                 values_to = "mean_exp",
                 names_to = "protein")
  result_scaled = aggregate(formula, data = met, FUN = mean, na.rm = TRUE)
  result_scaled[all_proteins] = lapply(result_scaled[all_proteins], range01)
  result_scaled = result_scaled %>% 
    pivot_longer(cols = 2:(length(all_proteins) + 1),
                 values_to = "sc_m_exp",
                 names_to = "protein")
  result_perc = aggregate(formula, data = met, FUN = percentage_nonzero) %>% 
    pivot_longer(cols = 2:(length(all_proteins) + 1),
                 values_to = "perc_exp",
                 names_to = "protein")
  result$perc_exp = result_perc$perc_exp
  result$sc_m_exp = result_scaled$sc_m_exp
  result$resolution = res
  colnames(result)[1] = 'cluster'
  dotplots[[s]] = rbind(dotplots[[s]], result)
  dotplots[[s]] = dotplots[[s]] %>% select(resolution, cluster, everything())
  
  # Add fold changes
  sobj = CreateSeuratObject(counts = t(met[all_proteins]), meta.data = met[! colnames(met) %in% all_proteins])
  sobj@assays$RNA@scale.data = as.matrix(sobj@assays$RNA@data)
  Idents(sobj) = sobj@meta.data[,paste0('RNA_snn_res.', res)]
  result = FindAllMarkers(sobj, min.pct = 0.001, logfc.threshold = 0.001) %>% filter(p_val_adj < 0.05) %>% mutate(gene = gsub("-", "_", gene))
  dotplots[[s]] = merge(dotplots[[s]], result[], by.x=c("cluster", "protein"), by.y = c("cluster", "gene"), all.x=T, all.y=F)
  
  dtp = as.data.table(dotplots[[s]])
  dtp = dtp[,.SD[which.max(mean_exp)],by=protein]
  protein_order = as.character(dtp[order(-dtp$cluster, dtp$mean_exp),]$protein)
  protein_order = unique(c("CD326","CD19","CD45","CD3","CD4","CD8", protein_order))
  dotplots[[s]]$protein = factor(dotplots[[s]]$protein, levels = protein_order)
  
  # stupid hack to try and cluster clusters together to make cell annotation easier - doesn't really do much
  # cluster_order = hclust(dist(dotplots[[s]] %>% 
  #                   filter(protein %in% c("CD326", "CD19", "CD3", "CD4", "CD8", "CD45")) %>% 
  #                   select(cluster, protein, mean_exp) %>% 
  #                   tidyr::pivot_wider(id_cols = 'cluster', names_from = 'protein', values_from = 'mean_exp') %>%
  #                   tibble::column_to_rownames("cluster")))$order - 1
  cluster_order = 0:100
  dotplots[[s]]$cluster = factor(as.character(dotplots[[s]]$cluster),
                                 levels = as.character(cluster_order))
  
  pdf(file.path(data_dir, s, 'clustering.pdf'), height = 4, width = 20)
  print(
    ggarrange(addSmallLegend(ggplot(met, aes_string('UMAP_1', 'UMAP_2', color=paste0("RNA_snn_res.", res))) + 
                               geom_point(size=1) + 
                               scale_color_manual(values = c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Set2"))) +
                               labs(color = 'cluster', title = paste0("Resolution: ", res))),
              ggplot(dotplots[[s]] %>% filter(protein %in% proteins_to_show & resolution == res), aes(x=cluster,y=protein,colour=mean_exp,size=perc_exp)) + 
                geom_point() + 
                labs(colour="Avg expr.", size="% Expr.") +
                scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                theme(axis.text.x = element_text(face='bold',colour = c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"))[cluster_order + 1])) + 
                scale_color_viridis_c(),
              ggplot(dotplots[[s]] %>% filter(protein %in% proteins_to_show & resolution == res), aes(x=cluster,y=protein,colour=sc_m_exp,size=perc_exp)) + 
                geom_point() + 
                labs(colour="Scaled\navg. expr.", size="% Expr.") +
                scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                theme(axis.text.x = element_text(face='bold',colour = c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"))[cluster_order + 1])) + 
                scale_color_viridis_c(),
              ggplot(dotplots[[s]] %>% filter(protein %in% proteins_to_show & resolution == res), aes(x=cluster,y=protein,colour=avg_log2FC,size=-log10(p_val_adj))) + 
                geom_point() + 
                labs(colour="LogFC", size="-Log(adj-P)") +
                scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                theme(axis.text.x = element_text(face='bold',colour = c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"))[cluster_order + 1])) + 
                scale_colour_gradientn(colors=c('darkblue','#f7f7f7','darkred'), limits=c(-max(abs(dotplots[[s]]['avg_log2FC']), na.rm = T), max(abs(dotplots[[s]]['avg_log2FC']), na.rm = T))),
              ncol = 4, nrow = 1)
  )
  dev.off()
  
  pdf(file.path(data_dir, s, 'clustering_all_markers.pdf'), height = length(all_proteins)/8, width = 20)
  print(
    ggarrange(addSmallLegend(ggplot(met, aes_string('UMAP_1', 'UMAP_2', color=paste0("RNA_snn_res.", res))) + 
                               geom_point(size=1) + 
                               scale_color_manual(values = c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Set2"))) +
                               labs(color = 'cluster', title = paste0("Resolution: ", res))),
              ggplot(dotplots[[s]] %>% filter(resolution == res), aes(x=cluster,y=protein,colour=mean_exp,size=perc_exp)) + 
                geom_point() + 
                labs(colour="Avg expr.", size="% Expr.") +
                scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                theme(axis.text.x = element_text(face='bold',colour = c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3")))) + 
                scale_color_viridis_c() +
                scale_size_continuous(range = c(1,3)),
              ggplot(dotplots[[s]] %>% filter(resolution == res), aes(x=cluster,y=protein,colour=sc_m_exp,size=perc_exp)) + 
                geom_point() + 
                labs(colour="Scaled\navg. expr.", size="% Expr.") +
                scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                theme(axis.text.x = element_text(face='bold',colour = c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3")))) + 
                scale_color_viridis_c() +
                scale_size_continuous(range = c(1,3)),
              ggplot(dotplots[[s]] %>% filter(resolution == res), aes(x=cluster,y=protein,colour=avg_log2FC,size=-log10(p_val_adj))) + 
                geom_point() + 
                labs(colour="LogFC", size="-Log(adj-P)") +
                scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                theme(axis.text.x = element_text(face='bold',colour = c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3")))) + 
                scale_colour_gradientn(colors=c('darkblue','#f7f7f7','darkred'), limits=c(-max(abs(dotplots[[s]]['avg_log2FC']), na.rm = T), max(abs(dotplots[[s]]['avg_log2FC']), na.rm = T))) +
                scale_size_continuous(range = c(1,3)),
              ncol = 4, nrow = 1)
  )
  dev.off()
}


# Manually annotate clusters
fct_orders = c("CD4 T","CD8 T","CD4-CD8- T","CD4+CD8+ T","CD3- T (CD7+)","CD3- T (CD7-)","Gamma Delta", "B cell","NK","Myeloid", "DC/Monocyte", "Fibroblast","Unknown","Doublet","Dead")

for (s in samples_to_run_unsupervised_clustering_on) {
  met = metadatas[[s]]
  all_proteins = gsub("[-/]", "_", colnames(assays[[s]]))
  colnames(assays[[s]]) = all_proteins
  all_proteins = all_proteins[all_proteins %like% "^[A-Za-z0-9_.]*$"]
  assays[[s]] = assays[[s]][,all_proteins]
  met = cbind(met, assays[[s]][rownames(met),])
  do_plot = F
  
  # Just take highest resolution by default (except in some few cases)
  res = 1
  if (s %like% "AP_21_RA") {
    res = 0.8 # or 0.5
  } else if (s == 'AP23-RA_-_2_lanes') {
    res = 0.8
  } else if (s == 'IBM-GMR-muscle') {
    res = 0.5 # or 0.35
  } else if (s == 'KASH-_synovium_fluid') {
    res = 0.5
  }
  
  if (s == "AP_21_RA") {
    cluster_names = c(
      "0"="NK", #?
      "1"="DC/Monocyte",
      "2"="Fibroblast", #?
      "3"="Dead",
      "4"="Unknown", #NK?
      "5"="Fibroblast",
      "6"="CD8 T",
      "7"="B cell",
      "8"="Myeloid"
    )
    do_plot = T
    
  } else if (s == 'AP23-RA_-_2_lanes') {
    cluster_names = c(
      "0"="Fibroblast", #?
      "1"="Myeloid",
      "2"="B cell",
      "3"="CD8 T",
      "4"="CD8 T",
      "5"="Myeloid",
      "6"="NK", #?
      "7"="Fibroblast", #?
      "8"="B cell",
      "9"="Unknown",
      "10"="NK", #?
      "11"="Unknown"
    )
    do_plot = T
    
  } else if (s == 'IBM-GMR-muscle') {
    cluster_names = c(
      "0"="DC/Monocyte", #?
      "1"="Dead", #? -- these are everywhere, even with higher resolutions
      "2"="CD8 T",
      "3"="CD8 T",
      "4"="Fibroblast", #?
      "5"="DC/Monocyte", #?
      "6"="CD4-CD8- T",
      "7"="DC/Monocyte", #?
      "8"="Fibroblast",
      "9"="CD8 T"
    )
    do_plot = T
    
  } else if (s == 'KASH-_synovium_fluid') {
    cluster_names = c(
      "0"="Myeloid", #?
      "1"="CD8 T",
      "2"="CD4 T",
      "3"="CD4 T",
      "4"="B cell",
      "5"="Unknown", #?
      "6"="DC/Monocyte", # Myeloid?
      "7"="NK", #?
      "8"="DC/Monocyte" # Myeloid?
    )
    do_plot = T
    
  }
  
  if (do_plot) {
    met$celltype_name1 = plyr::revalue(met[[paste0("RNA_snn_res.", res)]], cluster_names)
    met$celltype_name2 = paste0(met[[paste0("RNA_snn_res.", res)]], ": ", met$celltype_name1)
    
    temp_df = data.frame(num = as.numeric(names(cluster_names)), name = unname(cluster_names))
    temp_df$name = factor(temp_df$name, fct_orders)
    temp_df$pasted = paste(temp_df$num, temp_df$name, sep=': ')
    cluster_order = temp_df[order(temp_df$name, temp_df$num),'pasted']
    
    
    # For Celltype1
    formula = as.formula(paste0('cbind(', paste(all_proteins, collapse = ","), ") ~ celltype_name1"))
    result = aggregate(formula, data = met, FUN = mean, na.rm = TRUE)
    result = result %>%
      pivot_longer(cols = 2:(length(all_proteins) + 1),
                   values_to = "mean_exp",
                   names_to = "protein")
    result_perc = aggregate(formula, data = met, FUN = percentage_nonzero) %>%
      pivot_longer(cols = 2:(length(all_proteins) + 1),
                   values_to = "perc_exp",
                   names_to = "protein")
    result_scaled = aggregate(formula, data = met, FUN = mean, na.rm = TRUE)
    result_scaled[all_proteins] = lapply(result_scaled[all_proteins], range01)
    result_scaled = result_scaled %>%
      pivot_longer(cols = 2:(length(all_proteins) + 1),
                   values_to = "sc_m_exp",
                   names_to = "protein")
    result$sc_m_exp = result_scaled$sc_m_exp
    result$perc_exp = result_perc$perc_exp
    result$resolution = res
    colnames(result)[1] = 'cluster'
    result_celltype_1 = result %>% select(resolution, cluster, everything())
    
    # Add fold changes
    sobj = CreateSeuratObject(counts = t(met[all_proteins]), meta.data = met[! colnames(met) %in% all_proteins])
    sobj@assays$RNA@scale.data = as.matrix(sobj@assays$RNA@data)
    Idents(sobj) = sobj@meta.data$celltype_name1
    fc_result = FindAllMarkers(sobj, min.pct = 0.001, logfc.threshold = 0.001) %>% filter(p_val_adj < 0.05) %>% mutate(gene = gsub("-", "_", gene))
    result_celltype_1 = merge(result_celltype_1, fc_result, by.x=c("cluster", "protein"), by.y = c("cluster", "gene"), all.x=T, all.y=F)
    
    dtp = as.data.table(result)
    dtp = dtp[,.SD[which.max(mean_exp)],by=protein]
    protein_order_celltype_1 = as.character(dtp[order(-dtp$cluster, dtp$mean_exp),]$protein)
    
    # For Celltype2
    formula = as.formula(paste0('cbind(', paste(all_proteins, collapse = ","), ") ~ celltype_name2"))
    result = aggregate(formula, data = met, FUN = mean, na.rm = TRUE) %>%
      pivot_longer(cols = 2:(length(all_proteins) + 1),
                   values_to = "mean_exp",
                   names_to = "protein")
    result_perc = aggregate(formula, data = met, FUN = percentage_nonzero) %>%
      pivot_longer(cols = 2:(length(all_proteins) + 1),
                   values_to = "perc_exp",
                   names_to = "protein")
    result_scaled = aggregate(formula, data = met, FUN = mean, na.rm = TRUE)
    result_scaled[all_proteins] = lapply(result_scaled[all_proteins], range01)
    result_scaled = result_scaled %>%
      pivot_longer(cols = 2:(length(all_proteins) + 1),
                   values_to = "sc_m_exp",
                   names_to = "protein")
    result$sc_m_exp = result_scaled$sc_m_exp
    result$perc_exp = result_perc$perc_exp
    result$resolution = res
    colnames(result)[1] = 'cluster'
    result_celltype_2 = result %>% select(resolution, cluster, everything())
    
    # Add fold changes
    sobj = CreateSeuratObject(counts = t(met[all_proteins]), meta.data = met[! colnames(met) %in% all_proteins])
    sobj@assays$RNA@scale.data = as.matrix(sobj@assays$RNA@data)
    Idents(sobj) = sobj@meta.data$celltype_name2
    fc_result = FindAllMarkers(sobj, min.pct = 0.001, logfc.threshold = 0.001) %>% filter(p_val_adj < 0.05) %>% mutate(gene = gsub("-", "_", gene))
    result_celltype_2 = merge(result_celltype_2, fc_result, by.x=c("cluster", "protein"), by.y = c("cluster", "gene"), all.x=T, all.y=F)
    
    
    
    dtp = as.data.table(result)
    dtp = dtp[,.SD[which.max(mean_exp)],by=protein]
    protein_order_celltype_2 = as.character(dtp[order(-dtp$cluster, dtp$mean_exp),]$protein)
    result_celltype_2$cluster = factor(result_celltype_2$cluster, levels = cluster_order)
    
    
    pdf(file.path(data_dir, s, 'clustering_final.pdf'), height = 4.85, width = 20)
    print(
      annotate_figure(top=s, ggarrange(addSmallLegend(ggplot(met, aes(UMAP_1, UMAP_2, color=celltype_name1)) +
                                                        geom_point(size=1) +
                                                        labs(color = 'Cluster') +
                                                        scale_color_manual(values = colours)),
                                       ggplot(result_celltype_1 %>% filter(protein %in% proteins_to_show), aes(x=cluster,y=factor(protein, levels = rev(proteins_to_show)),colour=mean_exp,size=perc_exp)) +
                                         geom_point() +
                                         labs(colour="Avg expr.", size="% Expr.") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_color_viridis_c(),
                                       ggplot(result_celltype_1 %>% filter(protein %in% proteins_to_show), aes(x=cluster,y=factor(protein, levels = rev(proteins_to_show)),colour=sc_m_exp,size=perc_exp)) +
                                         geom_point() +
                                         labs(colour="Scaled\nAvg expr.", size="% Expr.") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_color_viridis_c(),
                                       ggplot(result_celltype_1 %>% filter(protein %in% proteins_to_show), aes(x=cluster,y=factor(protein, levels = rev(proteins_to_show)),colour=avg_log2FC,size=-log10(p_val_adj))) +
                                         geom_point() +
                                         labs(colour="LogFC", size="-Log(adj-P)") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_colour_gradientn(colors=c('darkblue','#f7f7f7','darkred'), limits=c(-max(abs(result_celltype_1['avg_log2FC']), na.rm = T), max(abs(result_celltype_1['avg_log2FC']), na.rm = T))),
                                       ncol = 4, nrow = 1))
    )
    
    print(
      annotate_figure(top=s, ggarrange(addSmallLegend(ggplot(met, aes(UMAP_1, UMAP_2, color=celltype_name2)) +
                                                        geom_point(size=1) +
                                                        labs(color = 'Cluster') +
                                                        scale_color_manual(values = c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")))),
                                       ggplot(result_celltype_2 %>% filter(protein %in% proteins_to_show), aes(x=cluster,y=factor(protein, levels = rev(proteins_to_show)),colour=mean_exp,size=perc_exp)) +
                                         geom_point() +
                                         labs(colour="Avg expr.", size="% Expr.") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_color_viridis_c(),
                                       ggplot(result_celltype_2 %>% filter(protein %in% proteins_to_show), aes(x=cluster,y=factor(protein, levels = rev(proteins_to_show)),colour=sc_m_exp,size=perc_exp)) +
                                         geom_point() +
                                         labs(colour="Scaled\nAvg expr.", size="% Expr.") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_color_viridis_c(),
                                       ggplot(result_celltype_2 %>% filter(protein %in% proteins_to_show), aes(x=cluster,y=factor(protein, levels = rev(proteins_to_show)),colour=avg_log2FC,size=-log10(p_val_adj))) +
                                         geom_point() +
                                         labs(colour="LogFC", size="-Log(adj-P)") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_colour_gradientn(colors=c('darkblue','#f7f7f7','darkred'), limits=c(-max(abs(result_celltype_2['avg_log2FC']), na.rm = T), max(abs(result_celltype_2['avg_log2FC']), na.rm = T))),
                                       ncol = 4, nrow = 1))
    )
    dev.off()
    
    
    pdf(file.path(data_dir, s, 'clustering_final_all_markers.pdf'), height = length(all_proteins)/8 + 0.4, width = 20)
    print(
      annotate_figure(top=s, ggarrange(addSmallLegend(ggplot(met, aes(UMAP_1, UMAP_2, color=celltype_name1)) +
                                                        geom_point(size=1) +
                                                        labs(color = 'Cluster') +
                                                        scale_color_manual(values = colours)),
                                       ggplot(result_celltype_1, aes(x=cluster,y=factor(protein, levels = protein_order_celltype_1),colour=mean_exp,size=perc_exp)) +
                                         geom_point() +
                                         labs(colour="Avg expr.", size="% Expr.") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_color_viridis_c() +
                                         scale_size_continuous(range = c(0.2,3)),
                                       ggplot(result_celltype_1, aes(x=cluster,y=factor(protein, levels = protein_order_celltype_1),colour=sc_m_exp,size=perc_exp)) +
                                         geom_point() +
                                         labs(colour="Scaled\navg. expr.", size="% Expr.") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_color_viridis_c() +
                                         scale_size_continuous(range = c(0.2,3)),
                                       ggplot(result_celltype_1, aes(x=cluster,y=factor(protein, levels = protein_order_celltype_1),colour=avg_log2FC,size=-log10(p_val_adj))) +
                                         geom_point() +
                                         labs(colour="LogFC", size="-Log(adj-P)") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_colour_gradientn(colors=c('darkblue','#f7f7f7','darkred'), limits=c(-max(abs(result_celltype_1['avg_log2FC']), na.rm = T), max(abs(result_celltype_1['avg_log2FC']), na.rm = T))) +
                                         scale_size_continuous(range = c(0.2,3)),
                                       ncol = 4, nrow = 1))
    )
    
    print(
      annotate_figure(top=s, ggarrange(addSmallLegend(ggplot(met, aes(UMAP_1, UMAP_2, color=celltype_name2)) +
                                                        geom_point(size=1) +
                                                        labs(color = 'Cluster') +
                                                        scale_color_manual(values = c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")))),
                                       ggplot(result_celltype_2, aes(x=cluster,y=factor(protein, levels = protein_order_celltype_2),colour=mean_exp,size=perc_exp)) +
                                         geom_point() +
                                         labs(colour="Avg expr.", size="% Expr.") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_color_viridis_c() +
                                         scale_size_continuous(range = c(0.2,3)),
                                       ggplot(result_celltype_2, aes(x=cluster,y=factor(protein, levels = protein_order_celltype_2),colour=sc_m_exp,size=perc_exp)) +
                                         geom_point() +
                                         labs(colour="Scaled\navg. expr.", size="% Expr.") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_color_viridis_c() +
                                         scale_size_continuous(range = c(0.2,3)),
                                       ggplot(result_celltype_2, aes(x=cluster,y=factor(protein, levels = protein_order_celltype_2),colour=avg_log2FC,size=-log10(p_val_adj))) +
                                         geom_point() +
                                         labs(colour="LogFC", size="-Log(adj-P)") +
                                         labs(y='Protein') +
                                         scale_size_continuous(limits = c(0,100), breaks=c(0,50,100)) +
                                         theme(axis.text.x = element_text(face='bold', angle = 90, hjust = 1, vjust = 0.5)) +
                                         scale_colour_gradientn(colors=c('darkblue','#f7f7f7','darkred'), limits=c(-max(abs(result_celltype_2['avg_log2FC']), na.rm = T), max(abs(result_celltype_2['avg_log2FC']), na.rm = T))) +
                                         scale_size_continuous(range = c(0.2,3)),
                                       ncol = 4, nrow = 1))
    )
    dev.off()
    
    write.csv(met, file.path(data_dir, s, 'unsupervised_meta_data.csv'))
  }
}


# Save new cluster data to Rdata
load('tapestri_data.Rdata')
for (s in all_samples) {
  met = metadatas[[s]]
  # Add unsupervised clustering information if it is available
  if (file.exists(file.path(data_dir, s, 'unsupervised_meta_data.csv'))) {
    unsupervised_clusters = read.csv(file.path(data_dir, s, 'unsupervised_meta_data.csv'), row.names = 1)
    if ('celltype_name1' %in% colnames(unsupervised_clusters)) {
      met$unsupervised_cluster_name = unsupervised_clusters[rownames(met), 'celltype_name1']
    } else {
      met$unsupervised_cluster_name = 'Unknown'
    }
  } else {
    met$unsupervised_cluster_name = 'Unknown'
  }
  met$short_unsupervised_cluster_name = gsub(":.*", "", met$unsupervised_cluster_name)
  metadatas[[s]] = met
  
  
  # For PDF of all unsupervised dotplots/umaps
  pdf_subset(file.path(data_dir, s, 'clustering_final_all_markers.pdf'), pages = 1)
  pdf_subset(file.path(data_dir, s, 'clustering_final.pdf'), pages = 1)
}
save(assays,metadatas,all_samples,files,file='tapestri_data.Rdata')

# Create PDF of all unsupervised dotplots/umaps
all_markers = paste(data_dir, all_samples, 'clustering_final_all_markers_output.pdf', sep='/')
selected_markers = paste(data_dir, all_samples, 'clustering_final_output.pdf', sep='/')
pdf_combine(all_markers, output = 'Unsupervised_analysis_dotplots_all_markers.pdf')
pdf_combine(selected_markers, output = 'Unsupervised_analysis_dotplots.pdf')
unlink(c(all_markers, selected_markers))

