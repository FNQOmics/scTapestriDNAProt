library(Matrix)
library(dplyr)
library(ggplot2)
library(rhdf5)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pdftools)
library(data.table)
library(gridExtra)
library(tidyr)
library(DEoptim)
source("~/Documents/Jerome/Kirby/Scripts/jeromefuncs.R")

theme_set(theme_classic())
colours = set1
saved_colours = colours
colours = c(colours, 'cyan', brewer.pal(8, 'Set3')[c(1,3,6,7)], "orange")
names(colours) = c("CD4-CD8- T","CD4 T","B cell","CD8 T","CD3- T (CD7-)","Myeloid","Doublet","Dead","CD4+CD8+ T", 'Unknown', "NK", "DC/Monocyte", "CD3- T (CD7+)", "Fibroblast", "CD3- T")
plot_density = function(protein, threshold=NULL, plot_dd=F, plot_conservatives=NULL, plot_mut=F) {
  g = ggplot(df, aes_string(protein))
  if (plot_dd) {
    g = g + geom_vline(mapping = aes_string(xintercept = protein), data = df %>% filter(supervised_cell_type == 'Dead'), color=colours['Dead'], linewidth = 0.1) +
      geom_vline(mapping = aes_string(xintercept = protein), data = df %>% filter(supervised_cell_type == 'Doublet'), color=colours['Doublet'], linewidth = 0.1)
  }
  if (plot_mut) {
    g = g + geom_vline(mapping = aes_string(xintercept = protein, color='Mutation'), data = df %>% filter(Mutation != "None"), linewidth = 0.1) +
      scale_color_manual(values=mut_colours)
  }
  g = g + geom_density()
  if (is.numeric(threshold)) {
    g = g +
      geom_vline(xintercept = threshold, color = 'red') +
      annotate("text", x=threshold + 0.1, y=0.5, label=threshold,hjust=0, color='red')
  }
  if (!is.null(plot_conservatives)) {
    for (i in plot_conservatives) {
      g = g +
        geom_vline(xintercept = i, color = 'red', linetype='dashed') +
        annotate("text", x=i + 0.1, y=0.6, label=i,hjust=0, color='red')
    }
  }
  return(g + theme(legend.position = 'none'))
}
plot_2d_density = function(protein1, protein2, threshold1 = NULL, threshold2 = NULL, legend = F, color = F, colour_by = 'supervised_cell_type') {
  if (color) {
    g = ggplot(df, aes_string(protein1, protein2)) +
      geom_point(aes_string(color = colour_by), size=0.5, alpha=0.8) +
      scale_color_manual(values=c(colours, mut_colours))
  } else {
    g = ggplot(df, aes_string(protein1, protein2)) +
      geom_point(color = 'grey', size=0.5)
  }
  g = g + geom_density_2d(contour_var = "ndensity", color = 'black')
  if (! legend) g = g + theme(legend.position='none')
  if (is.numeric(threshold1)) g = g + geom_vline(xintercept = threshold1, color = 'red')
  if (is.numeric(threshold2)) g = g + geom_hline(yintercept = threshold2, color = 'red')
  return(g)
}

classify_with_thresholds <- function(thresholds, df) {
  # Apply rules to classify data points
  tmp = df
  #########################
  ######### RULES #########
  #########################
  # Note that the rules of annotation will change based on what cells are being annotated, and what proteins are in the panel.
  tmp$predicted = 'Unknown'
  tmp[tmp['CD45'] >= thresholds[1] & tmp['CD19'] < thresholds[2] & tmp['CD3'] >= thresholds[3] & tmp['CD4'] < thresholds[4] & tmp['CD8'] >= thresholds[5] & tmp['CD90'] < thresholds[6], "predicted"] = 'CD8 T'
  tmp[tmp['CD45'] >= thresholds[1] & tmp['CD19'] < thresholds[2] & tmp['CD3'] >= thresholds[3] & tmp['CD4'] >= thresholds[4] & tmp['CD8'] < thresholds[5] & tmp['CD90'] < thresholds[6], "predicted"] = 'CD4 T'
  tmp[tmp['CD45'] >= thresholds[1] & tmp['CD19'] < thresholds[2] & tmp['CD3'] >= thresholds[3] & tmp['CD4'] < thresholds[4] & tmp['CD8'] < thresholds[5] & tmp['CD90'] < thresholds[6], "predicted"] = 'CD4-CD8- T'
  tmp[tmp['CD45'] >= thresholds[1] & tmp['CD19'] < thresholds[2] & tmp['CD3'] < thresholds[3] & tmp['CD4'] < thresholds[4] & tmp['CD8'] < thresholds[5] & tmp['CD90'] < thresholds[6], "predicted"] = 'CD3- T'
  tmp[tmp['CD45'] >= thresholds[1] & tmp['CD19'] >= thresholds[2] & tmp['CD3'] < thresholds[3] & tmp['CD4'] < thresholds[4] & tmp['CD8'] < thresholds[5] & tmp['CD90'] < thresholds[6], "predicted"] = 'B cell'
  tmp[tmp['CD45'] < thresholds[1] & tmp['CD19'] < thresholds[2] & tmp['CD3'] < thresholds[3] & tmp['CD4'] < thresholds[4] & tmp['CD8'] < thresholds[5] & tmp['CD90'] >= thresholds[6], "predicted"] = 'Fibroblast'
  tmp[tmp['CD45'] >= thresholds[1] & tmp['CD19'] < thresholds[2] & tmp['CD3'] >= thresholds[3] & tmp['CD4'] >= thresholds[4] & tmp['CD8'] >= thresholds[5] & tmp['CD90'] < thresholds[6], "predicted"] = 'CD4+CD8+ T'
  tmp[(tmp['CD3'] >= thresholds[3] | tmp['CD19'] >= thresholds[2]) & tmp['CD45'] < thresholds[1], "predicted"] = 'Dead'
  tmp[tmp['CD45'] >= thresholds[1] & (tmp['CD4'] >= thresholds[4] | tmp['CD8'] >= thresholds[5]) & tmp['CD3'] < thresholds[3], "predicted"] = 'Dead'
  tmp[(tmp['CD8'] >= thresholds[5] | tmp['CD4'] >= thresholds[4] | tmp['CD3'] >= thresholds[3]) & tmp['CD19'] >= thresholds[2], "predicted"] = 'Doublet'
  tmp[tmp['CD90'] >= thresholds[6] & (tmp['CD8'] >= thresholds[5] | tmp['CD4'] >= thresholds[4] | tmp['CD3'] >= thresholds[3] | tmp['CD19'] >= thresholds[2] | tmp['CD45'] >= thresholds[1]), "predicted"] = 'Doublet'
  #########################
  # Calculate accuracy
  accuracy <- mean(tmp$predicted == df$short_unsupervised_cluster_name_renamed)
  
  return(-accuracy) # Minimize negative accuracy
}

find_optimal_thresholds = function(df, proteins_of_interest, lower_quantile = 0.05, upper_quantile = 0.95, max_iterations = 500) {
  # Set bounds for thresholds (assuming features are within [0, 1])
  lower_bounds <- sapply(df[proteins_of_interest], function(feature) return(quantile(feature, lower_quantile)))
  upper_bounds <- sapply(df[proteins_of_interest], function(feature) return(quantile(feature, upper_quantile)))
  names(lower_bounds) = names(upper_bounds) = proteins_of_interest
  
  # Run Differential Evolution to find optimal thresholds
  de_result <- DEoptim(fn = classify_with_thresholds, 
                       lower = lower_bounds, upper = upper_bounds, 
                       df = df, 
                       DEoptim.control(itermax = max_iterations,
                                       strategy = 6,
                                       VTR = -1))
  return(de_result)
}

# Collect mutation data
if (file.exists("../Variant_data_from_matt/Variant_barcodes.csv")) {
  mutation_data_short = read.csv("../Variant_data_from_matt/Variant_barcodes.csv")
  mutation_data = mutation_data_short %>%
    pivot_longer(cols = which(colnames(mutation_data_short) == 'Barcodes'):ncol(mutation_data_short),
                 values_to = "Barcodes",
                 names_to = NULL) %>%
    filter(!is.null(Barcodes) & Barcodes != '')
  mutation_data$Sample = ""
} else {
  mutation_data = data.frame("Sample"=c(), "Barcodes"=c(), "Patient"=c())
}

for (p in unique(mutation_data$Patient)) {
  for (s in all_samples) {
    if (s %like% p) {
      mutation_data[mutation_data$Patient == p, "Sample"] = s
      break
    }
  }
}
mutation_data$Barcodes = gsub("-1$", "", mutation_data$Barcodes)
mut_colours = c(saved_colours, brewer.pal(10, "Set3")[c(1:8,10:12)], brewer.pal(8,"Set2"))[0:length(unique(mutation_data$Mutation))]
names(mut_colours) = sort(unique(mutation_data$Mutation))
mut_colours = c(mut_colours[1:length(sort(unique(mutation_data$Mutation)))], "None" = 'lightgrey', "Multiple"='cyan')
mut_colours = mut_colours[!is.na(mut_colours)]

load('tapestri_data.Rdata')
thresholds = annotation_accuracy = list()
for (s in all_samples) {
  print(s)
  df = as.data.frame(assays[[s]])
  met = metadatas[[s]]
  met = cbind(met[,! colnames(met) %in% colnames(df)], df[rownames(met),])
  
  proteins_of_interest = c("CD19","CD3","CD4","CD45","CD8", "CD90")
  proteins_of_interest = proteins_of_interest[proteins_of_interest %in% colnames(df)]
  
  met$short_unsupervised_cluster_name_renamed = plyr::revalue(met$short_unsupervised_cluster_name, replace = c(
    "CD8+CD4+" = "CD8+CD4+ T",
    "B Cell/Plasma" = "B cell",
    "B cells" = "B cell",
    "unsure" = "Unknown",
    "Gamma delta" = "Unknown",
    "Gamma Delta" = "Unknown",
    "Plasma" = "B cell",
    "plasma" = "B cell",
    "NK" = "Unknown",
    "Rogue" = "CD4-CD8- T",
    "rogue" = "CD4-CD8- T",
    "DC/Monocytes"="Unknown",
    "DC/Monocyte"="Unknown",
    "CD3- T (CD7-)"="CD3- T",
    "CD3- T (CD7+)"="CD3- T",
    "mix" = "Unknown",
    "Mix" = "Unknown",
    "CD4 T: CD3 high" = "CD4 T",
    "CD4 T: CD3 int" = "CD4 T",
    "CD8 T: CD3 high" = "CD8 T",
    "CD8 T: CD3 int" = "CD8 T",
    "CD4 T: CD7−"="CD4 T",
    "CD4 T: CD7+"="CD4 T",
    "CD8 T: CD3 high"="CD8 T",
    "CD8 T: CD3 intermeidate"="CD8 T",
    "Myeloid"="Unknown"
  ))
  
  # Find optimal thresholds which match unsupervised annotation
  de_results = find_optimal_thresholds(met, proteins_of_interest)
  thresholds[[s]] = optimal_thresholds = round(de_results$optim$bestmem,3)
  annotation_accuracy[[s]] = threshold_annotation_accuracy = de_results$optim$bestval
  
  df[paste0(proteins_of_interest,'_cat')] = "-"
  for (p in proteins_of_interest) {
    df[df[[p]] >= optimal_thresholds[p], paste0(p, '_cat')] = "+"
  }

  df$supervised_cell_type = "Unknown"
  df[df$CD45_cat == '+' & df$CD19_cat == '-' & df$CD3_cat == '+' & df$CD4_cat == '-' & df$CD8_cat == '+' & df$CD90_cat == '-', "supervised_cell_type"] = 'CD8 T'
  df[df$CD45_cat == '+' & df$CD19_cat == '-' & df$CD3_cat == '+' & df$CD4_cat == '+' & df$CD8_cat == '-' & df$CD90_cat == '-', "supervised_cell_type"] = 'CD4 T'
  df[df$CD45_cat == '+' & df$CD19_cat == '-' & df$CD3_cat == '+' & df$CD4_cat == '-' & df$CD8_cat == '-' & df$CD90_cat == '-', "supervised_cell_type"] = 'CD4-CD8- T'
  df[df$CD45_cat == '+' & df$CD19_cat == '-' & df$CD3_cat == '-' & df$CD4_cat == '-' & df$CD8_cat == '-' & df$CD90_cat == '-', "supervised_cell_type"] = 'CD3- T'
  df[df$CD45_cat == '+' & df$CD19_cat == '+' & df$CD3_cat == '-' & df$CD4_cat == '-' & df$CD8_cat == '-' & df$CD90_cat == '-', "supervised_cell_type"] = 'B cell'
  df[df$CD45_cat == '-' & df$CD19_cat == '-' & df$CD3_cat == '-' & df$CD4_cat == '-' & df$CD8_cat == '-' & df$CD90_cat == '+', "supervised_cell_type"] = 'Fibroblast'
  df[df$CD45_cat == '+' & df$CD19_cat == '-' & df$CD3_cat == '+' & df$CD4_cat == '+' & df$CD8_cat == '+' & df$CD90_cat == '-', "supervised_cell_type"] = 'CD4+CD8+ T'
  df[(df$CD3_cat == '+' | df$CD19_cat == '+') & df$CD45_cat == '-', "supervised_cell_type"] = 'Dead'
  df[df$CD45_cat == '+' & (df$CD4_cat == '+' | df$CD8_cat == '+') & df$CD3_cat == '-', "supervised_cell_type"] = 'Dead'
  df[(df$CD8_cat == '+' | df$CD4_cat == '+' | df$CD3_cat == '+') & df$CD19_cat == '+', "supervised_cell_type"] = 'Doublet'
  df[df$CD90_cat == '+' & (df$CD8_cat == '+' | df$CD4_cat == '+' | df$CD3_cat == '+' | df$CD19_cat == '+' | df$CD45_cat == '+'), "supervised_cell_type"] = 'Doublet'
  
  cat_columns = paste0(proteins_of_interest,'_cat')
  df$markers = do.call(paste0, as.data.frame(t(apply(as.matrix(df[,cat_columns]), 1, function(x) x=paste0(gsub("_cat", "", cat_columns), x)))))
  plots = list()
  for (x1 in proteins_of_interest) {
    plots[[x1]] = plot_density(x1, threshold = optimal_thresholds[x1])
    for (x2 in proteins_of_interest) {
      if (x1 != x2) {
        plots[[paste(x1,x2)]] = plot_2d_density(x1,x2,optimal_thresholds[x1],optimal_thresholds[x2],color=T)
      }
    }
    plots[[paste0(x1,'dd')]] = plot_density(x1, threshold = optimal_thresholds[x1], plot_dd = T)
  }
  
  # Add any mutation data to metadata and create plots
  plots_mut = c()
  if (s %in% mutation_data$Sample) {
    sample_mutation_data = mutation_data %>% filter(Sample == s)
    # 
    # # Special sample where only Matt's pipeline picks up the variant and not the export NGT.csv
    # if (s == "2088_NA_Gut_LymphEPCAM_Tapestri_r1_bNA") {
    #   variants = read.delim(file.path("Data", s, 'data/I659L_barcode_alleles.tsv'), header = F)
    #   rownames(variants) = gsub("-1", "", variants$V1)
    #   df$Variant_STAT3_I659L = plyr::revalue(variants[rownames(df),"V2"], replace = c('0/0'="WT", '0/1'="HET", '0/2'="WT", './.'="No_call"))
    #   table(df$Variant_STAT3_I659L)
    # }
    
    if (file.exists(file.path("Data", s, 'data/export/NGT.csv'))) {
      mutations_to_look_for = unique(mutation_data %>% filter(Sample == s) %>% pull(Mutation))
      mut_alleles = read.csv(file.path("Data", s, 'data/export/NGT.csv'))
      mut_alleles$Sample = s
      mut_alleles$CellID = paste(mut_alleles$Sample, gsub("-1", "", mut_alleles$Cell), sep = '.')
      variants = read.csv(file.path("Data", s, 'data/export/Variants.csv'))
      mutations_modified = gsub(" ", ":p.", mutations_to_look_for)
      names(mutations_modified) = mutations_to_look_for
      for (m in names(mutations_modified)) {
        chrom_location = make.names(variants[variants$Protein %like% mutations_modified[[m]] | variants$Protein == mutations_modified[[m]], 'Variant'])
        if (length(chrom_location) > 0) {
          df[unlist(lapply(strsplit(mut_alleles$CellID, split='[.]'), function(x) {x[2]})), make.names(paste0("Variant_", gsub(" ", "_", m)))] = as.character(mut_alleles[,make.names(chrom_location)])
          df[,make.names(paste0("Variant_", gsub(" ", "_", m)))] = plyr::revalue(df[,make.names(paste0("Variant_", gsub(" ", "_", m)))], replace = c('0'="WT", '1'="HET", '2'="HOM", '3'="No_call"))
        } else {
          print(paste("Problem with",s,"-- mutation:",m))
        }
      }
    }
    
    df = merge(df, sample_mutation_data[c("Barcodes", "Mutation")], by.x='row.names', by.y = 'Barcodes', all.x = T, all.y=F)
    df[df$Row.names %in% df[duplicated(df$Row.names),"Row.names"],"Mutation"] = "Multiple"
    df = df[!duplicated(df$Row.names),]
    rownames(df) = df$Row.names
    df$Row.names = NULL
    df[is.na(df$Mutation) | df$Mutation == '', "Mutation"] = "None"
    
    for (x1 in proteins_of_interest) {
      plots_mut[[x1]] = plot_density(x1, threshold = optimal_thresholds[x1])
      for (x2 in proteins_of_interest) {
        if (x1 != x2) {
          plots_mut[[paste(x1,x2)]] = plot_2d_density(x1,x2,optimal_thresholds[x1],optimal_thresholds[x2],color=T, colour_by = 'Mutation')
        }
      }
      plots_mut[[paste0(x1,'dd')]] = plot_density(x1, threshold = optimal_thresholds[x1], plot_mut = T)
    }
  } else {
    df$Mutation = 'None'
  }
  
  met[colnames(df)] = NULL
  met = cbind(met, df[rownames(met),])
  
  pdf("~/tmp1.pdf", height = 4, width = 3)
  grid.table(data.frame(table(df$supervised_cell_type)) %>% dplyr::rename(supervised_cell_type='Var1'))
  dev.off()
  
  pdf("~/tmp2.pdf", height = 12, width = 5)
  num_dead_doublet_unknown = length(which(df$supervised_cell_type %in% c('Dead','Doublet','Unknown')))
  grid.table(rbind(data.frame(head(sort(table(df[df$supervised_cell_type %in% c('Dead', 'Doublet', 'Unknown'),'markers']),decreasing = T), 40)), data.frame('Var1'='SUM', 'Freq'=num_dead_doublet_unknown)) %>% dplyr::rename('dead/unknown/doublet markers'='Var1'))
  dev.off()
  
  pdf("~/tmp3.pdf", height = 4, width = 10)
  grid.table(as.data.frame(Matrix(xtabs(~supervised_cell_type + Mutation, df))))
  dev.off()
  
  pdf("~/tmp4.pdf", height = length(proteins_of_interest) * 3, width = length(proteins_of_interest) * 3)
  print(ggarrange(plotlist = plots,
                  nrow = length(proteins_of_interest),
                  ncol = length(proteins_of_interest) + 1,
                  common.legend = T,
                  label.y = 1.05,
                  legend.grob = get_legend(ggplot(df, aes(supervised_cell_type, supervised_cell_type, color=supervised_cell_type)) +
                                             geom_point() +
                                             theme(legend.position = 'top') +
                                             scale_colour_manual(values=colours) +
                                             labs(colour="Cell type"))))
  if (length(plots_mut) > 0) {
    print(ggarrange(plotlist = plots_mut,
                    nrow = length(proteins_of_interest),
                    ncol = length(proteins_of_interest) + 1,
                    common.legend = T,
                    label.y = 1.05,
                    legend.grob = get_legend(ggplot(df, aes(Mutation, Mutation, color=Mutation)) +
                                               geom_point() +
                                               theme(legend.position = 'top') +
                                               scale_colour_manual(values=mut_colours) +
                                               labs(colour="Mutation"))))
  }
  
  plots2 = c()
  plots3 = c()
  for (p in names(optimal_thresholds)) {
    plots2[[p]] = ggplot() +
      geom_point(data = met[met[[paste0(p, "_cat")]] == '-',], mapping = aes_string('UMAP_1', 'UMAP_2', color=p)) +
      geom_point(data = met[met[[paste0(p, "_cat")]] == '+',], mapping = aes_string('UMAP_1', 'UMAP_2', fill=p), pch=21, color=brewer.pal(8, "Purples")[7]) +
      theme(axis.ticks = element_blank(), axis.text = element_blank()) +
      scale_fill_gradientn(colours=brewer.pal(8, "Purples"), limits=c(min(met[[p]], na.rm = T), max(met[[p]], na.rm = T))) +
      scale_color_gradientn(colours=brewer.pal(8, "Purples"), limits=c(min(met[[p]], na.rm = T), max(met[[p]], na.rm = T)))
    
    
    plots3[[p]] = ggplot(met, aes_string('UMAP_1', 'UMAP_2', color=p)) +
      geom_point(size=1) +
      theme(axis.ticks = element_blank(), axis.text = element_blank()) +
      scale_color_gradientn(colours=brewer.pal(8, "Purples"))
  }
  egg::ggarrange(plots = plots3, ncol = 2, nrow = 3)
  egg::ggarrange(plots = plots2, ncol = 2, nrow = 3)
  dev.off()
  
  if (length(plots_mut) > 0) {
    pdf("~/tmp5.pdf", height = 14, width = 17, onefile = F)
    egg::ggarrange(
      ggplot(met, aes(UMAP_1, UMAP_2, color=supervised_cell_type)) + geom_point() + scale_color_manual(values=colours) + theme(axis.ticks = element_blank(), axis.text = element_blank()),
      ggplot(met, aes(UMAP_1, UMAP_2, color=Mutation)) + geom_point() + scale_color_manual(values=mut_colours) + theme(axis.ticks = element_blank(), axis.text = element_blank()),
      ggplot(met, aes(UMAP_1, UMAP_2, color=short_unsupervised_cluster_name)) + geom_point() + scale_color_manual(values=colours) + theme(axis.ticks = element_blank(), axis.text = element_blank()),
      ncol = 2, nrow = 2)
    dev.off()
  } else {
    pdf("~/tmp5.pdf", height = 14, width = 20, onefile = F)
    egg::ggarrange(
      ggplot(met, aes(UMAP_1, UMAP_2, color=supervised_cell_type)) + geom_point() + scale_color_manual(values=colours) + theme(axis.ticks = element_blank(), axis.text = element_blank()),
      ggplot(met, aes(UMAP_1, UMAP_2, color=markers)) + geom_point() + theme(axis.ticks = element_blank(), axis.text = element_blank()),
      ggplot(met, aes(UMAP_1, UMAP_2, color=short_unsupervised_cluster_name)) + geom_point() + scale_color_manual(values=colours) + theme(axis.ticks = element_blank(), axis.text = element_blank()),
      ncol = 2, nrow = 2)
    dev.off()
  }
  
  pdf("~/tmp6.pdf", height = 6, width = 12)
  if (length(unique(met$Mutation)) > 1) {
    long_met = rbind(data.frame('CellType' = met$supervised_cell_type, 'Mutation'=met$Mutation, 'ClusterType'="Supervised"), data.frame('CellType' = met$short_unsupervised_cluster_name, 'Mutation'=met$Mutation, 'ClusterType'="Unsupervised"))
    
    print(
      ggpubr::ggarrange(
        ggpubr::ggarrange(
          ggplot(long_met, aes(x=CellType, fill=Mutation)) + geom_bar() + facet_wrap(~ClusterType, scales = 'free_x') + scale_fill_manual(values=c(mut_colours, colours)) + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)),
          ggplot(long_met %>% filter(Mutation != 'None'), aes(x=CellType, fill=Mutation)) + geom_bar() + facet_wrap(~ClusterType, scales = 'free_x') + scale_fill_manual(values=c(mut_colours, colours)) + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)),
          ncol = 2, nrow = 1,
          common.legend = T,
          legend.grob = get_legend(ggplot(met, aes(x=Mutation,fill=Mutation)) +
                                     geom_bar() +
                                     theme(legend.position = 'top') +
                                     scale_fill_manual(values=mut_colours) +
                                     labs(colour="Mutation"))),
        ggpubr::ggarrange(
          ggplot(long_met, aes(x=Mutation, fill=CellType)) + geom_bar() + facet_wrap(~ClusterType, scales = 'free_x') + scale_fill_manual(values=c(mut_colours, colours)) + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)),
          ggplot(long_met %>% filter(Mutation != 'None'), aes(x=Mutation, fill=CellType)) + geom_bar() + facet_wrap(~ClusterType, scales = 'free_x') + scale_fill_manual(values=c(mut_colours, colours)) + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)),
          ncol = 2, nrow = 1,
          common.legend = T,
          legend = 'right',
          legend.grob = get_legend(addSmallLegend(ggplot(long_met, aes(x=CellType,fill=CellType)) +
                                                    geom_bar() +
                                                    theme(legend.position = 'right') +
                                                    scale_fill_manual(values=colours) +
                                                    labs(colour="CellType")))),
        ncol = 2, nrow = 1)
    )
  }
  dev.off()
  
  pdf("~/tmp7.pdf", height = 6, width = 7)
  print(ggtexttable(as.data.frame(Matrix(xtabs(~supervised_cell_type + short_unsupervised_cluster_name_renamed, met)))))
  dev.off()
  
  # SAVE OUTPUT
  pdfs = c("~/tmp1.pdf", "~/tmp2.pdf","~/tmp3.pdf", "~/tmp4.pdf", "~/tmp5.pdf", "~/tmp6.pdf", "~/tmp7.pdf")
  pdftools::pdf_combine(input = pdfs, output=file.path('Data', s, paste0(s,'_protein_thresholds.pdf')))
  unlink(pdfs)
  met$orig.ident = NULL
  write.csv(met %>% select(supervised_cell_type,markers,Mutation,!!!syms(proteins_of_interest), !!!syms(paste0(proteins_of_interest,'_cat')), everything()), file.path('Data', s, paste0(s,'_tapestri_data.csv')))
}


