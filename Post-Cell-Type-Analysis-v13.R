
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ##      ###     ##       ##        ##        ##      ##   ########    ##     #######                               
##      ####     ####    ##      ####       ##         ##    ##   ##           ##    ##                                               
##     ##  ##    ## ##   ##     ##  ##      ##          ##  ##    ##           ##    ##           
##    ##    ##   ##  ##  ##    ##    ##     ##           ####      #######     ##     #######                                           
##    ########   ##   ## ##   ##########    ##            ##             ##    ##           ##                                      
##    ##    ##   ##    ####   ##      ##    ##            ##             ##    ##           ##                                          
##    ##    ##   ##     ###   ##      ##    ########      ##       #######     ##     #######                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#GOAL:  TCC processed data

{

setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 scRNAseq - canine")

library("Seurat")
library("ggplot2")

########################################################################################################################################
#TEST SAVED FILE
########################################################################################################################################


load("snRNAseq_Tux_cellTypes.RData")

colors<-c("blue","skyblue","darkcyan")

DimPlot(snRNAseq_organd_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid bldder cancer")+
  scale_color_manual(values=colors)




colors<-c("darkblue","blue","skyblue","darkcyan","darkred","black","darkgreen","grey40","skyblue")
# Visualize the assigned cell types in the UMAP plot
DimPlot(snRNAseq_invivo_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo bladder cancer")+
  scale_color_manual(values=colors)

########################################################################################################################################
#TEST SAVED FILE
########################################################################################################################################


library("Seurat")
library("ggplot2")


# Plot the clusters using UMAP
DimPlot(snRNAseq_invivo_cellTypes, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("in vivo Tux snRNAseq")



colors<-c("darkblue","blue","skyblue","darkcyan","darkred","black","darkgreen","grey40","skyblue")
# Visualize the assigned cell types in the UMAP plot
DimPlot(snRNAseq_invivo_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo Tux lung")+
  scale_color_manual(values=colors)





############################################################
#STANDARD SEURAT FIGS
############################################################

library(Seurat)
library(tidyverse)

seurat_obj<-snRNAseq_organd_cellTypes

# Set "cell_type" as the grouping identity
Idents(seurat_obj) <- seurat_obj@meta.data$cell_type

# Exclude genes containing "ENSGAL"
valid_genes <- rownames(seurat_obj)[!grepl("ENSCAF", rownames(seurat_obj))]

# Subset Seurat object to exclude these genes
seurat_obj <- subset(seurat_obj, features = valid_genes)

# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(top_genes_per_cell_type, "organoid_TopGenes.csv", row.names = FALSE)

# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("organoid_DotPlot.pdf", plot = dotplot)

# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("Expression of", gene)) +
    theme_minimal()
})

pdf("organoid_FeaturePlots.pdf", width = 10, height = 10)
for (plot in feature_plots) {
  print(plot)
}
dev.off()

# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = top_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("organoid_Heatmap.pdf", plot = heatmap)


#############################
# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
#############################
library(patchwork) # For combining plots into a grid


top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

# Create individual feature plots for each gene
feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(gene) +
    theme_minimal()
})

# Combine all feature plots into a 5x2 grid
combined_plot <- wrap_plots(feature_plots, ncol = 2)

# Save the combined plot as a single PDF
ggsave("organoid_FeaturePlots.pdf", plot = combined_plot, width = 15, height = 8)

}



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ##      ###     ##       ##        ##        ##      ##   ########    ##     #######                               
##      ####     ####    ##      ####       ##         ##    ##   ##           ##    ##                                               
##     ##  ##    ## ##   ##     ##  ##      ##          ##  ##    ##           ##    ##           
##    ##    ##   ##  ##  ##    ##    ##     ##           ####      #######     ##     #######                                           
##    ########   ##   ## ##   ##########    ##            ##             ##    ##           ##                                      
##    ##    ##   ##    ####   ##      ##    ##            ##             ##    ##           ##                                          
##    ##    ##   ##     ###   ##      ##    ########      ##       #######     ##     #######                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#GOAL:  standard seurat figures for B816 Kidney scRNAseq:  in vivo and organoid

{

setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 scRNAseq - canine")

library("Seurat")
library("ggplot2")

########################################################################################################################################
#TEST SAVED FILE
########################################################################################################################################


load("B816_Kidney_cellTypes.RData")

colors<-c("skyblue","blue","darkblue")

# Visualize the assigned cell types in the UMAP plot
DimPlot(B816_Kidney_org_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid Kidney cancer")+
  scale_color_manual(values=colors)



colors<-c("darkred","royalblue", "dodgerblue4", "blue3", "blue", "navy", "royalblue3","black","grey")
# Visualize the assigned cell types in the UMAP plot
DimPlot(B816_Kidney_invivo_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo Kidney cancer")+
  scale_color_manual(values=colors)


############################################################
#STANDARD SEURAT FIGS
############################################################

library(Seurat)
library(tidyverse)

seurat_obj<-B816_Kidney_org_cellTypes

# Set "cell_type" as the grouping identity
Idents(seurat_obj) <- seurat_obj@meta.data$cell_type

# Exclude genes containing "ENSGAL"
valid_genes <- rownames(seurat_obj)[!grepl("ENSCAF", rownames(seurat_obj))]

# Subset Seurat object to exclude these genes
seurat_obj <- subset(seurat_obj, features = valid_genes)

# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(top_genes_per_cell_type, "B816_Kidney-organoid_02 TopGenes.csv", row.names = FALSE)

# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Kidney-organoid_03 DotPlot.pdf", plot = dotplot)

# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("Expression of", gene)) +
    theme_minimal()
})

pdf("B816_Kidney-organoid_04 Feature Plot.pdf", width = 10, height = 10)
for (plot in feature_plots) {
  print(plot)
}
dev.off()

# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Kidney-organoid_05 Heatmap.pdf", plot = heatmap)


#############################
# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
#############################
library(patchwork) # For combining plots into a grid


top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

# Create individual feature plots for each gene
feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(gene) +
    theme_minimal()
})

# Combine all feature plots into a 5x2 grid
combined_plot <- wrap_plots(feature_plots, ncol = 2)

# Save the combined plot as a single PDF
ggsave("B816_Kidney-organoid_04 Feature Plot.pdf", plot = combined_plot, width = 15, height = 8)

}


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ##      ###     ##       ##        ##        ##      ##   ########    ##     #######                               
##      ####     ####    ##      ####       ##         ##    ##   ##           ##    ##                                               
##     ##  ##    ## ##   ##     ##  ##      ##          ##  ##    ##           ##    ##           
##    ##    ##   ##  ##  ##    ##    ##     ##           ####      #######     ##     #######                                           
##    ########   ##   ## ##   ##########    ##            ##             ##    ##           ##                                      
##    ##    ##   ##    ####   ##      ##    ##            ##             ##    ##           ##                                          
##    ##    ##   ##     ###   ##      ##    ########      ##       #######     ##     #######                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#GOAL:  standard seurat figures for other B816 tissues

{

setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 scRNAseq - canine")

library("Seurat")
library("ggplot2")

####################################################################################################################################################################################################################################################################
#TEST SAVED FILE
####################################################################################################################################################################################################################################################################


load("B816_other_cellTypes.RData")



colors<-c("blue","skyblue","darkblue","purple","grey")

# Visualize the assigned cell types in the UMAP plot
DimPlot(B816_Pancreas_org_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid pancreatic cancer")+
  scale_color_manual(values=colors)



############################################################
#STANDARD SEURAT FIGS
############################################################

library(Seurat)
library(tidyverse)

seurat_obj<-B816_Pancreas_org_cellTypes

# Set "cell_type" as the grouping identity
Idents(seurat_obj) <- seurat_obj@meta.data$cell_type

# Exclude genes containing "ENSGAL"
valid_genes <- rownames(seurat_obj)[!grepl("ENSCAF", rownames(seurat_obj))]

# Subset Seurat object to exclude these genes
seurat_obj <- subset(seurat_obj, features = valid_genes)

# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(top_genes_per_cell_type, "B816_Pancreas-organoid_02 TopGenes.csv", row.names = FALSE)

# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Pancreas-organoid_03 DotPlot.pdf", plot = dotplot)

# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("Expression of", gene)) +
    theme_minimal()
})

pdf("B816_Pancreas-organoid_04 Feature Plot.pdf", width = 10, height = 10)
for (plot in feature_plots) {
  print(plot)
}
dev.off()

# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Pancreas-organoid_05 Heatmap.pdf", plot = heatmap)


#############################
# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
#############################
library(patchwork) # For combining plots into a grid


top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

# Create individual feature plots for each gene
feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(gene) +
    theme_minimal()
})

# Combine all feature plots into a 5x2 grid
combined_plot <- wrap_plots(feature_plots, ncol = 2)

# Save the combined plot as a single PDF
ggsave("B816_Pancreas-organoid_04 Feature Plot.pdf", plot = combined_plot, width = 15, height = 8)




####################################################################################################################################################################################################################################################################
#TEST SAVED FILE
####################################################################################################################################################################################################################################################################


colors<-c("blue","purple","blue3","skyblue","darkblue","darkcyan")

# Visualize the assigned cell types in the UMAP plot
DimPlot(B816_lung_org_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid lung cancer")+
  scale_color_manual(values=colors)




############################################################
#STANDARD SEURAT FIGS
############################################################

library(Seurat)
library(tidyverse)

seurat_obj<-B816_lung_org_cellTypes

# Set "cell_type" as the grouping identity
Idents(seurat_obj) <- seurat_obj@meta.data$cell_type

# Exclude genes containing "ENSGAL"
valid_genes <- rownames(seurat_obj)[!grepl("ENSCAF", rownames(seurat_obj))]

# Subset Seurat object to exclude these genes
seurat_obj <- subset(seurat_obj, features = valid_genes)

# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(top_genes_per_cell_type, "B816_Lung-organoid_02 TopGenes.csv", row.names = FALSE)

# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Lung-organoid_03 DotPlot.pdf", plot = dotplot)

# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("Expression of", gene)) +
    theme_minimal()
})

pdf("B816_Lung-organoid_04 Feature Plot.pdf", width = 10, height = 10)
for (plot in feature_plots) {
  print(plot)
}
dev.off()

# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Lung-organoid_05 Heatmap.pdf", plot = heatmap)


#############################
# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
#############################
library(patchwork) # For combining plots into a grid


top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

# Create individual feature plots for each gene
feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(gene) +
    theme_minimal()
})

# Combine all feature plots into a 5x2 grid
combined_plot <- wrap_plots(feature_plots, ncol = 2)

# Save the combined plot as a single PDF
ggsave("B816_Lung-organoid_04 Feature Plot.pdf", plot = combined_plot, width = 15, height = 8)













####################################################################################################################################################################################################################################################################
#TEST SAVED FILE
####################################################################################################################################################################################################################################################################



colors<-c("blue","skyblue","darkblue","purple")

# Visualize the assigned cell types in the UMAP plot
DimPlot(B816_Bladder_org_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid bldder cancer")+
  scale_color_manual(values=colors)


############################################################
#STANDARD SEURAT FIGS
############################################################

library(Seurat)
library(tidyverse)

seurat_obj<-B816_Bladder_org_cellTypes

# Set "cell_type" as the grouping identity
Idents(seurat_obj) <- seurat_obj@meta.data$cell_type

# Exclude genes containing "ENSGAL"
valid_genes <- rownames(seurat_obj)[!grepl("ENSCAF", rownames(seurat_obj))]

# Subset Seurat object to exclude these genes
seurat_obj <- subset(seurat_obj, features = valid_genes)

# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(top_genes_per_cell_type, "B816_Bladder-organoid_02 TopGenes.csv", row.names = FALSE)

# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Bladder-organoid_03 DotPlot.pdf", plot = dotplot)

# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("Expression of", gene)) +
    theme_minimal()
})

pdf("B816_Bladder-organoid_04 Feature Plot.pdf", width = 10, height = 10)
for (plot in feature_plots) {
  print(plot)
}
dev.off()

# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Bladder-organoid_05 Heatmap.pdf", plot = heatmap)


#############################
# (3) Grid of feature plots in UMAP space showing expression of top gene for each cell type
#############################
library(patchwork) # For combining plots into a grid


top_gene_by_cell_type <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  pull(gene)

# Create individual feature plots for each gene
feature_plots <- lapply(top_gene_by_cell_type, function(gene) {
  FeaturePlot(seurat_obj, features = gene) +
    ggtitle(gene) +
    theme_minimal()
})

# Combine all feature plots into a 5x2 grid
combined_plot <- wrap_plots(feature_plots, ncol = 2)

# Save the combined plot as a single PDF
ggsave("B816_Bladder-organoid_04 Feature Plot.pdf", plot = combined_plot, width = 15, height = 8)

}

####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   ##         ########        ##       ###     ##                                  
##     ##          ##         ##             ####      ####    ##                             
##    ##           ##         ##            ##  ##     ## ##   ##                                             
##    ##           ##         #######      ##    ##    ##  ##  ##                                  
##    ##           ##         ##          ##########   ##   ## ##                                   
##     ##          ##         ##          ##      ##   ##    ####                                   
##      ########   ########   ########    ##      ##   ##     ###                                   
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ##       ##  ###      ###        ##        #########                                                                                                     
##       ##       ##  ####    ####       ####       ##      ##                                                   
##       ##       ##  ## ##  ## ##      ##  ##      ##      ##                                                   
##       ##       ##  ##  ####  ##     ##    ##     ##      ##                                                                       
##       ##       ##  ##   ##   ##    ##########    #########                                                                  
##        ##     ##   ##        ##   ##        ##   ##                                                         
##         #######    ##        ##  ##          ##  ##                                                                                                               
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  
#NOTE:   make multi-tissue umap



setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 scRNAseq - canine")

library("Seurat")
library("ggplot2")
library(dplyr)

####################################################################################################################################################################################################################################################################
#CLEANUP ALL DATASETS
####################################################################################################################################################################################################################################################################


load("B816_other_cellTypes.RData")

load("B816_Kidney_cellTypes.RData")

################################
#CLEAN PANCREAS:
################################
B816_Pancreas_org_cellTypes$Tissue<-"Pancreas"

DimPlot(B816_Pancreas_org_cellTypes, reduction = "umap", group.by = "cell_type")
seurat1<-subset(B816_Pancreas_org_cellTypes, subset = cell_type != "Low quality")
seurat1<-subset(seurat1, subset = cell_type != "Epithel (beta)")
seurat1@meta.data$UMAP_1<-seurat1@reductions$umap@cell.embeddings[,1]
seurat1<-subset(seurat1, subset = UMAP_1 >(-4.5))

DimPlot(seurat1, reduction = "umap", group.by = "cell_type")



################################
#CLEAN PANCREAS:
################################
B816_Bladder_org_cellTypes$Tissue<-"Bladder"
DimPlot(B816_Bladder_org_cellTypes, reduction = "umap", group.by = "cell_type")
seurat2<-subset(B816_Bladder_org_cellTypes, subset = cell_type != "Epithel (inflamm)")
DimPlot(seurat2, reduction = "umap", group.by = "cell_type")




################################
#CLEAN Kidney:
################################
B816_Kidney_org_cellTypes$Tissue<-"Kidney"
DimPlot(B816_Kidney_org_cellTypes, reduction = "umap", group.by = "cell_type")
seurat3<-B816_Kidney_org_cellTypes




################################
#CLEAN Lung:
################################
B816_lung_org_cellTypes$Tissue<-"Lung"
DimPlot(B816_lung_org_cellTypes, reduction = "umap", group.by = "cell_type")
seurat4<-subset(B816_lung_org_cellTypes, subset = cell_type != "Epithel (neuro?)")
seurat4<-subset(seurat4, subset = cell_type != "Epithel (basal)")

DimPlot(seurat4, reduction = "umap", group.by = "cell_type")




####################################################################################################################################################################################################################################################################
#SIMPLE MERGE DATA:   Run integration on merged object using Seurat5
####################################################################################################################################################################################################################################################################


# Merge the Seurat objects while keeping raw normalized counts
seurat_merged <- merge(seurat1, y = list(seurat2, seurat3, seurat4), add.cell.ids = c("T1", "T2", "T3", "T4"))

# Normalize data (if not already normalized)
seurat_merged <- NormalizeData(seurat_merged)

# Find variable features
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)

# Scale the data and run PCA
seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
seurat_merged <- RunPCA(seurat_merged, npcs = 30, verbose = FALSE)

# Perform clustering and run UMAP
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:30)

# UMAP visualization
DimPlot(seurat_merged, reduction = "umap", group.by = "seurat_clusters")  # By cluster
DimPlot(seurat_merged, reduction = "umap", group.by = "Tissue") +
  ggtitle("B816 organoids: tissue") 

merge_types<-paste0(seurat_merged$Tissue,gsub("Epithel","",seurat_merged$cell_type))
seurat_merged$type2<-merge_types


DimPlot(seurat_merged, reduction = "umap", group.by = "type2")+
  ggtitle("B816 organoids: sub-types")




####################################################################################################################################################################################################################################################################
#FINALIZE FOR FIGURES AND APPS:
####################################################################################################################################################################################################################################################################
#confirm can subset by tissue
tmp<-subset(seurat_merged, subset = Tissue == "Bladder")

DimPlot(tmp, reduction = "umap", group.by = "type2")+
  ggtitle("B816 organoids: bladder")

seurat_filtered <- seurat_merged[!grepl("ENSCAFG", rownames(seurat_merged)), ]

saveRDS(seurat_filtered, file = "FINAL_B816-organoids.rds")

library(qs)
qsave(seurat_filtered, "FINAL_B816-organoids.qs")
#save(seurat_merged,file="FINAL_B816-organoids.RData")


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   ##         ########        ##       ###     ##                                  
##     ##          ##         ##             ####      ####    ##                             
##    ##           ##         ##            ##  ##     ## ##   ##                                             
##    ##           ##         #######      ##    ##    ##  ##  ##                                  
##    ##           ##         ##          ##########   ##   ## ##                                   
##     ##          ##         ##          ##      ##   ##    ####                                   
##      ########   ########   ########    ##      ##   ##     ###                                   
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     #######    #######   ###     ##    ##########                                 
##     ##     ##     ##     ##     ##          ##          ##        ####    ##        ##                      
##     ##     ##     ##    ##      ##          ##          ##        ## ##   ##        ##                       
##     ########      #######       ########     #######    #######   ##  ##  ##        ##                               
##     ##            ##    ##      ##                 ##   ##        ##   ## ##        ##                       
##     ##            ##     ##     ##                 ##   ##        ##    ####        ##                       
##     ##            ##      ##    ########     #######    #######   ##     ###        ##                                
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#FIGURES FOR PAPER (Seurat generic)


setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 scRNAseq - canine")

library("Seurat")
library("ggplot2")
library(dplyr)

#SAVE IN FORMS THAT LOAD FASTER
#load("FINAL_B816-organoids.RData")
#saveRDS(seurat_merged, file = "FINAL_B816-organoids.rds")
#qsave(seurat_merged, "FINAL_B816-organoids.qs")
seurat_merged <- qs::qread("FINAL_B816-organoids.qs")


all_cols<-scales::hue_pal()(13)
bladder_cols<-c("#F8766D","#E18A00","#BE9C00")
kidney_cols<-c("#8CAB00","#24B700","#00BE70")
lung_cols<-c("#00C1AB","#00BBDA","#8B93FF")
pancreas_cols<-c("#D575FE","#F962DD","#FF65AC")



umap_plot <- function(seurat_merged,color){
    DimPlot(seurat_merged, reduction = "umap", group.by = "type2") +
      ggtitle("B816 organoids: sub-types") + 
      theme(legend.position = "bottom") +
      scale_color_manual(values = color)
}

generate_feature_plot <- function(seurat_merged,gene="TOP2A") {
  #req(seurat_merged, gene)
  FeaturePlot(seurat_merged, features = gene) +
    scale_color_gradient(low = "grey", high = "red", name = paste0(gene, " (logCounts)")) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 18)
    )
}

####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#MERGED CELL FIGS:
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################

umap_plot(seurat_merged,color=c(bladder_cols,kidney_cols,lung_cols,pancreas_cols))
generate_feature_plot (seurat_merged,gene="TOP2A")

DimPlot(seurat_merged, reduction = "umap", group.by = "type2") +
  ggtitle("B816 organoids: sub-types")


DimPlot(seurat_merged, reduction = "umap", group.by = "Tissue") +
  ggtitle("B816 organoids: sub-types")


tmp<-subset(seurat_merged, subset = Tissue == "Pancreas")
umap_plot(tmp,color=c(pancreas_cols))
generate_feature_plot (tmp,gene="TOP2A")

##################################################
#GENE-list provided by Chris
##################################################


DimPlot(tmp, reduction = "umap", group.by = "type2")+
  ggtitle("B816 organoids: pancreas")+
  scale_color_manual(values=c(pancreas_cols))

FeaturePlot(tmp, features = "TOP2A")


marker_genes_tf<-read.csv("05 CLEAN B816 - umap and Dotplot/scRNAseq genes.csv",as.is=T,header=T)

top_genes<-marker_genes_tf$Gene


dotplot <- DotPlot(seurat_merged, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6,group.by="type2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("DotPlot - All-Cells - Markers.pdf", plot = dotplot)


# (4) Heatmap of the top 10 genes per cell type


seurat_merged <- ScaleData(seurat_merged, features = top_genes)


heatmap <- DoHeatmap(seurat_merged, features = unique(top_genes), size = 3,group.by="type2") +
  ggtitle("Marker Genes per Cell Type")

ggsave("Heatmap - All-Cells - Markers.pdf", plot = heatmap)

####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#INDIVIDUAL SAMPLE SEURAT BOILER PLATE:
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################


######################################
#MAKE INDIVIDUAL FILES
######################################

seurat_merged <- qs::qread("FINAL_B816-organoids.qs")
tmp_bladder <- subset(seurat_merged, subset = Tissue == "Bladder")
tmp_kidney <- subset(seurat_merged, subset = Tissue == "Kidney")
tmp_lung <- subset(seurat_merged, subset = Tissue == "Lung")
tmp_pancreas <- subset(seurat_merged, subset = Tissue == "Pancreas")


# Feature list for gene selection



######################################
#ALL FILES:
######################################


# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_merged, reduction = "umap", group.by = "type2") +
  ggtitle("Cell-Type:  organoid pancreatic cancer")+
  scale_color_manual(values=all_cols)

seurat_obj<-seurat_merged

Idents(seurat_obj) <- seurat_obj@meta.data$type2

seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_obj, reduction = "umap", group.by = "type2") +
  ggtitle("Cell-Type:  organoid pancreatic cancer")+
  scale_color_manual(values=all_cols)


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("ALL_B816_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("All_B816_Heatmap.pdf", plot = heatmap)




######################################
#PANCREAS
######################################


# Visualize the assigned cell types in the UMAP plot
DimPlot(tmp_pancreas, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid pancreatic cancer")+
  scale_color_manual(values=pancreas_cols)

seurat_obj<-tmp_pancreas

Idents(seurat_obj) <- seurat_obj@meta.data$cell_type


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Pancreas_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Pancreas_Heatmap.pdf", plot = heatmap)




######################################
#LUNG
######################################


# Visualize the assigned cell types in the UMAP plot
DimPlot(tmp_lung, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid lung cancer")+coord_cartesian(xlim = c(-12, -4),ylim=c(-1.5,6))+
  scale_color_manual(values=lung_cols)

seurat_obj<-tmp_lung

Idents(seurat_obj) <- seurat_obj@meta.data$cell_type


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Lung_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Lung_Heatmap.pdf", plot = heatmap)










######################################
#BLADDER
######################################



# Visualize the assigned cell types in the UMAP plot
DimPlot(tmp_bladder, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid bladder cancer")+coord_cartesian(xlim = c(1, 8))+
  scale_color_manual(values=bladder_cols)

seurat_obj<-tmp_bladder

Idents(seurat_obj) <- seurat_obj@meta.data$cell_type


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Bladder_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Bladder_Heatmap.pdf", plot = heatmap)




######################################
#KIDNEY
######################################


# Visualize the assigned cell types in the UMAP plot
DimPlot(tmp_kidney, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid kidney cancer")+
  scale_color_manual(values=kidney_cols)

seurat_obj<-tmp_kidney

Idents(seurat_obj) <- seurat_obj@meta.data$cell_type


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Kidney_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Kidney_Heatmap.pdf", plot = heatmap)



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   ##         ########        ##       ###     ##                                  
##     ##          ##         ##             ####      ####    ##                             
##    ##           ##         ##            ##  ##     ## ##   ##                                             
##    ##           ##         #######      ##    ##    ##  ##  ##                                  
##    ##           ##         ##          ##########   ##   ## ##                                   
##     ##          ##         ##          ##      ##   ##    ####                                   
##      ########   ########   ########    ##      ##   ##     ###                                   
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ##########   ##      ##   ##    ##                 
##       ##        ##      ##    ##  ##     
##       ##        ##      ##     ####                       
##       ##        ##      ##      ##          
##       ##        ##      ##     ####           
##       ##         ##    ##     ##  ##          
##       ##          ######     ##    ##           
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#DATE:  2-26-2025
#GOAL:  clean Tux datasets for shiny apps to aid POH collabs


setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 scRNAseq - canine")

library("Seurat")
library("ggplot2")
library(dplyr)
library(qs)

umap_plot <- function(seurat_merged){
  DimPlot(seurat_merged, reduction = "umap", group.by = "type2") +
    ggtitle("B816 organoids: sub-types") + 
    theme(legend.position = "bottom") 
    #+scale_color_manual(values = color)
}

generate_feature_plot <- function(seurat_merged,gene="TOP2A") {
  req(seurat_merged, gene)
  FeaturePlot(seurat_merged, features = gene) +
    scale_color_gradient(low = "grey", high = "red", name = paste0(gene, " (logCounts)")) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 18)
    )
}

load("snRNAseq_Tux_cellTypes_CLEAN.RData")



####################################################################################################################################################################################################################################################################
#0 Clean low read depth cluster:
####################################################################################################################################################################################################################################################################

{
load("snRNAseq_Tux_cellTypes.RData")



colors<-c("blue","skyblue","darkcyan")

DimPlot(snRNAseq_organd_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid bldder cancer")+
  scale_color_manual(values=colors)




colors<-c("darkblue","blue","skyblue","darkcyan","darkred","black","darkgreen","grey40","skyblue")
# Visualize the assigned cell types in the UMAP plot
DimPlot(snRNAseq_invivo_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo bladder cancer")+
  scale_color_manual(values=colors)



################################
#CLEAN In vivo:
################################
DimPlot(snRNAseq_invivo_cellTypes, reduction = "umap", group.by = "cell_type")
FeaturePlot(snRNAseq_invivo_cellTypes, features = "read_depth", cols = c("lightblue", "darkred"))
seurat1<-subset(snRNAseq_invivo_cellTypes, subset = cell_type != "Cancer")


DimPlot(seurat1, reduction = "umap", group.by = "cell_type")
FeaturePlot(seurat1, features = "read_depth", cols = c("lightblue", "darkred"))

colors<-c("darkblue","skyblue","darkcyan","darkred","black","darkgreen","grey40","skyblue")
# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat1, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo bladder cancer")+
  scale_color_manual(values=colors)

snRNAseq_invivo_cellTypes<-seurat1

save(list=c("snRNAseq_organd_cellTypes","snRNAseq_invivo_cellTypes"),file="snRNAseq_Tux_cellTypes_CLEAN.RData")

shiny_invivo<-snRNAseq_invivo_cellTypes[!grepl("ENSCAFG", rownames(snRNAseq_invivo_cellTypes)), ]
shiny_organd<-snRNAseq_organd_cellTypes[!grepl("ENSCAFG", rownames(snRNAseq_organd_cellTypes)), ]


qsave(shiny_invivo, "Tux_blca_invivo.qs")
qsave(shiny_organd, "Tux_blca_organd.qs")

}

####################################################################################################################################################################################################################################################################
#Tux Bladder App:  In vitro vs in vivo
####################################################################################################################################################################################################################################################################

colors<-c("blue","skyblue","darkcyan")

DimPlot(snRNAseq_organd_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid bldder cancer")+
  scale_color_manual(values=colors)




colors<-c("blue","skyblue","darkcyan","darkred","black","darkgreen","grey40","skyblue")
# Visualize the assigned cell types in the UMAP plot
DimPlot(snRNAseq_invivo_cellTypes, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  in vivo bladder cancer")+
  scale_color_manual(values=colors)





