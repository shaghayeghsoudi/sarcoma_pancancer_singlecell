library(Seurat)
library(stringr)
library(Matrix)
library(ggplot2)
library(DoubletFinder)
#library(glmGamPoi)
library(patchwork)
library(harmony)
library(dplyr)
library(qs)



# merge samples and save doublet-removed object
object_dir<-"/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed"
object_files <- list.files(object_dir, pattern = "\\.rds$", full.names = TRUE)

seurat_list <- lapply(object_files, readRDS)

# Extract names from filenames (without extension)
#object_names <- tools::file_path_sans_ext(basename(object_files))

# Assign sample_group metadata and names
#for (i in seq_along(seurat_list)) {
#  seurat_list[[i]]$sample_group <- object_names[i]
#}

#names(seurat_list) <- object_names

## merge all objects
merged_obj <- Reduce(function(x, y) merge(x, y), seurat_list)
class(merged_obj)

saveRDS(merged_obj , file = paste("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed_merged_objects/singlets_seurat_objects_all_three_merged.rds",sep = ""))
#merged_obj<-readRDS(file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed_merged_objects/singlets_seurat_objects_all_three_merged.rds")

merged_obj<-JoinLayers(merged_obj)   ## join layers


merged_obj <- NormalizeData(merged_obj) %>%
              FindVariableFeatures(nfeatures = 3000) %>%
              ScaleData() %>%
              RunPCA(npcs = 30)  ## run PCA


merged_obj <- RunHarmony(merged_obj, group.by.vars = "Sample_ID") ### run harmony
head(Embeddings(merged_obj, "harmony"))

merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

### save harmonized data ###
saveRDS(merged_obj,file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed_merged_objects/singlets_seurat_objects_all_three_merged_harmony_SampleID.rds")
qsave(merged_obj,file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed_merged_objects/singlets_seurat_objects_all_three_merged_harmony_SampleID.qs")
#merged_obj <- qread("merged_seurat_object.qs")

### do some plotting
pdf("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/plots/all_datasets_harmonized/dimplot_harmony.pdf", width = 16, height = 10)
dim_harmony_sampleid<-DimPlot(merged_obj, reduction= "umap", group.by = "Sample_ID", label = T,pt.size = 0.8)
dim_harmony_dataset<-DimPlot(merged_obj, reduction= "umap", group.by = "orig.ident", label = T,pt.size = 0.8)
plot_both<-(dim_harmony_sampleid|dim_harmony_dataset)
print(plot_both)
dev.off()



resolutions <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)  
pdf("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/plots/all_datasets_harmonized/clustered_resolurions.pdf", width = 20, height = 20)
umap_plots <- list()

for (res in resolutions) {
  umap_plots[[as.character(res)]] <- DimPlot(merged_obj, group.by = paste0("RNA_snn_res.", res), label = TRUE) + 
    ggtitle(paste("Resolution =", res))
}

combined_plot <- wrap_plots(umap_plots, ncol = 2) 
print(combined_plot)
dev.off()

pdf("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/plots/all_datasets_harmonized/feature_plots_harminized_sampleID.pdf", width = 20, height = 20)
plot_feature<-FeaturePlot(merged_obj, c("CD68","COL1A1","FLT1","NKG7","CCL4","SAT1","FN1","IGFBP3"))  ## TO DO: add some features
print(plot_feature)
dev.off()

#cd4_seurat$seurat_clusters <- cd4_seurat$RNA_snn_res.0.7 ## decided to go with res 0.7
#Idents(cd4_seurat)<-"RNA_snn_res.0.7"

pdf("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/plots/all_datasets_harmonized/clustered_resolurions0.5.pdf", width = 20, height = 20)
plot1_harmony_resolution<-DimPlot(merged_obj, group.by = "RNA_snn_res.0.5", label = TRUE)
print(plot1_harmony_resolution)
dev.off()

Idents(merged_obj) <- "RNA_snn_res.0.5" ## choose which you want 
# Assign resolution 0.5 clusters to seurat_clusters
merged_obj@meta.data$seurat_clusters <- merged_obj@meta.data$RNA_snn_res.0.5

# Verify consistency
identical(
  as.character(Idents(merged_obj)), 
  as.character(merged_obj@meta.data$seurat_clusters)
)  # Should return TRUE


#DimPlot(merged_obj, label = TRUE) 


