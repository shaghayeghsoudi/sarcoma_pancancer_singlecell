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
library(rjags) 
library(infercnv)
library(HDF5Array)
library(SummarizedExperiment)
library(rhdf5)


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

######## save finalized data (harmonized) ############

#saveRDS(merged_obj,file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed_merged_objects/singlets_seurat_objects_all_three_merged_harmony_SampleID.rds")
qsave(merged_obj,file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed_merged_objects/singlets_seurat_objects_all_three_merged_harmony_SampleID.qs")

### qsave(merged_obj,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_before_after_radiation/test_scripts/singlets_seurat_objects_only_ecotyper_merged_harmony_SampleID.qs")

#merged_obj <- qread("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/singlets_seurat_objects_only_ecotyper_merged_harmony_SampleID.qs")

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

pdf("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/dimplot_clustered_resolurions0.5.pdf", width = 20, height = 20)
plot1_harmony_resolution<-DimPlot(merged_obj, 
                                  group.by = "RNA_snn_res.0.5", 
                                  label = TRUE,
                                  pt.size = 0.9)
print(plot1_harmony_resolution)
dev.off()

Idents(merged_obj) <- "RNA_snn_res.0.5" ## choose which you want 
merged_obj@meta.data$seurat_clusters <- merged_obj@meta.data$RNA_snn_res.0.5 # Assign resolution 0.5 clusters to seurat_clusters

# Verify consistency
identical(
  as.character(Idents(merged_obj)), 
  as.character(merged_obj@meta.data$seurat_clusters)
)  # Should return TRUE

#DimPlot(merged_obj, label = TRUE) 
##Annotate Major Lineages (Use markers to label broad cell types (immune, stromal, potential tumor)

###########################################
####### From here work on inferCNV ########
###########################################


#1. Pan-Immune Markers (All Immune Cells)
#General leukocyte marker:
#  PTPRC (CD45)
#T cells:
#  CD3D, CD3E, CD3G, CD8A, CD4, TRAC, TRBC1/2
#B cells:
#  CD19, CD79A, MS4A1 (CD20), IGHG1, IGKC
#NK cells:
#  NKG7, GNLY, FCGR3A (CD16)
#Myeloid cells (Macrophages/Dendritic cells):
#  CD14, CD68, FCGR3A, CST3, LYZ
#Endothelial cells: PECAM1 (CD31), VWF
#Fibroblasts/Stromal cells: COL1A1, DCN, LUM


# Immune cells (diploid reference candidates)
pdf("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/feature_plots_resolution0.5_refence_high_confidence_immunecells.pdf", width = 22, height = 22)
feature_immune<-FeaturePlot(merged_obj, features = c("PTPRC","CD3D", "CD79A", "CD68","NKG7","CD14" ),raster=FALSE)  # T cells, B cells, Macrophages
print(feature_immune)
dev.off()

# Stromal cells (caution: may harbor CNVs)
pdf("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/feature_plots_resolution0.5_stromalcells.pdf", width = 22, height = 22)
feature_stromal<-FeaturePlot(merged_obj, features = c("COL1A1", "PECAM1", "LUM","DCN", "FAP", "VWF", "CLDN5"),raster=FALSE)  # Fibroblasts, Endothelial
print(feature_stromal)
dev.off()

# Liposarcoma cells (likely CNV-high)
FeaturePlot(merged_obj, features = c("MDM2", "CDK4", "HMGA2"))  # Hallmark liposarcoma amplific  


## option A #################################
## Option A: Immune Cells (Best Default) ####
immune_cells <- subset(merged_obj, idents = c(10,7,3,14,4))  # Clusters expressing CD3D/CD79A/CD68 
#immune_cells <- subset(merged_obj, idents = c(1,3,5,11))  # only ecotype for testing Clusters expressing CD3D/CD79A/CD68 


qsave(immune_cells,file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed_merged_objects/immune_cells_all_three_merged_harmony_0.5.qs")
DefaultAssay(merged_obj) <- "RNA". ## amek sure it is RNA

################################
### prepare to run inferCNV ####

# Step1: Validate Reference Cells
## Check CNV Baselines (Quick Sanity Check) ## (temporary inferCNV using candidate references)

# Create annotation file marking immune/stromal as "normal"
cell_annotations <- data.frame(
  cell_id = colnames(merged_obj),
  type = ifelse(
    colnames(merged_obj) %in% colnames(immune_cells),
    "normal",
    "tumor"
  )
)



#head(cell_annotations,20)
#cell_id   type
#1  UPS2_AAACCCAAGACAAGCC-1  tumor
#2  UPS2_AAACCCACAAGGACAC-1  tumor
#3  UPS2_AAACCCACAATAGGGC-1  tumor
write.table(cell_annotations, "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/cell_annotations_ecotyper.txt", sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE)


## get counts matrix
counts_matrix <- LayerData(merged_obj, assay = "RNA", layer = "counts")
colnames(counts_matrix)[1:5]

# 1. Validate inputs
stopifnot(all(cell_annotations$cell_id %in% colnames(counts_matrix)))
stopifnot("normal" %in% cell_annotations$type)

# Run inferCNV (quick test) 
## Key Note: If you passed the data frame itself to annotations_file, that’s incorrect — it needs to be a filename or path to a file.

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  gene_order_file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/gencode.v38.gene_positions.txt",
  annotations_file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/cell_annotations_ecotyper.txt",
  ref_group_names = "normal"
)



#infercnv_obj <- run(
#  infercnv_obj,
#  cutoff = 0.1,              # Adjust based on noise
#  smooth_method = "pyramidinal",  # Smooth focal CNVs
#  cluster_by_groups = FALSE , # Group by individual cells
#  out_dir = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts" 
#)

### for subclonal articheture detection
infercnv::run(
  infercnv_obj = infercnv_obj,
  cutoff = 0.1,
  cluster_by_groups = FALSE,
  k_obs_groups = 3,
  denoise = TRUE,
  HMM = TRUE,
  BayesMaxPNormal = 0.1,   ## You get too many CNV calls (try 0.2 or 0.3) or You’re missing known CNV regions (try 0.05 or 0.01)
  analysis_mode = "subclusters",
  tumor_subcluster_partition_method = "leiden",
  output_format = "hdf5",
  num_threads = 8,
  out_dir = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/inferCNV_outputs" ,
  
  # ---- Phylogeny-specific additions ----
  reassignCNVs = TRUE,
  cluster_cells = TRUE,
  plot_steps = FALSE  # optional, if you want fewer intermediate plots

)


# Check if reference cells show flat CNV profiles
plot_cnv(infercnv_obj, cluster_by_groups = TRUE)

  
################################################## 
##################################################
#### OptionB : very smart strategy using copyKAT
library(copykat)
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
library(rjags) 
library(rjags)
library(infercnv)

seurat_obj <- qread("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/singlets_seurat_objects_only_ecotyper_merged_harmony_SampleID.qs")
table(seurat_obj$Sample_ID) 
counts_matrix <- as.matrix(GetAssayData(seurat_obj, layer = "counts"))


copykat_results <- copykat(
  rawmat = counts_matrix,
  id.type = "S",  # "S" = gene symbols # Human genes ("S" for symbol, "E" for Ensembl)
  ngene.chr = 5,  # Min genes per chr for CNV estimation (default=5)
  win.size = 25,
  KS.cut = 0.1,   # Sensitivity threshold (lower = more sensitive)
  sam.name = "merged_samples",
  distance = "euclidean",  # Distance metric for clustering
  output.seg = TRUE,
  n.cores = 10
)

## Add predictions to Seurat metadata
seurat_obj$copykat_prediction <- copykat_results$prediction$prediction
table(seurat_obj$copykat_prediction)

## Use this to create the inferCNV annotations
cell_annotations <- data.frame(
  cell_id = colnames(seurat_obj),
  cell_type = ifelse(seurat_obj$copykat_prediction == "diploid", "normal", "tumor")
)

write.table(cell_annotations, "cell_annotations_infercnv.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


## Run inferCNV using CopyKAT-based normal cells as reference
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = GetAssayData(seurat_obj, slot = "counts"),
  annotations_file = "cell_annotations_infercnv.txt",
  gene_order_file = "gene_order_cleaned.txt",
  ref_group_names = "normal"
)

infercnv_obj <- run(infercnv_obj,
                    cutoff = 0.1,
                    cluster_by_groups = TRUE,
                    denoise = TRUE,
                    HMM = TRUE,
                    out_dir = "infercnv_from_copykat/"
)


#############################
#############################

### NOTE: It eats a ton of RAM (may crash or freeze your session)

## solution : Run CopyKAT on 1 sample at a time

library(Seurat)
library(copykat)
library(dplyr)

# 🧬 Load merged Seurat object
seurat_obj <- readRDS("merged_25_samples_seurat.rds")

# 🔁 Split by Sample_ID (adjust if your column is named differently)
seurat_list <- SplitObject(seurat_obj, split.by = "Sample_ID")

# 🗃 Initialize result collector
all_copykat_preds <- list()

# 🔁 Loop over each sample
for (sample_name in names(seurat_list)) {
  message("Running CopyKAT on sample: ", sample_name)
  
  obj <- seurat_list[[sample_name]]
  
  # Optional: filter for quality (adjust as needed)
  obj <- subset(obj, subset = nFeature_RNA > 200 & nCount_RNA > 500)
  
  # Convert to dense matrix
  counts <- as.matrix(GetAssayData(obj, layer = "counts"))
  
  # Run CopyKAT
  ck <- copykat(
    rawmat = counts,
    id.type = "S",
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    sam.name = sample_name,
    distance = "euclidean",
    output.seg = FALSE,
    n.cores = 4
  )
  
  # Collect predictions
  preds <- ck$prediction
  preds$cell_id <- rownames(preds)
  preds$sample <- sample_name
  all_copykat_preds[[sample_name]] <- preds
}

### Combine All Results
# 🧩 Merge predictions from all samples
copykat_combined <- bind_rows(all_copykat_preds)

# 🧠 Add to metadata in original object
copykat_df <- copykat_combined %>%
  select(cell_id, prediction) %>%
  column_to_rownames("cell_id")

# Match and assign
seurat_obj$copykat_prediction <- NA
common_cells <- intersect(colnames(seurat_obj), rownames(copykat_df))
seurat_obj$copykat_prediction[common_cells] <- copykat_df[common_cells, "prediction"]

# Check result
table(seurat_obj$copykat_prediction)


### Now Proceed with inferCNV ###
cell_annotations <- data.frame(
  cell_id = colnames(seurat_obj),
  cell_type = ifelse(seurat_obj$copykat_prediction == "diploid", "normal", "tumor")
)

write.table(cell_annotations, "cell_annotations_copykat.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)







  
