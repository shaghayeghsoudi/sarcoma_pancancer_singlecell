#### Rscript to run infrCNV #####

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


merged_obj <- qread("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/singlets_seurat_objects_only_ecotyper_merged_harmony_SampleID.qs")

Idents(merged_obj) <- "RNA_snn_res.0.5" ## choose which you want 
merged_obj@meta.data$seurat_clusters <- merged_obj@meta.data$RNA_snn_res.0.5 # Assign resolution 0.5 clusters to seurat_clusters

# Verify consistency
identical(
  as.character(Idents(merged_obj)), 
  as.character(merged_obj@meta.data$seurat_clusters)
)  # Should return TRUE


#immune_cells <- subset(merged_obj, idents = c(10,7,3,14,4))
immune_cells <- subset(merged_obj, idents = c(1,3,5,11))  # only ecotype for testing Clusters expressing CD3D/CD79A/CD68 
DefaultAssay(merged_obj) <- "RNA"   ## amek sure it is RNA

out_dir<-"/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/inferCNV_outputs/"
### prepare to run inferCNV ####
##Sample-Specific Analysis (Gold Standard)

# Pseudocode for per-sample workflow
for (sample in unique(merged_obj$Sample_ID)) {
  # Subset cells for this sample
  sample_obj <- subset(merged_obj, subset = Sample_ID == sample)
  
  ## Identify sample-specific normal cells (immune/stromal)
  ## Build annotation table
  cell_annotations <- data.frame(
    cell_id = colnames(sample_obj),
    group = ifelse(
      colnames(sample_obj) %in% colnames(immune_cells),
      "normal",
      paste0("tumor_", sample)
    )
  )
    
  # make and save annotation file
  sample<-gsub(" ", "_", sample)
  annot_file <- paste0(out_dir,"infercnv_", sample, "_annotations.txt", sep = "")

  write.table(cell_annotations, annot_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Save count matrix
  counts_matrix <- LayerData(sample_obj, assay = "RNA", layer = "counts")
  #counts_file <- paste0(out_dir,"infercnv_", sample, "_counts.txt", sep = "")
  #write.table(as.matrix(counts_matrix), counts_file, sep = "\t", quote = FALSE)

  # 2. PROPER MATRIX CONVERSION
  if (inherits(counts_matrix, "dgCMatrix")) {
  # For large datasets: Efficient sparse-to-dense conversion
  counts_matrix <- as.matrix(counts_matrix)
  } else if (is(counts_matrix, "DelayedArray")) {
  # For HDF5-backed data
  counts_matrix <- as.matrix(counts_matrix)
  }
   
# 3. VALIDATE MATRIX
stopifnot(
  is.matrix(counts_matrix),          # Must be TRUE now
  !inherits(counts_matrix, "dgCMatrix"), 
  all(is.finite(counts_matrix))     # No NA/NaN/Inf
)



  
  # Run inferCNV per sample
  # Create inferCNV object
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = annot_file,
    delim = "\t",
    gene_order_file = "/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/gencode.v38.gene_positions.txt",
    ref_group_names = c("normal")
  )

  
  # Run inferCNV
  infercnv::run(
    infercnv_obj = infercnv_obj,
    cutoff = 0.1,
    cluster_by_groups = FALSE,
    k_obs_groups = 3,
    denoise = TRUE,
    HMM = TRUE,
    BayesMaxPNormal = 0.1,
    analysis_mode = "subclusters",
    tumor_subcluster_partition_method = "leiden",
    output_format = "hdf5",
    num_threads = 8,
    out_dir = paste0(out_dir,"infercnv_output_", sample,sep = ""),
    
    # ---- Phylogeny-specific additions ----
    plot_steps = TRUE,
    resume_mode = FALSE,
    no_plot = FALSE,
    png_res = 300,
    plot_tree = TRUE,                   # Generates Newick-format tree
    plot_tree_scale = 0.5,              # Adjust branch length scaling
    tree_method = "NJ",                 # Neighbor-joining (default)
    dynamic_resize = 0.8, 

     # ---- Subclone refinement ----
    #min_cells_per_subcluster  = 10,          # Avoid tiny subclones
    max_centered_threshold = 3,        
    reassignCNVs = TRUE,
    
    # ---- HMM tuning ----
    HMM_type = "i6",                        # More CNV states (default: i3)
    HMM_report_by = "subcluster"      
  )
}




