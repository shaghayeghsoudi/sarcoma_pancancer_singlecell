#!/usr/bin/env Rscript

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

#merged_obj <- qread("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/test_scripts/singlets_seurat_objects_only_ecotyper_merged_harmony_SampleID.qs")
merged_obj <- qread("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/doublet_removed_merged_objects/singlets_seurat_objects_all_three_merged_harmony_SampleID.qs")


Idents(merged_obj) <- "RNA_snn_res.0.5" ## choose which you want 
merged_obj@meta.data$seurat_clusters <- merged_obj@meta.data$RNA_snn_res.0.5 # Assign resolution 0.5 clusters to seurat_clusters

# Verify consistency
identical(
  as.character(Idents(merged_obj)), 
  as.character(merged_obj@meta.data$seurat_clusters)
)  # Should return TRUE


immune_cells <- subset(merged_obj, idents = c(3,5,7,10))
#immune_cells <- subset(merged_obj, idents = c(1,3,5,11))  # only ecotype for testing Clusters expressing CD3D/CD79A/CD68 
DefaultAssay(merged_obj) <- "RNA"   ## amek sure it is RNA

out_dir<-"/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/inferCNV/"

### prepare to run inferCNV ####
## making small subsets
#samples<-c("UPS2", "LMS2","UPS1 pt" , "UPS3" , "LMS1 Met pt" ,   "LMS3_T1 Met pt", "LMS3_T2 Met pt" )
samples<-c("MFS3","WDLPS1", "DDLPS1","UPS4", "UPS5","UPS6", "UPS7", "UPS8" , "DDLPS2", "UPS9")
#samples<-c("MFS3 pt" ,"DDLPS1pt","UPS4pt", "UPS5pt" ,"UPS6pt","UPS7pt" , "UDSCS1pt" "DDLPS2pt")


# Pseudocode for per-sample workflow
#for (sample in unique(merged_obj$Sample_ID)) {

for (sample in unique(samples)) {
  
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

  ## PROPER MATRIX CONVERSION
  if (inherits(counts_matrix, "dgCMatrix")) {
  # For large datasets: Efficient sparse-to-dense conversion
  counts_matrix <- as.matrix(counts_matrix)
  } else if (is(counts_matrix, "DelayedArray")) {
  # For HDF5-backed data
  counts_matrix <- as.matrix(counts_matrix)
  }
   
  ## VALIDATE MATRIX
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
  cutoff = 0.1,                              # Keep this for sensitivity 
  cluster_by_groups = FALSE,               # Critical for subclustering (optimal:0.1)
  k_obs_groups = 3,                        # *** Adjust based on expected subclones,If k_obs_groups = 3 gives too coarse of a split, try k_obs_groups = 4 or 5
  denoise = TRUE,                          # Reduces noise
  HMM = TRUE,                              # Essential for CNV calls
  BayesMaxPNormal = 0.1,                   # You get too many CNV calls (try 0.2 or 0.3) or You’re missing known CNV regions (try 0.05 or 0.01) (optimal:0.1)
  analysis_mode = "subclusters",           # Enables subclone detection
  tumor_subcluster_partition_method = "leiden",
  ## tumor_subcluster_partition_method = "random_trees",
  ## tumor_subcluster_partition_method = "kmeans"   # Then make sure this line stays:k_obs_groups = 3  # Try 3–6 depending on desired resolution
  ### output_format = "hdf5", does not support
  num_threads = 30,
  out_dir = paste0(out_dir, "infercnv_output_", sample),

  # Plotting and output
  ## write_phylo=T
  plot_steps = TRUE,
  resume_mode = FALSE,
  no_plot = FALSE,
  png_res = 300,
  #plot_tree = TRUE,           # did not support        # Generates Newick-format tree
  #plot_tree_scale = 0.5,      # did not support       # Adjust branch length scaling
  # tree_method = "NJ",        # did not support        # Neighbor-joining (default)
  #dynamic_resize = 0.8, 

  # Subclone refinement (only supported part)
  max_centered_threshold = 3,
  reassignCNVs = TRUE,

  # HMM tuning
  HMM_type = "i6",
  HMM_report_by = "subcluster"
)

#message("inferCNV done with k=3 Check your output folders!")

}


########################
########## END #########
########################

### suggested by UPylogplot 
#infercnv_obj = infercnv::run(infercnv_obj,
#  cutoff=1,out_dir="output_dir",
#  cluster_by_groups=FALSE,
#  plot_steps=T,
#  scale_data=T,
#  denoise=T,
#  noise_filter=0.12,
#  analysis_mode='subclusters',
#  HMM_type='i6')


##### NEXT STEPS ######
## Link the CNV profile to gene expression in Seurat
#### Extract subcluster labels and map to Seurat
#subclusters <- read.table("infercnv_output/infercnv_subclusters.observation_groupings.txt")
#colnames(subclusters) <- c("cell", "subclone")

### Add to Seurat metadata
#seurat_obj$cnv_subclone <- subclusters$subclone[match(Cells(seurat_obj), subclusters$cell)]

### Now you can visualize subclones in UMAP:
#DimPlot(seurat_obj, group.by = "cnv_subclone", label = TRUE)

## Compute CNV burden per cell
#cnv_matrix <- infercnv_obj@expr.data
#cnv_score <- apply(cnv_matrix, 2, function(x) mean(abs(x - 1)))
#seurat_obj$cnv_score <- cnv_score
#FeaturePlot(seurat_obj, "cnv_score")

## 3. Identify sample- or clone-specific CNV patterns
#Check which clones are enriched in:

#Specific patient samples
#Diagnosis vs relapse
#Spatial zones, etc.

########################################################
############## Creating phylogentics trees #############

# Step 1: Extract CNV profiles per subclone (use final infercnv_obj)

# Load inferCNV final object
#load("infercnv_output/run.final.infercnv_obj")

# Expression matrix of inferred CNVs
#cnv_matrix <- infercnv_obj@expr.data  # genes x cells

# Read subcluster labels
#subclusters <- read.table("infercnv_output/infercnv_subclusters.observation_groupings.txt")
#colnames(subclusters) <- c("cell", "subclone")

# Average CNV profile per subclone
#library(Matrix)
#library(dplyr)

#subclone_profiles <- do.call(cbind, lapply(unique(subclusters$subclone), function(clone) {
#  cells <- subclusters$cell[subclusters$subclone == clone]
#  if (length(cells) > 1) {
#    rowMeans(cnv_matrix[, cells, drop = FALSE])
#  } else {
#    cnv_matrix[, cells]
#  }
#}))
#colnames(subclone_profiles) <- unique(subclusters$subclone)


# Step 2: Compute pairwise distances and build a tree
#library(ape)

# Compute distance matrix between subclones
#dist_mat <- dist(t(subclone_profiles))  # Transpose: subclones as rows

# Build phylogenetic tree (Neighbor-Joining)
#tree <- nj(dist_mat)

# Plot
#plot(tree, main = "Phylogenetic Tree of CNV Subclones")

#Step 3: Annotate the tree (optional but powerful)

#library(ggtree)

# Create a data.frame for annotations
#tree_anno <- data.frame(
#  subclone = tree$tip.label,
#  color = "black"  # or assign based on sample, timepoint, etc.
#)

# Plot with ggtree
#p <- ggtree(tree) %<+% tree_anno + 
#  geom_tiplab(aes(color = color), size = 4) +
#  theme_tree2()
#
#print(p)
##########################
# Run inferCNV
#  infercnv::run(
#    infercnv_obj = infercnv_obj,
#    cutoff = 0.1,
#    cluster_by_groups = FALSE,
#    k_obs_groups = 3,
#    denoise = TRUE,
#    HMM = TRUE,
#    BayesMaxPNormal = 0.1,
#    analysis_mode = "subclusters",
#    tumor_subcluster_partition_method = "leiden",
#    output_format = "hdf5",
#    num_threads = 8,
#    out_dir = paste0(out_dir,"infercnv_output_", sample,sep = ""),
#    
#    # ---- Phylogeny-specific additions ----
#    plot_steps = TRUE,
#    resume_mode = FALSE,
#    no_plot = FALSE,
#    png_res = 300,
#    plot_tree = TRUE,                   # Generates Newick-format tree
#    plot_tree_scale = 0.5,              # Adjust branch length scaling
#    tree_method = "NJ",                 # Neighbor-joining (default)
#    dynamic_resize = 0.8, 
#
#     # ---- Subclone refinement ----
#    #min_cells_per_subcluster  = 10,          # Avoid tiny subclones
#    max_centered_threshold = 3,        
#    reassignCNVs = TRUE,
#    
#    # ---- HMM tuning ----
#    HMM_type = "i6",                        # More CNV states (default: i3)
#    HMM_report_by = "subcluster"      
#  )


