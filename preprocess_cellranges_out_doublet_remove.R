rm(list = ls())

## Step 1: Load and Create Seurat Objects for Each Dataset
# Begin by loading each dataset into a Seurat object. Even if the datasets are aligned to different reference genomes, this initial step is the same.


library(Seurat)
library(stringr)
library(Matrix)
library(ggplot2)
library(DoubletFinder)
library(glmGamPoi)
library(patchwork)

#library(DoubletFinder)
### read 10X cell ranger output, create Seurat object, do filtering and save the object in RDS format

full_meta<-read.csv(file="/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/metadata/sarcoma_paired_before_after_radiation_metadata.csv", header = TRUE, sep = ",")
root<-("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/")

data_dir <- paste(root,"kalbasi_surgery/", sep = "")  ### Adjust -> all files should be in the same folder ###
files<-dir(data_dir)  ### all files are in the same folder (not seperated by patient)
# Print folder names for debugging
print(files)

# Filter files to include only those with '_matrix.mtx.gz'
files <- files[grep('matrix.mtx.gz', files)]

# Extract file stems (sample identifiers) from the file names
file_stems <- str_replace_all(basename(files), '_matrix.mtx.gz', '')



# construct a list of seurat objects for each sample by iteratively loading each file
seurat_list <- lapply(file_stems, function(file){
     print(file)
     mat <- readMM(paste0(data_dir,file,'_matrix.mtx.gz'))
     genes <- read.csv(file=paste0(data_dir,file,'_features.tsv.gz'), sep='\t', header=FALSE)
     #genes<-genes[!duplicated(genes$V2),]
     barcodes <- read.csv(file=paste0(data_dir,file,'_barcodes.tsv.gz'), sep='\t', header=FALSE)
     colnames(mat) <- barcodes$V1
     rownames(mat) <- genes$V2


    ### removed dulicated genes
    if(sum(duplicated(rownames(mat))) >0){   

         rownames(mat) <- make.unique(rownames(mat))
         mat<- mat[!duplicated(rownames(mat)), ]
         rownames(mat) <- make.names(rownames(mat), unique = TRUE)
    }
      
    cur_seurat <- CreateSeuratObject(
    counts = mat,
    project = basename(data_dir)
    )

    ### add mt percentage 
    cur_seurat[["percent.mt"]] <- PercentageFeatureSet(cur_seurat, pattern = "^MT[-\\.]")

    ### add related information to the meta data
    full_meta_foc<-full_meta[full_meta$Patient_ID == file,]
    cur_seurat$Sample_ID <- full_meta_foc$Sample_ID 
    cur_seurat$Patient_ID<- file
    cur_seurat$Dataset<- full_meta_foc$Dataset
    cur_seurat$Library_type<- full_meta_foc$Sequencing_Plaftorm
    cur_seurat$Sex<- full_meta_foc$Sex
    cur_seurat$Age<- full_meta_foc$Age
    cur_seurat$Tumour_type<- full_meta_foc$Tumour_Type
    cur_seurat$genetic<- full_meta_foc$Genetic
    cur_seurat$Tumor_site<- full_meta_foc$Tumour_Site
    cur_seurat$Tumour_location<- full_meta_foc$Tumour_Location
    cur_seurat$Treatment_timepoint <- full_meta_foc$Sample_Treatment_Timepoint
    cur_seurat$Tumour_treatment<-full_meta_foc$Tumour_treatment
    cur_seurat$Tumour_grade<-full_meta_foc$Tumour_Grade

    # Create new cell names by combining metadata with existing cell names
    new_names <- paste(cur_seurat@meta.data$Sample_ID,colnames(cur_seurat), sep = "_")

    # Rename cells in the Seurat object
    cur_seurat <- RenameCells(object = cur_seurat, new.names = new_names)
    
    cur_seurat_filt<- subset(cur_seurat,
    subset = nFeature_RNA > 200 &
    nFeature_RNA < 7000 & 
    percent.mt < 20)

    # plot distributions of QC metrics, grouped by SampleID
    #png('figures/basic_qc.png', width=10, height=10, res=200, units='in')
    #VlnPlot(
    #    cur_seurat_filt, 
    #    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
    #    ncol = 3,
    #    pt.size=0.1, 
    #    group.by = "SampleID",
    #    cols = c("blue", "red","green")   )        


    
}) ### function


### merge all objects into one and do some plotting
seurat_obj_study <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])


# Define the folder path to save plots
folder_path <- paste("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/plots/",basename(data_dir), sep = "")

# Create the folder
if (!dir.exists(folder_path)) {  # Check if the folder already exists
  dir.create(folder_path)
  cat("Folder created at:", folder_path, "\n")
} else {
  cat("Folder already exists at:", folder_path, "\n")
}


png(paste(folder_path,"scatter_basic_qc.png", sep = "/"), width=8, height=10, res=200, units='in')
plot_scatter<-FeatureScatter(
    seurat_obj_study, 
    feature1 = "nCount_RNA",
    feature2="nFeature_RNA",
    group.by="Sample_ID",
    pt.size = 1,
    raster=FALSE)
print(plot_scatter)    
dev.off()     
#FeatureScatter(seurat_obj_study, feature1 = "nCount_RNA",feature2="percent.mt",group.by='SampleID')


# plot distributions of QC metrics, grouped by SampleID
png(paste(folder_path,"VlnPlot_basic_qc.png", sep = "/"), width=10, height=11, res=200, units='in')
plot_feature<-VlnPlot(
    seurat_obj_study,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by="Sample_ID",
    pt.size=0,
    ncol = 1
    )
print(plot_feature)    
dev.off()       


# plot the number of cells in each sample post filtering
df_cell_number <- as.data.frame(rev(table(seurat_obj_study$Sample_ID)))
colnames(df_cell_number) <- c('Sample_ID', 'n_cells')
png(paste(folder_path,"basic_cells_per_sample_filtered.png", sep = "/"),width=7, height=3, res=200, units='in')
p <- ggplot(df_cell_number, aes(y=n_cells, x=reorder(Sample_ID, -n_cells), fill=Sample_ID)) +
    geom_bar(stat='identity') + 
    scale_y_continuous(expand = c(0,0)) +
    NoLegend() + RotatedAxis() +
    ylab(expression(italic(N)[cells])) + xlab('Sample ID') +
    ggtitle(paste('Total cells post-filtering:', sum(df_cell_number$n_cells))) +
    theme(
    panel.grid.minor=element_blank()
    #panel.grid.major.y=linewidth(colour="lightgray", size=0.5),
)
print(p)
dev.off()


### remove suspicious or samples with low number of cell
## double check is any sampelID is missing
#sum(is.na(seurat_obj_study$Sample_ID))
#table(seurat_obj_study$Sample_ID, useNA = "ifany")
# Show metadata for cells with missing Sample_ID
#seurat_obj_study@meta.data[is.na(seurat_obj_study$Sample_ID), ] |> head()
table(seurat_obj_study@meta.data$Sample_ID)


seurat_obj_study <- subset(seurat_obj_study, subset = Sample_ID != "UPS9pt")
table(seurat_obj_study@meta.data$Sample_ID)

saveRDS(seurat_obj_study , file = paste("/home/shsoudi/sarcoma_pancancer_single-cell/sarcoma_before_after_radiation/objects/",basename(data_dir),"_seurat_object.rds",sep = ""))
# sc_tumor<-readRDS(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/ecotyper_all_seurat_object.rds")

####################################
#### find and remove doublets ######
####################################

## NOTE: DoubletFinder should not be run on aggregated data and should be run on a per sample basis
sc_tumor<- seurat_obj_study
table(sc_tumor$Sample_ID)

### join layers into ONE layer before doublet finder function
sc_tumor_joined<-JoinLayers(sc_tumor)

#doubletpercent = 0.075
#doubletpc = seq(30)
#doubletresolution = 1
sc_tumor_joined_split <- SplitObject(sc_tumor_joined, split.by = "Sample_ID") 


# loop through samples to find doublets
for (i in 1:length(sc_tumor_joined_split)) {   ### loop through easch sample in the seurat object
  # print the sample we are on
  #print(paste0("Sample_ID ",i))
  
     # Pre-process seurat object with standard seurat workflow
     tumor_sample <- NormalizeData(sc_tumor_joined_split[[i]],
     normalization.method = "LogNormalize"  ### lognormalize is a simpler method
     )


     #### Alternative *normalization method* with more power :scTransform (it does normalization, scaling and finding variable feature in one step)
     #sc_tumor <- scTransform(
      #    sc_tumor,
      #    vars.to.regress="percent.mt",
      #    variable.features.n = 3000
      #)    ### The normalized data are under $SCT in the assay slot


     tumor_sample <- FindVariableFeatures(tumor_sample,
     selection.method = "vst",
     nfeatures = 2000) 
     tumor_sample <- ScaleData(tumor_sample)
     tumor_sample <- RunPCA(tumor_sample)
     ElbowPlot(tumor_sample)


     # Find significant PCs
     stdv <- tumor_sample[["pca"]]@stdev
     sum.stdv <- sum(tumor_sample[["pca"]]@stdev)
     percent.stdv <- (stdv / sum.stdv) * 100
     cumulative <- cumsum(percent.stdv)
     co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
     co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
     min.pc <- min(co1, co2)
     min.pc


     # finish pre-processing
     tumor_sample <- FindNeighbors(object = tumor_sample, dims = 1:min.pc)     ### alternative: dims = doubletpc     
     tumor_sample <- FindClusters(object = tumor_sample)     ### alternative: resolution = doubletresolution or 0.1
     tumor_sample <- RunUMAP(tumor_sample, dims = 1:min.pc) ### alternative: dims = doubletpc
     # pK identification (no ground-truth)
     #sweep.list <- paramSweep(tumor_sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
     sweep.list <- paramSweep(tumor_sample, PCs = 1:min.pc, sct = FALSE)
     sweep.stats <- summarizeSweep(sweep.list)
     bcmvn <- find.pK(sweep.stats)

    ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line()

     # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
     bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
     optimal.pk <- bcmvn.max$pK
     optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
    ## Homotypic doublet proportion estimate
    annotations <- tumor_sample@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations) 
    nExp.poi <- round(optimal.pk * nrow(tumor_sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

    # run DoubletFinder
    tumor_sample <- doubletFinder(seu = tumor_sample, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
   metadata_sample <- tumor_sample@meta.data
   colnames(metadata_sample)[21] <- "doublet_finder"
   tumor_sample@meta.data <- metadata_sample 
  
  ### plot doublets
  plot_singlet_doublet<-DimPlot(tumor_sample, reduction = 'umap', group.by = "doublet_finder")
  table(tumor_sample@meta.data$doublet_finder)

  # subset and save
  tumor_singlets <- subset(tumor_sample, doublet_finder == "Singlet")
  sc_tumor_joined_split[[i]] <- tumor_singlets
  plot_doublet<-DimPlot(sc_tumor_joined_split[[i]], reduction = 'umap', group.by = "doublet_finder")
  
  ### plot 
  sample<-unique(sc_tumor_joined_split[[i]]@meta.data$Patient_ID)
  png(paste(folder_path,"/","doublets_",sample,".png", sep = ""),width=9, height=5, res=200, units='in')
  plot_both<-plot_singlet_doublet + plot_doublet
  print(plot_both)
  dev.off()
  

  remove(tumor_singlets)


} #### end of i loop



