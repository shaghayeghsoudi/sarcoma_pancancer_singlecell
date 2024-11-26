

rm(list = ls())

## Step 1: Load and Create Seurat Objects for Each Dataset
# Begin by loading each dataset into a Seurat object. Even if the datasets are aligned to different reference genomes, this initial step is the same.


library(Seurat)
library(stringr)
library(Matrix)
library(ggplot2)
library(DoubletFinder)
library(glmGamPoi)
#library(DoubletFinder)
### read 10X cell ranger output, create Seurat object, do filtering and save the object in RDS format

full_meta<-read.csv(file="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/metadata/sarcoma_pancancer_metadata_updated.csv", header = TRUE, sep = ",")
root<-("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/")

data_dir <- paste(root,"ecotyper/ecotyper_all/", sep = "")  ### Adjust -> all files should be in the same folder ###
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
    cur_seurat[["percent.mt"]] <- PercentageFeatureSet(cur_seurat, pattern = "^MT-")

    ### add related information to the meta data
    full_meta_foc<-full_meta[full_meta$Patient == file,]
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
folder_path <- paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/",basename(data_dir), sep = "")

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
    group.by="Sample_ID")
print(plot_scatter)    
dev.off()     
#FeatureScatter(seurat_obj_study, feature1 = "nCount_RNA",feature2="percent.mt",group.by='SampleID')


# plot distributions of QC metrics, grouped by SampleID
png(paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/",basename(data_dir),"_basic_qc.png", sep = ""), width=6, height=10, res=200, units='in')
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

png(paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/",basename(data_dir),"basic_cells_per_sample_filtered.png", sep = ""),width=7, height=3, res=200, units='in')
print(p)
dev.off()

saveRDS(seurat_obj_study , file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/",basename(data_dir),"_seurat_object.rds",sep = ""))
# sc_tumor<-readRDS(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/Liu_primary_all_seurat_object.rds")

####################################
#### find and remove doublets ######
####################################

## NOTE: DoubletFinder should not be run on aggregated data and should be run on a per sample basis
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
     normalization.method = "LogNormalize",   ### lognormalize is a simpler method
     scale.factor = 10000
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
     tumor_sample <- RunUMAP(tumor_sample, dims = 1:min.pc)    ### alternative: dims = doubletpc
     tumor_sample <- FindNeighbors(object = tumor_sample, dims = 1:min.pc)     ### alternative: dims = doubletpc     
     tumor_sample <- FindClusters(object = tumor_sample, resolution = 1)     ### alternative: resolution = doubletresolution or 0.1
  
     # pK identification (no ground-truth)
     #sweep.list <- paramSweep(tumor_sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
     sweep.list <- paramSweep(tumor_sample, PCs = 1:min.pc, sct = FALSE)
     sweep.stats <- summarizeSweep(sweep.list)
     bcmvn <- find.pK(sweep.stats)



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
  
  # subset and save
  tumor_singlets <- subset(tumor_sample, doublet_finder == "Singlet")
  sc_tumor_joined_split[[i]] <- tumor_singlets
  remove(tumor_singlets)


} #### end of i loop



# converge mouse.split
seurat_obj_study_singlets <- merge(x=sc_tumor_joined_split[[1]], y=sc_tumor_joined_split[2:length(sc_tumor_joined_split)])

saveRDS(seurat_obj_study_singlets , file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/",basename(data_dir),"singlets_seurat_object.rds",sep = ""))
seurat_obj_study_singlets_joined<-JoinLayers(seurat_obj_study_singlets)   ## join layers



## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#seu_kidney <- doubletFinder(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_kidney <- doubletFinder(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)



# To load the saved Seurat object later, use:
# loaded_seurat <- readRDS("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/ecotyper_primary_seurat_object.rds")


### END ####




########################################################################################
########## Lindy's paper (Plate base Single-cell - processed by the authors) ###########
#########################################################################################

rm(list = ls())

## Step 1: Load and Create Seurat Objects for Each Dataset
# Begin by loading each dataset into a Seurat object. Even if the datasets are aligned to different reference genomes, this initial step is the same.


library(Seurat)
library(stringr)
library(Matrix)
library(ggplot2)
library(DoubletFinder)
library(glmGamPoi)


ew <- readRDS("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/ewing/ewing_data.RDS")
ew_meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/ewing/ewing_meta_forSharing.txt", header = T, sep = " ")
ew_object <- CreateSeuratObject(counts = ew,meta.data = ew_meta)


#sc_tumor_joined_split <- SplitObject(sc_tumor_joined, split.by = "Sample_ID") 

#ew_object@meta.data$orig.ident<-paste(ew_object@meta.data$orig.ident, ew_object@meta.data$sample_id , sep = "_")
table(rowSums(ew_object@assays$RNA) > 0)
ew_object[["percent.mt"]] <- PercentageFeatureSet(ew_object, pattern = "^MT-")

head(ew_object@meta.data)


no_unk <- !grepl("UNK", rownames(ew_object@meta.data))

# Logical condition: rows that have 0 in the second column
# Replace 'Column2Name' with the actual name of the second column
zero_in_ncell <- ew_object@meta.data$n_cell != 0


# Combine both conditions using the AND operator
combined_condition <- no_unk  & zero_in_ncell 

# Subset the Seurat object using the combined condition
subset_ew_object <- subset(ew_object, cells = colnames(ew_object)[combined_condition])
subset_ew_object [["percent.mt"]] <- PercentageFeatureSet(subset_ew_object , pattern = "^MT-")

subset_ew_object_primary<-subset(subset_ew_object, subset = sample_id != "ES-016-meta")
subset_ew_object_primary<-subset(subset_ew_object_primary, subset = sample_id != "ES-036")
subset_ew_object_primary<-subset(subset_ew_object_primary, subset = sample_id != "ES-049")


## This did not work
#samples_to_exclude <- c("ES-016-meta", "ES-036","ES-049")  # Replace with your actual sample names
#seurat_obj_filtered <- subset(subset_ew_object, subset = !(sample_id %in% `samples_to_exclude`))


# Define the sample_IDs you want to exclude
### add required missing informtion in the metadata to make it consistent 
subset_ew_object_primary$Dataset<-"Visser et al-Ewing Sarcoma"
subset_ew_object_primary$Library_type<-"plate-based"
subset_ew_object_primary$Sex<-"Unknown"
subset_ew_object_primary$Age<-"Unknown"
subset_ew_object_primary$Tumour_treatment<-"Unknown"
subset_ew_object_primary$Tumour_grade<-"Unknown"

### select matched columns 
subset_ew_object_primary@meta.data <- subset_ew_object_primary@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","sample_id","patient_id","Dataset","Library_type","Sex","Age","tumor_type","fusion_type","sample_type","sample_body_site","sample_treatment_timepoint" ,"Tumour_treatment", "Tumour_grade")]
colnames(subset_ew_object_primary@meta.data)<-c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Sample_ID","Patient_ID","Dataset","Library_type", "Sex","Age" ,"Tumour_type" ,"genetic","Tumor_site" ,"Tumour_location", "Treatment_timepoint", "Tumour_treatment", "Tumour_grade")


subset_ew_object_primary_filt<- subset(subset_ew_object_primary,
    subset = nFeature_RNA > 200 &
    nFeature_RNA < 7000 & 
    percent.mt < 20)


# Verify the result by checking the metadata



png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Ewing_Visser_scatter_basic_qc.png", width=8, height=10, res=200, units='in')
FeatureScatter(
    subset_ew_object_primary_filt, 
    feature1 = "nCount_RNA",
    feature2="nFeature_RNA",
    group.by="Sample_ID")
dev.off()     
#FeatureScatter(seurat_obj_study, feature1 = "nCount_RNA",feature2="percent.mt",group.by='SampleID')



# plot distributions of QC metrics, grouped by SampleID
png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Ewing_Visser_basic_qc.png", width=6, height=10, res=200, units='in')
VlnPlot(
    subset_ew_object_primary_filt,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by="Sample_ID",
    pt.size=0,
    ncol = 1
    )
dev.off()       


# plot the number of cells in each sample post filtering
df_cell_number <- as.data.frame(rev(table(subset_ew_object_primary_filt$Sample_ID)))
colnames(df_cell_number) <- c('Sample_ID', 'n_cells')
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

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Ewing_Visser_basic_cells_per_sample_filtered.png",width=7, height=3, res=200, units='in')
print(p)
dev.off()

saveRDS(subset_ew_object_primary_filt , file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/Ewing_Visser_seurat_object.rds")
# sc_tumor<

###################################
################# End #############
###################################

### Jerby paper (synovial sarcoma)
rm(list = ls())

## Step 1: Load and Create Seurat Objects for Each Dataset
# Begin by loading each dataset into a Seurat object. Even if the datasets are aligned to different reference genomes, this initial step is the same.
#### Jerby dataset

library(Seurat)
library(stringr)
library(Matrix)
library(ggplot2)
library(DoubletFinder)
library(glmGamPoi)


#load expression matrix (Jerby papper)
tumor_counts <- read.csv(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/Jerby/GSM3770931/GSM3770931_SyS.tumors_counts.csv", row.names = 1, header = T)

#remove  mets
tumor_exp_primary<-tumor_counts[,(!grepl("met",colnames(tumor_counts)) & !grepl("SS10",colnames(tumor_counts)) & !grepl("SS5",colnames(tumor_counts)) & !grepl("SS13",colnames(tumor_counts))),]

#load metadata df 
#tumor_meta <- read.csv(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/Jerby/GSM3770931/GSM3770931_SyS.tumors_tpm.csv", row.names = 1, header = T)
## Create Seurat object
SS_obj_primary <- CreateSeuratObject(counts = as.matrix(tumor_exp_primary))
                                
full_meta<-read.csv(file="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/meta_for_testing.csv", header = TRUE, sep = ",")
full_meta_sys<-full_meta[full_meta$Tumour_Type=="Synovial Sarcoma",]

#extract_first_three_or_met_pt <- function(x) {
#  if (grepl("Met|pt", x, ignore.case = TRUE)) {
#    # If "Met" or "pt" is found, extract the first three characters along with "Met" or "pt"
#    return(sub("^(.{3}.*?(Met|pt).*).*$", "\\1", x, ignore.case = TRUE))
#  } else {
#    # Otherwise, just extract the first three characters
#    return(substr(x, 1, 3))
#  }
#}

#SS_obj@meta.data$new_column <- sapply(SS_obj@meta.data$orig.ident, extract_first_three_or_met_pt)




extract_before_third_alpha <- function(x) {
  str_extract(x, "^[^a-zA-Z]*[a-zA-Z][^a-zA-Z]*[a-zA-Z][^a-zA-Z]*")
}

# Create a new column with the extracted part
SS_obj_primary$Patient_ID <- sapply(SS_obj_primary$orig.ident, extract_before_third_alpha)
Patient_ID<-unique(SS_obj_primary$Patient_ID )

SS_tumor_joined_split <- SplitObject(SS_obj_primary, split.by = "Patient_ID") 


for(i in 1:length(SS_tumor_joined_split)){

   full_meta_focal<-full_meta_sys[full_meta_sys$Patient_ID== unique(SS_tumor_joined_split[[i]]$Patient_ID),]

   SS_tumor_joined_split[[i]][["percent.mt"]] <- PercentageFeatureSet(SS_tumor_joined_split[[i]], pattern = "^MT-")

    SS_tumor_joined_split[[i]]$Sample_ID<-full_meta_focal$Sample_ID
    SS_tumor_joined_split[[i]]$Dataset<-full_meta_focal$Dataset
    SS_tumor_joined_split[[i]]$Library_type<-full_meta_focal$Sequencing_Plaftorm
    SS_tumor_joined_split[[i]]$Sex<-full_meta_focal$Sex
    SS_tumor_joined_split[[i]]$Age<-full_meta_focal$Age
    SS_tumor_joined_split[[i]]$Tumour_type<- full_meta_focal$Tumour_Type
    SS_tumor_joined_split[[i]]$genetic<- full_meta_focal$Genetic
    SS_tumor_joined_split[[i]]$Tumor_site<- full_meta_focal$Tumour_Site
    SS_tumor_joined_split[[i]]$Tumour_location<- full_meta_focal$Tumour_Location
    SS_tumor_joined_split[[i]]$Treatment_timepoint <- full_meta_focal$Sample_Treatment_Timepoint
    SS_tumor_joined_split[[i]]$Tumour_treatment<-full_meta_focal$Tumour_treatment
    SS_tumor_joined_split[[i]]$Tumour_grade<-full_meta_focal$Tumour_Grade


    SS_tumor_joined_split[[i]]<- subset(SS_tumor_joined_split[[i]],
    subset = nFeature_RNA > 200 &
    nFeature_RNA < 7000 & 
    percent.mt < 20)

    metadata <- SS_tumor_joined_split[[i]]@meta.data

    # Switch the two columns by reordering them
    metadata <- metadata[, c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Sample_ID","Patient_ID","Dataset","Library_type","Sex" ,"Age", "Tumour_type", "genetic" ,"Tumor_site", "Tumour_location" ,"Treatment_timepoint" ,"Tumour_treatment" ,"Tumour_grade")]

# Assign the updated metadata back to the Seurat object
     SS_tumor_joined_split[[i]]@meta.data <- metadata



}


seurat_obj_study_SS  <- merge(x= SS_tumor_joined_split[[1]], y= SS_tumor_joined_split[2:length( SS_tumor_joined_split)])
#seurat_obj_study_singlets_joined<-JoinLayers(seurat_obj_study_SS)   ## join layers


png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby31_Synovial_scatter_basic_qc.png", width=8, height=10, res=200, units='in')
FeatureScatter(
    seurat_obj_study_SS, 
    feature1 = "nCount_RNA",
    feature2="nFeature_RNA",
    group.by="Sample_ID")
dev.off()     
#FeatureScatter(seurat_obj_study, feature1 = "nCount_RNA",feature2="percent.mt",group.by='SampleID')



# plot distributions of QC metrics, grouped by SampleID
png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby31_Synovial_basic_qc.png", width=6, height=10, res=200, units='in')
VlnPlot(
    seurat_obj_study_SS,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by="Sample_ID",
    pt.size=0,
    ncol = 1
    )
dev.off()       


# plot the number of cells in each sample post filtering
df_cell_number <- as.data.frame(rev(table(seurat_obj_study_SS$Sample_ID)))
colnames(df_cell_number) <- c('Sample_ID', 'n_cells')
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

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby31_Synovial_basic_cells_per_sample_filtered.png",width=7, height=3, res=200, units='in')
print(p)
dev.off()



saveRDS(seurat_obj_study_SS , file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/Synovial31_Jerby_seurat_object.rds")
SS_singlets_joined<-JoinLayers(seurat_obj_study_SS)   ## join layers



####################################
#### find and remove doublets ######
####################################

## NOTE: DoubletFinder should not be run on aggregated data and should be run on a per sample basis

table(SS_singlets_joined$Sample_ID)

#doubletpercent = 0.075
#doubletpc = seq(30)
#doubletresolution = 1

tumor_joined_split <- SplitObject(SS_singlets_joined, split.by = "Sample_ID") 


# loop through samples to find doublets
for (i in 1:length(tumor_joined_split)) {   ### loop through easch sample in the seurat object
  # print the sample we are on
  #print(paste0("Sample_ID ",i))
  
     # Pre-process seurat object with standard seurat workflow
     tumor_sample <- NormalizeData(tumor_joined_split[[i]],
     normalization.method = "LogNormalize",   ### lognormalize is a simpler method
     scale.factor = 10000
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
     tumor_sample <- RunUMAP(tumor_sample, dims = 1:min.pc)    ### alternative: dims = doubletpc
     tumor_sample <- FindNeighbors(object = tumor_sample, dims = 1:min.pc)     ### alternative: dims = doubletpc     
     tumor_sample <- FindClusters(object = tumor_sample, resolution = 1)     ### alternative: resolution = doubletresolution or 0.1
  
     # pK identification (no ground-truth)
     #sweep.list <- paramSweep(tumor_sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
     sweep.list <- paramSweep(tumor_sample, PCs = 1:min.pc, sct = FALSE)
     sweep.stats <- summarizeSweep(sweep.list)
     bcmvn <- find.pK(sweep.stats)



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
  
  # subset and save
  tumor_singlets <- subset(tumor_sample, doublet_finder == "Singlet")
  tumor_joined_split[[i]] <- tumor_singlets
  remove(tumor_singlets)


} #### end of i loop



# converge mouse.split
seurat_obj_study_singlets <- merge(x=tumor_joined_split[[1]], y=tumor_joined_split[2:length(tumor_joined_split)])

table(seurat_obj_study_singlets$Sample_ID)

#### plot singlets #####

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby31_Synovial_singlets_scatter_basic_qc.png", width=8, height=10, res=200, units='in')
FeatureScatter(
    seurat_obj_study_singlets, 
    feature1 = "nCount_RNA",
    feature2="nFeature_RNA",
    group.by="Sample_ID")
dev.off()     
#FeatureScatter(seurat_obj_study, feature1 = "nCount_RNA",feature2="percent.mt",group.by='SampleID')



# plot distributions of QC metrics, grouped by SampleID
png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby31_Synovial_singlets_basic_qc.png", width=6, height=10, res=200, units='in')
VlnPlot(
    seurat_obj_study_singlets,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by="Sample_ID",
    pt.size=0,
    ncol = 1
    )
dev.off()       


# plot the number of cells in each sample post filtering
df_cell_number <- as.data.frame(rev(table(seurat_obj_study_singlets$Sample_ID)))
colnames(df_cell_number) <- c('Sample_ID', 'n_cells')
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

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby31_Synovial_Singlets_basic_cells_per_sample_filtered.png",width=7, height=3, res=200, units='in')
print(p)
dev.off()



saveRDS(seurat_obj_study_singlets , file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/Synovial31_Jerby_seurat_object_Singlets.rds")
seurat_obj_study_singlets_joined<-JoinLayers(seurat_obj_study_singlets)   ## join layers



## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#seu_kidney <- doubletFinder(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_kidney <- doubletFinder(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)



# To load the saved Seurat object later, use:
# loaded_seurat <- readRDS("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/ecotyper_primary_seurat_object.rds")


#####################################
############# Jerby32-10X ###########


rm(list = ls())

library(Seurat)
library(stringr)
library(Matrix)
library(ggplot2)
library(DoubletFinder)
library(glmGamPoi)

#load expression matrix (Jerby papper)
tumor_counts_10X <- read.csv(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/Jerby/GSM3770932/GSM3770932_SyS.tumors10x_counts.csv.gz", row.names = 1, header = T)

#remove  mets
#tumor_exp_primary_10X<-tumor_counts_10X[,(!grepl("met",colnames(tumor_counts_10X))),]

#load metadata df 
#tumor_meta <- read.csv(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/Jerby/GSM3770931/GSM3770931_SyS.tumors_tpm.csv", row.names = 1, header = T)
## Create Seurat object
SS_obj_primary_10X <- CreateSeuratObject(counts = as.matrix(tumor_counts_10X))


SS_obj_primary_10X$Patient_ID<-sub("^([^\\.]*\\.[^\\.]*)\\..*", "\\1",SS_obj_primary_10X@meta.data$orig.ident)
                                
full_meta<-read.csv(file="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/meta_for_testing.csv", header = TRUE, sep = ",")
full_meta_sys<-full_meta[full_meta$Tumour_Type=="Synovial Sarcoma",]



Patient_ID<-unique(SS_obj_primary_10X$Patient_ID )

SS_tumor_split_10X <- SplitObject(SS_obj_primary_10X, split.by = "Patient_ID") 


for(i in 1:length(SS_tumor_split_10X)){

   full_meta_focal<-full_meta_sys[full_meta_sys$Patient_ID== unique(SS_tumor_split_10X[[i]]$Patient_ID),]

   SS_tumor_split_10X[[i]][["percent.mt"]] <- PercentageFeatureSet(SS_tumor_split_10X[[i]], pattern = "^MT-")

    SS_tumor_split_10X[[i]]$Sample_ID<-full_meta_focal$Sample_ID
    SS_tumor_split_10X[[i]]$Dataset<-full_meta_focal$Dataset
    SS_tumor_split_10X[[i]]$Library_type<-full_meta_focal$Sequencing_Plaftorm
    SS_tumor_split_10X[[i]]$Sex<-full_meta_focal$Sex
    SS_tumor_split_10X[[i]]$Age<-full_meta_focal$Age
    SS_tumor_split_10X[[i]]$Tumour_type<- full_meta_focal$Tumour_Type
    SS_tumor_split_10X[[i]]$genetic<- full_meta_focal$Genetic
    SS_tumor_split_10X[[i]]$Tumor_site<- full_meta_focal$Tumour_Site
    SS_tumor_split_10X[[i]]$Tumour_location<- full_meta_focal$Tumour_Location
    SS_tumor_split_10X[[i]]$Treatment_timepoint <- full_meta_focal$Sample_Treatment_Timepoint
    SS_tumor_split_10X[[i]]$Tumour_treatment<-full_meta_focal$Tumour_treatment
    SS_tumor_split_10X[[i]]$Tumour_grade<-full_meta_focal$Tumour_Grade


    # Create new cell names by combining metadata with existing cell names
    new_names <- paste(SS_tumor_split_10X[[i]]@meta.data$Sample_ID,colnames(SS_tumor_split_10X[[i]]), sep = "_")

    # Rename cells in the Seurat object
    SS_tumor_split_10X[[i]] <- RenameCells(object = SS_tumor_split_10X[[i]], new.names = new_names)


    SS_tumor_split_10X[[i]]<- subset(SS_tumor_split_10X[[i]],
    subset = nFeature_RNA > 200 &
    nFeature_RNA < 7000 & 
    percent.mt < 20)

     meta10X<- SS_tumor_split_10X[[i]]@meta.data

    # Switch the two columns by reordering them
    meta10X <- meta10X[, c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Sample_ID","Patient_ID","Dataset","Library_type","Sex" ,"Age", "Tumour_type", "genetic" ,"Tumor_site", "Tumour_location" ,"Treatment_timepoint" ,"Tumour_treatment" ,"Tumour_grade")]

# Assign the updated metadata back to the Seurat object
     SS_tumor_split_10X[[i]]@meta.data <- meta10X



}


seurat_obj_study_SS10X  <- merge(x= SS_tumor_split_10X[[1]], y= SS_tumor_split_10X[2:length( SS_tumor_split_10X)])
#seurat_obj_study_singlets_joined<-JoinLayers(seurat_obj_study_SS10X )   ## join layers

### exclude met sample 

seurat_obj_study_SS10X_primary<-subset(seurat_obj_study_SS10X, subset = Sample_ID != "SyS13_10x Met pt")




png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby32_10X_Synovial_primary_scatter_basic_qc.png", width=8, height=10, res=200, units='in')
FeatureScatter(
    seurat_obj_study_SS10X_primary , 
    feature1 = "nCount_RNA",
    feature2="nFeature_RNA",
    group.by="Sample_ID")
dev.off()     
#FeatureScatter(seurat_obj_study, feature1 = "nCount_RNA",feature2="percent.mt",group.by='SampleID')



# plot distributions of QC metrics, grouped by SampleID
png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby32_10X_Synovial_primary_Synovial_basic_qc.png", width=6, height=10, res=200, units='in')
VlnPlot(
    seurat_obj_study_SS10X_primary,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by="Sample_ID",
    pt.size=0,
    ncol = 1
    )
dev.off()       


# plot the number of cells in each sample post filtering
df_cell_number <- as.data.frame(rev(table(seurat_obj_study_SS10X_primary$Sample_ID)))
colnames(df_cell_number) <- c('Sample_ID', 'n_cells')
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

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby32_10X_Synovial_primary_basic_cells_per_sample_filtered.png",width=7, height=3, res=200, units='in')
print(p)
dev.off()



saveRDS(seurat_obj_study_SS10X_primary , file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/Synovial32_10X_Jerby_seurat_object.rds")
SS_singlets_joined<-JoinLayers(seurat_obj_study_SS)   ## join layers



###############################################
#### find and remove doublets (Jerby 10X)######
###############################################

## NOTE: DoubletFinder should not be run on aggregated data and should be run on a per sample basis

table(seurat_obj_study_SS10X_primary$Sample_ID)

#doubletpercent = 0.075
#doubletpc = seq(30)
#doubletresolution = 1

tumor_joined_split <- SplitObject(seurat_obj_study_SS10X_primary, split.by = "Sample_ID") 


# loop through samples to find doublets
for (i in 1:length(tumor_joined_split)) {   ### loop through easch sample in the seurat object
  # print the sample we are on
  #print(paste0("Sample_ID ",i))
  
     # Pre-process seurat object with standard seurat workflow
     tumor_sample <- NormalizeData(tumor_joined_split[[i]],
     normalization.method = "LogNormalize",   ### lognormalize is a simpler method
     scale.factor = 10000
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
     tumor_sample <- RunUMAP(tumor_sample, dims = 1:min.pc)    ### alternative: dims = doubletpc
     tumor_sample <- FindNeighbors(object = tumor_sample, dims = 1:min.pc)     ### alternative: dims = doubletpc     
     tumor_sample <- FindClusters(object = tumor_sample, resolution = 1)     ### alternative: resolution = doubletresolution or 0.1
  
     # pK identification (no ground-truth)
     #sweep.list <- paramSweep(tumor_sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
     sweep.list <- paramSweep(tumor_sample, PCs = 1:min.pc, sct = FALSE)
     sweep.stats <- summarizeSweep(sweep.list)
     bcmvn <- find.pK(sweep.stats)



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
  
  # subset and save
  tumor_singlets <- subset(tumor_sample, doublet_finder == "Singlet")
  tumor_joined_split[[i]] <- tumor_singlets
  remove(tumor_singlets)


} #### end of i loop



# converge mouse.split
seurat_obj_study_singlets <- merge(x=tumor_joined_split[[1]], y=tumor_joined_split[2:length(tumor_joined_split)])

table(seurat_obj_study_singlets$Sample_ID)

#### plot singlets #####

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby32_10X_Synovial_singlets_scatter_basic_qc.png", width=8, height=10, res=200, units='in')
FeatureScatter(
    seurat_obj_study_singlets, 
    feature1 = "nCount_RNA",
    feature2="nFeature_RNA",
    group.by="Sample_ID")
dev.off()     
#FeatureScatter(seurat_obj_study, feature1 = "nCount_RNA",feature2="percent.mt",group.by='SampleID')



# plot distributions of QC metrics, grouped by SampleID
png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby32_10X_Synovial_singlets_basic_qc.png", width=6, height=10, res=200, units='in')
VlnPlot(
    seurat_obj_study_singlets,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by="Sample_ID",
    pt.size=0,
    ncol = 1
    )
dev.off()       


# plot the number of cells in each sample post filtering
df_cell_number <- as.data.frame(rev(table(seurat_obj_study_singlets$Sample_ID)))
colnames(df_cell_number) <- c('Sample_ID', 'n_cells')
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

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/plots/Jerby32_10X_Synovial_Singlets_basic_cells_per_sample_filtered.png",width=7, height=3, res=200, units='in')
print(p)
dev.off()



saveRDS(seurat_obj_study_singlets , file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/Synovial32_10X_Jerby_seurat_object_Singlets.rds")
seurat_obj_study_singlets_joined<-JoinLayers(seurat_obj_study_singlets)   ## join layers


#####################################################
#####################################################

### merge various objects

ecotyper_singlet<-readRDS(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/ecotyper_primary_allsinglets_seurat_object.rds")
ecotyper_singlet_joined<-JoinLayers(ecotyper_singlet)   ## join layers

Jerby10X_singlet<-readRDS(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_pancancer/objects/Synovial32_10X_Jerby_seurat_object_Singlets.rds")
Jerby10X_singlet_joined<-JoinLayers(Jerby10X_singlet)


aa<-merge(ecotyper_singlet_joined,Jerby10X_singlet_joined) %>% 
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 50) %>%
    RunUMAP(dims = 1:20)


plot1 <- DimPlot(aa, group.by="Sample_ID")
plot1 <- DimPlot(aa, group.by="Dataset")


plot2 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size =
   0.1)
plot1 + plot2 + plot_layout(widths = c(1.5, 2))


#######################################
#######################################

### make consitent hg19 and hg38 genome assembly ###


# Load hg19 and hg38 datasets
hg19_data <- Read10X(data.dir = "path_to_hg19_data/")
hg38_data <- Read10X(data.dir = "path_to_hg38_data/")

# Create Seurat objects
hg19 <- CreateSeuratObject(counts = hg19_data, project = "hg19")
hg38 <- CreateSeuratObject(counts = hg38_data, project = "hg38")

# Normalize each dataset
hg19 <- NormalizeData(hg19) %>% FindVariableFeatures()
hg38 <- NormalizeData(hg38) %>% FindVariableFeatures()


# Harmonize gene names (ensure they match between datasets)
# You can use biomaRt or custom mappings to convert gene IDs
# Ensure that both datasets have the same set of genes after this step
common_genes <- intersect(rownames(hg19), rownames(hg38))
hg19 <- hg19[common_genes, ]
hg38 <- hg38[common_genes, ]

# Integration using Seurat's integration method
anchors <- FindIntegrationAnchors(object.list = list(hg19, hg38))
combined <- IntegrateData(anchorset = anchors)

# Scale data and run PCA
combined <- ScaleData(combined) %>% RunPCA()

# Run clustering and UMAP for visualization
combined <- FindNeighbors(combined, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)


# Visualize the integrated data
DimPlot(combined, reduction = "umap", group.by = "orig.ident")



##### other way ####
library("dbplyr")
library(biomaRt)

# Select mart for hg19 and hg38
hg19_mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
hg38_mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")




# Get gene symbols or IDs from hg19 and map to hg38
hg19_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = hg19_mart)
hg38_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = hg38_mart)



# Find common genes based on gene symbols or IDs
common_genes <- intersect(hg19_genes$external_gene_name, hg38_genes$external_gene_name)

library(Seurat)

# Create Seurat objects for hg19 and hg38 datasets
hg19 <- CreateSeuratObject(counts = hg19_data, project = "hg19")
hg38 <- CreateSeuratObject(counts = hg38_data, project = "hg38")


# Normalize and identify variable features for each dataset
hg19 <- NormalizeData(hg19) %>% FindVariableFeatures()
hg38 <- NormalizeData(hg38) %>% FindVariableFeatures()

# Find integration anchors between hg19 and hg38 datasets
anchors <- FindIntegrationAnchors(object.list = list(hg19, hg38), dims = 1:30)

# Integrate the datasets
combined <- IntegrateData(anchorset = anchors, dims = 1:30)



# Scale and perform PCA
combined <- ScaleData(combined) %>% RunPCA()

# Perform UMAP for visualization
combined <- RunUMAP(combined, dims = 1:30)

# Plot UMAP colored by dataset origin
DimPlot(combined, reduction = "umap", group.by = "orig.ident")


library(harmony)
combined <- RunHarmony(combined, group.by.vars = "orig.ident")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)