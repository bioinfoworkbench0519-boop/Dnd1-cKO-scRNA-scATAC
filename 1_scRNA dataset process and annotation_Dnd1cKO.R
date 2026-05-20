# 1.scRNA Rawdata process

# 1.1 library packages

library(Seurat)
library(dplyr)
set.seed(1234)
library(reticulate)
library(scales)
library(DoubletFinder)
library(patchwork)
library(velocyto.R)
library(SeuratWrappers)
library(stringr)


# 1.2 load 10X genomics datasets

setwd("~/")

process_sample <- function(data.dir, project, 
                           dims_use = 1:30,
                           PCs_df = 1:10,
                           doublet_rate = 0.1,
                           min_features = 1000) {
  
  # Read data
  obj <- Read10X(data.dir = data.dir) %>%
    CreateSeuratObject(project = project, min.cells = 3, min.features = 200)
  
  # QC
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj <- subset(obj, subset = nFeature_RNA > min_features & 
                          nCount_RNA < 50000 & 
                          percent.mt < 25)
  
  # Preprocessing
  obj <- obj %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA()
  
  # Clusterung
  obj <- FindNeighbors(obj, dims = dims_use) %>%
         FindClusters(resolution = 0.5) %>%
         RunUMAP(dims = dims_use)
  
  # DoubletFinder
  sweep.res <- paramSweep_v3(obj, PCs = PCs_df, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  homotypic.prop <- modelHomotypic(obj$seurat_clusters)
  nExp <- round(doublet_rate * ncol(obj))
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  
  obj <- doubletFinder_v3(obj,
                         PCs = PCs_df,
                         pN = 0.25,
                         pK = pK,
                         nExp = nExp.adj,
                         reuse.pANN = FALSE,
                         sct = FALSE)
  

  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  
  obj_singlet <- subset(obj, subset = get(df_col) == "Singlet")
  
  return(obj_singlet)
}

# Set parameters 
params <- list(
  E17.5_cKO_doubletfinder = list(dir="~/E17.5_cKO", project="E17.5_cKO", dims=1:40, doublet=0.1, min_features=1000),
  E16.5_cKO_doubletfinder = list(dir="~/E16.5_cKO", project="E16.5_cKO", dims=1:30, doublet=0.15, min_features=0),
  E15.5_cKO_doubletfinder = list(dir="~/E15.5_cKO", project="E15.5_cKO", dims=1:24, doublet=0.1, min_features=1000),
  E15.5_ctrl_doubletfinder = list(dir="~/E15.5_ctrl", project="E15.5_ctrl", dims=1:28, doublet=0.1, min_features=1000),
  E14.5_cKO_doubletfinder = list(dir="~/E14.5_cKO", project="E14.5_cKO", dims=1:28, doublet=0.12, min_features=1000),
  E14.5_ctrl_doubletfinder = list(dir="~/E14.5_ctrl", project="E14.5_ctrl", dims=1:28, doublet=0.1, min_features=1000)
)

# Analysis loop
results <- list()

for (name in names(params)) {
  p <- params[[name]]
  
  results[[name]] <- process_sample(
    data.dir = p$dir,
    project = p$project,
    dims_use = p$dims,
    doublet_rate = p$doublet,
    min_features = p$min_features
  )
  
  saveRDS(results[[name]], paste0(name, ".rds"))
}

# 1.3 Perform SCTransform on multiple samples

############################################################
# Mouse fetal testis scRNA-seq analysis pipeline
# Seurat integration / annotation / germ cell subclustering
# RNA velocity preprocessing
############################################################


############################################################
# 1. Load Seurat objects
############################################################

sample_paths <- list(
  E14.5_cKO  = "E14.5_cKO_doubletfinder.rds",
  E14.5_ctrl = "E14.5_ctrl_doubletfinder.rds",
  E15.5_cKO  = "E15.5_cKO_doubletfinder.rds",
  E15.5_ctrl = "E15.5_ctrl_doubletfinder.rds",
  E16.5_cKO  = "E16.5_doubletfinder.rds",
  E17.5_cKO  = "E17.5_doubletfinder.rds"
)

seu.list <- lapply(sample_paths, readRDS)

############################################################
# 2. Remove Dnd1-positive cells from cKO samples
############################################################

filter_samples <- c(
  "E14.5_cKO",
  "E15.5_cKO",
  "E16.5_cKO",
  "E17.5_cKO"
)

for(i in filter_samples){
  seu.list[[i]] <- subset(seu.list[[i]], subset = Dnd1 < 0.5)
}

############################################################
# 3. Whole dataset integration
############################################################

merged <- merge(
  seu.list[[1]],
  y = seu.list[-1]
)

obj.list <- SplitObject(merged, split.by = "orig.ident")

obj.list <- lapply(obj.list, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = 2000
  )
})

features <- SelectIntegrationFeatures(object.list = obj.list)

obj.list <- lapply(obj.list, function(x){
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:30
)

combined <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

combined <- RunUMAP(
  combined,
  reduction = "pca",
  dims = 1:30
)

combined <- FindNeighbors(
  combined,
  reduction = "pca",
  dims = 1:30
)

combined <- FindClusters(
  combined,
  resolution = 0.1
)


############################################################
# 4. UMAP visualization
############################################################

png(
  "ALL/cKO_merge_umap.png",
  width = 1000,
  height = 1000
)

DimPlot(
  combined,
  reduction = "umap",
  label = TRUE,
  label.size = 9,
  repel = TRUE
) +
  NoLegend() +
  NoAxes()

dev.off()

############################################################
# 5. FeaturePlot helper
############################################################

DefaultAssay(combined) <- "RNA"

Featureforanno <- function(obj, gene, celltype){

  FeaturePlot(
    obj,
    features = gene,
    cols = c("lightgrey", "red"),
    max.cutoff = 2.5
  ) +
    NoAxes() +
    theme(
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 20,
        face = "italic"
      ),
      legend.position = c(0,1),
      legend.justification = c(0,1),
      legend.text = element_text(size = 15)
    ) +
    labs(
      title = celltype,
      subtitle = gene
    )
}

############################################################
# 6. Marker visualization
############################################################

png(
  "ALL/cKO_merge_featureplots.png",
  width = 1000,
  height = 1400
)

p1  <- Featureforanno(combined, "Pou5f1",  "Germ Cells")
p2  <- Featureforanno(combined, "Dazl",    "Germ Cells")

p3  <- Featureforanno(combined, "Amh",     "Sertoli Cells")
p4  <- Featureforanno(combined, "Sox9",    "Sertoli Cells")

p5  <- Featureforanno(combined, "Cyp11a1", "Leydig Cells")
p6  <- Featureforanno(combined, "Hsd3b1",  "Leydig Cells")

p7  <- Featureforanno(combined, "Igf1",    "Interstitial/Stromal")
p8  <- Featureforanno(combined, "Acta2",   "Interstitial/Stromal")

p9  <- Featureforanno(combined, "Pecam1",  "Endothelial Cells")
p10 <- Featureforanno(combined, "Esam",    "Endothelial Cells")

p11 <- Featureforanno(combined, "Cd68",    "Macrophage")
p12 <- Featureforanno(combined, "Lyz2",    "Macrophage")

(p1|p3|p5) /
(p2|p4|p6) /
(p7|p9|p11) /
(p8|p10|p12)

dev.off()

############################################################
# 7. Cell type annotation
############################################################

new.cluster.ids <- c(
  "Interstitial/Stromal Cells",
  "Sertoli Cells",
  "PGCs/ECCs",
  "Interstitial/Stromal Cells",
  "Sertoli Cells",
  "Leydig Cells",
  "Interstitial/Stromal Cells",
  "Interstitial/Stromal Cells",
  "PGCs/ECCs",
  "Macrophages",
  "Endothelial Cells"
)

names(new.cluster.ids) <- levels(combined)

combined <- RenameIdents(
  combined,
  new.cluster.ids
)

combined$orig.ident.v2 <- combined$orig.ident

current.labels <- c(
  "E14.5_cKO",
  "E14.5_ctrl",
  "E15.5_cKO",
  "E15.5_ctrl",
  "E16.5_cKO",
  "E17.5_cKO"
)

new.labels <- c(
  "E14.5_Dnd1-cKO",
  "E14.5_Dnd1-control",
  "E15.5_Dnd1-cKO",
  "E15.5_Dnd1-control",
  "E16.5_Dnd1-cKO",
  "E17.5_Dnd1-cKO"
)

combined@meta.data$orig.ident.v2 <- plyr::mapvalues(
  x    = combined@meta.data$orig.ident,
  from = current.labels,
  to   = new.labels
)

saveRDS(
  combined,
  "all_merge_cellanno.rds"
)

############################################################
# 8. Germ cell subclustering
############################################################

germ <- subset(
  combined,
  idents = "PGCs/ECCs",
  subset =
    orig.ident %in% c(
      "E14.5_cKO",
      "E15.5_cKO",
      "E16.5_cKO",
      "E17.5_cKO"
    )
)

DefaultAssay(germ) <- "RNA"

germ.list <- SplitObject(
  germ,
  split.by = "orig.ident"
)

germ.list <- lapply(germ.list, function(x){

  x <- NormalizeData(x)

  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = 2000
  )

})

germ.features <- SelectIntegrationFeatures(
  object.list = germ.list
)

germ.list <- lapply(germ.list, function(x){

  x <- ScaleData(
    x,
    features = germ.features,
    verbose = FALSE
  )

  x <- RunPCA(
    x,
    features = germ.features,
    verbose = FALSE
  )

})

germ.anchors <- FindIntegrationAnchors(
  object.list = germ.list,
  anchor.features = germ.features,
  reduction = "rpca",
  dims = 1:25
)

germ.combined <- IntegrateData(
  anchorset = germ.anchors
)

DefaultAssay(germ.combined) <- "integrated"

germ.combined <- ScaleData(
  germ.combined,
  verbose = FALSE
)

germ.combined <- RunPCA(
  germ.combined,
  npcs = 30,
  verbose = FALSE
)

germ.combined <- RunUMAP(
  germ.combined,
  reduction = "pca",
  dims = 1:28
)

germ.combined <- FindNeighbors(
  germ.combined,
  reduction = "pca",
  dims = 1:28
)

germ.combined <- FindClusters(
  germ.combined,
  resolution = 0.5
)

saveRDS(
  germ.combined,
  "cKO_germ_merge.rds"
)

############################################################
# 9. Remove somatic contaminants
############################################################

germ.clean <- subset(
  germ.combined,
  idents = c(
    "0","1","2","3",
    "4","5","6","8"
  )
)

new.cluster.ids <- as.character(0:7)

names(new.cluster.ids) <- levels(germ.clean)

germ.clean <- RenameIdents(
  germ.clean,
  new.cluster.ids
)

saveRDS(
  germ.clean,
  "cKO_germ_merge_remsoma.rds"
)

############################################################
# 10. RNA velocity preprocessing
############################################################

setwd("~")

############################################################
# Load loom files
############################################################

loom.files <- list(
  E14.5_ctrl = "E14.5_ctrl.loom",
  E14.5_cKO  = "E14.5_cKO.loom",
  E15.5_ctrl = "E15.5_ctrl.loom",
  E15.5_cKO  = "E15.5_cKO.loom",
  E16.5_cKO  = "E16.5_cKO.loom",
  E17.5_cKO  = "E17.5_cKO.loom"
)

loom.list <- lapply(
  loom.files,
  ReadVelocity
)

velo.list <- lapply(
  loom.list,
  as.Seurat
)

############################################################
# Rename velocity barcodes
############################################################

rename_velocity_cells <- function(
  seu,
  start,
  end
){

  newnames <- str_sub(
    colnames(seu),
    start,
    end
  )

  newnames <- paste0(
    newnames,
    "-1"
  )

  RenameCells(
    seu,
    new.names = newnames
  )
}

velo.list$E14.5_ctrl <- rename_velocity_cells(
  velo.list$E14.5_ctrl,
  8, 23
)

velo.list$E14.5_cKO <- rename_velocity_cells(
  velo.list$E14.5_cKO,
  8, 23
)

velo.list$E15.5_ctrl <- rename_velocity_cells(
  velo.list$E15.5_ctrl,
  12, 27
)

velo.list$E15.5_cKO <- rename_velocity_cells(
  velo.list$E15.5_cKO,
  8, 23
)

velo.list$E16.5_cKO <- rename_velocity_cells(
  velo.list$E16.5_cKO,
  8, 23
)

velo.list$E17.5_cKO <- rename_velocity_cells(
  velo.list$E17.5_cKO,
  8, 23
)

############################################################
# Match cells between Seurat and velocity objects
############################################################

for(i in names(velo.list)){

  common.cells <- intersect(
    colnames(seu.list[[i]]),
    colnames(velo.list[[i]])
  )

  seu.list[[i]] <- subset(
    seu.list[[i]],
    cells = common.cells
  )

  velo.list[[i]] <- subset(
    velo.list[[i]],
    cells = common.cells
  )

  velo.list[[i]]@assays$RNA <-
    seu.list[[i]]@assays$RNA

  velo.list[[i]]$orig.ident <- i
}

############################################################
# Save velocity-ready Seurat objects
############################################################

for(i in names(velo.list)){

  saveRDS(
    velo.list[[i]],
    paste0(i, "_velo.rds")
  )

}


#--------------------------------------------------------
# Load data
#--------------------------------------------------------

sample_paths <- list(
  E14.5_cKO  = "E14.5_cKO_velo.rds",
  E14.5_ctrl = "E14.5_ctrl_velo.rds",
  E15.5_cKO  = "E15.5_cKO_velo.rds",
  E15.5_ctrl = "E15.5_ctrl_velo.rds",
  E16.5      = "E16.5_cKO_velo.rds",
  E17.5      = "E17.5_cKO_velo.rds"
)

samples <- lapply(sample_paths, readRDS)

#--------------------------------------------------------
# Set RNA assay
#--------------------------------------------------------

samples <- lapply(samples, function(x) {
  DefaultAssay(x) <- "RNA"
  x
})

#--------------------------------------------------------
# Dnd1 filtering for cKO samples
#--------------------------------------------------------

cko_samples <- c("E14.5_cKO", "E15.5_cKO", "E16.5", "E17.5")

samples[cko_samples] <- lapply(samples[cko_samples], function(x) {
  subset(x, subset = Dnd1 < 0.5)
})

#--------------------------------------------------------
# Merge objects
#--------------------------------------------------------

combined <- merge(
  x = samples[[1]],
  y = samples[-1]
)

#========================================================
# Whole testis integration
#========================================================

obj.list <- SplitObject(combined, split.by = "orig.ident")

obj.list <- lapply(obj.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = 2000
  )
})

features <- SelectIntegrationFeatures(object.list = obj.list)

obj.list <- lapply(obj.list, function(x) {
  x <- ScaleData(
    x,
    features = features,
    verbose = FALSE
  )

  x <- RunPCA(
    x,
    features = features,
    verbose = FALSE
  )

  x
})

anchors <- FindIntegrationAnchors(
  object.list      = obj.list,
  anchor.features  = features,
  reduction        = "rpca",
  dims             = 1:30
)

combined.integrated <- IntegrateData(
  anchorset = anchors,
  dims      = 1:30
)

#--------------------------------------------------------
# Clustering
#--------------------------------------------------------

DefaultAssay(combined.integrated) <- "integrated"

combined.integrated <- ScaleData(
  combined.integrated,
  verbose = FALSE
)

combined.integrated <- RunPCA(
  combined.integrated,
  npcs    = 30,
  verbose = FALSE
)

combined.integrated <- RunUMAP(
  combined.integrated,
  reduction = "pca",
  dims      = 1:30
)

combined.integrated <- FindNeighbors(
  combined.integrated,
  reduction = "pca",
  dims      = 1:30
)

combined.integrated <- FindClusters(
  combined.integrated,
  resolution = 0.1
)

saveRDS(
  combined.integrated,
  "all_merge_velo.rds"
)

#--------------------------------------------------------
# Cell type annotation
#--------------------------------------------------------

new.cluster.ids <- c(
  "Interstitial/Stromal Cells",
  "Sertoli Cells",
  "PGCs/ECCs",
  "Interstitial/Stromal Cells",
  "Sertoli Cells",
  "Leydig Cells",
  "Interstitial/Stromal Cells",
  "Interstitial/Stromal Cells",
  "PGCs/ECCs",
  "Macrophages",
  "Endothelial Cells"
)

names(new.cluster.ids) <- levels(combined.integrated)

cKO_cell_anno <- RenameIdents(
  combined.integrated,
  new.cluster.ids
)

saveRDS(
  cKO_cell_anno,
  "cKO_merge_Dnd0.5_anno_velo.rds"
)


#========================================================
# Germ cell subclustering
#========================================================

combined.annotated <- readRDS(
  "cKO_merge_Dnd0.5_anno_velo.rds"
)

germ_merge <- subset(
  combined.annotated,
  idents = "Germ Cells/ECCs",
  subset =
    orig.ident %in% c(
      "E14.5_ctrl",
      "E14.5_cKO",
      "E15.5_ctrl",
      "E15.5_cKO"
    )
)

DefaultAssay(germ_merge) <- "spliced"

#--------------------------------------------------------
# Integration
#--------------------------------------------------------

germ.list <- SplitObject(
  germ_merge,
  split.by = "orig.ident"
)

germ.list <- lapply(germ.list, function(x) {

  x <- NormalizeData(x)

  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = 2000
  )

  x
})

germ.features <- SelectIntegrationFeatures(
  object.list = germ.list
)

germ.list <- lapply(germ.list, function(x) {

  x <- ScaleData(
    x,
    features = germ.features,
    verbose = FALSE
  )

  x <- RunPCA(
    x,
    features = germ.features,
    verbose = FALSE
  )

  x
})

germ.anchors <- FindIntegrationAnchors(
  object.list     = germ.list,
  anchor.features = germ.features,
  reduction       = "rpca",
  dims            = 1:20
)

germ.integrated <- IntegrateData(
  anchorset = germ.anchors
)

#--------------------------------------------------------
# Clustering
#--------------------------------------------------------

DefaultAssay(germ.integrated) <- "integrated"

germ.integrated <- ScaleData(
  germ.integrated,
  verbose = FALSE
)

germ.integrated <- RunPCA(
  germ.integrated,
  npcs    = 30,
  verbose = FALSE
)

germ.integrated <- RunUMAP(
  germ.integrated,
  reduction = "pca",
  dims      = 1:30
)

germ.integrated <- FindNeighbors(
  germ.integrated,
  reduction = "pca",
  dims      = 1:30
)

germ.integrated <- FindClusters(
  germ.integrated,
  resolution = 0.85
)

#--------------------------------------------------------
# Remove unwanted clusters
#--------------------------------------------------------

germ.integrated <- subset(
  germ.integrated,
  idents = c(
    "0", "1", "2", "3", "4",
    "5", "6", "7", "8", "9", "12"
  )
)

new.cluster.ids <- c(
  "0", "1", "2", "3", "4",
  "5", "6", "7", "8", "9", "10"
)

names(new.cluster.ids) <- levels(germ.integrated)

germ.integrated <- RenameIdents(
  germ.integrated,
  new.cluster.ids
)

#--------------------------------------------------------
# Save object
#--------------------------------------------------------

saveRDS(
  germ.integrated,
  "E14.5_E15.5_germ_merge_velo.rds"
)