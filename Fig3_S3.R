# =========================================================
# Figure 3 & S3 analysis
# =========================================================

#--------------------------------------------------------
# Libraries
#--------------------------------------------------------

library(Seurat)
library(dplyr)
library(reticulate)
library(scales)
library(ggplot2)
library(reshape2)
library(patchwork)

set.seed(1234)

#========================================================
# Load data
#========================================================

cKO_germ <- readRDS(
  "cKO_germ_merge_remsoma.rds"
)

EC <- subset(
  cKO_germ,
  idents = c("4", "7")
)

DefaultAssay(EC) <- "RNA"

FSC <- readRDS(
  "FS_merge_before_scale.rds" #Kinoshita M et.al. 2021
)

ES <- readRDS(
  "merge_ES2i_Epi_blast_old.rds" #Chen et.al. 2016
)

E6.75 <- readRDS(
  "merge_6.75_sub_before_scale.rds" #Mohammed H et.al. 2017
)

E7.0 <- readRDS(
  "E7.0_2.rds" #Pijuan-Sala et.al. 2019
)

#========================================================
# Metadata formatting
#========================================================

cl.number <- paste(
  "cl",
  Idents(EC),
  sep = ""
)

EC$celltype.orig.ident <- paste(
  EC$orig.ident.v2,
  cl.number,
  sep = "_"
)

Idents(EC) <- "celltype.orig.ident"

FSC$celltype.orig.ident   <- FSC$cell_type
ES$celltype.orig.ident    <- ES$cell_type
E6.75$celltype.orig.ident <- E6.75$cell_type
E7.0$celltype.orig.ident  <- E7.0$celltype

#--------------------------------------------------------
# Rename E7.0 labels
#--------------------------------------------------------

current.labels <- c(
  "Epiblast",
  "Primitive Streak",
  "Anterior Primitive Streak"
)

new.labels <- c(
  "E7.0 ~ E7.5",
  "Primitive Streak",
  "Anterior Primitive Streak"
)

E7.0@meta.data$celltype.orig.ident <- plyr::mapvalues(
  x    = E7.0@meta.data$celltype.orig.ident,
  from = current.labels,
  to   = new.labels
)

#--------------------------------------------------------
# Downsample E7.0 cells
#--------------------------------------------------------

E7.0 <- subset(
  E7.0,
  cells = sample(
    Cells(E7.0),
    300
  )
)

unique(E7.0@meta.data$celltype.orig.ident)

#========================================================
# Merge datasets
#========================================================

obj.list <- merge(
  x = EC,
  y = c(
    E7.0,
    FSC,
    ES,
    E6.75
  )
)

obj.list <- SplitObject(
  obj.list,
  split.by = "orig.ident.v2"
)

#========================================================
# Add metadata
#========================================================

#--------------------------------------------------------
# Current study
#--------------------------------------------------------

obj.list$E14.5_cKO$study <- "current study"
obj.list$E15.5_cKO$study <- "current study"
obj.list$E16.5_cKO$study <- "current study"
obj.list$E17.5_cKO$study <- "current study"

obj.list$E14.5_cKO$cell_type <- "E14.5_cKO_cl7"
obj.list$E15.5_cKO$cell_type <- "E15.5_cKO_cl7"
obj.list$E16.5_cKO$cell_type <- "E16.5_cKO_cl7"
obj.list$E17.5_cKO$cell_type <- "E17.5_cKO_cl4,7"

#--------------------------------------------------------
# Reference datasets
#--------------------------------------------------------

obj.list$E7.0$cell_type <- E7.0$celltype

obj.list$FSCs$study  <- "Kinoshita M, et al."
obj.list$ES$study    <- "Chen G, et al."
obj.list$E6.75$study <- "Mohammed H, et al."
obj.list$E7.0$study  <- "Pijuan-Sala B, et.al."

#========================================================
# Data integration
#========================================================

obj.list <- lapply(
  X = obj.list,
  FUN = function(x){

    x <- NormalizeData(x)

    x <- FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = 2000
    )

    x
  }
)

anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  dims        = 1:8,
  k.filter    = 31
)

combined <- IntegrateData(
  anchorset = anchors,
  dims      = 1:20,
  k.weight  = 31
)

#========================================================
# Clustering
#========================================================

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(
  combined,
  verbose = FALSE
)

combined <- RunPCA(
  combined,
  npcs    = 30,
  verbose = FALSE
)

combined0 <- RunUMAP(
  combined,
  reduction = "pca",
  dims      = 1:30
)

combined0 <- FindNeighbors(
  combined0,
  reduction = "pca",
  dims      = 1:30
)

combined0 <- FindClusters(
  combined0,
  resolution = 0.1
)

#--------------------------------------------------------
# Rename celltype labels
#--------------------------------------------------------

combined0$celltype.orig.ident2 <-
  combined0$celltype.orig.ident

current.labels <- c(
  "E14.5_Dnd1-cKO_cl7",
  "E15.5_Dnd1-cKO_cl7",
  "E16.5_Dnd1-cKO_cl7",
  "E17.5_Dnd1-cKO_cl7",
  "E17.5_Dnd1-cKO_cl4",
  "E3.5",
  "E4.5",
  "E5.5",
  "E6.75",
  "Primitive Streak",
  "Anterior Primitive Streak",
  "ES_2iL",
  "FSC",
  "EpiSC"
)

new.labels2 <- c(
  "Cluster7",
  "Cluster7",
  "Cluster7",
  "Cluster7",
  "Cluster4",
  "E3.5",
  "E4.5",
  "E5.5",
  "E6.75",
  "PS",
  "APS",
  "ES_2iL",
  "FSC",
  "EpiSC"
)

combined1@meta.data$celltype.orig.ident2 <-
  plyr::mapvalues(
    x    = combined1@meta.data$celltype.orig.ident2,
    from = current.labels,
    to   = new.labels2
  )

combined1$celltype.orig.ident2 <- factor(
  combined1$celltype.orig.ident2,
  levels = c(
    "Cluster7",
    "Cluster4",
    "E3.5",
    "E4.5",
    "E5.5",
    "E6.75",
    "Primitive Streak",
    "Anterior Primitive Streak",
    "ES_2iL",
    "FSC",
    "EpiSC"
  )
)

#========================================================
# Figure 3a
# DotPlot
#========================================================

DefaultAssay(combined1) <- "RNA"

DotPlot(
  combined1,
  features = c(
    "Tbx3", "Tfcp2l1", "Esrrb",
    "Sox4", "Pou3f1", "Fgf5",
    "Hhex", "Foxa2", "Mixl1"
  ),
  group.by = "celltype.orig.ident3",
  cols      = c("blue", "red"),
  dot.scale = 10
) +
  RotatedAxis() +
  guides(
    size = guide_legend(
      title = "Pct.Exp",
      order = 1
    ),

    color = guide_colorbar(
      title = "Avg.Exp",
      order = 2
    )
  ) +
  theme(
    legend.text = element_text(
      family = "Arial",
      size   = 35
    ),

    axis.title = element_blank(),

    axis.text = element_text(
      family = "Arial",
      size   = 30
    ),

    legend.title = element_text(
      family = "Arial",
      size   = 35
    ),

    legend.key.size = unit(1.2, "cm")
  ) +
  scale_y_discrete(
    limits = rev(c(
      "Cluster7",
      "Cluster4",
      "E3.5",
      "E4.5",
      "E5.5",
      "E6.75",
      "PS",
      "APS",
      "ES_2iL",
      "FSC",
      "EpiSC"
    ))
  )

#========================================================
# Figure 3b-e
# UMAP visualization
#========================================================

combined1$celltype.orig.ident3 <- factor(
  combined1$celltype.orig.ident3,
  levels = c(
    "Cluster7",
    "Cluster4",
    "E3.5",
    "E4.5",
    "E5.5",
    "E6.75",
    "PS",
    "APS",
    "ES_2iL",
    "FSC",
    "EpiSC"
  )
)


base_theme <- theme(
  plot.title = element_text(family = "Arial", size = 40, hjust = 0),
  legend.position = "left",
  legend.text = element_text(family = "Arial", size = 38)
)

make_dimplot <- function(cols, title, shuffle = FALSE) {
  DimPlot(
    combined1,
    reduction = "umap",
    cols = cols,
    group.by = "celltype.orig.ident3",
    pt.size = 2.5,
    shuffle = shuffle
  ) +
    NoLegend() +
    NoAxes() +
    base_theme +
    labs(title = title) +
    guides(color = guide_legend(
      override.aes = list(size = 10),
      ncol = 1
    ))
}


cols_embryo_c7 <- c(
  "#FFCC00","grey85","#FF0000","#FFFF00",
  "#00FFFF","#1478FF","#ec44ec","#25f073",
  "grey85","grey85","grey85"
)

cols_embryo_c4 <- c(
  "grey85","#e96417","#FF0000","#FFFF00",
  "#00FFFF","#1478FF","#ec44ec","#25f073",
  "grey85","grey85","grey85"
)

cols_culture_c7 <- c(
  "#FFCC00","grey85","grey85","grey85",
  "grey85","grey85","grey85","grey85",
  "#00FF00","#a052e0","#62b962"
)

cols_culture_c4 <- c(
  "grey85","#e96417","grey85","grey85",
  "grey85","grey85","grey85","grey85",
  "#00FF00","#a052e0","#62b962"
)

p1 <- make_dimplot(cols_embryo_c7, "Embryo+Cluster7", shuffle = TRUE)
p2 <- make_dimplot(cols_embryo_c4, "Embryo+Cluster4", shuffle = TRUE)

p3 <- make_dimplot(cols_culture_c7, "Culture+Cluster7", shuffle = TRUE)
p4 <- make_dimplot(cols_culture_c4, "Culture+Cluster4", shuffle = TRUE)

(p1 | p3) /
(p2 | p4) 

#========================================================
# Figure 3f
#========================================================

DimPlot(
  combined1,
  reduction  = "umap",
  label      = TRUE,
  label.size = 15,
  repel      = TRUE
) +
  NoLegend() +
  NoAxes() +
  theme(
    text = element_text(
      family = "Arial",
      size   = 20,
      vjust  = -8
    )
  )

#========================================================
# Figure 3g
# Cluster composition
#========================================================

table(combined1$seurat_clusters[
  grepl("cl7", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("cl4", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("E3.5", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("E4.5", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("E5.5", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("E6.75", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("Primitive", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("Anterior", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("ES_2iL", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("FSC", combined1$celltype.orig.ident2)
])

table(combined1$seurat_clusters[
  grepl("EpiSC", combined1$celltype.orig.ident2)
])

#========================================================
# Relative proportion plot
#========================================================

proportion <- read.csv(
  "embryo_merge_proportion.csv",
  row.names = 1,
  header = TRUE
)

proportion <- proportion[c(10, 11, 1:9), ]

proportion$stage <- rownames(proportion)

proportion.long <- melt(
  proportion,
  id.vars      = "stage",
  variable.name = "Cluster",
  value.name    = "Count"
)

proportion.long <- proportion.long %>%
  group_by(stage) %>%
  mutate(
    Relative = Count / sum(Count)
  )

proportion.long$stage <- factor(
  proportion.long$stage,
  levels = c(
    "Cluster7",
    "Cluster4",
    "E3.5",
    "E4.5",
    "E5.5",
    "E6.75",
    "PS",
    "APS",
    "ES_2i/LIF",
    "FSC",
    "EpiSC"
  )
)

ggplot(
  proportion.long,
  aes(
    x    = stage,
    y    = Relative,
    fill = Cluster
  )
) +
  geom_bar(
    stat     = "identity",
    position = "fill",
    color    = "black"
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = NULL,
    x     = NULL,
    y     = "Relative Proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),

    axis.text = element_text(size = 14),

    axis.title.x = element_text(size = 16),

    axis.title.y = element_text(size = 16),

    plot.title = element_text(
      size = 16,
      face = "bold"
    ),

    legend.text = element_text(size = 14),

    legend.title = element_text(size = 16)
  ) +
  scale_fill_manual(
    values = c(
      "cluster0" = "tomato",
      "cluster1" = "lightgreen",
      "cluster2" = "skyblue",
      "cluster3" = "#b172d8"
    )
  )

#========================================================
# Supplementary Figure 3a
#========================================================

DimPlot(
  combined1,
  reduction = "umap",
  group.by  = "study"
) +
  NoAxes() +
  theme(
    plot.title = element_blank(),

    legend.text = element_text(size = 10)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 7),
      ncol = 1
    )
  )


#========================================================
# Supplementary Figure 3b-d
#========================================================

base_theme <- theme(
  plot.title = element_text(size = 30),
  legend.position = "left",
  legend.text = element_text(size = 30)
)


make_dimplot <- function(cols, title) {
  DimPlot(
    combined1,
    reduction = "umap",
    cols = cols,
    group.by = "celltype.orig.ident3",
    pt.size = 2.5
  ) +
    NoLegend() +
    NoAxes() +
    base_theme +
    labs(title = title) +
    guides(color = guide_legend(
      override.aes = list(size = 10),
      ncol = 1
    ))
}


cols_es2il <- c(
  "grey85","grey85","#FF0000","#FFFF00",
  "#00FFFF","#1478FF","#ec44ec","#25f073",
  "#00FF00","grey85","grey85"
)

cols_fsc <- c(
  "grey85","grey85","#FF0000","#FFFF00",
  "#00FFFF","#1478FF","#ec44ec","#25f073",
  "grey85","#a052e0","grey85"
)

cols_episc <- c(
  "grey85","grey85","#FF0000","#FFFF00",
  "#00FFFF","#1478FF","#ec44ec","#25f073",
  "grey85","grey85","#62b962"
)

p1 <- make_dimplot(cols_es2il, "Embryo+ES_2iL")
p2 <- make_dimplot(cols_fsc, "Embryo+FSC")
p3 <- make_dimplot(cols_episc, "Embryo+EpiSC")

p1 | p2 | p3