# =========================================================
# Figure 5 / Figure S5 analysis
# scRNA-seq analysis of E13.5–E15.5 PGCs
# =========================================================

# =========================================================
# Load libraries
# =========================================================
library(Seurat)
library(dplyr)
library(reticulate)
library(scales)
library(ggplot2)
library(velocyto.R)
library(SeuratWrappers)

set.seed(1234)

# =========================================================
# Load datasets
# =========================================================
combined_germ <- readRDS("E14.5_E15.5_germ_merge_velo.rds")
cKO_germ <- readRDS("cKO_germ_merge_remsoma.rds")

E13.5 <- readRDS("E13.5_scRNA_sub.rds") #Nguyen D.H. et.al. 2021

# =========================================================
# Figure 4a
# Marker gene DotPlot
# =========================================================
DefaultAssay(cKO_germ) <- "RNA"

marker.genes <- c(
  "Trp53", "Lefty1", "Lefty2", "Nodal",
  "Nanos2", "Mael", "Rhox6", "Rhox9",
  "Kit", "Pecam1", "Car4"
)

DotPlot(
  cKO_germ,
  features = marker.genes,
  cols = c("blue", "red"),
  dot.scale = 8
) +
  RotatedAxis() +
  guides(
    size = guide_legend(title = "Pct.Exp"),
    color = guide_colorbar(title = "Avg.Exp")
  ) +
  theme(
    legend.text = element_text(size = 20),
    axis.title = element_blank(),
    axis.text = element_text(size = 17)
  )

# =========================================================
# Figure 4b-e
# Integration of E13.5 and E14.5/E15.5 PGCs
# =========================================================

# Define AP populations
E13.5_AP <- WhichCells(E13.5, idents = "5")
Fig5_AP <- WhichCells(combined_germ, idents = "8")

E13.5$orig.ident.v2 <- "E13.5_WT"

DefaultAssay(combined_germ) <- "RNA"

# Merge datasets
combined.list <- merge(
  x = combined_germ,
  y = E13.5
)

# Split by sample
combined.list <- SplitObject(
  combined.list,
  split.by = "orig.ident"
)

# Normalize and identify variable features
combined.list <- lapply(
  X = combined.list,
  FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = 2000
    )
  }
)

# Integration
anchors <- FindIntegrationAnchors(
  object.list = combined.list,
  dims = 1:20
)

combined <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

# Dimensional reduction and clustering
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)

combined <- RunPCA(
  combined,
  npcs = 30,
  verbose = FALSE
)

combined1 <- RunUMAP(
  combined,
  reduction = "pca",
  dims = 1:20
)

combined1 <- FindNeighbors(
  combined1,
  reduction = "pca",
  dims = 1:20
)

combined1 <- FindClusters(
  combined1,
  resolution = 0.7
)

# =========================================================
# Figure 4b
# UMAP colored by sample identity
# =========================================================
DimPlot(
  combined_germ,
  reduction = "umap",
  group.by = "orig.ident.v2"
) +
  NoAxes() +
  theme(
    plot.title = element_blank(),
    legend.text = element_text(size = 25),
    legend.position = c(0, 0.9)
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 10),
    ncol = 1
  ))

# =========================================================
# Figure 4c
# Cluster UMAP
# =========================================================

DimPlot(
  combined_germ,
  reduction = "umap",
  label = TRUE,
  label.size = 7
) +
  NoLegend() +
  NoAxes() +
  theme(
    plot.title = element_text(size = 20, vjust = -8)
  )

# =========================================================
# Figure 4d
# Highlight E14.5/E15.5 AP PGCs
# =========================================================

DimPlot(
  combined1,
  reduction = "umap",
  cells.highlight = Fig5_AP
) +
  NoLegend() +
  NoAxes() +
  theme(
    plot.title = element_text(size = 30, vjust = -8)
  )

# =========================================================
# Figure 4e
# Highlight E13.5 AP PGCs
# =========================================================

DimPlot(
  combined1,
  reduction = "umap",
  cells.highlight = E13.5_AP
) +
  NoLegend() +
  NoAxes() +
  theme(
    plot.title = element_text(size = 30, vjust = -8)
  )

# =========================================================
# Figure 4f
# Cluster composition tables
# =========================================================
table(
  combined_germ$seurat_clusters[
    grepl("E14.5_ctrl", combined_germ$orig.ident)
  ]
)

table(
  combined_germ$seurat_clusters[
    grepl("E14.5_cKO", combined_germ$orig.ident)
  ]
)

table(
  combined_germ$seurat_clusters[
    grepl("E15.5_ctrl", combined_germ$orig.ident)
  ]
)

table(
  combined_germ$seurat_clusters[
    grepl("E15.5_cKO", combined_germ$orig.ident)
  ]
)

table(combined_germ$orig.ident)

# =========================================================
# Figure 4g
# RNA velocity analysis
# =========================================================
DefaultAssay(combined_germ) <- "spliced"

combined_germ <- RunVelocity(
  object = combined_germ,
  deltaT = 1,
  kCells = 25,
  fit.quantile = 0.02
)

# Cell colors
ident.colors <- scales::hue_pal()(
  n = length(levels(combined_germ))
)

names(ident.colors) <- levels(combined_germ)

cell.colors <- ident.colors[
  Idents(combined_germ)
]

names(cell.colors) <- colnames(combined_germ)

# Velocity plot
show.velocity.on.embedding.cor(
  emb = Embeddings(
    combined_germ,
    reduction = "umap"
  ),
  vel = Tool(
    combined_germ,
    slot = "RunVelocity"
  ),
  n = 300,
  scale = "sqrt",
  cell.colors = ac(
    x = cell.colors,
    alpha = 0.5
  ),
  cex = 0.8,
  arrow.scale = 2,
  show.grid.flow = TRUE,
  min.grid.cell.mass = 0.5,
  grid.n = 40,
  arrow.lwd = 1,
  do.par = FALSE
)

# =========================================================
# Figure S5a
# E13.5 germ cell UMAP
# =========================================================
DimPlot(
  E13.5,
  reduction = "umap",
  label = TRUE,
  label.size = 9
) +
  NoLegend() +
  NoAxes() +
  theme(
    plot.title = element_text(size = 20, vjust = -8)
  )

# =========================================================
# Figure S5b
# E13.5 marker DotPlot
# =========================================================
DotPlot(
  E13.5,
  features = c(
    "Trp53", "Lefty1", "Lefty2", "Nodal",
    "Nanos2", "Mael", "Rhox6", "Rhox9",
    "Kit", "Pecam1"
  ),
  cols = c("blue", "red"),
  dot.scale = 8
) +
  RotatedAxis() +
  guides(
    size = guide_legend(title = "Pct.Exp"),
    color = guide_colorbar(title = "Avg.Exp")
  ) +
  theme(
    legend.text = element_text(size = 16),
    axis.title = element_blank(),
    axis.text = element_text(size = 16)
  )

# =========================================================
# Figure S5c
# E13.5 marker DotPlot
# =========================================================

DimPlot(combined_germ, reduction="umap", group.by = "orig.ident") +
NoAxes() +
theme(plot.title = element_blank(),
      legend.text = element_text(size = 20), 
      legend.position = c(0, 1),
      legend.justification = c(0, 1)) +
guides(color = guide_legend(override.aes = list(size=10), nrow=4))

# =========================================================
# Figure S5d
# =========================================================

DimPlot(combined_germ, 
        reduction="umap",
        label=T,
        label.size = 10,
        repel = T) + 
 NoLegend() + 
 NoAxes() + 
 theme(plot.title = element_text(size = 20,vjust = -8))

# =========================================================
# Figure S5e
# DotPlot for cKO germ cells
# =========================================================
DotPlot(
  combined_germ,
  features = c(
    "Trp53", "Lefty1", "Lefty2", "Nodal",
    "Nanos2", "Mael", "Rhox6", "Rhox9",
    "Kit", "Pecam1"
  ),
  cols = c("blue", "red"),
  dot.scale = 8
) +
  RotatedAxis() +
  guides(
    size = guide_legend(title = "Pct.Exp"),
    color = guide_colorbar(title = "Avg.Exp")
  ) +
  theme(
    legend.text = element_text(size = 16),
    axis.title = element_blank(),
    axis.text = element_text(size = 16)
  )