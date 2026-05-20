#========================================================
# Figure 1 analysis
#========================================================

#--------------------------------------------------------
# Libraries
#--------------------------------------------------------

library(Seurat)
library(dplyr)
library(reticulate)
library(scales)
library(ggplot2)
library(ggrepel)
library(patchwork)

set.seed(1234)

#========================================================
# Load data
#========================================================

setwd("~")

all_cell_anno <- readRDS("all_merge_cellanno.rds")

#========================================================
# Figure 1b
# UMAP by sample
#========================================================

DimPlot(
  all_cell_anno,
  reduction  = "umap",
  group.by   = "orig.ident.v2",
  label.size = 7,
  repel      = TRUE
) +
  NoAxes() +
  theme(
    plot.title = element_blank(),
    legend.text = element_text(size = 17)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 10),
      ncol = 1
    )
  )


#========================================================
# FeaturePlot function
#========================================================

Featureforanno <- function(gene, celltype){

  FeaturePlot(
    all_cell_anno,
    features   = gene,
    cols       = c("lightgrey", "red"),
    max.cutoff = 2.5
  ) +
    NoAxes() +
    theme(
      plot.title = element_text(size = 24),

      plot.subtitle = element_text(
        hjust = 0.5,
        size  = 24,
        face  = "italic"
      ),

      legend.position      = c(0, 1),
      legend.justification = c(0, 1.05),

      legend.text = element_text(size = 15)
    ) +
    labs(
      title    = celltype,
      subtitle = gene
    )
}

#========================================================
# Figure 1c
# Marker gene FeaturePlots
#========================================================

p1  <- Featureforanno("Pou5f1",  "PGCs/ECCs")
p2  <- Featureforanno("Dazl",    "PGCs")
p3  <- Featureforanno("Amh",     "Sertoli Cells")
p4  <- Featureforanno("Sox9",    "Sertoli Cells")
p5  <- Featureforanno("Cyp11a1", "Leydig Cells")
p6  <- Featureforanno("Hsd3b1",  "Leydig Cells")
p7  <- Featureforanno("Igf1",    "Interstitial/Stromal Cells")
p8  <- Featureforanno("Acta2",   "Interstitial/Stromal Cells")
p9  <- Featureforanno("Pecam1",  "Endothelial Cells")
p10 <- Featureforanno("Esam",    "Endothelial Cells")
p11 <- Featureforanno("Cd68",    "Macrophages")
p12 <- Featureforanno("Lyz2",    "Macrophages")

(p1 | p3 | p5) /
(p2 | p4 | p6) /
(p7 | p9 | p11) /
(p8 | p10 | p12)

#========================================================
# Figure 1d
# Annotated UMAP
#========================================================

DimPlot(
  all_cell_anno,
  reduction  = "umap",
  label      = TRUE,
  label.size = 14,
  repel      = TRUE
) +
  NoLegend() +
  NoAxes() +
  theme(
    plot.title = element_text(
      size  = 20,
      vjust = -8
    )
  )

#========================================================
# Figure 1e
# DEG heatmap
#========================================================

DefaultAssay(all_cell_anno) <- "RNA"

markers <- FindAllMarkers(
  all_cell_anno,
  only.pos        = TRUE,
  min.pct         = 0.3,
  logfc.threshold = 0.25
)

top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

all_cell_anno <- ScaleData(all_cell_anno)

DoHeatmap(
  subset(all_cell_anno, downsample = 50),
  features = top10$gene
)