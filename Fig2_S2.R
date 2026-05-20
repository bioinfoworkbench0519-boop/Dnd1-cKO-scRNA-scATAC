#========================================================
# Figure 2 & S2 analysis
#========================================================

#--------------------------------------------------------
# Libraries
#--------------------------------------------------------

library(Seurat)
library(dplyr)
library(reticulate)
library(scales)
library(ggplot2)
library(rtracklayer)
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(patchwork)

set.seed(1234)

#========================================================
# Load data
#========================================================

cKO_germ <- readRDS("cKO_germ_merge_remsoma.rds")

#========================================================
# Figure 2a
# UMAP grouped by developmental stage
#========================================================

cKO_germ$orig.ident.v3 <- cKO_germ$orig.ident.v2

current.labels <- c(
  "E14.5_Dnd1-cKO",
  "E15.5_Dnd1-cKO",
  "E16.5_Dnd1-cKO",
  "E17.5_Dnd1-cKO"
)

new.labels <- c(
  "E14.5",
  "E15.5",
  "E16.5",
  "E17.5"
)

cKO_germ@meta.data$orig.ident.v3 <- plyr::mapvalues(
  x    = cKO_germ@meta.data$orig.ident.v3,
  from = current.labels,
  to   = new.labels
)

cKO_germ$orig.ident.v3 <- factor(
  cKO_germ$orig.ident.v3,
  levels = c("E14.5", "E15.5", "E16.5", "E17.5")
)

DimPlot(
  cKO_germ,
  reduction = "umap",
  group.by  = "orig.ident.v3"
) +
  NoAxes() +
  theme(
    plot.title = element_blank(),

    legend.text = element_text(size = 25),

    legend.position = "top"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 10)
    )
  )

#========================================================
# Figure 2b
# Annotated UMAP
#========================================================

DimPlot(
  cKO_germ,
  reduction  = "umap",
  label      = TRUE,
  label.size = 13,
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
# Figure 2c
# Cell number summary
#========================================================

E17.5 <- subset(
  cKO_germ,
  subset = orig.ident == "E17.5_cKO"
)

E16.5 <- subset(
  cKO_germ,
  subset = orig.ident == "E16.5_cKO"
)

E15.5 <- subset(
  cKO_germ,
  subset = orig.ident == "E15.5_cKO"
)

E14.5 <- subset(
  cKO_germ,
  subset = orig.ident == "E14.5_cKO"
)

table(
  cKO_germ$seurat_clusters[
    grepl("E17.5", cKO_germ$orig.ident)
  ]
)

table(
  cKO_germ$seurat_clusters[
    grepl("E16.5", cKO_germ$orig.ident)
  ]
)

table(
  cKO_germ$seurat_clusters[
    grepl("E15.5", cKO_germ$orig.ident)
  ]
)

table(
  cKO_germ$seurat_clusters[
    grepl("E14.5", cKO_germ$orig.ident)
  ]
)

length(E17.5@meta.data$orig.ident)
length(E16.5@meta.data$orig.ident)
length(E15.5@meta.data$orig.ident)
length(E14.5@meta.data$orig.ident)

#========================================================
# Figure 2d
# DotPlot across all stages
#========================================================

marker.genes.all <- c(
  "Dazl", "Ddx4", "Piwil2",
  "Nanog", "Lin28a", "Sox2",
  "Tbx3", "Tfcp2l1", "Esrrb",
  "Nodal", "Fgf5", "Otx2",
  "Cer1", "T", "Mixl1"
)

DotPlot(
  cKO_germ,
  features  = marker.genes.all,
  cols      = c("blue", "red"),
  dot.scale = 8,
  scale.max = 100
) +
  RotatedAxis() +
  guides(
    size = guide_legend(
      title = "Pct.Exp",
      order = 2
    ),

    color = guide_colorbar(
      title = "Avg.Exp",
      order = 1
    )
  ) +
  theme(
    legend.text = element_text(size = 20),

    axis.title = element_blank(),

    axis.text = element_text(size = 20)
  )

#========================================================
# Figure 2e-h
# Stage-specific DotPlots
#========================================================

marker.genes.early <- c(
  "Dazl", "Ddx4", "Piwil2",
  "Nanog", "Lin28a", "Sox2",
  "Tbx3", "Tfcp2l1", "Esrrb",
  "Nodal", "Fgf5", "Otx2"
)

marker.genes.late <- c(
  marker.genes.early,
  "Cer1", "T", "Mixl1"
)

#--------------------------------------------------------
# E14.5
#--------------------------------------------------------

DotPlot(
  E14.5,
  features  = marker.genes.early,
  cols      = c("blue", "red"),
  dot.scale = 8
) +
  RotatedAxis() +
  guides(
    size = guide_legend(
      title = "Pct.Exp",
      order = 2
    ),

    color = guide_colorbar(
      title = "Avg.Exp",
      order = 1
    )
  ) +
  theme(
    legend.text = element_text(size = 20),

    axis.title = element_blank(),

    axis.text = element_text(size = 20)
  )

#--------------------------------------------------------
# E15.5
#--------------------------------------------------------

DotPlot(
  E15.5,
  features  = marker.genes.early,
  cols      = c("blue", "red"),
  dot.scale = 8
) +
  RotatedAxis() +
  guides(
    size = guide_legend(
      title = "Pct.Exp",
      order = 2
    ),

    color = guide_colorbar(
      title = "Avg.Exp",
      order = 1
    )
  ) +
  theme(
    legend.text = element_text(size = 20),

    axis.title = element_blank(),

    axis.text = element_text(size = 20)
  )

#--------------------------------------------------------
# E16.5
#--------------------------------------------------------

DotPlot(
  E16.5,
  features  = marker.genes.early,
  cols      = c("blue", "red"),
  dot.scale = 8,
  scale.max = 100
) +
  RotatedAxis() +
  guides(
    size = guide_legend(
      title = "Pct.Exp",
      order = 2
    ),

    color = guide_colorbar(
      title = "Avg.Exp",
      order = 1
    )
  ) +
  theme(
    legend.text = element_text(size = 20),

    axis.title = element_blank(),

    axis.text = element_text(size = 20)
  )

#--------------------------------------------------------
# E17.5
#--------------------------------------------------------

DotPlot(
  E17.5,
  features  = marker.genes.late,
  cols      = c("blue", "red"),
  dot.scale = 8
) +
  RotatedAxis() +
  guides(
    size = guide_legend(
      title = "Pct.Exp",
      order = 2
    ),

    color = guide_colorbar(
      title = "Avg.Exp",
      order = 1
    )
  ) +
  theme(
    legend.text = element_text(size = 20),

    axis.title = element_blank(),

    axis.text = element_text(size = 20)
  )