# =========================================================
# Figure 6 / Figure S6 analysis
# scATAC-seq analysis using ArchR and Monocle3
# =========================================================

# =========================================================
# Load libraries
# =========================================================
library(ArchR)
library(monocle3)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(BSgenome.Mmusculus.UCSC.mm10)

# =========================================================
# Initial settings
# =========================================================
set.seed(1234)

# Mouse genome
addArchRGenome("mm10")

# =========================================================
# Load datasets
# =========================================================
all_cell_anno <- readRDS("all_merge_cellanno.rds")

proj <- loadArchRProject(path = "Merge")

proj_EG <- loadArchRProject(path = "Merge_EG2")

# =========================================================
# Figure 6a-f
# =========================================================

proj_EG$Sample1 <- mapLabels(
  proj_EG$Sample,
  newLabels = c("E14.5", "E15.5", "E16.5", "E17.5"),
  oldLabels = c(
    "E14.5_Dnd1-cKO",
    "E15.5_Dnd1-cKO",
    "E16.5_Dnd1-cKO",
    "E17.5_Dnd1-cKO"
  )
)

# ---------------------------------------------------------
# Figure 6a
# ---------------------------------------------------------

Monocle_ATP <- getMonocleTrajectories(
  ArchRProj = proj_EG,
  name = "Monocle_ATP",
  useGroups = paste0("C", 1:21),
  principalGroup = "C20",
  groupBy = "Clusters_EG2",
  embedding = "UMAP_EG1",
  clusterParams = list(k = 50),
  seed = 1
)

proj_EG <- addMonocleTrajectory(
  ArchRProj = proj_EG,
  name = "Monocle_ATP",
  useGroups = paste0("C", 1:21),
  groupBy = "Clusters_EG2",
  monocleCDS = Monocle_ATP,
  force = TRUE
)

Monocle_ATP <- monocle3::learn_graph(
  Monocle_ATP,
  learn_graph_control = list(
    minimal_branch_len = 5
  )
)

Monocle_ATP <- monocle3::order_cells(
  Monocle_ATP,
  reduction_method = "UMAP"
)


monocle3::plot_cells(
  Monocle_ATP,
  color_cells_by = "Sample1",
  group_cells_by = "Sample1",
  label_branch_points = FALSE,
  label_cell_groups = FALSE,
  cell_size = 0.6,
  label_roots = FALSE,
  label_leaves = FALSE,
  show_trajectory_graph = FALSE,
  trajectory_graph_segment_size = 1,
  graph_label_size = 4,
  rasterize = TRUE
) +
  theme(
    panel.border = element_blank(),
    plot.title = element_blank(),
    legend.text = element_text(
      family = "Arial",
      size = 15
    ),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  NoAxes() +
  guides(
    color = guide_legend(
      nrow = 2,
      override.aes = list(size = 5)
    )
  ) +
  scale_color_manual(
    values = c("red", "blue", "green", "purple")
  )

# =========================================================
# Figure 6b
# Cluster visualization
# =========================================================

proj_EG$Clusters_EG2_R <- mapLabels(
  proj_EG$Clusters_EG2,
  newLabels = c(
    "C1","C2","C3","C4","C5","C6","C7",
    "C8","C9","C10","C11","C12","C13",
    "C14","C15","C16","C17","C18","C19",
    "C20","C21"
  ),
  oldLabels = c(
    "C17","C16","C15","C21","C19","C20",
    "C18","C12","C11","C10","C13","C14",
    "C2","C3","C4","C6","C5","C1",
    "C9","C8","C7"
  )
)

proj_EG@cellColData$Clusters_EG2_R <- gsub(
  "C",
  "",
  proj_EG@cellColData$Clusters_EG2_R
)

pal <- c(
  "1"="red",
  "2"="#9595f7e2",
  "3"="#2a8d2a",
  "4"="#8a14d3",
  "5"="#e96d08",
  "6"="#e0d100",
  "7"="#9fbe9f",
  "8"="#957aa5",
  "9"="#e26757",
  "10"="#9b9bd1",
  "11"="#7deb7d",
  "12"="#a8a86e",
  "13"="#f797a7",
  "14"="orange",
  "15"="#9e5c27",
  "16"="skyblue",
  "17"="#7aeb7a",
  "18"="#5757e2",
  "19"="brown",
  "20"="#8f8f00",
  "21"="black"
)

monocle3::plot_cells(
  Monocle_ATP,
  color_cells_by = "Clusters_EG2_R",
  group_cells_by = "Clusters_EG2_R",
  label_branch_points = FALSE,
  group_label_size = 10,
  cell_size = 0.6,
  label_roots = FALSE,
  label_leaves = FALSE,
  show_trajectory_graph = FALSE,
  trajectory_graph_segment_size = 1,
  graph_label_size = 4
) +
  theme(
    panel.border = element_blank(),
    plot.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_blank()
  ) +
  NoAxes() +
  guides(
    color = guide_legend(
      nrow = 2,
      override.aes = list(size = 5)
    )
  ) +
  scale_color_manual(values = pal)

# ---------------------------------------------------------
# Figure 6c
# ---------------------------------------------------------

monocle3::plot_cells(
  Monocle_ATP,
  color_cells_by = "pseudotime",
  label_branch_points = FALSE,
  label_roots = TRUE,
  label_leaves = FALSE,
  trajectory_graph_color = "blue",
  trajectory_graph_segment_size = 1,
  graph_label_size = 4
) +
  theme(
    legend.title = element_text(size = 23),
    legend.text = element_text(size = 20),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1, "cm")
  ) +
  NoAxes() +
  scale_color_viridis_c(
    limits = c(0, 6.5),
    oob = scales::squish,
    option = "plasma"
  ) +
  labs(color = "Pseudotime")

# =========================================================
# Figure 6d
# =========================================================

monocle3::plot_cells(Monocle_ATP,
                     color_cells_by = "Clusters_EG2_R", 
                     group_cells_by = "Clusters_EG2_R", 
                     label_branch_points = FALSE, 
                     group_label_size = 10, 
                     cell_size = 0.6,
                     label_roots = TRUE, 
                     label_leaves = FALSE, 
                     trajectory_graph_color = "blue",
                     trajectory_graph_segment_size = 1, 
                     graph_label_size = 4) +
          theme(panel.border = element_blank(), 
                plot.title = element_blank(),
                legend.text = element_text(size = 15),
                legend.title = element_blank()) + 
          NoAxes() + 
          guides(color = guide_legend(nrow = 2, 
                                      override.aes = list(size = 5))) +
          scale_color_manual(values=pal)

# =========================================================
# Figure 6e
# =========================================================

cKO_germ$ident <- Idents(cKO_germ)

proj_EG <- addGeneIntegrationMatrix(
  ArchRProj = proj_EG,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_EG1",
  reducedDims = "IterativeLSI_EG1",
  seRNA = cKO_germ,
  addToArrow = FALSE,
  groupRNA = "ident",
  dimsToUse = 1:25,
  groupATAC = "Clusters_EG2",
  nameCell = "predictedCell_Un_EG",
  nameGroup = "predictedGroup_Un_EG",
  nameScore = "predictedScore_Un_EG"
)

monocle3::plot_cells(Monocle_ATP,
                     color_cells_by = "predictedGroup_Un_EG", 
                     group_cells_by = "predictedGroup_Un_EG", 
                     label_branch_points = FALSE, 
                     group_label_size = 10, 
                     cell_size = 0.6,
                     label_roots = TRUE, 
                     label_leaves = FALSE, 
                     trajectory_graph_color = "blue",
                     trajectory_graph_segment_size = 1, 
                     graph_label_size = 4) +
          theme(panel.border = element_blank(), 
                plot.title = element_blank(),
                legend.text = element_text(size = 15),
                legend.title = element_blank()) + 
          NoAxes() + 
          guides(color = guide_legend(nrow = 2, 
                                      override.aes = list(size = 5))) +
          scale_color_manual(values=c("#F8766D",
                                      "#CD9600",
                                      "#7CAE00",
                                      "#00BE67",
                                      "#00BFC4",
                                      "#00A9FF", 
                                      "#C77CFF",
                                      "#FF61CC"))

# =========================================================
# Figure 6f
# =========================================================

proj_EG <- addImputeWeights(
  proj_EG,
  reducedDims = "IterativeLSI_EG1"
)

plotEmbedding(ArchRProj = proj_EG,
           colorBy = "GeneScoreMatrix",
           size = 0.6,
           name = "Car4", 
           embedding = "UMAP_EG1", 
           imputeWeights = getImputeWeights(proj_EG),
           plotAs = "points") +
 theme(panel.border = element_blank(),
       plot.title = element_blank(),
       legend.text = element_text(family = "Arial", size = 15),
       legend.title = element_text(family = "Arial",size = 20)) + 
 NoAxes()

# =========================================================
# Figure 6g
# Gene score visualization
# =========================================================

plotGroups(
  ArchRProj = proj_EG,
  groupBy = "Clusters_EG2_R",
  colorBy = "GeneScoreMatrix",
  name = "Car4",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  imputeWeights = getImputeWeights(proj_EG)
) +
  theme(
    plot.title = element_blank(),
    legend.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 18)
  ) +
  labs(
    x = "Cluster",
    y = "Log2(NormCounts + 1)"
  )

# =========================================================
# Figure 6h
# =========================================================

proj_EG <- addGroupCoverages(
  ArchRProj = proj_EG,
  groupBy = "Clusters_EG2"
)

pathToMacs2 <- findMacs2()

proj_EG <- addReproduciblePeakSet(
  ArchRProj = proj_EG,
  groupBy = "Clusters_EG2",
  pathToMacs2 = pathToMacs2
)

proj_EG <- addPeakMatrix(proj_EG)

proj_EG <- addMotifAnnotations(
  ArchRProj = proj_EG,
  motifSet = "cisbp",
  name = "Motif"
)

markerTest <- getMarkerFeatures(
  ArchRProj = proj_EG,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters_EG2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = proj_EG,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

mlogp_mat <- assay(
  enrichMotifs,
  "mlog10Padj"
)

mlogp_mat_C4 <- mlogp_mat[C4g, ]

rownames(mlogp_mat_C4) <- gsub(
  "_.*",
  "",
  rownames(mlogp_mat_C4)
)

rownames(mlogp_mat_C4) <- toupper(
  rownames(mlogp_mat_C4)
)

normalize_row <- function(x) {
  (x - min(x)) /
    (max(x) - min(x)) * 100
}

mlogp_mat_C4_NM <- t(
  t(
    apply(
      mlogp_mat_C4,
      1,
      normalize_row
    )
  )
)

col_fun <- colorRamp2(
  c(0, 50, 100),
  c("White", "blue", "black")
)

ht <- Heatmap(
  mlogp_mat_C4_NM,
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 12),
  cluster_rows = FALSE
)


ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot"
)

# =========================================================
# Figure 6i
# =========================================================

DefaultAssay(cKO_germ) <- "RNA"

DotPlot(
  cKO_germ,
  features = c(
    "Pitx3","Pitx2","Pitx1","Dmbx1",
    "Obox3","Pnp","Otx1","Otx2",
    "Crx","Gsc","Obox6","Zfx"
  ),
  cols = c("blue", "red"),
  dot.scale = 8,
  scale.min = 0,
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
    axis.text = element_text(size = 17)
  )

# =========================================================
# Figure S6a
# =========================================================

plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP",
  labelAsFactors = FALSE,
  labelMeans = FALSE
) +
  theme(
    panel.border = element_blank(),
    plot.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_blank()
  ) +
  NoAxes() +
  guides(
    color = guide_legend(
      nrow = 2,
      override.aes = list(size = 5)
    )
  )

# =========================================================
# Figure S6b
# =========================================================

new.cluster.ids <- c(
  "Interstitial/Stromal Cells",
  "Sertoli Cells",
  "PGCs/ECCs",
  "Leydig Cells",
  "Macrophages",
  "Endothelial Cells"
)

names(new.cluster.ids) <- levels(all_cell_anno)

all_cell_anno <- RenameIdents(
  all_cell_anno,
  new.cluster.ids
)

all_cell_anno$ident <- Idents(all_cell_anno)

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = all_cell_anno,
  addToArrow = FALSE,
  groupRNA = "ident",
  dimsToUse = 1:30,
  groupATAC = "Clusters",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

plotEmbedding(proj,
              colorBy = "cellColData",
              name = "predictedGroup_Un",
              labelMeans = FALSE,
              labelAsFactors = FALSE,
              rastr = TRUE,
              labelSize = 7,
              pal = c("#F564E3",
                      "#F8766D",
                      "#00BFC4",
                      "#619CFF",
                      "#00BA38",
                      "#879F00")) +
      theme(panel.border = element_blank(),
            plot.title = element_blank(),
            legend.text = element_text(size = 15),
            legend.title = element_blank()) + 
      NoAxes() + 
      guides(color = guide_legend(nrow = 3, 
             override.aes = list(size = 5)))


# =========================================================
# Figure S6c
# =========================================================

cell_data <- getCellColData(
  ArchRProj = proj,
  select = c("Sample", "predictedGroup_Un")
)

table(
  cell_data[
    cell_data$Sample == "E14.5_Dnd1-cKO",
  ]$predictedGroup_Un
)

table(
  cell_data[
    cell_data$Sample == "E15.5_Dnd1-cKO",
  ]$predictedGroup_Un
)

table(
  cell_data[
    cell_data$Sample == "E16.5_Dnd1-cKO",
  ]$predictedGroup_Un
)

table(
  cell_data[
    cell_data$Sample == "E17.5_Dnd1-cKO",
  ]$predictedGroup_Un
)

# ---------------------------------------------------------
# Fragment counts
# ---------------------------------------------------------
fragment_counts <- getCellColData(
  proj,
  select = "nFrags"
)

median(
  fragment_counts[
    grepl("E17.5", rownames(fragment_counts)),
  ]
)

median(
  fragment_counts[
    grepl("E16.5", rownames(fragment_counts)),
  ]
)

median(
  fragment_counts[
    grepl("E15.5", rownames(fragment_counts)),
  ]
)

median(
  fragment_counts[
    grepl("E14.5", rownames(fragment_counts)),
  ]
)

# ---------------------------------------------------------
# TSS enrichment
# ---------------------------------------------------------
tss_enrichment <- getCellColData(
  proj,
  select = "TSSEnrichment"
)

median(
  tss_enrichment[
    grepl("E17.5", rownames(tss_enrichment)),
  ]
)

median(
  tss_enrichment[
    grepl("E16.5", rownames(tss_enrichment)),
  ]
)

median(
  tss_enrichment[
    grepl("E15.5", rownames(tss_enrichment)),
  ]
)

median(
  tss_enrichment[
    grepl("E14.5", rownames(tss_enrichment)),
  ]
)

# =========================================================
# Figure S6d
# =========================================================

markerGenes  <- c("Pou5f1", "Dazl", #PGCs/ECCs
                   "Car4", "Nanog", #Sertoli Cells
                   "Ddx4", "Sox2", #Leidig cells
                   "Nodal", "Otx2", #Interstitial/Stromal
                   "Fgf5", "Cer1", #Macrophage
                   "T", "Mixl1") #Endothelial Cells

theme_umap_clean <- function(x) {
  ggplot2::theme_void() +
  ggplot2::theme(
      legend.position = c(0.2, 0.3),
      legend.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.key.size = unit(0.9, "cm"),
      plot.title = element_text(size = 30,hjust = 0.5,face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 30, face = "italic")
    )
}

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Pou5f1", embedding = "UMAP", quantCut = c(0, 0.8),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "PGCs/ECCs", subtitle = "Pou5f1")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Ddx4", embedding = "UMAP", quantCut = c(0, 0.88),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "PGCs", subtitle = "Ddx4")

p3 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Amh", embedding = "UMAP", quantCut = c(0, 0.8),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Sertoli Cells", subtitle = "Amh")

p4 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Ptgds", embedding = "UMAP", quantCut = c(0, 0.88),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Sertoli Cells", subtitle = "Ptgds")

p5 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Cyp17a1", embedding = "UMAP", quantCut = c(0, 0.88),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Leydig Cells", subtitle = "Cyp17a1")

p6 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Star", embedding = "UMAP", quantCut = c(0, 0.88),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Leydig Cells", subtitle = "Star")

p7 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Col1a1", embedding = "UMAP", quantCut = c(0, 0.88),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Interstitial/Stromal Cells", subtitle = "Col1a1")

p8 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Acta2", embedding = "UMAP", quantCut = c(0, 0.88),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Interstitial/Stromal Cells", subtitle = "Acta2")

p9 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Pecam1", embedding = "UMAP", quantCut = c(0, 1),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Endothelial Cells", subtitle = "Pecam1")

p10 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Cdh5", embedding = "UMAP", quantCut = c(0, 0.98),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Endothelial Cells", subtitle = "Cdh5")

p11 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Adgre1", embedding = "UMAP", quantCut = c(0.1, 0.88),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() +
                   ggplot2::labs(title = "Macrophages", subtitle = "Adgre1")

p12 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", 
                   name = "Cd68", embedding = "UMAP", quantCut = c(0.1, 0.95),
                   imputeWeights = NULL) +
                   #scale_fill_gradient(low = "lightgrey", high = "red") +
                   scale_fill_viridis_c(option = "magma") +
                   theme_umap_clean() + 
                   ggplot2::labs(title = "Macrophages", subtitle = "Cd68")

(p1|p3|p5)/(p2|p4|p6)/(p7|p9|p11)/(p8|p10|p12)