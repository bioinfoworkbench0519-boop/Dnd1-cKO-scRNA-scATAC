# =========================================================
# Comprehensive scATAC-seq analysis pipeline
# ArchR / Monocle3 workflow for Dnd1-cKO testes
# =========================================================

# =========================================================
# Load libraries
# =========================================================
library(ArchR)
library(monocle3)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ArchRtoSignac)
library(Signac)
library(stringr)

# =========================================================
# Initial settings
# =========================================================
set.seed(1)

addArchRGenome("mm10")

setwd("~")

# =========================================================
# Create Arrow files
# =========================================================

inputFiles <- c(
  "./inputFiles/E14.5_cKO/fragments.tsv.gz",
  "./inputFiles/E15.5_cKO/fragments.tsv.gz",
  "./inputFiles/E16.5_cKO/fragments.tsv.gz",
  "./inputFiles/E17.5_cKO/fragments.tsv.gz"
)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = c(
    "E14.5_Dnd1-cKO",
    "E15.5_Dnd1-cKO",
    "E16.5_Dnd1-cKO",
    "E17.5_Dnd1-cKO"
  ),
  filterTSS = 4,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles <- loadArchRProject(
  path = "inputFiles"
)

# =========================================================
# Doublet detection
# =========================================================

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# =========================================================
# Create ArchR project
# =========================================================

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Merge",
  copyArrows = TRUE
)

# Remove doublets
proj <- filterDoublets(proj)

# Save and reload
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "Merge",
  overwrite = TRUE,
  load = TRUE,
  dropCells = FALSE
)

proj <- loadArchRProject(
  path = "Merge",
  force = FALSE,
  showLogo = TRUE
)

# =========================================================
# Dimensional reduction and clustering
# =========================================================

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  name = "Clusters",
  resolution = 0.9,
  force = TRUE
)

# =========================================================
# Cluster heatmap
# =========================================================

cM_LSI <- confusionMatrix(
  paste0(proj$Clusters),
  paste0(proj$Sample)
)

cM_LSI <- cM_LSI / Matrix::rowSums(cM_LSI)

# =========================================================
# UMAP visualization
# =========================================================

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.15,
  force = TRUE
)

# =========================================================
# Marker accessibility analysis
# =========================================================

markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

proj <- addImputeWeights(proj)

# =========================================================
# Save project
# =========================================================

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "Merge",
  overwrite = TRUE,
  load = TRUE,
  dropCells = FALSE
)

# =========================================================
# Germ cell extraction
# =========================================================

proj_EG <- subset(
  proj,
  proj$Clusters == "C5" |
    proj$Clusters == "C6" |
    proj$Clusters == "C7" |
    proj$Clusters == "C8" |
    proj$Clusters == "C9" |
    proj$Clusters == "C10"
)

saveArchRProject(
  ArchRProj = proj_EG,
  outputDirectory = "Merge_EG2",
  overwrite = TRUE,
  load = TRUE,
  dropCells = FALSE
)

# =========================================================
# Re-clustering of germ cells
# =========================================================

proj_EG <- loadArchRProject(
  path = "Merge_EG2",
  force = FALSE,
  showLogo = TRUE
)

proj_EG <- addIterativeLSI(
  ArchRProj = proj_EG,
  useMatrix = "TileMatrix",
  name = "IterativeLSI_EG1",
  iterations = 3,
  clusterParams = list(
    resolution = c(1.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  force = TRUE
)

proj_EG <- addClusters(
  input = proj_EG,
  reducedDims = "IterativeLSI_EG1",
  name = "Clusters_EG1",
  resolution = 0.9,
  force = TRUE
)

proj_EG <- addClusters(
  input = proj_EG,
  reducedDims = "IterativeLSI_EG1",
  name = "Clusters_EG2",
  resolution = 2.5,
  force = TRUE
)

proj_EG <- addUMAP(
  ArchRProj = proj_EG,
  reducedDims = "IterativeLSI_EG1",
  name = "UMAP_EG1",
  nNeighbors = 30,
  minDist = 0.15,
  force = TRUE,
  seed = 4
)

# =========================================================
# Harmony batch correction
# =========================================================

proj_EG <- addHarmony(
  ArchRProj = proj_EG,
  reducedDims = "IterativeLSI_EG1",
  scaleDims = NULL,
  corCutOff = 0.18,
  name = "Harmony_EG1",
  groupBy = "Sample",
  force = TRUE
)

proj_EG <- addClusters(
  input = proj_EG,
  reducedDims = "Harmony_EG1",
  method = "Seurat",
  name = "Clusters_Harmony_EG1",
  resolution = 1.5,
  force = TRUE
)

proj_EG <- addUMAP(
  ArchRProj = proj_EG,
  reducedDims = "Harmony_EG1",
  name = "UMAP_Harmony_EG1",
  nNeighbors = 30,
  minDist = 0.15,
  force = TRUE
)


# =========================================================
# RNA integration
# =========================================================

seRNA <- readRDS(
  "cKO_germ_merge_remsoma.rds"
)

seRNA$ident <- Idents(seRNA)

proj_EG <- addGeneIntegrationMatrix(
  ArchRProj = proj_EG,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_EG1",
  reducedDims = "IterativeLSI_EG1",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "ident",
  dimsToUse = 1:20,
  groupATAC = "Clusters_EG2",
  nameCell = "predictedCell_Un_EG",
  nameGroup = "predictedGroup_Un_EG",
  nameScore = "predictedScore_Un_EG"
)

# =========================================================
# Peak calling
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

# =========================================================
# Motif enrichment analysis
# =========================================================

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