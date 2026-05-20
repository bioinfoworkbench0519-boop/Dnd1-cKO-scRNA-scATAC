# =========================================================
# MouseGastrulationData preprocessing
# =========================================================

# ---------------------------
# Load libraries
# ---------------------------
library(Seurat)
library(MouseGastrulationData)
library(SingleCellExperiment)

# ---------------------------
# Set working directory
# ---------------------------
setwd("~")

# =========================================================
# Load embryo datasets
# =========================================================

embryo_samples <- list(
  
  E7.0_1  = 14,
  E7.0_2  = 15,
  E7.25_1 = 20,
  E7.5_1  = 6,
)

embryo_data <- lapply(embryo_samples, EmbryoAtlasData)

# =========================================================
# Function for Seurat conversion
# =========================================================

process_embryo_sample <- function(
  sce_object,
  sample_name,
  output_name,
  celltypes_keep = c(
    "Primitive Streak",
    "Anterior Primitive Streak",
    "Epiblast"
  )
) {

  # Set gene symbols as rownames
  rownames(sce_object) <- rowData(sce_object)$SYMBOL

  # Remove doublets and stripped cells
  singlets <- sce_object[
    ,
    !(colData(sce_object)$doublet |
      colData(sce_object)$stripped)
  ]

  # Convert to Seurat
  seu <- as.Seurat(singlets, data = NULL)

  # Subset selected celltypes
  seu <- subset(
    seu,
    celltype %in% celltypes_keep
  )

  # Add metadata
  seu$orig.ident <- sample_name

  # Replace RNA assay with original counts
  seu@assays$RNA <- seu@assays$originalexp

  # Save object
  saveRDS(seu, paste0(output_name, ".rds"))

  return(seu)
}

# =========================================================
# Process selected datasets
# =========================================================

# ---------------------------
# E7.0 sample 1
# ---------------------------
E7.0_1_seu <- process_embryo_sample(
  sce_object = embryo_data$E7.0_3,
  sample_name = "E7.0",
  output_name = "E7.0_3"
)

# ---------------------------
# E7.0 sample 2
# ---------------------------
E7.0_2_seu <- process_embryo_sample(
  sce_object = embryo_data$E7.0_2,
  sample_name = "E7.0",
  output_name = "E7.0_2"
)

# ---------------------------
# E7.25 sample 2
# ---------------------------
E7.25_1_seu <- process_embryo_sample(
  sce_object = embryo_data$E7.25_2,
  sample_name = "E7.25",
  output_name = "E7.25_2"
)

# ---------------------------
# E7.5 sample 4
# ---------------------------
E7.5_1_seu <- process_embryo_sample(
  sce_object = embryo_data$E7.5_4,
  sample_name = "E7.5",
  output_name = "E7.5_4"
)