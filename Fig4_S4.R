# =========================================================
# Figure 4 and S4 analysis
# =========================================================

# ---------------------------
# Load libraries
# ---------------------------
library(Seurat)
library(dplyr)
library(reticulate)
library(scales)
library(ggplot2)
library(rtracklayer)
library(ggpubr)

set.seed(1234)

# =========================================================
# Load data
# =========================================================

cKO_germ <- readRDS("cKO_germ_merge_remsoma.rds")

DefaultAssay(cKO_germ) <- "RNA"

# Split by developmental stage
E14.5 <- subset(cKO_germ, subset = orig.ident == "E14.5_cKO")
E15.5 <- subset(cKO_germ, subset = orig.ident == "E15.5_cKO")
E16.5 <- subset(cKO_germ, subset = orig.ident == "E16.5_cKO")
E17.5 <- subset(cKO_germ, subset = orig.ident == "E17.5_cKO")

# =========================================================
# Table S1
# Marker gene analysis
# =========================================================

markers <- FindAllMarkers(
  E16.5,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top2000 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 2000, wt = avg_log2FC)

cl7_top2000 <- top2000 %>%
  filter(cluster == "7")

# =========================================================
# Figure 4 a-d
# DotPlot of Car4-related genes
# =========================================================

plot_car4_dotplot <- function(
  seu_object,
  stage_name,
  output_file
) {

  tiff(
    output_file,
    width  = 1200,
    height = 1650,
    res    = 400
  )

  p <- DotPlot(
    seu_object,
    features = c("Car4", "Nanog", "Otx2"),
    cols = c("blue", "red"),
    dot.scale = 8,
    scale.max = 100,
    col.min = -0.5,
    col.max = 1.8
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
      legend.text = element_text(size = 15),
      axis.title  = element_blank(),
      axis.text   = element_text(size = 17)
    ) +
    labs(title = stage_name)

  print(p)
  dev.off()
}

# Generate plots
plot_car4_dotplot(E14.5, "E14.5_cKO", "E14.5_Car4_Dot.tif")
plot_car4_dotplot(E15.5, "E15.5_cKO", "E15.5_Car4_Dot.tif")
plot_car4_dotplot(E16.5, "E16.5_cKO", "E16.5_Car4_Dot.tif")
plot_car4_dotplot(E17.5, "E17.5_cKO", "E17.5_Car4_Dot.tif")

# =========================================================
# Figure 4v
# Statistical analysis and boxplot
# =========================================================

# ---------------------------
# Two-way ANOVA
# ---------------------------

Car4q <- read.csv("Car4_quant.csv", header = TRUE)

anova_result <- aov(
  formation ~ cell * condition,
  data = Car4q
)

summary_res <- summary(anova_result)

write.csv(
  as.data.frame(summary_res[[1]]),
  "two-way-annovar.csv"
)

# ---------------------------
# Tukey-Kramer post hoc test
# ---------------------------

tukey_res <- TukeyHSD(anova_result)

write.csv(
  as.data.frame(tukey_res$`cell:condition`),
  "Tukey-Krammer.csv"
)

# ---------------------------
# Boxplot visualization
# ---------------------------

Car4q <- read.csv(
  "Car4_quant.csv",
  header = TRUE
)

Car4q2$cell.condition <- factor(
  Car4q2$cell.condition,
  levels = c(
    "AFX_CAR4-",
    "AFX_CAR4+",
    "2i/LIF_CAR4-",
    "2i/LIF_CAR4+"
  )
)


p <- ggplot(
  data = Car4q2,
  mapping = aes(
    x = cell.condition,
    y = formation,
    fill = cell.condition
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    coef = 2.5
  ) +
  geom_jitter(
    shape = 19,
    size = 1,
    width = 0.1
  ) +
  scale_fill_manual(
    values = c(
      "AFX_CAR4-"     = "orange",
      "AFX_CAR4+"     = "#ec4242",
      "2i/LIF_CAR4-"  = "skyblue",
      "2i/LIF_CAR4+"  = "#3939e4"
    )
  ) +
  scale_y_continuous(
    limits = c(-0.1, 17.5)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  theme(
    axis.text.x = element_text(
      size  = 12,
      angle = 45,
      vjust = 0.9,
      hjust = 0.9
    ),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  geom_signif(
    comparisons = list(
      c("AFX_CAR4-", "AFX_CAR4+")
    ),
    map_signif_level = TRUE,
    annotations = "***",
    textsize = 7,
    y = 14.6
  ) +
  geom_signif(
    comparisons = list(
      c("AFX_CAR4+", "2i/LIF_CAR4+")
    ),
    map_signif_level = TRUE,
    annotations = "***",
    textsize = 7,
    y = 15.5
  )

print(p)