############################################################
## 03_PCA_allSamples_topVarGenes.R
## PCA on VST using top N most variable genes (replicate-level).
## - No labels inside plot; no ellipse/circle
## - Distinct colors per group; HSFA2 italic in legend
## - No auto-saving
############################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

source("config.R")
source("99_helpers.R")

TOP_VAR <- 3000

dds_int <- readRDS(CFG$dds_rds)
dds_g   <- build_dds_g(dds_int)
vsd <- vst(dds_g, blind=FALSE)
mat <- assay(vsd)

rv <- rowVars_fast(mat)
top <- order(rv, decreasing=TRUE)[seq_len(min(TOP_VAR, length(rv)))]
mat_top <- t(mat[top, , drop=FALSE])

pca <- prcomp(mat_top, center=TRUE, scale.=FALSE)
pct <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  group = colData(vsd)$group,
  stringsAsFactors=FALSE
)
df$group <- factor(as.character(df$group), levels=levels(colData(dds_g)$group))

COLS <- c(
  Col0_NS   = "#4D4D4D",
  Col0_T2H  = "#2CA02C",  # green
  Col0_HS   = "#D62728",  # red
  hsfa2_NS  = "#E377C2",  # pink
  hsfa2_T2H = "#7B2CBF",  # violet
  hsfa2_HS  = "#1F77B4"   # blue
)

breaks <- names(COLS)
lab_parse <- c(
  "Col-0~NS", "Col-0~T2H", "Col-0~HS",
  "italic(HSFA2)~NS", "italic(HSFA2)~T2H", "italic(HSFA2)~HS"
)

base_size <- 18
g <- ggplot(df, aes(PC1, PC2, color=group)) +
  geom_point(size=4.2, alpha=0.95) +
  scale_color_manual(values=COLS, breaks=breaks, labels=parse(text=lab_parse)) +
  labs(
    x = paste0("PC1 (", sprintf("%.1f", pct[1]), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", pct[2]), "%)"),
    color = NULL
  ) +
  theme_classic(base_size=base_size) +
  theme(
    axis.title = element_text(size=base_size+2, face="plain"),
    axis.text  = element_text(size=base_size, face="plain"),
    legend.text = element_text(size=base_size),
    legend.key.height = unit(0.75, "cm"),
    legend.key.width  = unit(0.75, "cm")
  ) +
  guides(color=guide_legend(override.aes=list(size=7)))

print(g)

## Optional manual save:
# ggsave(file.path(CFG$figures_dir, "PCA_topVar3000.pdf"), g, width=7.5, height=6)
# ggsave(file.path(CFG$figures_dir, "PCA_topVar3000.png"), g, width=7.5, height=6, dpi=300)
