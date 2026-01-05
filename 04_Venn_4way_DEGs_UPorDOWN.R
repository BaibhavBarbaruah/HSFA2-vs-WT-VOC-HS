############################################################
## 04_Venn_4way_DEGs_UPorDOWN.R
## 4-way overlaps across within-genotype stress contrasts:
##   Col0_T2H vs Col0_NS
##   Col0_HS  vs Col0_NS
##   hsfa2_T2H vs hsfa2_NS
##   hsfa2_HS  vs hsfa2_NS
## Choose UP or DOWN sets by thresholds.
## If ggVennDiagram is installed, it will draw a Venn; otherwise prints overlap table.
############################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

source("config.R")
source("99_helpers.R")

PADJ_CUT <- 0.05
LFC_CUT  <- 1
MODE <- "UP"   # "UP" or "DOWN"

dds_int <- readRDS(CFG$dds_rds)
dds_g   <- build_dds_g(dds_int)

contrasts <- list(
  Col0_T2H = c("group","Col0_T2H","Col0_NS"),
  Col0_HS  = c("group","Col0_HS","Col0_NS"),
  hsfa2_T2H = c("group","hsfa2_T2H","hsfa2_NS"),
  hsfa2_HS  = c("group","hsfa2_HS","hsfa2_NS")
)

get_set <- function(con) {
  rr <- as.data.frame(results(dds_g, contrast=con, alpha=PADJ_CUT)) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange), padj < PADJ_CUT)
  if (MODE == "UP") {
    rr <- rr %>% dplyr::filter(log2FoldChange >= LFC_CUT)
  } else {
    rr <- rr %>% dplyr::filter(log2FoldChange <= -LFC_CUT)
  }
  clean_tair(rr$gene_id)
}

sets <- lapply(contrasts, get_set)
print(lapply(sets, length))

## overlap counts table
nm <- names(sets)
for (i in 1:(length(nm)-1)) {
  for (j in (i+1):length(nm)) {
    cat(nm[i], "∩", nm[j], "=", length(intersect(sets[[i]], sets[[j]])), "\n")
  }
}

if (requireNamespace("ggVennDiagram", quietly=TRUE)) {
  suppressPackageStartupMessages({ library(ggplot2) })
  p <- ggVennDiagram::ggVennDiagram(sets, label_alpha=0) +
    ggplot2::theme_void()
  print(p)
} else {
  message("Install ggVennDiagram if you want the plotted Venn: install.packages('ggVennDiagram')")
}
