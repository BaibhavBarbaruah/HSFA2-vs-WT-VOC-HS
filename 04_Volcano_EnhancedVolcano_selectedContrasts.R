############################################################
## 04_Volcano_EnhancedVolcano_selectedContrasts.R
## Volcano plots for selected contrasts using EnhancedVolcano.
## (Requires EnhancedVolcano installed.)
## No auto-saving.
############################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

if (!requireNamespace("EnhancedVolcano", quietly=TRUE)) {
  stop("EnhancedVolcano not installed. Install with: BiocManager::install('EnhancedVolcano')")
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(EnhancedVolcano)
})

source("config.R")
source("99_helpers.R")

PADJ_CUT <- 0.05
LFC_CUT  <- 1
TOP_LAB  <- 10

dds_int <- readRDS(CFG$dds_rds)
dds_g   <- build_dds_g(dds_int)

CONTRASTS <- list(
  `Col-0 HS vs NS`  = c("group","Col0_HS","Col0_NS"),
  `Col-0 T2H vs NS` = c("group","Col0_T2H","Col0_NS"),
  `HSFA2 HS vs NS`  = c("group","hsfa2_HS","hsfa2_NS"),
  `HSFA2 T2H vs NS` = c("group","hsfa2_T2H","hsfa2_NS")
)

for (nm in names(CONTRASTS)) {
  res <- as.data.frame(results(dds_g, contrast=CONTRASTS[[nm]], alpha=PADJ_CUT)) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange))

  ## label top genes by padj
  lab <- res %>% dplyr::arrange(padj) %>% head(TOP_LAB) %>% dplyr::pull(gene_id)

  p <- EnhancedVolcano::EnhancedVolcano(
    res,
    lab = res$gene_id,
    selectLab = lab,
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = PADJ_CUT,
    FCcutoff = LFC_CUT,
    title = nm,
    subtitle = paste0('DESeq2; BH padj<', PADJ_CUT, ', |LFC|>=', LFC_CUT),
    caption = '',
    labSize = 5.0,
    pointSize = 2.5,
    legendLabSize = 12,
    legendIconSize = 4.5
  )
  print(p)
}
