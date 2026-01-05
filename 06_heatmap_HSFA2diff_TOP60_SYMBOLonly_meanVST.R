############################################################
## 06_heatmap_HSFA2diff_TOP60_SYMBOLonly_meanVST.R
## Mean VST heatmap (6 groups) for genes differing in HSFA2 vs Col-0
## within each treatment (NS, T2H, HS).
## - Select TOP 60 by best padj then best |LFC| across those 3 contrasts
## - Row names: SYMBOL-only (drops genes without SYMBOL)
## - No auto-saving
############################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

ht_opt$message <- FALSE

source("config.R")
source("99_helpers.R")

TOP_N <- 60
PADJ_CUT <- 0.05
ABS_LFC_CUT <- 1

BASE_SIZE    <- 16
COLNAME_SIZE <- 16
ROWNAME_SIZE <- 12
REF_COLS <- c("#3e6aa0", "#979293", "#cfb479", "#f7c853", "#df7a3d")

dds_int <- readRDS(CFG$dds_rds)
dds_g   <- build_dds_g(dds_int)
meanVST <- get_meanVST_6groups(dds_g)

contrasts <- list(
  HSFA2_vs_Col0_NS  = c("group","hsfa2_NS","Col0_NS"),
  HSFA2_vs_Col0_T2H = c("group","hsfa2_T2H","Col0_T2H"),
  HSFA2_vs_Col0_HS  = c("group","hsfa2_HS","Col0_HS")
)

get_lfc_padj <- function(con, nm){
  as.data.frame(results(dds_g, contrast=con, alpha=PADJ_CUT)) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::transmute(
      gene_id,
      !!paste0("lfc_", nm)  := log2FoldChange,
      !!paste0("padj_", nm) := padj
    )
}

tabs <- lapply(names(contrasts), function(nm) get_lfc_padj(contrasts[[nm]], nm))
df <- Reduce(function(a,b) dplyr::full_join(a,b, by="gene_id"), tabs)
df <- df %>% dplyr::filter(gene_id %in% rownames(meanVST))

lfc_cols  <- grep("^lfc_",  colnames(df), value=TRUE)
padj_cols <- grep("^padj_", colnames(df), value=TRUE)

best_padj <- rep(NA_real_, nrow(df))
best_abs_lfc <- rep(NA_real_, nrow(df))
for (i in seq_len(nrow(df))) {
  lfcs  <- as.numeric(df[i, lfc_cols,  drop=TRUE])
  padjs <- as.numeric(df[i, padj_cols, drop=TRUE])
  ok <- (!is.na(lfcs) & !is.na(padjs))
  if (any(ok)) {
    best_padj[i]    <- min(padjs[ok], na.rm=TRUE)
    best_abs_lfc[i] <- max(abs(lfcs[ok]), na.rm=TRUE)
  }
}
df$best_padj <- best_padj
df$best_abs_lfc <- best_abs_lfc

rank_tbl <- df %>% dplyr::filter(!is.na(best_padj) & !is.na(best_abs_lfc)) %>%
  dplyr::arrange(best_padj, dplyr::desc(best_abs_lfc))

cand <- rank_tbl %>% dplyr::filter(best_padj < PADJ_CUT, best_abs_lfc >= ABS_LFC_CUT)
if (nrow(cand) < TOP_N) cand <- rank_tbl %>% dplyr::filter(best_padj < PADJ_CUT)
if (nrow(cand) < TOP_N) cand <- rank_tbl

## Build SYMBOL-only list (skip genes w/o SYMBOL; avoid duplicate SYMBOL labels)
want <- cand$gene_id
sym_map <- map_TAIR_to_SYMBOL_only(want)
picked_genes <- character(0)
picked_labels <- character(0)

for (g in want) {
  sym <- sym_map[[g]]
  if (is.na(sym) || sym == "") next
  if (sym %in% picked_labels) next
  picked_genes <- c(picked_genes, g)
  picked_labels <- c(picked_labels, sym)
  if (length(picked_genes) >= TOP_N) break
}

if (length(picked_genes) == 0) stop("No SYMBOL-mapped genes found. If you want TAIR IDs as fallback, tell me.")
if (length(picked_genes) < TOP_N) message("Only ", length(picked_genes), " genes had unique SYMBOL names (requested ", TOP_N, ").")

mat_sel <- meanVST[picked_genes, , drop=FALSE]
rownames(mat_sel) <- picked_labels

vals <- as.numeric(mat_sel)
brks <- quantile(vals, probs=c(0.02, 0.35, 0.60, 0.80, 0.98), na.rm=TRUE)
col_fun <- colorRamp2(as.numeric(brks), REF_COLS)

ht <- Heatmap(
  mat_sel,
  name = "mean VST",
  col  = col_fun,
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  show_row_names  = TRUE,
  row_names_side  = "right",
  row_names_gp    = gpar(fontsize=ROWNAME_SIZE, fontface="plain", col="black"),
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp  = gpar(fontsize=COLNAME_SIZE, fontface="plain", col="black"),
  column_title = paste0("TOP ", nrow(mat_sel), " genes differing in HSFA2 vs Col-0 within treatments (mean VST)"),
  column_title_gp = gpar(fontsize=BASE_SIZE+4, fontface="plain", col="black")
)

draw(ht)
