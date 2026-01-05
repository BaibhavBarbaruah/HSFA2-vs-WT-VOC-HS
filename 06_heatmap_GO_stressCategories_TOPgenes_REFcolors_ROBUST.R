############################################################
## 06_heatmap_GO_stressCategories_TOPgenes_REFcolors_ROBUST.R
## - Genes UP in any within-genotype stress contrast
## - Assign to stress categories via GO BP offspring sets
## - Mean VST (6 groups), ref color scheme, no row names
## - Robust GO keytype selection (no SYMBOL crash)
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
  library(AnnotationDbi)
  library(org.At.tair.db)
  library(GO.db)
})

ht_opt$message <- FALSE

source("config.R")
source("99_helpers.R")

PADJ_CUT      <- 0.05
MIN_LFC_UP    <- 0
N_PER_CAT     <- 35
MAX_TOTAL_GENES <- 1200

BASE_SIZE     <- 16
COLNAME_SIZE  <- 16
SPLIT_TITLE   <- 14

REF_COLS <- c("#3e6aa0", "#979293", "#cfb479", "#f7c853", "#df7a3d")

safe_mapIds <- function(db, keys, keytype, column, multiVals="list") {
  tryCatch({
    AnnotationDbi::mapIds(db, keys=keys, keytype=keytype, column=column, multiVals=multiVals)
  }, error=function(e) {
    setNames(vector("list", length(keys)), keys)
  })
}

pick_best_keytype <- function(keys_vec, candidates) {
  keys_vec <- unique(keys_vec)
  best <- candidates[1]
  best_n <- -1
  for (kt in candidates) {
    x <- safe_mapIds(org.At.tair.db, keys=keys_vec, keytype=kt, column="GO", multiVals="list")
    n_ok <- sum(!vapply(x, function(z) is.null(z) || length(z)==0 || all(is.na(z)), logical(1)))
    if (n_ok > best_n) {
      best_n <- n_ok
      best <- kt
    }
  }
  list(keytype=best, n_mapped=best_n)
}

dds_int <- readRDS(CFG$dds_rds)
dds_g   <- build_dds_g(dds_int)
meanVST <- get_meanVST_6groups(dds_g)

stress_contrasts <- list(
  Col0_T2H_vs_Col0_NS   = c("group","Col0_T2H","Col0_NS"),
  Col0_HS_vs_Col0_NS    = c("group","Col0_HS","Col0_NS"),
  HSFA2_T2H_vs_HSFA2_NS = c("group","hsfa2_T2H","hsfa2_NS"),
  HSFA2_HS_vs_HSFA2_NS  = c("group","hsfa2_HS","hsfa2_NS")
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

tabs <- lapply(names(stress_contrasts), function(nm) get_lfc_padj(stress_contrasts[[nm]], nm))
df <- Reduce(function(a,b) dplyr::full_join(a,b, by="gene_id"), tabs)
df <- df %>% dplyr::filter(gene_id %in% rownames(meanVST))

lfc_cols  <- grep("^lfc_",  colnames(df), value=TRUE)
padj_cols <- grep("^padj_", colnames(df), value=TRUE)

best_padj_up <- rep(NA_real_, nrow(df))
best_lfc_up  <- rep(NA_real_, nrow(df))
for (i in seq_len(nrow(df))) {
  lfcs  <- as.numeric(df[i, lfc_cols,  drop=TRUE])
  padjs <- as.numeric(df[i, padj_cols, drop=TRUE])
  ok <- (!is.na(lfcs) & lfcs > MIN_LFC_UP & !is.na(padjs))
  if (any(ok)) {
    best_padj_up[i] <- min(padjs[ok], na.rm=TRUE)
    best_lfc_up[i]  <- max(lfcs[ok],  na.rm=TRUE)
  }
}
df$best_padj_up <- best_padj_up
df$best_lfc_up  <- best_lfc_up

up_genes <- df %>%
  dplyr::filter(!is.na(best_padj_up) & best_padj_up < PADJ_CUT & !is.na(best_lfc_up) & best_lfc_up > MIN_LFC_UP) %>%
  dplyr::arrange(best_padj_up) %>%
  dplyr::pull(gene_id)

cat("Upregulated DE genes (any stress contrast):", length(up_genes), "\n")
if (length(up_genes) == 0) stop("No UP genes found.")

up_genes_clean <- clean_tair(up_genes)
looks_tair <- mean(grepl("^AT[1-5MC]G\\d{5}$", up_genes_clean)) > 0.6

available_kt <- keytypes(org.At.tair.db)
preferred <- if (looks_tair) {
  c("TAIR", "TAIRLocus", "TAIRLOCUS", "GENENAME", "SYMBOL")
} else {
  c("SYMBOL", "GENENAME", "TAIR", "TAIRLocus", "TAIRLOCUS")
}
try_kt <- unique(intersect(preferred, available_kt))
if (length(try_kt) == 0) try_kt <- available_kt[1]

kt_pick <- pick_best_keytype(up_genes_clean, try_kt)
keytype_use <- kt_pick$keytype
cat("Using keytype for GO mapping:", keytype_use, "| genes with GO mapped:", kt_pick$n_mapped, "\n")
if (kt_pick$n_mapped <= 0) stop("GO mapping failed for all candidate keytypes.")

go_list_clean <- safe_mapIds(org.At.tair.db, keys=unique(up_genes_clean), keytype=keytype_use, column="GO", multiVals="list")
all_go <- unique(unlist(go_list_clean, use.names=FALSE))
all_go <- all_go[!is.na(all_go)]
if (length(all_go) == 0) stop("GO mapping returned 0 GO IDs.")

go_anno <- AnnotationDbi::select(GO.db, keys=all_go, columns=c("TERM","ONTOLOGY"), keytype="GOID") %>%
  dplyr::distinct(GOID, .keep_all=TRUE)

bp_goids <- go_anno %>% dplyr::filter(ONTOLOGY=="BP") %>% dplyr::pull(GOID) %>% unique()
bp_goids <- bp_goids[!is.na(bp_goids)]
if (length(bp_goids) == 0) stop("No BP GOIDs found.")

gene_bp_goids <- lapply(seq_along(up_genes), function(i){
  g_clean <- up_genes_clean[i]
  gids <- go_list_clean[[g_clean]]
  if (is.null(gids) || length(gids)==0 || all(is.na(gids))) return(character(0))
  gids <- as.character(gids)
  gids <- gids[!is.na(gids)]
  intersect(unique(gids), bp_goids)
})
names(gene_bp_goids) <- up_genes

category_roots <- list(
  "Heat / Temperature" = c("GO:0009408"),
  "Oxidative / ROS"    = c("GO:0006979"),
  "Hypoxia / Low O2"   = c("GO:0001666"),
  "Defense / Immune"   = c("GO:0006952"),
  "Detox / Xenobiotic" = c("GO:0009636"),
  "ER / Protein Fold"  = c("GO:0006457", "GO:0034976"),
  "Wound"              = c("GO:0009611"),
  "Osmotic / Salt"     = c("GO:0006970", "GO:0009651"),
  "Cold"               = c("GO:0009409"),
  "UV / Light Stress"  = c("GO:0009411", "GO:0009644")
)

get_bp_offspring_set <- function(root_goids){
  out <- character(0)
  for (r in root_goids) {
    off <- AnnotationDbi::mget(r, GOBPOFFSPRING, ifnotfound=NA)[[1]]
    if (length(off)==1 && is.na(off)) off <- character(0)
    out <- c(out, r, off)
  }
  unique(out)
}
category_sets <- lapply(category_roots, get_bp_offspring_set)

assigned_cat <- rep(NA_character_, length(up_genes))
names(assigned_cat) <- up_genes
for (cat_name in names(category_sets)) {
  go_set <- category_sets[[cat_name]]
  hit <- vapply(up_genes, function(g){ length(intersect(gene_bp_goids[[g]], go_set)) > 0 }, logical(1))
  pick <- names(hit)[hit & is.na(assigned_cat[names(hit)])]
  if (length(pick) > 0) assigned_cat[pick] <- cat_name
}

cat_genes <- names(assigned_cat)[!is.na(assigned_cat)]
cat("Genes assigned to stress categories:", length(cat_genes), "\n")
if (length(cat_genes) == 0) stop("0 assignments.")

score_tbl <- df %>%
  dplyr::select(gene_id, best_padj_up, best_lfc_up) %>%
  dplyr::filter(gene_id %in% cat_genes) %>%
  dplyr::mutate(category = assigned_cat[gene_id]) %>%
  dplyr::arrange(best_padj_up, dplyr::desc(best_lfc_up))

picked <- c(); picked_cat <- c()
for (cat_name in names(category_roots)) {
  g <- score_tbl %>% dplyr::filter(category == cat_name)
  if (nrow(g) == 0) next
  take <- head(g$gene_id, N_PER_CAT)
  picked <- c(picked, take)
  picked_cat <- c(picked_cat, rep(cat_name, length(take)))
}
if (length(picked) > MAX_TOTAL_GENES) {
  picked <- picked[1:MAX_TOTAL_GENES]
  picked_cat <- picked_cat[1:MAX_TOTAL_GENES]
}
dup <- duplicated(picked)
picked <- picked[!dup]
picked_cat <- picked_cat[!dup]
picked_cat <- factor(picked_cat, levels=names(category_roots))

mat_sel <- meanVST[picked, , drop=FALSE]

vals <- as.numeric(mat_sel)
brks <- quantile(vals, probs=c(0.02, 0.35, 0.60, 0.80, 0.98), na.rm=TRUE)
col_fun <- colorRamp2(as.numeric(brks), REF_COLS)

ht <- Heatmap(
  mat_sel,
  name = "mean VST",
  col  = col_fun,
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  show_row_names  = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize=COLNAME_SIZE, fontface="plain", col="black"),
  row_split = picked_cat,
  row_title_gp = gpar(fontsize=SPLIT_TITLE, fontface="plain", col="black"),
  row_title_rot = 0,
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize=BASE_SIZE, fontface="plain"),
    labels_gp = gpar(fontsize=BASE_SIZE-2, fontface="plain")
  )
)

draw(ht)
