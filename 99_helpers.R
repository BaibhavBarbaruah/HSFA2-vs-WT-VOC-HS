## 01_scripts/99_helpers.R
## Shared helpers (no auto-saving).

normalize_cd <- function(cd) {
  cd <- as.data.frame(cd)
  cd$genotype <- as.character(cd$genotype)
  cd$genotype[cd$genotype %in% c("Col-0","Col0","COL0","col0","Col")] <- "Col0"
  cd$genotype[cd$genotype %in% c("hsfa2","HSFA2","Hsfa2")]           <- "hsfa2"

  cd$treatment <- as.character(cd$treatment)
  cd$treatment[cd$treatment %in% c("Control","control","NS","ns")] <- "NS"
  cd$treatment[cd$treatment %in% c("T2H","t2h")]                     <- "T2H"
  cd$treatment[cd$treatment %in% c("HS","hs","Heat","heat")]     <- "HS"

  cd$genotype  <- factor(cd$genotype,  levels=c("Col0","hsfa2"))
  cd$treatment <- factor(cd$treatment, levels=c("NS","T2H","HS"))

  if (any(is.na(cd$genotype)) || any(is.na(cd$treatment))) {
    stop("Found NA in genotype/treatment after normalization. Check colData.")
  }

  cd$group <- factor(
    paste0(cd$genotype, "_", cd$treatment),
    levels=c("Col0_NS","Col0_T2H","Col0_HS","hsfa2_NS","hsfa2_T2H","hsfa2_HS")
  )
  cd
}

pretty_group <- function(x) {
  x <- gsub("^Col0_",  "Col-0 ", x)
  x <- gsub("^hsfa2_", "HSFA2 ", x)
  gsub("_", " ", x)
}

clean_tair <- function(x) toupper(gsub("\\.\\d+$", "", x))

rowVars_fast <- function(m) {
  if (requireNamespace("matrixStats", quietly=TRUE)) {
    matrixStats::rowVars(m)
  } else {
    apply(m, 1, var)
  }
}

build_dds_g <- function(dds_int) {
  suppressPackageStartupMessages({
    library(DESeq2)
  })
  cd <- normalize_cd(colData(dds_int))
  dds_g <- DESeqDataSetFromMatrix(countData=counts(dds_int), colData=cd, design=~group)
  keep <- rowSums(counts(dds_g) >= 10) >= 3
  dds_g <- dds_g[keep, ]
  dds_g <- DESeq(dds_g, quiet=TRUE)
  dds_g
}

get_meanVST_6groups <- function(dds_g) {
  suppressPackageStartupMessages({
    library(DESeq2)
  })
  vsd <- vst(dds_g, blind=FALSE)
  mat <- assay(vsd)
  grp <- colData(dds_g)$group
  meanVST <- sapply(levels(grp), function(g){ rowMeans(mat[, grp==g, drop=FALSE]) })
  colnames(meanVST) <- vapply(colnames(meanVST), pretty_group, character(1))
  meanVST[, c("Col-0 NS","Col-0 T2H","Col-0 HS","HSFA2 NS","HSFA2 T2H","HSFA2 HS"), drop=FALSE]
}

map_TAIR_to_SYMBOL_only <- function(tair_ids) {
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.At.tair.db)
  })
  tair_ids2 <- clean_tair(tair_ids)
  sym <- tryCatch({
    AnnotationDbi::mapIds(org.At.tair.db, keys=tair_ids2, keytype="TAIR", column="SYMBOL", multiVals="first")
  }, error=function(e) {
    setNames(rep(NA_character_, length(tair_ids2)), tair_ids2)
  })
  sym_vec <- as.character(sym[tair_ids2])
  names(sym_vec) <- tair_ids
  sym_vec
}
