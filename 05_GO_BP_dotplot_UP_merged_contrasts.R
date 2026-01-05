############################################################
## 05_GO_BP_dotplot_UP_merged_contrasts.R
## GO BP enrichment for UP genes across key contrasts (merged dotplot).
## - Redundancy reduced via simplify() if available, else heuristic
## - Points (not bars); legend points enlarged
## - No auto-saving
############################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.At.tair.db)
  library(clusterProfiler)
})

source("config.R")
source("99_helpers.R")

PADJ_CUT <- 0.05
LFC_CUT  <- 1
TOP_TERMS_PER <- 10
MAX_TERMS_TOTAL <- 25

dds_int <- readRDS(CFG$dds_rds)
dds_g   <- build_dds_g(dds_int)

CONTRASTS <- list(
  `Col-0 T2H vs NS`      = c("group","Col0_T2H","Col0_NS"),
  `Col-0 HS vs NS`       = c("group","Col0_HS","Col0_NS"),
  `HSFA2 T2H vs NS`      = c("group","hsfa2_T2H","hsfa2_NS"),
  `HSFA2 HS vs NS`       = c("group","hsfa2_HS","hsfa2_NS"),
  `HSFA2 vs Col-0 (NS)`  = c("group","hsfa2_NS","Col0_NS"),
  `HSFA2 vs Col-0 (T2H)` = c("group","hsfa2_T2H","Col0_T2H"),
  `HSFA2 vs Col-0 (HS)`  = c("group","hsfa2_HS","Col0_HS")
)

make_term_key <- function(desc) {
  x <- tolower(desc)
  x <- gsub("cellular|process|response|regulation|positive|negative|of|to|in|via|and|the|a|an", " ", x)
  x <- gsub("[^a-z ]", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  x
}

run_enrich_up <- function(res_df, contrast_name) {
  res_df <- res_df %>% dplyr::filter(!is.na(padj), !is.na(log2FoldChange))
  up <- res_df %>% dplyr::filter(padj < PADJ_CUT, log2FoldChange >= LFC_CUT) %>% dplyr::pull(gene_id)
  if (length(up) < 10) return(NULL)
  up <- clean_tair(up)

  eg <- suppressMessages(clusterProfiler::enrichGO(
    gene = up,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = PADJ_CUT,
    qvalueCutoff = PADJ_CUT,
    readable = TRUE
  ))
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)

  if (requireNamespace("GOSemSim", quietly=TRUE)) {
    eg2 <- tryCatch({
      clusterProfiler::simplify(eg, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang")
    }, error=function(e) eg)
    eg <- eg2
  } else {
    tab <- as.data.frame(eg)
    tab$key <- make_term_key(tab$Description)
    tab <- tab %>% dplyr::arrange(p.adjust) %>% dplyr::group_by(key) %>% dplyr::slice(1) %>% dplyr::ungroup()
    eg@result <- tab
  }

  tab <- as.data.frame(eg)
  tab$Contrast <- contrast_name
  tab$GeneRatioNum <- sapply(tab$GeneRatio, function(gr){
    sp <- strsplit(gr, "/")[[1]]
    as.numeric(sp[1]) / as.numeric(sp[2])
  })

  tab %>% dplyr::arrange(p.adjust) %>% head(TOP_TERMS_PER)
}

all_tabs <- list()
for (nm in names(CONTRASTS)) {
  rr <- as.data.frame(results(dds_g, contrast=CONTRASTS[[nm]], alpha=PADJ_CUT)) %>%
    tibble::rownames_to_column("gene_id")
  out <- run_enrich_up(rr, nm)
  if (!is.null(out)) all_tabs[[nm]] <- out
}

df <- dplyr::bind_rows(all_tabs)
if (nrow(df) == 0) stop("No GO results produced. Try lowering LFC_CUT or PADJ_CUT.")

df2 <- df %>%
  dplyr::group_by(Description) %>%
  dplyr::summarize(best_p = min(p.adjust, na.rm=TRUE), .groups="drop") %>%
  dplyr::arrange(best_p) %>%
  head(MAX_TERMS_TOTAL)

keep_terms <- df2$Description
df <- df %>% dplyr::filter(Description %in% keep_terms)

ord <- df %>% dplyr::group_by(Description) %>% dplyr::summarize(best_p=min(p.adjust), .groups="drop") %>% dplyr::arrange(best_p)
df$Description <- factor(df$Description, levels=rev(ord$Description))

base_size <- 16
p <- ggplot(df, aes(x=GeneRatioNum, y=Description, color=Contrast, size=Count)) +
  geom_point(alpha=0.95) +
  theme_classic(base_size=base_size) +
  labs(x="Gene ratio", y=NULL, color=NULL, size="Gene count") +
  theme(
    axis.text.y  = element_text(size=base_size),
    axis.text.x  = element_text(size=base_size),
    legend.text  = element_text(size=base_size),
    legend.title = element_text(size=base_size)
  ) +
  guides(color=guide_legend(override.aes=list(size=7)),
         size =guide_legend(override.aes=list(alpha=1)))

print(p)

## Optional manual save:
# ggsave(file.path(CFG$figures_dir, "GO_BP_UP_merged_dotplot.pdf"), p, width=10, height=8)
# ggsave(file.path(CFG$figures_dir, "GO_BP_UP_merged_dotplot.png"), p, width=10, height=8, dpi=300)
