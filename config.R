## config.R
## Repository-relative paths for the downstream RNA-seq scripts.

data_dir <- Sys.getenv("HSFA2_DATA_DIR", unset="data")

CFG <- list(
  dds_rds     = file.path(data_dir, "dds_Col0_vs_hsfa2.rds"),
  out_dir     = "results",
  figures_dir = file.path("results", "figures"),
  tables_dir  = file.path("results", "tables")
)

dir.create(CFG$figures_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(CFG$tables_dir, recursive=TRUE, showWarnings=FALSE)

message("Using dds_rds: ", CFG$dds_rds)
if (!file.exists(CFG$dds_rds)) {
  stop(
    "DESeq2 input object not found at ", CFG$dds_rds, ". ",
    "Place dds_Col0_vs_hsfa2.rds in the data directory or set HSFA2_DATA_DIR."
  )
}
