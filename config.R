## 01_scripts/config.R
## Edit paths here.

CFG <- list(
  dds_rds = "F:/RNA seq/hsfa2_vs_Col0_combined_20251224/dds_Col0_vs_hsfa2.rds",
  out_dir     = "../02_results",
  figures_dir = "../02_results/figures",
  tables_dir  = "../02_results/tables"
)

message("Using dds_rds: ", CFG$dds_rds)
if (!file.exists(CFG$dds_rds)) stop("dds RDS not found. Edit config.R and point to your local file.")
