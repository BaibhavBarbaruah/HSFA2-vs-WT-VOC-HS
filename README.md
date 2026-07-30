# HSFA2 and Col-0 responses to T2H and heat stress

This repository contains R scripts used for the downstream RNA-seq analyses associated with the manuscript:

*Trans-2-hexenal priming supports HSFA2-independent thermotolerance in Arabidopsis thaliana*

The analysis compares Col-0 and *hsfa2* seedlings under non-stressed (NS), trans-2-hexenal-primed (T2H), and heat-stressed (HS) conditions.

## Data availability

The raw RNA-seq data are deposited in DDBJ under BioProject **PRJDB39904**.

- Col-0: SAMD01789795–SAMD01789803; DRX872982–DRX872990; DRR895156–DRR895164
- *hsfa2*: SAMD01915725–SAMD01915733; DRX1025185–DRX1025193; DRR1049487–DRR1049495

The scripts begin with a prepared DESeq2 object named `dds_Col0_vs_hsfa2.rds`. Raw-read processing and construction of this object are not included in this repository.

## Repository contents

- `config.R`: input and output paths
- `99_helpers.R`: shared functions for sample annotation, DESeq2 analysis, and gene annotation
- `03_PCA_allSamples_topVarGenes.R`: PCA using the 3,000 most variable genes
- `04_Venn_4way_DEGs_UPorDOWN.R`: four-way overlap analysis of up- or downregulated genes
- `04_Volcano_EnhancedVolcano_selectedContrasts.R`: volcano plots for the main treatment contrasts
- `05_GO_BP_dotplot_UP_merged_contrasts.R`: Gene Ontology enrichment of upregulated genes
- `06_heatmap_*.R`: heatmaps of treatment-responsive and HSFA2-dependent expression patterns

## Requirements

The scripts require R and the following packages:

- DESeq2
- ggplot2
- dplyr
- tibble
- EnhancedVolcano
- ggVennDiagram
- clusterProfiler
- AnnotationDbi
- org.At.tair.db
- ComplexHeatmap
- circlize
- GO.db

`matrixStats` and `GOSemSim` are used when available.

## Running the scripts

1. Clone or download the repository.
2. Create a folder named `data` in the repository root.
3. Place `dds_Col0_vs_hsfa2.rds` in the `data` folder.
4. Start R with the repository root as the working directory.
5. Run the required script with `source()`, for example:

```r
source("03_PCA_allSamples_topVarGenes.R")
```

To keep the input file elsewhere, set the `HSFA2_DATA_DIR` environment variable before running a script:

```r
Sys.setenv(HSFA2_DATA_DIR="path/to/input_directory")
source("03_PCA_allSamples_topVarGenes.R")
```

Figures are displayed in the active R graphics device. Scripts containing optional `ggsave()` lines can be edited to save figures in the `results/figures` directory.
