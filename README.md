# schard <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/sergiolitwiniuk85/schard_plus/actions/workflows/check.yaml/badge.svg)](https://github.com/sergiolitwiniuk85/schard_plus/actions/workflows/check.yaml)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**schard** is a pure R package to read, write, validate, and interactively QC
[scanpy](https://scanpy.readthedocs.io/) H5AD (AnnData) files.
No Python, no reticulate — just `rhdf5` under the hood.

---

## Quick Start

```r
# Install
remotes::install_github("sergiolitwiniuk85/schard_plus")

# Load
library(schard)

# Validate before loading
h5ad_validate("data.h5ad")

# Interactive QC dashboard
launch_qc_dashboard("data.h5ad")
```

Opens a Shiny app. Adjust sliders, explore UMAP, click **Save Filtered H5AD**
to write the filtered dataset — all without leaving R.

---

## Features

### Read
| Function | Output | Use case |
|---|---|---|
| `h5ad2list()` | `list` | Low-level access to raw data |
| `h5ad2sce()` | `SingleCellExperiment` | Bioconductor workflows |
| `h5ad2seurat()` | `Seurat` | Seurat v4/v5 workflows |
| `h5ad2seurat_spatial()` | `Seurat` or list | Visium spatial data |
| `h5ad2data.frame()` | `data.frame` | Quick metadata access |
| `h5ad2Matrix()` | `dgCMatrix` | Pull expression or embedding matrix |
| `h5ad2images()` | list of images | Visium spatial images |

### Validate
| Function | What it checks |
|---|---|
| `h5ad_validate()` | H5AD structure, required groups, encoding types, index uniqueness — without loading full data |

### Write
| Function | Input | Writes to H5AD |
|---|---|---|
| `write_h5ad()` | `SingleCellExperiment`, `Seurat`, or `list` | `X`, `obs`, `var`, `obsm`, `uns` → scanpy-compatible `.h5ad` |

### Interactive QC Dashboard
| Function | What it does |
|---|---|
| `launch_qc_dashboard()` | Shiny app with real-time QC sliders, UMAP exploration, per-cluster impact, replicate tracking, and **Save Filtered H5AD** button |

### Batch QC Reports
| Function | What it does |
|---|---|
| `qc_report()` | Generates static HTML report + CSV + plots — no browser needed |

### QC Utilities
| Function | Use case |
|---|---|
| `detect_qc_cols()` | Auto-detect QC columns in cell metadata |
| `compute_qc_metrics()` | Compute QC metrics from obs |
| `apply_filter_thresholds()` | Apply threshold filters programmatically |
| `impact_by_replicate()` | Per-replicate cell loss analysis |

---

## Installation

```r
# From GitHub (one package, everything included)
remotes::install_github("sergiolitwiniuk85/schard_plus")

# Load
library(schard)
```

---

## Reading H5AD files

```r
# Load a public dataset
download.file(
  "https://datasets.cellxgene.cziscience.com/8cc521c8-c4ff-4cba-a07b-cae67a9dcba9.h5ad",
  "sn.heart.h5ad"
)

# As SingleCellExperiment
sce <- h5ad2sce("sn.heart.h5ad")

# As Seurat
seu <- h5ad2seurat("sn.heart.h5ad")

# Quick metadata
obs <- h5ad2data.frame("sn.heart.h5ad", "obs")
head(obs)

# Pull a reduction
umap <- t(h5ad2Matrix("sn.heart.h5ad", "/obsm/X_umap"))
plot(umap[, 1:2], pch = 16, cex = 0.4, col = factor(obs$cell_state))
```

### Spatial (Visium)

```r
seu <- h5ad2seurat_spatial("visium_sample.h5ad")
Seurat::SpatialPlot(seu, features = "total_counts")
```

---

## Writing H5AD files

Write R objects back to scanpy-compatible H5AD:

```r
# Round-trip a SingleCellExperiment
sce <- h5ad2sce("input.h5ad")
write_h5ad(sce, "output.h5ad")

# scanpy can now read it:
# >>> import scanpy as ad
# >>> adata = ad.read_h5ad("output.h5ad")
# >>> adata
# AnnData object with n_obs × n_vars = 1000 × 20000
#     obs: 'cell_type', 'batch'
#     var: 'gene_symbol'
#     obsm: 'X_umap', 'X_pca'

# From a Seurat object
seu <- h5ad2seurat("input.h5ad")
write_h5ad(seu, "output.h5ad")

# From a list (manual construction)
write_h5ad(
  list(
    X   = Matrix::rsparsematrix(100, 50, 0.1),
    obs = data.frame(cell_type = sample(c("A", "B"), 100, TRUE)),
    var = data.frame(gene_symbol = paste0("gene", 1:50)),
    obsm = list(X_umap = matrix(runif(200), 100, 2))
  ),
  "output.h5ad"
)
```

### What gets written

| H5AD slot | SCE source | Seurat source |
|---|---|---|
| `/X` | First assay (counts) | `GetAssayData(slot = "counts")` |
| `/obs` | `colData()` | `meta.data` |
| `/var` | `rowData()` | Feature metadata + `rownames` |
| `/obsm` | `reducedDims()` | `Reductions()` |
| `/uns` | `metadata()` | _(not mapped yet)_ |

The output is compatible with `scanpy.read_h5ad()` — same encoding-type conventions, CSR sparse matrix format, and factor→categorical mapping.

---

## Interactive QC Dashboard

Launch a full-featured Shiny app from any supported input:

```r
# From an H5AD file
launch_qc_dashboard("data.h5ad")

# From any loaded object
launch_qc_dashboard(sce)     # SingleCellExperiment
launch_qc_dashboard(seu)     # Seurat
launch_qc_dashboard(data)    # list from h5ad2list()
```

### Dashboard features

| Feature | Description |
|---|---|
| **QC sliders** | Real-time thresholds for Mito%, gene counts, UMI counts, doublet scores |
| **Interactive UMAP** | Plotly with dynamic coloring (pass/fail, cluster, replicate, any metadata column) |
| **Per-cluster impact** | DT table showing cells lost per cluster at current thresholds |
| **Biological replicate tracking** | Per-replicate summary with configurable loss warnings |
| **Save Filtered H5AD** | Writes filtered dataset to a new H5AD file with progress bar |
| **Export HTML report** | Self-contained summary with thresholds, plots, and replicate impact |

---

## CLI / Batch Mode

Run everything from the terminal without interactive R or Shiny:

```bash
# Full QC report — no browser needed
Rscript -e "schard::qc_report('data.h5ad', 'qc_output/')"

# With custom thresholds
Rscript -e "
  schard::qc_report('data.h5ad', 'qc_output/',
    thresholds = list(pct_mito = 10, n_genes = c(500, 6000)))
"
```

From within R:

```r
library(schard)
qc_report("data.h5ad", "qc_output/",
  thresholds = list(pct_mito = 10, n_genes = c(500, 6000)))
```

### What you get

| File | Content |
|---|---|
| `qc_report.html` | Self-contained summary with thresholds, plots, and replicate impact |
| `filtered_cells.csv` | Cell barcodes with pass/fail and fail reason |
| `umap_qc.png` | UMAP colored by pass/fail |
| `qc_{metric}.png` | QC metric distributions (one per metric) |
| `impact_by_replicate.png` | Cell loss per biological replicate |

### Pipeline example

```bash
# Validate → QC report → save filtered H5AD, all from the terminal
Rscript -e "
  library(schard)
  h5ad_validate('data.h5ad')
  qc_report('data.h5ad', 'qc_output/')
"
```

---

## H5AD Validation

Check file integrity before loading large datasets:

```r
h5ad_validate("large_dataset.h5ad")
# → TRUE / FALSE with detailed messages

# Loading from a validated file
if (h5ad_validate("data.h5ad")) {
  sce <- h5ad2sce("data.h5ad")
}
```

Validates: file exists, required groups present, encoding types correct,
index uniqueness, shape consistency — without loading the full expression matrix.

---

## Quality & Metrics

| Metric | Value |
|---|---|
| Cyclomatic complexity (avg) | **4.9** |
| Functions | **63** |
| Test files | **10** (170+ expectations) |
| CI | **GitHub Actions** (ubuntu + macOS) |

Complexity reports are in [`metrics/`](metrics/).

---

## Comparison with alternatives

| Feature | **schard** | sceasy | SeuratDisk |
|---|---|---|---|
| **Read H5AD** | ✅ Pure R, rhdf5 | ✅ Uses reticulate | ✅ Uses h5Seurat intermediate |
| **Write H5AD** | ✅ v1.1.0+ | ✅ | ✅ (via h5Seurat) |
| **QC Dashboard** | ✅ Built-in | ❌ | ❌ |
| **H5AD Validation** | ✅ Built-in | ❌ | ❌ |
| **Spatial (Visium)** | ✅ | Partial | ❌ |
| **Python dependency** | ❌ None | ✅ Requires | ❌ None |
| **Performance** | Fast (direct HDF5) | Moderate (Python bridge) | Moderate (format conversion) |

---

## License

GPL-3. Original author: Pavel Mazin ([cellgeni/schard](https://github.com/cellgeni/schard)).
Dashboard and extended features: Sergio Litwiniuk.
