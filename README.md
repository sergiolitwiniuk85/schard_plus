# schard <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/sergiolitwiniuk85/schard_plus/actions/workflows/check.yaml/badge.svg)](https://github.com/sergiolitwiniuk85/schard_plus/actions/workflows/check.yaml)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**schard** is a pure R package to read and write [scanpy](https://scanpy.readthedocs.io/) H5AD (AnnData) files.
No Python, no reticulate — just `rhdf5` under the hood.

---

## Features

### Read (original)
| Function | Output | Use case |
|---|---|---|
| `h5ad2list()` | `list` | Low-level access to raw data |
| `h5ad2sce()` | `SingleCellExperiment` | Bioconductor workflows |
| `h5ad2seurat()` | `Seurat` | Seurat v4/v5 workflows |
| `h5ad2seurat_spatial()` | `Seurat` or list | Visium spatial data |
| `h5ad2data.frame()` | `data.frame` | Quick metadata access |
| `h5ad2Matrix()` | `dgCMatrix` | Pull expression or embedding matrix |
| `h5ad2images()` | list of images | Visium spatial images |

### Write (new in v1.1.0)
| Function | Input | Writes to H5AD |
|---|---|---|
| `write_h5ad()` | `SingleCellExperiment`, `Seurat`, or `list` | `X`, `obs`, `var`, `obsm`, `uns` → scanpy-compatible `.h5ad` |

---

## Installation

```r
# From GitHub
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
sce <- schard::h5ad2sce("sn.heart.h5ad")

# As Seurat
seu <- schard::h5ad2seurat("sn.heart.h5ad")

# Quick metadata
obs <- schard::h5ad2data.frame("sn.heart.h5ad", "obs")
head(obs)

# Pull a reduction
umap <- t(schard::h5ad2Matrix("sn.heart.h5ad", "/obsm/X_umap"))
plot(umap[, 1:2], pch = 16, cex = 0.4, col = factor(obs$cell_state))
```

### Spatial (Visium)

```r
seu <- schard::h5ad2seurat_spatial("visium_sample.h5ad")
Seurat::SpatialPlot(seu, features = "total_counts")
```

---

## Writing H5AD files

Write R objects back to scanpy-compatible H5AD:

```r
# Round-trip a SingleCellExperiment
sce <- schard::h5ad2sce("input.h5ad")
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
seu <- schard::h5ad2seurat("input.h5ad")
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

## Quality & Metrics

| Metric | Before | After |
|---|---|---|
| Cyclomatic complexity (avg) | 9.8 | **6.6** 🔽 |
| Functions | 9 | **21** |
| Tests | 0 | **13** (49 expectations) |
| CI | None | **GitHub Actions** (ubuntu + macOS) |

Complexity baseline and post-refactor reports are in [`metrics/`](metrics/).

---

---

## schardExplorer — Interactive QC Dashboard

[schardExplorer](schardExplorer/) is a companion R package that provides interactive and batch QC for single-cell data.

### Visual mode

```r
library(schardExplorer)

# Launch interactive dashboard from any loaded object
launch_qc_dashboard("data.h5ad")
launch_qc_dashboard(sce)     # SingleCellExperiment
launch_qc_dashboard(seu)     # Seurat
launch_qc_dashboard(data)    # list from h5ad2list()
```

Opens a Shiny app with:
- **QC sliders** — real-time thresholds for Mito%, gene counts, UMI counts, doublet scores
- **Interactive UMAP** — plotly with dynamic coloring (pass/fail, cluster, replicate, any metadata)
- **Per-cluster impact** — DT table showing cells lost per cluster at current thresholds
- **Biological replicate tracking** — per-replicate summary with configurable loss warnings
- **CSV + HTML export** — filtered cells, thresholds, and per-replicate stats

### Batch mode

```r
# Generate static report — no browser needed, runs anywhere
qc_report("data.h5ad", "qc_output/",
  thresholds = list(pct_mito = 10, n_genes = c(500, 6000)))
```

Produces: `filtered_cells.csv`, UMAP/QC distribution plots, impact barplots, and self-contained HTML summary.

### Install

```r
remotes::install_github("sergiolitwiniuk85/schard_plus", subdir = "schardExplorer")
```

See the [schardExplorer README](schardExplorer/) for full documentation.

---

## Comparison with alternatives

| Feature | **schard** | sceasy | SeuratDisk |
|---|---|---|---|
| **Read H5AD** | ✅ Pure R, rhdf5 | ✅ Uses reticulate | ✅ Uses h5Seurat intermediate |
| **Write H5AD** | ✅ **v1.1.0** | ✅ | ✅ (via h5Seurat) |
| **QC Dashboard** | ✅ **schardExplorer** | ❌ | ❌ |
| **Spatial (Visium)** | ✅ | Partial | ❌ |
| **Python dependency** | ❌ None | ✅ Requires | ❌ None |
| **Performance** | Fast (direct HDF5) | Moderate (Python bridge) | Moderate (format conversion) |
| **Maintenance** | Active | Low | Low |

---

## License

GPL-3. Original author: Pavel Mazin ([cellgeni/schard](https://github.com/cellgeni/schard)).
