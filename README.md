## SCOUT: Ornstein–Uhlenbeck modelling of gene expression evolution on single-cell lineage trees

Given a single-cell lineage tree, fit gene expression to evolutionary models to profile selection and drift.

<img src="https://github.com/hrstuart/SCOUT/blob/main/data/F1_GraphicalAbstract.svg" alt="Graphical abstract of SCOUT" width="450" />

### Overview

SCOUT is an R package that applies Ornstein-Uhlenbeck (OU) modeling to analyze gene expression evolution along single-cell lineage trees. The package enables researchers to:

- Fit multiple evolutionary models (Brownian Motion, Ornstein-Uhlenbeck) to gene expression data
- Compare models to identify genes under selection vs. neutral drift
- Analyze regime-specific evolutionary dynamics across lineage trajectories
- Perform model selection using information criteria (AICc)

### Features

- **Multiple evolutionary models**: Support for BM1 (single-rate Brownian Motion), OU1 (single-optimum OU), and OUM (multi-optimum OU) models
- **Parallel processing**: Efficient model fitting across multiple genes using parallel computation
- **Flexible data preparation**: Built-in normalization, smoothing, and ancestral state inference
- **Comprehensive model selection**: Automated calculation of AICc weights, delta AICc, and other fit metrics
- **Regime annotation**: Support for custom regime classifications across the lineage tree

### Installation

```r
# Install dependencies
devtools::install_github("Galaxeee/TedSim")

# Install SCOUT from GitHub
devtools::install_github('https://github.com/hrstuart/SCOUT/')
```

### Dependencies

SCOUT requires the following R packages:
- `dplyr` - Data manipulation
- `ape` - Phylogenetic analyses
- `OUwie` - OU model fitting
- `foreach`, `doParallel`, `parallel` - Parallel processing
- `castor` - Tree manipulation and ancestral state reconstruction
- `reshape2` - Data reshaping
- `caret` - Model evaluation
- `TedSim` - Simulation framework (install from GitHub: `Galaxeee/TedSim`)

### Quick Start

#### Step 1: Prepare your data

SCOUT requires two main inputs:

1. **Phylogenetic tree**: A Newick format file containing the single-cell lineage tree
2. **Metadata**: A CSV file with cell barcodes, gene expression values, and regime annotations

#### Step 2: Prepare inputs for model fitting

```r
# Prepare data with normalization and smoothing
inputs <- prepare_data(
  tree_path = tree,           # path to newick file
  metadata_path = meta,       # path to metadata file
  outpath = outdir,           # output directory
  regimes = c("BM1", "OU1", "OUM"),  # substitute OUM for column names in metadata with regime annotation
  normalize = TRUE,           # if using counts (applies log normalization)
  smoothing = 10              # smoothing parameter for count data
)
```

**Example metadata format:** 
| cellBC             | OU4 | HES4      | ISG15     | AGRN      | SDF4      |
|--------------------|-----|-----------|-----------|-----------|-----------|
| ACGGAGAGTAAGCACG-1 | LL  | 0         | 0         | 0         | 0.577392  |
| ACTTACTAGGAATGGA-1 | LL  | 0.7890783 | 0.1583906 | 0.5226525 | 0.1583906 |
| AGCGGTCCAACACGCC-1 | RL  | 0.197258  | 0.3619425 | 0         | 0.6271318 |
| AGTGTCATCGGAGCAA-1 | M   | 0.3356066 | 0.3356066 | 0.5864398 | 0.3356066 |
| ATAGACCCAGCCAGAA-1 | LL  | 0         | 0.7329661 | 0.3077897 | 0         |
| CATCGGGCATGACATC-1 | Liv | 0.3847468 | 0.3847468 | 0         | 0.6619065 |

- The first column contains cell barcodes
- `OU4` is a regime column (containing regime labels like LL, RL, M, Liv)
- Remaining columns contain gene expression values (can be raw counts or normalized)

#### Step 3: Fit models in parallel

```r
# Fit evolutionary models to gene expression data
results <- fitModel(
  inputs,
  cores = 4,                          # number of parallel cores
  write = TRUE,                       # save results to disk
  outpath = outdir,                   # output directory
  testgenes = c('HES4', 'ISG15', 'AGRN', 'SDF4')  # genes to test (leave empty to test all)
)
```

#### Step 4: Calculate model fit metrics

```r
# Summarize and calculate fit metrics (AICc weight, delta AICc, etc.)
results <- get_fit_metrics(results, write = FALSE)
```

#### Step 5: Perform model selection

```r
# Returns dataframe with final model per gene
final_results <- process_results(results)
```

### Understanding the Models

- **BM1** (Brownian Motion): Models trait evolution as a random walk with constant variance. Assumes no selection.
- **OU1** (Single Ornstein-Uhlenbeck): Models evolution toward a single optimum with selection strength α and trait optimum θ.
- **OUM** (Multiple Optima OU): Models evolution with different trait optima for different regimes (e.g., different cell types or lineages).

### Additional Parameters

The `prepare_data()` function supports many optional parameters:

- `species_key`: Column name for species/cell ID (default: looks for 'species' column)
- `anc_infer`: Ancestral state inference method ('ape', 'castor', or 'skip')
- `resolve`: Tree multifurcation resolution method
- `algorithm`: OUwie algorithm to use (default: 'three.point')
- `floor`: Numeric floor for zero counts
- `blacklist`: Vector of column names to exclude from analysis

### Examples and Vignettes

See detailed examples in the `examples/` directory:
  - **SimulationVignette.ipynb**: Walkthrough using simulated data to demonstrate SCOUT workflow
  - **simulate_data/**: Scripts for generating simulation datasets

### Output

SCOUT produces several outputs:

- Model fit objects for each gene-regime combination
- AICc values and weights for model comparison
- Parameter estimates (α, σ², θ) for OU models
- Final model selection results with best-fit model per gene

### Citation

If you use SCOUT in your research, please cite:

[Citation information to be added]

### License

MIT License

### Issues and Support

For bug reports and feature requests, please open an issue on the [GitHub repository](https://github.com/hrstuart/SCOUT/issues). 

