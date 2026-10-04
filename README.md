<img src="docs/assets/sparq-hero.png" align="right" width="150" alt="SPARQ spatial-transcriptomics perturbation illustration">

# SPARQ

**Spatial Perturbation Analysis for Reliability Quantification**

[![R-CMD-check](https://github.com/rahhulse100/SPARQ1/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rahhulse100/SPARQ1/actions/workflows/R-CMD-check.yaml)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.1-276DC3)](https://www.r-project.org/)

<br clear="right">

SPARQ asks whether an analysis result remains consistent when plausible spot loss is introduced into a spatial-transcriptomic section. It does not denoise the data, replace an analysis method, or establish biological truth. It quantifies how sensitive a specified result is to a defined perturbation model.

SPARQ supports four built-in output types:

- **scalar**: one numeric estimate, score, or statistic;
- **ranked**: an ordered feature, pathway, or interaction list;
- **partition**: spot-level domain or cluster labels;
- **graph**: canonical edge IDs describing a network or neighborhood graph.

## Installation

Install the development package from GitHub:

```r
install.packages("remotes")
remotes::install_github("rahhulse100/SPARQ1")
library(SPARQ)
```

## Standard workflow

Give SPARQ a spot-level data frame and a function that reruns the analysis after rows have been removed. For a scalar result:

```r
set.seed(1)
spot_data <- data.frame(
  signature_score = rnorm(100),
  x = runif(100),
  y = runif(100)
)

assessment <- sparq_assess_workflow(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  output_type = "scalar",
  stress_model = "uniform_random",
  retention = 0.75,
  n_iterations = 100,
  result_id = "section_001",
  reference_scale = 1,
  seed = 1
)

assessment$summary
assessment$comparisons
```

SPARQ prints a compact execution log by default: the reference analysis,
perturbation progress, and the final support result. Printing the assessment
itself gives the same decision-relevant summary without requiring users to
search through a wide table. Set `verbose = FALSE` for notebooks, batch jobs,
or programmatic runs where console output is not wanted.

```r
assessment

# <SPARQ workflow assessment: section_001>
#   Output type: scalar
#   Support quality: 0.91 (adequate support)
#   Completed perturbations: 100/100; failed: 0
```

For domain analysis, the supplied analysis function must return a named vector of spot IDs and domain labels. For ranked and graph analyses, it must return an ordered unique feature vector or canonical edge-ID vector, respectively.

## Reusing an existing R script

`sparq_assess_script()` allows a user to reuse a script rather than writing an R function. SPARQ injects the current perturbed data under `input_object`; the script must create `result_object`.

```r
assessment <- sparq_assess_script(
  data = spot_data,
  script = "banksy_analysis.R",
  input_object = "seurat_obj",
  result_object = "banksy_domains",
  output_type = "partition",
  retention = 0.75,
  n_iterations = 100
)
```

The script must use the injected object rather than repeatedly reading an unchanged input file.

## Spatial-object adapters

`sparq_assess_spatial()` accepts a Seurat object, SpatialExperiment object,
Giotto object, or a coordinate-bearing data frame. SPARQ retains the original
object and subsets it for every perturbation, so the supplied analysis function
receives the object type it already expects.

```r
assessment <- sparq_assess_spatial(
  object = seurat_spatial_object,
  analysis_function = function(x) {
    # Rerun BANKSY or another domain method on x.
    # Return named spot labels, e.g. setNames(new_domains, colnames(x)).
  },
  output_type = "partition",
  retention = 0.75,
  n_iterations = 100
)
```

For custom spatial classes, provide an adapter with an `extract()` function
that returns spot IDs and named coordinates, plus a `subset()` function that
returns the object restricted to retained spot IDs. See
`?sparq_prepare_spatial` for the exact contract.

## Adaptive stopping

For expensive analyses, `sparq_assess_adaptive()` evaluates the support score in batches and stops when a bootstrap interval lies entirely above or below a user-specified decision threshold. It does not infer reliability from a single baseline analysis.

```r
adaptive_assessment <- sparq_assess_adaptive(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  comparator = sparq_compare_scalar,
  retention = 0.75,
  support_threshold = 0.80,
  batch_size = 10,
  min_iterations = 50,
  max_iterations = 100
)

adaptive_assessment$adaptive_history
adaptive_assessment$adaptive_decision
```

## Reports and audit bundles

`sparq_report()` returns an inspectable report object. With an `output_dir`, it
writes a self-contained audit bundle: the support summary, perturbation table,
failure ledger, settings, reference result, any retained perturbation results,
fragility tables, and a PDF overview. SPARQ never creates an output directory
unless one is explicitly supplied.

```r
report <- sparq_report(
  assessment,
  output_dir = "sparq_section_001_report"
)

report
report$output_files
```

## Precomputed outputs

If reference and perturbed outputs have already been generated elsewhere, use `sparq_assess_precomputed()` or `sparq_assess_precomputed_files()`. A single baseline output alone cannot be assessed for perturbation robustness.

## Interpretation

SPARQ returns a continuous support-quality score, instability summaries, failed-iteration records, and—when perturbed outputs are retained—fragility localization. A high score supports consistency under the specified perturbation model; it does not prove a result is biologically correct or unchanged under every possible source of variation.

## Development status

SPARQ is an installable R source package. The current development focus is broader workflow adapters, adaptive stopping validation for every output type, documentation, and external biological validation.
