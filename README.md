# SPARQ

**Spatial Perturbation Analysis for Reliability Quantification**

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

The script must use the injected object rather than repeatedly reading an unchanged input file. Current generic workflow functions accept a data frame; direct Seurat and Visium object adapters are not yet implemented.

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

## Precomputed outputs

If reference and perturbed outputs have already been generated elsewhere, use `sparq_assess_precomputed()` or `sparq_assess_precomputed_files()`. A single baseline output alone cannot be assessed for perturbation robustness.

## Interpretation

SPARQ returns a continuous support-quality score, instability summaries, failed-iteration records, and—when perturbed outputs are retained—fragility localization. A high score supports consistency under the specified perturbation model; it does not prove a result is biologically correct or unchanged under every possible source of variation.

## Development status

SPARQ is an installable R source package. The current development focus is broader workflow adapters, adaptive stopping validation for every output type, documentation, and external biological validation.
