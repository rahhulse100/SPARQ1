\name{sparq_known_truth_ablation}
\alias{sparq_ablation_score}
\alias{sparq_known_truth_ablation}
\alias{sparq_export_ablation}
\title{Ablate support-quality components}
\description{
Compare the complete support score with versions omitting uncertainty or
reproducibility support and with prespecified fixed-threshold baselines.
}
\usage{
sparq_ablation_score(rows,
  mode = c("full", "no_uncertainty", "no_reproducibility"))
sparq_known_truth_ablation(benchmark_result,
  fixed_thresholds = c(0.25, 0.5, 0.75), ...)
sparq_export_ablation(ablation_result, output_dir)
}
\arguments{
  \item{rows}{Benchmark rows containing recorded support components.}
  \item{mode}{Complete score or one component-deletion variant.}
  \item{benchmark_result}{A Figure 3C--E benchmark result or its result table.}
  \item{fixed_thresholds}{Distinct thresholds chosen before inspecting ablation
  performance.}
  \item{ablation_result}{Output from \code{sparq_known_truth_ablation}.}
  \item{output_dir}{A new directory for exported tables.}
  \item{...}{Bootstrap, confidence, and seed controls.}
}
\value{
An ablation result containing replicate metrics, threshold values, curves, and
summary tables. Exported artifacts retain the score definition and threshold
rule used for every comparison.
}
\examples{
rows <- data.frame(support_quality = 0.8, iteration_support = 1,
  bias_support = 0.9, uncertainty_support = 0.7,
  reproducibility_support = 0.95)
sparq_ablation_score(rows, "no_uncertainty")
}
