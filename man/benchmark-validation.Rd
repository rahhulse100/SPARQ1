\name{sparq_known_truth_benchmark}
\alias{sparq_known_truth_benchmark}
\alias{sparq_known_truth_benchmark_parallel}
\alias{sparq_methods_design}
\alias{sparq_read_design}
\alias{sparq_write_design}
\alias{sparq_export_benchmark}
\alias{sparq_known_truth_stress}
\alias{sparq_stress_response}
\alias{sparq_export_stress}
\alias{sparq_binary_metrics}
\alias{sparq_curves}
\alias{sparq_metric_summary}
\title{Known-truth benchmarking and stress-response analysis}
\description{
Build, validate, summarize, and export the versioned synthetic benchmarks used
to test whether support quality tracks defined instability conditions.
}
\section{Calling conventions}{\preformatted{
sparq_methods_design()
sparq_read_design(path)
sparq_write_design(path, design = sparq_methods_design())
sparq_known_truth_benchmark(benchmark_data, ...)
sparq_known_truth_benchmark_parallel(benchmark_data, n_workers = 1, ...)
sparq_known_truth_stress(benchmark_result, ...)
sparq_binary_metrics(truth, score, threshold)
sparq_curves(truth, score)
sparq_metric_summary(metrics, ...)
}}
\arguments{
  \item{path}{Path to a versioned SPARQ design file.}
  \item{design}{A design produced by \code{sparq_methods_design}.}
  \item{benchmark_data}{Figure 3B-style benchmark rows with IDs, difficulty,
  and latent instability.}
  \item{n_workers}{Number of workers used for independent benchmark replicates.}
  \item{benchmark_result}{A completed known-truth benchmark.}
  \item{truth}{Known binary labels.}
  \item{score}{Continuous support-quality scores.}
  \item{threshold}{Classification threshold.}
  \item{metrics}{Replicate-level metric table.}
  \item{...}{Versioned design, export path, and bootstrap settings.}
}
\value{
Benchmark functions return per-sample results, replicate metrics, threshold
summaries, reliability summaries, and failure ledgers. Export helpers write
named TSV and RDS artifacts without overwriting a nonempty output directory.
}
\examples{
design <- sparq_methods_design()
design$protocol_version
sparq_binary_metrics(c(1, 1, 0, 0), c(0.9, 0.8, 0.2, 0.1), threshold = 0.5)
}
