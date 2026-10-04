\name{sparq_calibrate_threshold}
\alias{sparq_calibrate_threshold}
\alias{sparq_calibration_estimators}
\alias{sparq_calibrate}
\alias{sparq_decide}
\alias{sparq_calibrate_scalar}
\alias{sparq_calibration_policy}
\alias{sparq_finalize_support}
\alias{sparq_reliability_curve}
\title{Calibrate and interpret SPARQ decisions}
\description{
These functions separate continuous support-quality assessment from optional
threshold calibration and scalar correction policy.
}
\section{Calling conventions}{\preformatted{
sparq_calibrate_threshold(training_score, training_truth, ...)
sparq_decide(support_quality, threshold)
sparq_calibrate_scalar(observed_value, perturbed_values,
  reference_type = c("none", "simulation", "replicate", "pathology"), ...)
sparq_reliability_curve(support_quality, true_error, ...)
}}
\arguments{
  \item{training_score}{Continuous support-quality values from training cases.}
  \item{training_truth}{Known binary stable/unstable labels for training cases.}
  \item{support_quality}{One or more support-quality values.}
  \item{threshold}{A calibrated or explicitly prespecified threshold.}
  \item{observed_value}{Reference scalar estimate.}
  \item{perturbed_values}{Scalar results from perturbation runs.}
  \item{reference_type}{Evidence source for whether scalar correction is allowed.}
  \item{true_error}{Known simulation error or latent-instability measure.}
  \item{...}{Bootstrap, confidence, reference, and calibration controls.}
}
\value{
Threshold calibration returns selected thresholds and uncertainty summaries.
\code{sparq_decide} validates and returns a scalar-calibration decision table.
Scalar calibration returns the original or conditionally corrected value with an
explicit status.
Reliability functions return binned association summaries; they do not claim
probability calibration unless separately validated.
}
\examples{
scores <- c(0.95, 0.90, 0.20, 0.10)
truth <- c(1, 1, 0, 0)
calibration <- sparq_calibrate_threshold(scores, truth, B = 20)
calibration$threshold
}
