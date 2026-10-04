\name{sparq_validate_results}
\alias{sparq_validate_results}
\alias{sparq_validate_columns}
\alias{sparq_safe_cor}
\alias{sparq_mode_sign_fraction}
\alias{sparq_topk_jaccard}
\alias{sparq_topk_replacement}
\title{Validate inputs and calculate small reproducibility utilities}
\description{
Validation and utility functions used by SPARQ workflows and available for
inspection in custom analyses.
}
\section{Calling conventions}{\preformatted{
sparq_validate_results(data, required_cols = c("result_id", "theta_full",
  "theta_pert", "scenario", "iteration"))
sparq_validate_columns(data, required_cols)
sparq_safe_cor(x, y, method = "spearman")
sparq_mode_sign_fraction(x)
sparq_topk_jaccard(full_ids, pert_ids, k = 100)
sparq_topk_replacement(full_ids, pert_ids, k = 100)
}}
\arguments{
  \item{data}{A data frame to validate.}
  \item{required_cols}{Required column names.}
  \item{x, y}{Numeric vectors for a safe correlation.}
  \item{method}{Correlation method passed to \code{stats::cor}.}
  \item{full_ids, pert_ids}{Ranked feature-ID vectors.}
  \item{k}{Top-ranked set size.}
}
\value{
Validation functions return invisibly on success or stop with an actionable
error. Other functions return a numeric statistic or \code{NA} when it is not
defined.
}
\examples{
sparq_topk_jaccard(c("A", "B", "C"), c("B", "A", "D"), k = 2)
sparq_safe_cor(c(1, 2, 3), c(3, 2, 1))
}
