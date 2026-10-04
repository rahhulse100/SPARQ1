\name{sparq_prepare_ranked_feature_input}
\alias{sparq_prepare_ranked_feature_input}
\alias{sparq_validate_ranked_feature}
\alias{sparq_rank_turnover}
\alias{sparq_bias_ranked_feature}
\alias{sparq_calibrate_ranked_feature}
\alias{sparq_sample_fragility_ranked_feature}
\alias{sparq_pipeline_ranked_feature}
\alias{sparq_bias}
\alias{sparq_support}
\alias{sparq_reproducibility_feature}
\alias{sparq_pipeline}
\title{Analyze reproducibility of ranked features}
\description{
Prepare and assess repeated ranked-feature outputs, including feature turnover,
bias, uncertainty, reproducibility, and optional scalar correction policy.
}
\usage{
sparq_prepare_ranked_feature_input(data, ...)
sparq_rank_turnover(data, top_k = 100, ...)
sparq_pipeline_ranked_feature(data, top_k = 100, ...)
sparq_pipeline(data, ...)
}
\arguments{
  \item{data}{A result table containing the columns required by the selected
  ranked-feature helper.}
  \item{top_k}{Number of top-ranked features to evaluate.}
  \item{...}{Column names, reference scale, and calibration controls.}
}
\value{
Functions return standardized feature-level tables or a combined pipeline
result. They require \pkg{data.table}, which is installed with SPARQ.
}
\examples{
\dontrun{
# See the ranked-feature vignette for a complete input table and pipeline.
ranked_result <- sparq_pipeline_ranked_feature(feature_results, top_k = 100)
}
}
