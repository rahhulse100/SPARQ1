\name{sparq_assess_frontier}
\alias{sparq_assess_frontier}
\alias{print.sparq_reliability_frontier}
\alias{summary.sparq_reliability_frontier}
\alias{plot.sparq_reliability_frontier}
\alias{sparq_infer_cohort_reliability}
\alias{print.sparq_cohort_reliability}
\alias{summary.sparq_cohort_reliability}
\alias{as.data.frame.sparq_cohort_reliability}
\title{Estimate reliability frontiers and cohort-level reliability}
\description{
Evaluates a conclusion over a fixed grid of retained spot fractions and combines
section-level preservation counts using an empirical-Bayes beta-binomial model.
}
\usage{
sparq_assess_frontier(
  data, analysis_function,
  output_type = c("scalar", "ranked", "partition", "graph"),
  retentions = c(0.5, 0.6, 0.7, 0.8, 0.9, 1),
  stress_model = c("uniform_random", "contiguous_hole", "none", "custom"),
  reliability_tolerance, required_reliability = 0.9, confidence = 0.95,
  top_k = 100, x_col = NULL, y_col = NULL, custom_function = NULL,
  loss_model = NULL, result_id = "result", reference_scale = 1,
  seed = 1, verbose = TRUE, ...
)

sparq_infer_cohort_reliability(
  cohort, confidence = 0.95,
  section_metadata = NULL, sample_id_col = "sample_id"
)
}
\arguments{
  \item{data}{Spot-level data frame for one section.}
  \item{analysis_function}{Function rerunning the analysis after spot loss.}
  \item{output_type}{Result type accepted by \code{sparq_assess_workflow}.}
  \item{retentions}{Unique evaluated retained fractions in \code{(0, 1]}.}
  \item{stress_model}{Built-in stress model used when \code{loss_model} is
  \code{NULL}.}
  \item{reliability_tolerance, required_reliability}{Reliability decision
  settings.}
  \item{confidence}{Family-wise confidence level across the prespecified
  retention grid. SPARQ uses a Bonferroni allocation before certifying a
  critical retention. For cohort inference, this is the posterior interval
  level.}
  \item{top_k, x_col, y_col, custom_function, loss_model, result_id,
  reference_scale, seed, verbose, ...}{Assessment settings forwarded to
  \code{sparq_assess_workflow}.}
  \item{cohort}{A \code{sparq_cohort_assessment} or a named list of assessments
  containing reliability certificates.}
  \item{section_metadata}{Optional one-row-per-section metadata.}
  \item{sample_id_col}{Section ID column in \code{section_metadata}.}
}
\value{
\code{sparq_assess_frontier} returns a frontier table, all individual
assessments, and a discrete critical retention. The critical retention is the
lowest evaluated fraction supported at every higher evaluated fraction. The
critical-retention statement has simultaneous coverage across the supplied grid
through a Bonferroni allocation of \code{confidence}; SPARQ does not
interpolate between the supplied grid points.

\code{sparq_infer_cohort_reliability} returns raw section certificates,
empirical-Bayes section posterior estimates, and a cohort beta-binomial summary.
Its intervals are posterior intervals, not independent biological replication
of spots.
}
\details{
The beta-binomial model represents uncertainty from the finite number of
perturbations within each section separately from observed heterogeneity in
preservation probability across sections. It does not treat spots as independent
biological samples.
}
\examples{
set.seed(1)
spots <- data.frame(value = rnorm(40), x = runif(40), y = runif(40))
frontier <- sparq_assess_frontier(
  spots,
  analysis_function = function(x) mean(x$value),
  output_type = "scalar",
  retentions = c(0.75, 1),
  reliability_tolerance = 0.2,
  required_reliability = 0.5,
  preset = "quick",
  n_iterations = 10,
  min_iterations = 5,
  verbose = FALSE
)
frontier$frontier
}
