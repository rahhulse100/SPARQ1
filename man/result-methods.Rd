\name{summary.sparq_assessment}
\alias{summary.sparq_assessment}
\alias{summary.sparq_workflow_assessment}
\alias{summary.sparq_script_assessment}
\alias{summary.sparq_spatial_assessment}
\alias{summary.sparq_adaptive_assessment}
\alias{as.data.frame.sparq_assessment}
\alias{as.data.frame.sparq_workflow_assessment}
\alias{as.data.frame.sparq_script_assessment}
\alias{as.data.frame.sparq_spatial_assessment}
\alias{as.data.frame.sparq_adaptive_assessment}
\alias{plot.sparq_assessment}
\alias{plot.sparq_workflow_assessment}
\alias{plot.sparq_script_assessment}
\alias{plot.sparq_spatial_assessment}
\alias{plot.sparq_adaptive_assessment}
\title{Summarize, tabulate, and plot SPARQ assessments}
\description{
Conventional R methods for completed SPARQ assessments.
}
\usage{
\method{summary}{sparq_assessment}(object, ...)
\method{as.data.frame}{sparq_assessment}(x, row.names = NULL, optional = FALSE, ...)
\method{plot}{sparq_assessment}(x, ...)
}
\arguments{
  \item{object, x}{A completed SPARQ assessment returned by a workflow, script,
  spatial-object, adaptive, or direct assessment interface.}
  \item{row.names, optional}{Passed to \code{as.data.frame}.}
  \item{...}{Graphical parameters for \code{plot}, or unused method arguments.}
}
\value{
\code{summary()} and \code{as.data.frame()} return the one-row support summary.
For adaptive assessments, \code{summary()} also includes the stopping decision.
\code{plot()} draws instability across perturbations and the component support
values, and returns the assessment invisibly.
}
\examples{
spots <- data.frame(value = 1:8, x = 1:8, y = 8:1)
fit <- sparq_assess_workflow(
  spots, function(x) mean(x$value), "scalar",
  preset = "quick", min_iterations = 3,
  instability_bootstrap_B = 20, verbose = FALSE
)
summary(fit)
as.data.frame(fit)
if (interactive()) plot(fit)
}
