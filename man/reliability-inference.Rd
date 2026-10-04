\name{sparq_reliability_confidence_sequence}
\alias{sparq_reliability_confidence_sequence}
\title{Anytime-valid sequential reliability inference}
\description{
Computes a beta-binomial mixture confidence sequence for the probability that a
perturbed analysis result remains within a prespecified normalized comparator
distance. The confidence sequence is valid across repeated looks at the same
accumulating perturbation stream and is used by adaptive reliability stopping.
}
\usage{
sparq_reliability_confidence_sequence(
  comparison_table, tolerance, required_reliability = 0.9,
  confidence = 0.95, reference_scale = 1,
  prior_success = 0.5, prior_failure = 0.5
)
}
\arguments{
  \item{comparison_table}{One row per attempted perturbation, with \code{status}
  and \code{instability} columns. Failed rows count as non-preserving attempts.}
  \item{tolerance}{Prespecified maximum normalized comparator distance that still
  preserves the result.}
  \item{required_reliability}{Minimum preservation probability needed for a
  \code{supported} decision.}
  \item{confidence}{Simultaneous confidence level across sequential looks.}
  \item{reference_scale}{Scale used to normalize comparator distance.}
  \item{prior_success, prior_failure}{Fixed positive beta-mixture parameters.
  The default is the Jeffreys mixture.}
}
\value{
A data frame with one row per perturbation count. It includes preservation
counts, an anytime-valid lower and upper confidence bound, and a
\code{supported}, \code{not_supported}, or \code{inconclusive} decision.
}
\details{
For each perturbation \eqn{k}, SPARQ forms a preservation event
\deqn{Z_k = 1\{d_k / c \leq \delta\},}
where \eqn{d_k} is comparator instability, \eqn{c} is \code{reference_scale},
and \eqn{\delta} is \code{tolerance}. The beta-binomial mixture confidence
sequence controls coverage under repeated interim looks, unlike repeatedly
applying a fixed-time interval and stopping when it crosses a boundary.
}
\examples{
comparisons <- data.frame(
  status = rep("ok", 20),
  instability = c(rep(0.02, 19), 0.20)
)
sequence <- sparq_reliability_confidence_sequence(
  comparisons,
  tolerance = 0.05,
  required_reliability = 0.70
)
tail(sequence, 1)
}
