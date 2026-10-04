\name{sparq_define_loss_model}
\alias{sparq_define_loss_model}
\alias{sparq_loss_probabilities}
\alias{sparq_apply_loss_model}
\title{Define and apply spatial spot-loss mechanisms}
\description{
Defines explicit fixed-retention spot-loss mechanisms for SPARQ perturbation
assessment. The mechanism is recorded with the assessment so that a reliability
claim is conditional on a transparent sampling-loss distribution.
}
\usage{
sparq_define_loss_model(
  type = c("uniform_random", "contiguous_hole", "quality_weighted",
    "edge_weighted", "custom"),
  quality_col = NULL,
  quality_direction = c("higher_quality_retained", "lower_quality_retained"),
  quality_strength = 1,
  edge_strength = 1,
  custom_function = NULL
)

sparq_loss_probabilities(data, loss_model, x_col = NULL, y_col = NULL)

sparq_apply_loss_model(data, loss_model, retention = 0.75,
  x_col = NULL, y_col = NULL, iteration = 1, seed = NULL)
}
\arguments{
  \item{type}{The loss mechanism. \code{quality_weighted} samples retained
  spots by ranked quality; \code{edge_weighted} preferentially removes spots
  farther from the coordinate-wise median tissue center.}
  \item{quality_col}{One finite numeric quality column for
  \code{quality_weighted} loss.}
  \item{quality_direction}{Whether larger quality values increase or decrease
  retention weight.}
  \item{quality_strength}{Nonnegative exponent controlling quality weighting.}
  \item{edge_strength}{Nonnegative edge-removal strength.}
  \item{custom_function}{Function returning exactly the retained integer row
  indices for custom loss.}
  \item{data}{Spot-level data frame.}
  \item{loss_model}{A \code{sparq_loss_model} created by
  \code{sparq_define_loss_model}.}
  \item{x_col, y_col}{Finite numeric coordinate columns required for spatial
  loss mechanisms.}
  \item{retention}{Fraction of spots retained; the retained count is fixed as
  \code{floor(n * retention)}, with a minimum of two spots.}
  \item{iteration, seed}{Reproducibility controls.}
}
\value{
\code{sparq_define_loss_model} returns a loss-model object.
\code{sparq_loss_probabilities} returns per-spot retention and relative loss
weights. \code{sparq_apply_loss_model} returns retained data, retained and
removed row indices, and a per-spot audit table.
}
\details{
These models specify a perturbation distribution; they do not establish that a
particular loss mechanism is the true experimental missingness process. Users
should choose a mechanism justified by the assay, tissue quality measurements,
or the scientific stress question.
}
\examples{
spots <- data.frame(
  quality = 1:20,
  x = runif(20),
  y = runif(20)
)
model <- sparq_define_loss_model(
  "quality_weighted",
  quality_col = "quality"
)
loss <- sparq_apply_loss_model(spots, model, retention = 0.75, seed = 1)
loss$audit
}
