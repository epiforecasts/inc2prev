% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/incidence.R
\name{incidence}
\alias{incidence}
\title{Incidence estimation model}
\usage{
incidence(
  prev,
  ab = NULL,
  vacc = NULL,
  init_ab = NULL,
  prob_detect,
  data_args = list(),
  model = i2p_model(),
  var_col = NULL,
  variables = NULL,
  quantiles = seq(0.05, 0.95, by = 0.05),
  samples = 100,
  keep_fit = FALSE,
  p,
  ...
)
}
\arguments{
\item{prev}{Observed positiviy prevalence data}

\item{ab}{Observed antibody prevalence data}

\item{vacc}{Observed vaccination data}

\item{init_ab}{Observed initial antibody data}

\item{prob_detect}{Time-varying probability of detection}

\item{data_args}{A list of arguments to pass to \code{i2p_data()}}

\item{model}{A stan model object as produced by \code{i2p_model()}
or a similar compiled stan model.}

\item{var_col}{Column of a common name in the data frames passed as data
that should be treated as indicator of multiple data sets that should be used
jointly for estimation}

\item{variables}{A vector of variables to summarise and return}

\item{quantiles}{A vector of quantiles to return in summarising}

\item{samples}{The number of samples to return}

\item{keep_fit}{Whether to retain the stan fit object}

\item{p}{A \code{progressr} function used when fitting multiple models
to track progress}

\item{...}{Additional arguments passed to \code{fit_fn}.}
}
\value{
A data.table including summarised estimates and parameter samples.
}
\description{
Fits an \code{inc2prev} model to recover incidence
from prevalence data using a probility of detection curve.
}
\seealso{
Other incidence: 
\code{\link{combine_incidence_fits}()},
\code{\link{incidence_lapply}()}
}
\concept{incidence}
