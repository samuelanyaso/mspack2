#' Simulated clustered current-status from a three-state tracking model
#'
#' Current status data from 200 clusters with informative cluster size from a three-state tracking model.
#' The cluster-correlated exit times from state 1 are simulated from a lognormal AFT model with two covariates:
#' a cluster-level exposure covariate denoted by \emph{Z1} and a subject-level continuous covariate denoted by \emph{Z2}.
#' The exit times from state 2 are simulated by a transformation which ensures that the exit times from state 2 are
#' greater than the exit times from state 1. The clustered inspection times are generated from a Weibull
#' distribution with parameters; shape=3 and scale=5. The cluster sizes are generated from a Poisson distribution where the mean
#' depends on \emph{Z1} and a cluster-specific random effects term, thus, inducing informative cluster size.
#' The data contains the IDs, inspection times states in which the subjects where occupying at their inspection time,
#' and information on the covariates..
#'
#' @format A data frame with 1389 rows and 7 variables:
#' \describe{
#'   \item{cID}{Cluster ID}
#'   \item{id}{Subject ID}
#'   \item{time}{The current status (inspection) time of subject j in cluster i.}
#'   \item{state}{State occupied when inspected. States 1, 2, 3 and 4 correspond to 'No', 'Mild', 'Moderate' and 'Severe' level of periodontal disease. State 4 is the absorbing state}
#'   \item{Z1}{Cluster-level exposure variable, 0: Unexposed, 1: Exposed}
#'   \item{Z2}{Subject-level continuous variable simulated from a normal distribution with mean=1 and sd=0.15}
#'   \item{csize}{Cluster size}
#'   ...
#' }
"CSdata"
