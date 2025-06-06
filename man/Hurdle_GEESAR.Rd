% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hurdle_GEESAR.R
\name{Hurdle_GEESAR}
\alias{Hurdle_GEESAR}
\title{Hurdle Model using GEESAR}
\usage{
Hurdle_GEESAR(
  formula,
  data,
  W,
  weights = NULL,
  toler = 1e-05,
  maxit = 200,
  trace = FALSE
)
}
\arguments{
\item{formula}{A formula specifying the model.}

\item{data}{The dataset.}

\item{W}{The spatial weight matrix.}

\item{weights}{Optional weights.}

\item{toler}{Convergence tolerance.}

\item{maxit}{Maximum number of iterations.}

\item{trace}{Logical. If TRUE, prints progress.}
}
\value{
A list containing the logit and Poisson-truncated models.
}
\description{
This function fits a hurdle model using GEESAR, consisting of:
(1) A logit model for zero vs. non-zero responses.
(2) A truncated Poisson model for positive counts.
}
\examples{
\donttest{
set.seed(123)
n <- 100
x <- rnorm(n)
y <- rpois(n, lambda = exp(0.5 * x))
y[rbinom(n, 1, 1/(1+exp(-0.5*x)))] <- 0  # Introduce zeros
W <- matrix(rbinom(n^2,1,0.2), n, n)  # Example spatial weight matrix
diag(W) <- 0
rtot <- rowSums(W)
W <- W/ifelse(rtot==0, 0.1, rtot)
model <- Hurdle_GEESAR(y ~ x, data = data.frame(y, x), W = W)
summary_SAR(model$logit_model)
summary_SAR(model$poisson_truncated_model)
}

}
