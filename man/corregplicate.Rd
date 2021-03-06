% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corregp.r
\name{corregplicate}
\alias{corregplicate}
\title{Repeated Correspondence Regression}
\usage{
corregplicate(formula, data, part = NULL, b = 100, r = 10, xep = TRUE,
  std = FALSE, rel = TRUE, phi = FALSE, chr = ".",
  b_scheme = "multinomial")
}
\arguments{
\item{formula}{A \code{\link[stats]{formula}} specification of which factors to cross with each other. The left-hand (\code{y}) side must be a single factor.
The right-hand side (\code{x}) can involve all the usual specifications of interactions and/or nested analyses.}

\item{data}{The data frame containing the variables specified in the \code{formula}.}

\item{part}{Character vector specifying the names of conditional factors (e.g. a factor partioning the levels of the left-hand side \code{y} into groups).
This argument is relevant for analyses in which one wants to remove between-item variation.}

\item{b}{Number of the bootstrap replications (simulations).}

\item{r}{Number of repeated calls to \code{\link{corregp}}.}

\item{xep}{Logical specifying whether to output the separate terms in the right-hand side (\code{x}) as components in a list.
If \code{FALSE}, then all \code{x} output is collected in a matrix.}

\item{std}{Logical specifying whether to output the standardized coordinates. Defaults to \code{FALSE}.}

\item{rel}{Logical specifying whether to divide the coordinates by the \code{sqrt} of their totals, so that one obtains coordinates for
the relative frequencies (as is customary in correspondence analysis). Defaults to \code{TRUE}.}

\item{phi}{Logical specifying whether to compute the output on the scale of the \emph{Chi-squared} value of the contingency table or of the \emph{Phi-squared} value
(which is \emph{Chi-squared} divided by \emph{N}). Reminiscent of \code{\link[MASS]{corresp}} in package \pkg{MASS}, defaults to \code{FALSE}.}

\item{chr}{Character specifying the separator string for constructing the interaction terms.}

\item{b_scheme}{Character specifying the sampling scheme for bootstrapping. Must match either \code{"multinomial"} (the default) or \code{"product-multinomial"}.}
}
\value{
An object of class "corregp" in which the bootstrap replications of all the repeated calls to \code{corregp} are put together.
}
\description{
A function for repeated correspondence regressions with bootstrapping in order to handle large data sets. This is essentially a wrapper
  \code{replicate(n = r, expr = corregp(...), simplify = FALSE)}, so it may dissappear in the future.
}
\seealso{
\code{\link{corregp}}.
}
