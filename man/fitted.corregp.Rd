% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corregp.r
\name{fitted.corregp}
\alias{fitted.corregp}
\alias{fitted.values.corregp}
\title{Extracting Fitted Values from Correspondence Regression}
\usage{
\method{fitted}{corregp}(object, parm = "all", nf = NULL, ...)

\method{fitted.values}{corregp}(object, parm = "all", nf = NULL, ...)
}
\arguments{
\item{object}{The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").}

\item{parm}{The parameter for which to extract the fitted values. Can be \code{"all"}, \code{"both"} (or abbreviations), \code{"y"} or \code{"x"} for
the fitted values of every cell in the data, but it can also be any vector of term names in X or level names in X. Defaults to \code{"all"}.}

\item{nf}{The number of dimensions to be retained in the reduced space. Defaults to all dimensions (no reduction).}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
A matrix or vector with the fitted values for the parameters of interest, based on the selected number of dimensions.
}
\description{
Method to extract the fitted values of a correspondence regression.
}
\details{
If all dimensions are retained, then the fitted values will only be equal to the observed counts if no conditioning factors were specified with the argument
  \code{"part"} in the \code{\link{corregp}} call. This is because the associations with the conditioning factors (in \code{"part"}) are not taken into
  account.
}
\examples{
\donttest{
data(HairEye)
haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
fitted(haireye.crg, parm = c("Hair", "Sex"), nf = 2)
fitted.values(haireye.crg, parm = c("Hair", "Sex"), nf = 2)
}
}
\seealso{
\code{\link{coef.corregp}}, \code{\link{residuals.corregp}}.
}
