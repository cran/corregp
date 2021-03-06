% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corregp.r
\name{confGet}
\alias{confGet}
\title{Getting \code{conf} Components from \code{corregp} Objects}
\usage{
confGet(crg, parm)
}
\arguments{
\item{crg}{The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").}

\item{parm}{The parameter for which to retrieve the \code{conf} components. Can be either \code{"y"}, \code{"x"}, or any vector of term names in X, level names
in X or level names in Y.}
}
\value{
A list of components selected with \code{parm}.
}
\description{
Internal function for retrieving the \code{conf} component(s) in a \code{corregp} object.
}
\details{
\code{confGet} is an internal function to be called by \code{\link{cint.corregp}}, \code{\link{cell.corregp}} or \code{\link{cell3d.corregp}}, but not by users.
}
