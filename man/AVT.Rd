% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corregp.r
\docType{data}
\name{AVT}
\alias{AVT}
\title{The Use of Linguistic Variants in Audio-Visual Translation (Subtitles)}
\format{A data frame with 3302 rows and 7 variables.
\itemize{
  \item{\code{Variant} }{The linguistic variant used in a set of alternatives (27 levels).}
  \item{\code{Variable} }{The linguistic variable specifying a set of alternatives (13 levels).}
  \item{\code{Variety} }{The dichotomization of \code{Variant} into standard and non-standard.}
  \item{\code{Speaker} }{The role of the speaker in the data (2 levels).}
  \item{\code{Language} }{The language (and source language) of the data (3 levels).}
  \item{\code{Language2} }{The same as \code{Language} but with the observations of level \code{intra.nl} set to \code{NA}.}
  \item{\code{Genre} }{The genre or register of the data (2 levels).}
}}
\source{
Prieels, L., I. Delaere, K. Plevoets and G. De Sutter (2015) A corpus-based multivariate analysis of linguistic norm-adherence in audiovisual and written
  translation. \emph{Across Languages and Cultures} \strong{16} (2), 209--231.
}
\description{
This data set was a follow-up study to the \code{\link{COMURE}} project and was conducted at the Department of Translation, Interpreting and Communication of
  Ghent University between 2014 and 2018.
}
\examples{
\donttest{
data(AVT)
# The execution of corregp may be slow, due to bootstrapping:  
avt.crg <- corregp(Variant ~ Speaker * Language * Genre, data = AVT, part = "Variable", b = 3000)
avt.crg
summary(avt.crg, parm = "b", add_ci = TRUE)
screeplot(avt.crg, add_ci = TRUE)
anova(avt.crg, nf = 2)
avt.col <- ifelse( xtabs(~ Variant + Variety, data = AVT)[, "Standard"] > 0, "blue", "red")
plot(avt.crg, x_ell = TRUE, xsub = c("Speaker", "Language", "Genre"), col_btm = avt.col, 
  col_top = "black")
}
}
