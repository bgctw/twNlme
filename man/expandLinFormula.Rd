\name{expandLinFormula}
\alias{expandLinFormula}
\title{expandLinFormula}
\description{Extends the linear formula to an expression involving coefficients}
\usage{
expandLinFormula(linForm, suffix = "", varNames = NULL)
}
\arguments{
  \item{linForm}{formula for fixed coefficients depending on linear term see \code{\link{.extractFixedList}}}
  \item{suffix}{base parameter name}
  \item{varNames}{variable names to use}
}
\details{parameter names will be lhs of the formula + suffix + .i wiht .0 for the intercept
if there is no intercept, starting from .1}
\value{String of formula with coefficients expanded to linear combinations of covariates and random effects.}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twNlme-package}}}
\examples{
data(modExampleStem)
fixForm <- twNlme:::.extractFixedList(modExampleStem)
expandLinFormula(fixForm[[1]])
}
