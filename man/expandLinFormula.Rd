\name{expandLinFormula}
\alias{expandLinFormula}
\title{expandLinFormula}
\description{Extends the linear formula to an expression involving coefficients}
\usage{expandLinFormula(linForm, suffix = "", varNames = NULL)}
\arguments{
  \item{linForm}{formula for fixed coefficients depending on linear term}
  \item{suffix}{base parameter name}
  \item{varNames}{variable names to use}
}
\details{parameter names will be lhs of the formula + suffix + .i wiht .0 for the intercept
if there is no intercept, starting from .1}


\author{Thomas Wutzler}



\seealso{\code{\link{twNlme-package}}}
\examples{
expandLinFormula(b0~1)
# note that si*age treated as multiplication instead of all interactions
expandLinFormula(linForm <- b0~si+log(age)+I(si+age)+si*age)	 
}
