\name{varResidPower}
\alias{varResidPower}
\title{varResidPower}
\description{Variance of residual of prediction for Power Variance without random effects e.g. \code{weights=varPower(form=~fitted(.))}}
\usage{
varResidPower(object, newdata, pred)
}
\arguments{
  \item{object}{the fitted object}
  \item{newdata}{the data frame with new predictors}
  \item{pred}{the prediction at population level}
}

\value{numeric vector of var(eps_i)}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twNlme-package}}}

