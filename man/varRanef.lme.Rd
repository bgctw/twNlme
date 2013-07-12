\name{varRanef.lme}
\alias{varRanef.lme}
\title{varRanef lme}
\description{Variance-Covariance of lme random effects (Psi)}
\usage{\method{varRanef}{lme}(object, varName = NULL, ...)}
\arguments{
  \item{object}{the fitted nlme or lme model}
  \item{varName}{for which parameter}
  \item{\dots}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{twNlme-package}}}
\examples{
fm1 <- lme(distance ~ age, data = Orthodont, random = ~age)
varRanef(fm1)
varFixef(fm1)
fm2 <- gls(distance ~ age, data = Orthodont)
varFixef(fm2)
varRanef(fm2)
}
