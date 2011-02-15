\name{twNlme-package}
\alias{twNlme-package}
\alias{twNlme}
\title{Prediction intervals for nlme and gnls models}
\description{Prediction intervals for nlme and gnls models}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}


\details{
The main purpose of the package is to provide standard errors for predictions of 
simple models fitted by \code{\link{nlme}} or \code{\link{gnls}}.

The main function is \code{\link{varPredictNlmeGnls}} using function \code{\link{attachVarPrep}}.

\describe{\item{Further functionality of package \code{twNlme} includes}{
\itemize{
\item Variance of sum of predictions \code{\link{varSumPredictNlmeGnls}}
\item treating gnls models as special random effects models with zero random effects, i.e accessing fixed (coefficients) and random (none) of gls models: \code{\link{fixef.gls}}, \code{\link{ranef.gls}}
\item extracting Variance Covariance matrix of fixed and random effects: \code{\link{varRanef.lme}}, \code{\link{varFixef.lme}} for both lme and gls models 
\item extracting Standard error of fixed effects: varRanef, \code{\link{seFixef.lme}} for both lme and gls models
}
}}

\describe{\item{Data}{
\itemize{
\item tree weights in relation to diameter and height: \code{\link{Wutzler08BeechStem}}
\item nlme model fit to the tree weight data: \code{\link{modExampleStem}}
}
}}
}


\keyword{ package }


