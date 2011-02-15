\name{attachVarPrep}
\alias{attachVarPrep}
\title{attachVarPrep}
\description{Attach partial derivative and residual variance functions to the nonlinear fitted object}
\usage{
attachVarPrep(object, form = object$call$model, fDerivFixef = NULL, fDerivRanef = NULL, fVarResidual = NULL)
}
\arguments{
  \item{object}{the fitted nlme object}
  \item{form}{the formula used to fit the object, either formula or string used for automated derivation}
  \item{fDerivFixef}{\code{function(nfit,newdata,pred)} of derivatives in respect to fixed effects at newdata}
  \item{fDerivRanef}{\code{function(nfit,newdata,pred)} of derivatives in respect to random effects at newdata}
  \item{fVarResidual}{\code{function(nfit,newdata,pred)} to calculate var(residual) at newdata}
}
\details{For usage with \code{\link{varPredictNlmeGnls}}, this function attaches \itemize{
\item derivative functions
\item residual variance function
\item Variance-Covariance methods to fitted object
}
\describe{ \item{Automatic derivation}{
If proper basic formula is given, \code{fDerivFixef} and \code{fDerivRanef} will be automatically derived from the model.
In this case fixed and random effects have to be specified as list of formulas for each coefficient instead of single formula 
e.g. \code{list(b0~1,b1~1)} instead of \code{list(b0+b1~1)}. 
Up until now, only single level random effects models are supported for automatic derivation.
}}

\describe{\item{Variance of Residuals}{
Providing no argument \code{fResidual} assumes iid residuals, i.e. \code{weights=NULL}. 
For other residual variance models. See e.g. \code{\link{varResidPower}} corresponding to \code{weights=varPower(form=~fitted(.))}
}}}
\value{nfit with additional entry \code{varPrep}, which is a list of 
\item{varFix}{variance-covariance matrix of fixed effects}
\item{varRan}{variance-covariance matrix of random effects}
\item{coefFix}{names of the fixed coefficients in gradiant function}
\item{coefRan}{names of the random coefficients in gradiant function}
\item{gradFix}{derivative function for fixed effects}
\item{gradRan}{derivative function for random effects}
\item{fVarResidual}{function to calculate residual variance}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twNlme-package}}}
\examples{
data(modExampleStem)
#mtrace(.covarMap)
nfit <- attachVarPrep( modExampleStem, form = "b0*dbh^b1*height^b2")

data(Wutzler08BeechStem)
newdata=data.frame(dbh=18.8, height=16.9, age=40, si=30, alt=470)
(uNew <- nfit$varPrep$gradFix(newdata=newdata))
(wNew <- nfit$varPrep$gradRan(newdata=newdata))

newdata=head(Wutzler08BeechStem)
(uNew <- nfit$varPrep$gradFix(newdata=newdata))
(wNew <- nfit$varPrep$gradRan(newdata=newdata))
}
