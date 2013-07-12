\name{attachVarPrep}
\alias{attachVarPrep}
\title{attachVarPrep}
\description{Attach partial derivative and residual variance functions to the nonlinear fitted object}
\usage{attachVarPrep(object, form = formula(object), fDerivFixef = NULL, 
    fDerivRanef = NULL, fVarResidual = NULL, fAddDummies = NULL)}
\arguments{
  \item{object}{the fitted nlme object}
  \item{form}{the formula used to fit the object, either formula or string used for automated derivation}
  \item{fDerivFixef}{\code{function(nfit,newdata,pred)} of derivatives in respect to fixed effects at newdata}
  \item{fDerivRanef}{\code{function(nfit,newdata,pred)} of derivatives in respect to random effects at newdata}
  \item{fVarResidual}{\code{function(nfit,newdata,pred)} to calculate var(residual) at newdata}
  \item{fAddDummies}{\code{function(newdata)} see "Handling categorial variables"}
}
\details{For usage with \code{\link{varPredictNlmeGnls}}, this function attaches to the fitted object \itemize{
\item derivative functions
\item residual variance function
\item Variance-Covariance methods 
}
\describe{ \item{Automatic derivation}{
If proper basic formula is given, \code{fDerivFixef} and \code{fDerivRanef} will be automatically derived from the model.
Up until now, only single level random effects models are supported for automatic derivation.
}}

\describe{ \item{Handling of categorial variables}{
item \code{fAddDummies=function(newdata)} of restult entry \code{varPrep} adds columns for dummy variables 
of categorial variables to newdata.
Default implementation supports only (and assumes) \code{contr.treatment} coding.
For other codings user must provide the function with argument \code{fAddDummies}. 
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
\item{fAddDummies}{function to add dummy columns for categorial variables to predictor data frame}
\item{fVarResidual}{function to calculate residual variance}
\item{fullFormula}{the extended formula as a string}}

\author{Thomas Wutzler}



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
