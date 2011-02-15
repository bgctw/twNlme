\name{varSumPredictNlmeGnls}
\alias{varSumPredictNlmeGnls}
\title{varSumPredictNlmeGnls}
\description{Variance of the sum of predictions taking care of covariances between single predictions.}
\usage{
varSumPredictNlmeGnls(object, newdata, pred = FALSE, retComponents = FALSE)
}
\arguments{
  \item{object}{the model fit object used for predictions, treated by \code{\link{attachVarPrep}}}
  \item{newdata}{dataframe of new predictors and covariates}
  \item{pred}{if TRUE, the predicted value (sum of predictions) is  returned in attribute pred}
  \item{retComponents}{if TRUE, the sum of the error components (fixed, random, noise) are returned in attributes "varFix","varRan","varResid"}
}
\details{Variance calculation is based on Taylor series expansion as described in appendix A2 by Wutzler08.

Performance of this function scales with n^2. So do not apply for too many records.}
\value{named vector
\item{pred}{sum of predictions}
\item{sdPred}{standard deviation of pred}
\item{varFix}{variance component due to uncertainty in fixed effects}
\item{varRan}{variance component due to random effects}
\item{varResid}{variance component due to residual variance}}
\references{Wutzler, T.; Wirth, C. & Schumacher, J. (2008) 
Generic biomass functions for Common beech (Fagus sylvatica L.) in Central Europe - predictions and components of uncertainty. 
Canadian Journal of Forest Research, 38, 1661-1675}
\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(modExampleStem)	# load the model, which has already been prepared for prediction
#-- prediction on with varying number of records
data(Wutzler08BeechStem)
(resNlme <- varSumPredictNlmeGnls(modExampleStem, head( Wutzler08BeechStem, n=10 )))	
(resNlme2 <- varSumPredictNlmeGnls(modExampleStem, head( Wutzler08BeechStem, n=180 )))
# plotting relative error components
barplot(c(sqrt(resNlme[-(1:2)])/resNlme[1], sqrt(resNlme2[-(1:2)])/resNlme2[1]) )
# note how the residual error declines with record number, 
# while the fixed and random error does does not decline
}
