\name{varPredictNlmeGnls}
\alias{varPredictNlmeGnls}
\title{varPredictNlmeGnls}
\description{Predictions including variance for basic nlme and gnls models.}
\usage{varPredictNlmeGnls(object, newdata, ...)}
\arguments{
  \item{object}{the model fit object used for predictions, treated by \code{\link{attachVarPrep}}}
  \item{newdata}{dataframe of new predictors and covariates}
  \item{\dots}{further arguments to \code{\link{predict.lme}} or \code{\link{predict.gls}}}
}
\details{Variance calculation is based on Taylor series expansion as described in appendix A1 by Wutzler08.

Fitted \code{object} needs to be prepared by function \code{\link{attachVarPrep}}.
If not done before, this function is called automatically within \code{varPredictNlmeGnls}. 
However, for finetuning or avoiding overhead in repeated calls, it
is recommended to explicitely call \code{\link{attachVarPrep}} before calling \code{varPredictNlmeGnls}.}
\value{numeric matrix with columns 
\item{fit}{predictions}
\item{varFix}{variance component due to uncertainty in fixed effects}
\item{varRan}{variance component due to uncertainty in random effects}
\item{varResid}{variance component due to residual error}
\item{sdPop}{standard deviation of prediction of a new population}
\item{sdInd}{standard deviation of prediction of a new individual}}
\references{Wutzler, T.; Wirth, C. & Schumacher, J. (2008) 
Generic biomass functions for Common beech (Fagus sylvatica L.) in Central Europe - predictions and components of uncertainty. 
Canadian Journal of Forest Research, 38, 1661-1675}
\author{Thomas Wutzler}



\seealso{\code{\link{varSumPredictNlmeGnls}}, \code{\link{twNlme-package}}}
\examples{
#----  fit a nlme and gnls model to data of stem weights 
data(Wutzler08BeechStem)

lmStart <- lm(log(stem) ~ log(dbh) + log(height), Wutzler08BeechStem )
nlmeFit <- nlme( stem~b0*dbh^b1*height^b2, data=Wutzler08BeechStem
			,fixed=list(b0 ~ si + age + alt, b1+b2 ~ 1)
			,random=  b0 ~ 1 | author
			,start=c( b0=c(as.numeric(exp(coef(lmStart)[1])),0,0,0), b1=as.numeric(coef(lmStart)[2]), b2=as.numeric(coef(lmStart)[3]) )
			,weights=varPower(form=~fitted(.))
			,method='REML'		# for unbiased error estimates
		)
summary(nlmeFit)
x3 <- update(nlmeFit, fixed=list(b0 ~ si * log(age), b1+b2 ~ 1))

gnlsFit <- gnls( stem~b0*dbh^b1*height^b2, data=Wutzler08BeechStem
	,params = list(b0 ~ si + age + alt, b1~1, b2 ~ 1)
	,start=c( b0=c(as.numeric(exp(coef(lmStart)[1])),0,0,0), b1=as.numeric(coef(lmStart)[2]), b2=as.numeric(coef(lmStart)[3]) )
	,weights=varPower(form=~fitted(.))
)
summary(gnlsFit)
fixef(gnlsFit)		# note the usage of fixef.gnls.
ranef(gnlsFit)		# note the usage of ranef.gnls.

#---- some artificial data for new prediction
nData <- data.frame( dbh=tmp <- seq(2,80,length.out=40), alt=median(Wutzler08BeechStem$alt), si=median(Wutzler08BeechStem$si) )
lmHeight <- lm( height ~ dbh, Wutzler08BeechStem)
nData$height <- predict(lmHeight, nData)
lmAge <- lm( age ~ dbh, Wutzler08BeechStem)
nData$age <- predict(lmAge, nData)

#---- do the prediction including variance calculation
# automatic derivation with accounting for residual variance model
nlmeFit <- attachVarPrep(nlmeFit, fVarResidual=varResidPower)
resNlme <- varPredictNlmeGnls(nlmeFit,nData)	

# plotting prediction and standard errors
plot( resNlme[,"fit"] ~ dbh, nData, type="l", xlim=c(40,80), lty="dashed")
lines( resNlme[,"fit"]+resNlme[,"sdPop"] ~ dbh, nData, col="maroon", lty="dashed" )
lines( resNlme[,"fit"]-resNlme[,"sdPop"] ~ dbh, nData, col="maroon", lty="dashed" )
lines( resNlme[,"fit"]+resNlme[,"sdInd"] ~ dbh, nData, col="orange", lty="dashed"  )
lines( resNlme[,"fit"]-resNlme[,"sdInd"] ~ dbh, nData, col="orange", lty="dashed" )

#---- handling special model of residual weights
# here we fit different power coefficients for authors
# for the prediction we take the mean, but because it appears in a nonlinear term
# we need a correction term
.tmp.f <- function(){	# takes long, so do not execute each test time 
	nlmeFitAuthor <- nlme( stem~b0*dbh^b1*height^b2, data=Wutzler08BeechStem
		,fixed=list(b0 ~ si + age + alt, b1+b2 ~ 1)
		,random=  b0 ~ 1 | author
		,start=c( b0=c(as.numeric(exp(coef(lmStart)[1])),0,0,0), b1=as.numeric(coef(lmStart)[2]), b2=as.numeric(coef(lmStart)[3]) )
		,weights=varPower(form=~fitted(.)|author)
	)
	#pred <- predict(modExampleStem, newdata, level=0)
	varResidPowerAuthor <- function(object,newdata,pred	){
		sigma <- object$sigma
		deltaAuthor <- coef(object$modelStruct$varStruct, allCoef = TRUE)
		delta <-  mean(deltaAuthor)
		varDelta <- var(deltaAuthor)
		sigma^2 * abs(pred)^(2*delta) * (1+2*log(abs(pred))^2*varDelta)
	}	
	#mtrace(varResidPowerAuthor)
	nlmeFitResidAuthor <- attachVarPrep(nlmeFitAuthor, fVarResidual=varResidPowerAuthor)
	resNlmeAuthor <- varPredictNlmeGnls(nlmeFitResidAuthor,nData)
	AIC(nlmeFitResidAuthor)
}

nlmeFit	#return for creating modExampleStem.RData
}
