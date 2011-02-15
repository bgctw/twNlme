varPredictNlmeGnls <- function(
	### Predictions including variance for basic nlme and gnls models.
	object			##<< the model fit object used for predictions, treated by \code{\link{attachVarPrep}}
	, newdata		##<< dataframe of new predictors and covariates
	, ...			##<< further arguments to \code{\link{predict.lme}} or \code{\link{predict.gls}}
){
	##details<<
	## Variance calculation is based on Taylor series expansion as described in appendix A1 by Wutzler08.
	##
	## Fitted \code{object} needs to be prepared by function \code{\link{attachVarPrep}}.
	## If not done before, this function is called automatically within \code{varPredictNlmeGnls}. 
	## However, for finetuning or avoiding overhead in repeated calls, it
	## is recommended to explicitely call \code{\link{attachVarPrep}} before calling \code{varPredictNlmeGnls}.
	
	##references<<
	## Wutzler, T.; Wirth, C. & Schumacher, J. (2008) 
	## Generic biomass functions for Common beech (Fagus sylvatica L.) in Central Europe - predictions and components of uncertainty. 
	## Canadian Journal of Forest Research, 38, 1661-1675

	##seealso<< \code{\link{varSumPredictNlmeGnls}}, \code{\link{twNlme-package}}
	
	pred <- predict( object, newdata, level=0, ...)	# the population level predictions at newdata
	if( !inherits(object,"nlmeVarPrep") )
		object <- attachVarPrep(object)
	#tmpf <- object$varPrep$gradFix; mtrace(tmpf); tmpf(newdata)
	uNew <- object$varPrep$gradFix(object, newdata, pred)
	varFix <- varFixef(object)
	vcFix <- rep(0,nrow(newdata))
	for( i in seq(along=vcFix) ){
		vcFix[i] <- t(uNew[i,]) %*% varFix %*% uNew[i,]
	}
	vcRan <- if( inherits(object,"lme")){
		wNew <- object$varPrep$gradRan(object,newdata,pred)
		varRan <- varRanef(object)
		vcRan <- rep(0,nrow(newdata))
		for( i in seq(along=vcRan) ){
			vcRan[i] <- t(wNew[i,]) %*% varRan %*% wNew[i,]
		}
		vcRan
	}else rep(0,nrow(newdata))
	vcResid <- object$varPrep$fVarResidual(object,newdata,pred)	
	
	#varSum <- vcFix + vcRan + vcResid
	sdPop <- sqrt(vcFix+vcRan)
	sdInd <- sqrt(vcFix+vcRan+vcResid)
	res <- cbind( pred, vcFix, vcRan, vcResid, sdPop, sdInd )
	##value<< numeric matrix with columns 
	colnames(res) <- c(
		fit="fit"				##<< predictions
		,varFix="varFix"		##<< variance component due to uncertainty in fixed effects	
		,varRan="varRan"		##<< variance component due to uncertainty in random effects	
		,varResid="varResid"	##<< variance component due to residual error
		,sdPop="sdPop"			##<< standard deviation of prediction of a new population
		,sdInd="sdInd"			##<< standard deviation of prediction of a new individual
	) 
	##end<<
	res
}
attr(varPredictNlmeGnls,"ex") <- function(){
	#----  fit a nlme and gnls model to data of stem weights 
	data(Wutzler08BeechStem)
	
	lmStart <- lm(log(stem) ~ log(dbh) + log(height), Wutzler08BeechStem )
	nlmeFit <- nlme( stem~b0*dbh^b1*height^b2, data=Wutzler08BeechStem
				,fixed=list(b0 ~ si + age + alt, b1~1, b2 ~ 1)
				,random=  b0 ~ 1 | author
				,start=c( b0=c(as.numeric(exp(coef(lmStart)[1])),0,0,0), b1=as.numeric(coef(lmStart)[2]), b2=as.numeric(coef(lmStart)[3]) )
				,weights=varPower(form=~fitted(.))
				,method='REML'		# for unbiased error estimates
			)
	summary(nlmeFit)
	
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
			,fixed=list(b0 ~ si + age + alt, b1~1, b2 ~ 1)
			,random=  b0 ~ 1 | author
			,start=c( b0=c(as.numeric(exp(coef(lmStart)[1])),0,0,0), b1=as.numeric(coef(lmStart)[2]), b2=as.numeric(coef(lmStart)[3]) )
			,weights=varPower(form=~fitted(.)|author)
		)
		#pred <- predict(modExampleStem, newdata, level=0)
		varResidPowerAuthor <- function(object,newdata,pred	){
			sigma <- object$sigma
			deltaAuthor <- coef(object$modelStruct$varStruct)
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

.tmp.f <- function(){
	#mtrace(varPredictNlmeGnls)
	gnlsFit <- attachVarPrep(gnlsFit, fVarResidual=varResidPowerFitted)	
	resGnls <- varPredictNlmeGnls(gnlsFit,nData)
	lines( resGnls[,"fit"] ~ dbh, nData, type="l", col="gray")
	lines( resGnls[,"fit"]+resGnls[,"sdPop"] ~ dbh, nData, col="blue" )
	lines( resGnls[,"fit"]-resGnls[,"sdPop"] ~ dbh, nData, col="blue" )
	lines( resGnls[,"fit"]+resGnls[,"sdInd"] ~ dbh, nData, col="orange" )
	lines( resGnls[,"fit"]-resGnls[,"sdInd"] ~ dbh, nData, col="orange" )
	
	fit2 <- update(gnlsFit, weights=varPower(fixed=0.6))
	tmp<-varPredictNlmeGnls(fit2,nData)
	
	form <- y ~ beta0 + beta1 * exp(beta2*x)
	xDat <- cbind(beta0=1, beta1=2, beta2=-0.7, x=1:10)
	newx <- as.data.frame(xDat)
	yTrue <- with(as.data.frame(xDat), beta0 + beta1 * exp(beta2*x))
	y <- yTrue + rnorm(length(yTrue), sd=min(yTrue)/5)
	mf <- as.data.frame(cbind( x=xDat[,"x"], y=y))
	plot(y~x, mf)
	lines(yTrue~x, mf )
	#m <- model.frame( delete.response(terms(form)), as.data.frame(xDat) )
	fm1 <- gnls(form,mf
		,params=beta0+beta1+beta2~1
		,start = xDat[1,1:3]+rnorm(3,sd=1e-4)
	)
	fDerivFixef <- derivFixef(fm1)
	with( as.data.frame(xDat), eval(dFix) )$gradient
	varRanef(fm1)
	varFixef(fm1)
	seFixef(fm1)
	fm2 <- gls(distance ~ age, data = Orthodont)
	varFixef(fm2)
	varRanef(fm2)
}

#, sdType = c(	##<< specify which kind of standard deviation is calculated (set single value to save minor computation time) 
#	##describe<<
#	,population="population"	##<< sd for prediction of a new population
#	,individual="individual"	##<< sd for prediction of a new individual
#)##end<<

.tmp.f <- function(){
	
	data(modExampleStem)
	nlmeFit2 <- update(modExampleStem, weights=NULL)
	
	
	#plot(nlmeFit)	
	#plot( stem ~ dbh, Wutzler08BeechStem, col=author )
	#tmp<-order(Wutzler08BeechStem$dbh); lines(predict(nlmeFit, newdata=Wutzler08BeechStem[tmp,], level=0)~dbh[tmp],Wutzler08BeechStem)
	#plot( fitted(gnlsFit) ~ I(fitted(gnlsFit)+resid(gnlsFit)) )
	#abline(0,1,col="gray")
	#points( fitted(nlmeFit) ~ I(fitted(nlmeFit)+resid(nlmeFit)), col="maroon" )
	
	data(modExampleStem)
	nfit <- attachVarPrep( modExampleStem, form = "b0*dbh^b1*height^b2")
	#mtrace(varPredictNlmeGls)
	varPredictNlmeGnls(nfit, newdata=data.frame(dbh=18.8, height=16.9, age=40, si=30, alt=470))
	
	data(Wutzler08BeechStem)
	varPredictNlmeGnls(nfit, newdata=head(Wutzler08BeechStem))
	
	
}

.tmp.f <- function(){
	data(modExampleStem)
	nfit <- modExampleStem
	form <- "b0*d^b1*h^b2"
}

varResidPower <- function(
	### Variance of residual of prediction for Power Variance without random effects e.g. \code{weights=varPower(form=~fitted(.))}
	object	##<< the fitted object
	,newdata	##<< the data frame with new predictors
	,pred	##<< the prediction at population level
){
	sigma <- object$sigma
	delta <- coef(object$modelStruct$varStruct)["power"]
	sigma^2 * abs(pred)^(2*delta)
	### numeric vector of var(eps_i)
	##seealso<< \code{\link{twNlme-package}}
}

#------------------- constructing the full formula 
.extractFixedList <- function( 
	### extracts the fixed effects in list form "list(b0=b0 ~ si + age + alt, b1=b1 ~ 1, b2=b2 ~ 1)"
	x					##<< the fitted nlme model
	,makeNames=TRUE 	##<< whether to attach an explicit name attribute
){
	if( inherits(x,"nlme") ){
		fixF <- x$call$fixed
	}else{
		if( !is.null(x$call$params) ){
			fixF <- x$call$params
		}else{
			fixF <- lapply( as.list(names(coefficients(x))), function(el){ eval(parse(text=paste("formula(",el,"~1)"))) } )
		}
	}
	if (inherits(fixF, "formula")){
		fixF = list( fixF )
	}
	if ( is.call(fixF) ) {
		fixF = eval(fixF)
	}
	#print("pred_beech_nl.formula.fixed:"); print(fixF);
	if( makeNames ){
		names( fixF ) <- as.character(lapply(fixF, function(el){el[[2]]} ))
	}
	#deparse(lapply(fixF, function(el) as.name(deparse(el))))
	fixF
	# list of formulas for each coeficient
}
attr(.extractFixedList,"ex") <- function(){
	data(modExampleStem)
	(tmp <- .extractFixedList(modExampleStem))
}


.extractRandomList <- function(
	### extracts the random effects in the form "list(b0=b0 ~ 1, b1=b1 ~ 1)"
	x
	, makeNames=TRUE
){
	if( inherits(x,"nlme")){
		ranF <- unlist(attributes(x$modelStruct$reStruct[[1]])$formula)
		if( makeNames ){
			names(ranF) <- as.character(lapply(ranF, function(elr){elr[[2]]} ))
		}
		#deparse(ranF, control=c("keepInteger") ) #no option showAttributes,"warnIncomplete" 
		ranF 
	}else{
		list()
	}
}
attr(.extractRandomList,"ex") <- function(){
	data(modExampleStem)
	(tmp <- .extractRandomList(modExampleStem))
}

expandLinFormula <- function( 
	### Extends the linear formula to an expression involving coefficients
	linForm		##<< formula for fixed coefficients depending on linear term see \code{\link{.extractFixedList}}
	, suffix="" ##<< base parameter name
	, varNames=NULL	##<< variable names to use
){
	##details<<
	## parameter names will be lhs of the formula + suffix + .i wiht .0 for the intercept
	## if there is no intercept, starting from .1
	value = NULL
	coef2 <- c("")[ FALSE ]
	parBaseName <- linForm[[2]]
	terms2 <- attr(terms( linForm ),"term.labels")
	if( length(terms2) > 0 ){
		if( 0 == length(attr(terms( linForm ),"intercept")) ){
			coef2 <- if( 0<length(varNames) ) varNames else paste( parBaseName, suffix, ".", seq(along=terms2), sep="" )
			value = paste( coef2, terms2, sep="*", collapse=" + " )
		}else{
			coef2 <- if( 0<length(varNames) ) varNames[-1] else paste( parBaseName, suffix, ".", seq(along=terms2), sep="" )
			value = paste( coef2, terms2, sep="*", collapse=" + " )
			coefb = if( 0<length(varNames) ) varNames[1] else paste( parBaseName,suffix,".0", sep="")
			coef2 = c( coefb, coef2 )
			value = paste( coefb, " + ", value, sep="")
		}
	}else{
		if( attr(terms( linForm ),"intercept")){
			value <- coef2 <- if( 0<length(varNames) ) varNames else paste( parBaseName, suffix, ".0", sep="" )
			#value = coef2 = paste( parBaseName,suffix,".0", value, sep="")
		}
	}
	# if names are given 
	
	attr(value,"coef") <- coef2
	attr(value,"parBaseName") <- parBaseName
	value
	### String of formula with coefficients expanded to linear combinations of covariates and random effects.
	##seealso<< \code{\link{twNlme-package}}
	}
attr(expandLinFormula,"ex") <- function(){
	data(modExampleStem)
	fixForm <- twNlme:::.extractFixedList(modExampleStem)
	expandLinFormula(fixForm[[1]])
}
#ex: fixF( string x$call$fixed, "b0" ) 
#results in 'b0.0 + b0.1*si + b0r.1'
#attributes coef and resp are attached listing all the fixed and random coefficients

#getTerms( fixF, 4 )

.covarMap <- function( 
	### extract the expanded linear formula with covariates and random effect for each parameter
	nfit = NULL		##<< the fitted nlme object 
	,fixedMap = .extractFixedList(nfit) 
	,randomMap = .extractRandomList(nfit)
){
	pMap <- list()
	pMapfc <- NULL #c("")[FALSE]
	pMaprc <- NULL #c("")[FALSE]
	#name list elements by RHS of the formulas
	#names(fixedMap) <- sapply( fixedMap, function(x) x[[2]] )
	#names(randomMap) <- sapply( randomMap, function(x) x[[2]] )
	params <- unique( c(names(fixedMap),names(randomMap) ))
	fixefNames <- names(fixef(nfit))
	ranefNames <- names(ranef(nfit))
	for( i in seq(along=params) ){
		val <- NULL
		valc <- NULL
		if( !is.null(fixedMap[[ params[i] ]]) ){
#			tmpNamesReplace <- c(
#				fixefNames[ grep( paste("^",params[i],"\\.",sep=""), fixefNames )]
#				,fixefNames[ grep( paste("^",params[i],"$",sep=""), fixefNames )]
#			)
#			tmp <- expandLinFormula(fixedMap[[i]], varNames=tmpNamesReplace)
			# deriv cant digest b0.(Intercept), so we have to create new names
			tmp <- expandLinFormula(fixedMap[[i]])
			val = c( val, tmp )
			pMapfc <- c( pMapfc, attr(tmp,"coef") ) 
		}
		if( !is.null(randomMap[[ params[i] ]]) ){
			#names are repeated from fixed effects, hence construct new names
			#tmpNamesReplace <- ranefNames[ grep( paste("^",params[i],".",sep=""), ranefNames )]
			tmp <- expandLinFormula(randomMap[[i]],suffix="r")
			val = c( val, tmp )
			pMaprc <- c( pMaprc, attr(tmp,"coef") ) 
		}
		if( is.null(val) ) val = params[i]
		val = paste( val, collapse="+" )
		pMap <- c( pMap, x = val	)	
		names(pMap)[ length(pMap) ] = as.character(params[i])
	}
	attr(pMap,"fixCoef") <- pMapfc
	attr(pMap,"ranCoef") <- pMaprc
	pMap
	### list with each entry representing a string of an expanded formula
	### additional all the coefficients are provided in attributes fixCoef and ranCoef
}
attr(.covarMap,"ex") <- function(){
	data(modExampleStem)
	#mtrace(.covarMap)
	(tmp <- .covarMap(modExampleStem))
}

attachVarPrep <- function(
	### Attach partial derivative and residual variance functions to the nonlinear fitted object
	object		##<< the fitted nlme object
	,form = object$call$model  	##<< the formula used to fit the object, either formula or string used for automated derivation 
	,fDerivFixef=NULL			##<< \code{function(nfit,newdata,pred)} of derivatives in respect to fixed effects at newdata 
	,fDerivRanef=NULL			##<< \code{function(nfit,newdata,pred)} of derivatives in respect to random effects at newdata 
	,fVarResidual=NULL			##<< \code{function(nfit,newdata,pred)} to calculate var(residual) at newdata
){
	##details<<
	## For usage with \code{\link{varPredictNlmeGnls}}, this function attaches \itemize{
	## \item derivative functions
	## \item residual variance function
	## \item Variance-Covariance methods to fitted object
	## }
	## \describe{ \item{Automatic derivation}{
	## If proper basic formula is given, \code{fDerivFixef} and \code{fDerivRanef} will be automatically derived from the model.
	## In this case fixed and random effects have to be specified as list of formulas for each coefficient instead of single formula 
	## e.g. \code{list(b0~1,b1~1)} instead of \code{list(b0+b1~1)}. 
	## Up until now, only single level random effects models are supported for automatic derivation.
	## }}
	
	##seealso<< \code{\link{twNlme-package}}
	

	if( is.null(fDerivFixef) | (0 < length(ranef(object) & is.null(fDerivRanef) )) ){
		formStr <- if( inherits(form,c("formula","call","language")) ) deparse(form[[length(form)]]) else form
		#formStr = "b0*dbh^b1*height^b2"
		#construct a mapping of each fixed/random parameter in original model formula onto its linear combination
		#formula.fixed and formula.random in pred_beech_nl.r
		tmp <- .covarMap( object )
		tmpf1str <- paste("~",formStr," ") #enclose by empty spaces
		#construct the full formula
		#replace parameter names by their extended linear combination of covariates
		#parametername anclosed by non-name chars
		#gsub("([^0-9A-Za-z_.])(b0)([^0-9A-Za-z_.])","\\1XXX\\3","b0.01~b0+b01")
		for( param in names(tmp) ){
			pat <- paste("([^0-9A-Za-z_.])(",param,")([^0-9A-Za-z_.])", sep="" )
			repl <- paste("\\1(",tmp[[param]],")\\3",sep="")
			tmpf1str <- gsub( pat, repl, tmpf1str )
		}
		tmpf2str <- gsub("I\\(","\\(",gsub(":", "*", tmpf1str))
		tmpf <- as.formula(tmpf2str)
		
		fixCoefNames <- attr(tmp,"fixCoef")
		if( is.null(fixCoefNames) ) fixCoefNames <- names(fixef(object))
		ranCoefNames <- attr(tmp,"ranCoef")
		if( is.null(ranCoefNames) ) ranCoefNames <- names(ranef(object))
		
		gradRanExp <- if( 0 < length(ranCoefNames) ){
				deriv( tmpf, ranCoefNames )
			}else{
				expression(numeric(0))
			}
	}else{
		if( is.null(fDerivRanef)) fDerivRanef <- expression(numeric(0))
	}
	
	coefGradFixed <- {tmp <- fixef(object); names(tmp) <- fixCoefNames; as.list(tmp) }
	coefGradRan <- as.list(structure( rep(0,length(ranCoefNames)), names=ranCoefNames))
	coefGrad <- c(coefGradFixed,coefGradRan)
	gradFix <- if( !is.null(fDerivFixef) ) fDerivFixef else {
		 gradFixExp <- deriv(tmpf, fixCoefNames )
		 function(object,newdata,pred){ attr( with( coefGrad, with(newdata, eval(gradFixExp) )),"gradient") }
	}
	gradRan <- if( !is.null(fDerivRanef) ) fDerivRanef else {
		function(object,newdata,pred){ attr( with( coefGrad, with(newdata, eval(gradRanExp) )),"gradient")}
	} 

	##details<< \describe{\item{Variance of Residuals}{
	## Providing no argument \code{fResidual} assumes iid residuals, i.e. \code{weights=NULL}. 
	## For other residual variance models. See e.g. \code{\link{varResidPower}} corresponding to \code{weights=varPower(form=~fitted(.))}
	##}}
	if( 0==length(fVarResidual) ) fVarResidual <-	function(object,...) object$sigma^2

	if( !inherits(object,"nlmeVarPrep"))	class(object) <- c("nlmeVarPrep", class(object))
	
	#calulate vector of derivative functions
	##value<< nfit with additional entry \code{varPrep}, which is a list of 
	object$varPrep <- list(
		varFix	= varFixef(object)	##<< variance-covariance matrix of fixed effects
		,varRan	= varRanef(object)	##<< variance-covariance matrix of random effects
		,coefFix = fixCoefNames		##<< names of the fixed coefficients in gradiant function
		,coefRan = ranCoefNames		##<< names of the random coefficients in gradiant function
		,gradFix = gradFix			##<< derivative function for fixed effects
		,gradRan = gradRan			##<< derivative function for random effects
		,fVarResidual = fVarResidual	##<< function to calculate residual variance
		)
	##end<<
	object		
	}
attr(attachVarPrep,"ex") <- function(){
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

#,formStr = formStr
#,fullFormula = tmpf


varSumPredictNlmeGnls <- function( 
	### Variance of the sum of predictions taking care of covariances between single predictions.
	object					##<< the model fit object used for predictions, treated by \code{\link{attachVarPrep}}
	, newdata				##<< dataframe of new predictors and covariates
	, pred=FALSE			##<< if TRUE, the predicted value (sum of predictions) is  returned in attribute pred
	, retComponents=FALSE	##<< if TRUE, the sum of the error components (fixed, random, noise) are returned in attributes "varFix","varRan","varResid" 
){
	##details<<
	## Variance calculation is based on Taylor series expansion as described in appendix A2 by Wutzler08.
	## 
	## Performance of this function scales with n^2. So do not apply for too many records.
	
	##references<<
	## Wutzler, T.; Wirth, C. & Schumacher, J. (2008) 
	## Generic biomass functions for Common beech (Fagus sylvatica L.) in Central Europe - predictions and components of uncertainty. 
	## Canadian Journal of Forest Research, 38, 1661-1675
	
	#reps specifies a vector of replicated observations (e.g. three identical rows are just represented by one row with reps=3)
	#reps are not working yet rows are hidden
	reps=rep(1,nrow(newdata))
	if( length(reps) < nrow(newdata) ){
		print("length of vector of replications does not match length of dataframe")
		return(NA)
	}
	if( !inherits(object,"nlmeVarPrep") )
		object <- attachVarPrep(object)
	#calculates Variance at position in new dataframe
	#pfit must be a nlme or gls fit, that has been extended by prepVarFunc
	n <- nrow(newdata)
	resFix <- resRan <- 0
	tmpFixVar <- object$varPrep$varFix
	uNew <- object$varPrep$gradFix(object, newdata, pred)
	ynew <- predict( object, newdata, level=0 )
	if( inherits(object, "nlme") ){
		tmpRanVar <- object$varPrep$varRan
		wNew <- object$varPrep$gradRan(object,newdata,pred)
		if( n > 1 ){
			for( i in 1:(n-1) ){
				#cat(i,", ")
				#+covar(i,i)
				resFix <- resFix + ( t(uNew[i,]) %*% tmpFixVar %*% uNew[i,] )*reps[i]
				resRan <- resRan + ( t(wNew[i,]) %*% tmpRanVar %*% wNew[i,] )*reps[i]
				#+covar(i,j)+covar(j,i)=2*covar(i,j)
				for( j in (i+1):n ){
					resFix <- resFix + ( t(uNew[i,]) %*% tmpFixVar %*% uNew[j,] )*reps[i]*reps[j]
					resRan <- resRan + ( t(wNew[i,]) %*% tmpRanVar %*% wNew[j,] )*reps[i]*reps[j]
				}#for j
			}#for i
		}#n>1
		i<-n
		#covar(n,n)
		resFix <- resFix + ( t(uNew[i,]) %*% tmpFixVar %*% uNew[i,] )*reps[i]
		resRan <- resRan + ( t(wNew[i,]) %*% tmpRanVar %*% wNew[i,] )*reps[i]
	}else{
		if( n > 1 ){
			for( i in 1:(n-1) ){
				#cat(i,", ")
				#+covar(i,i)
				resFix <- resFix + ( t(uNew[i,]) %*% tmpFixVar %*% uNew[i,] )*reps[i]
				#resRan <- resRan + ( t(tmpRanGrad[i,]) %*% tmpRanVar %*% tmpRanGrad[i,] )*reps[i]
				#+covar(i,j)+covar(j,i)=2*covar(i,j)
				for( j in (i+1):n ){
					resFix <- resFix + 2*( t(uNew[i,]) %*% tmpFixVar %*% uNew[j,] )*reps[i]*reps[j]
					#resRan <- resRan + 2*( t(tmpRanGrad[i,]) %*% tmpRanVar %*% tmpRanGrad[j,] )*reps[i]*reps[j]
				}#for j
			}#for i
		}#n>1
		i<-n
		#covar(n,n)
		resFix <- resFix + ( t(uNew[i,]) %*% tmpFixVar %*% uNew[i,] )*reps[i]
		#		resRan <- resRan + ( t(tmpRanGrad[i,]) %*% tmpRanVar %*% tmpRanGrad[i,] )*reps[i]
	}
	resResid <- sum( object$varPrep$fVarResidual(object,newdata,ynew)*reps )
	##value<< named vector
	res <- c(	
		pred = sum(ynew*reps)	##<< sum of predictions  
		,sdPred = sqrt(resFix + resRan + resResid)	##<< standard deviation of pred
		,varFix = resFix		##<< variance component due to uncertainty in fixed effects	
		,varRan = resRan		##<< variance component due to random effects
		,varResid = resResid	##<< variance component due to residual variance
		)
	##end<<
}
attr(varSumPredictNlmeGnls,"ex") <- function(){
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



