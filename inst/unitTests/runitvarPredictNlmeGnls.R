.setUp <-function () {
	data(Wutzler08BeechStem)
	data(modExampleStem)
	vars <- list(nData = head(Wutzler08BeechStem))
	attach(vars)
}

.tearDown <- function () {
	detach()
}

test..gsubFormTerm <- function(){
	formTerms <- c("I(si + alt):log(age)","age +I(si + alt)","SI(si + alt):log(age)") 
	#gsub("\\bI\\(","\\1\\(",formTerms)
	res <- try( .gsubFormTerm(formTerms) )
	if( inherits(res,"try-error"))
		res <- twNlme:::.gsubFormTerm(formTerms)
	checkEquals( c("(si+alt)*log(age)","age+(si+alt)","SI(si+alt)*log(age)"), res)
}

test.basemodel <- function(){
	nFit <- modExampleStem
	pnFit <- attachVarPrep(nFit, fVarResidual=varResidPower)
	resNlme <- varPredictNlmeGnls(pnFit,nData)	
	checkTrue( all(resNlme > 0))
}

test.IinFormula <- function(){
	nFit <- update(modExampleStem, fixed=list(b0 ~ I(si+alt) * log(age), b1~1, b2 ~ 1), start=fixef(modExampleStem))
	#str(nFit$plist)
	comp <- compb <- nFit$plist[[1]]
	try({
			(termVec <- .extractFixedList(nFit, excludeIntercept=TRUE)[[1]])
			(formRes <- .constructLinFormulaString(termVec, "b0")) 
			grad <- deriv(parse(text=formRes$form), formRes$coefNames)
			tmp <- with( tail(Wutzler08BeechStem), with( as.list({tmp<-fixef(nFit);names(tmp)<-formRes$coefNames;tmp}), eval( grad) ))
		})
	pnFit <- attachVarPrep(nFit, fVarResidual=varResidPower)
	resNlme <- varPredictNlmeGnls(pnFit,nData)	
	checkTrue( all(resNlme > 0))
}

test.sumFormula <- function(){
	nFit <- update(modExampleStem, fixed=list(b0+b1+b2 ~ 1), start=fixef(modExampleStem)[c(1,5,6)])
	#str(nFit$plist)
	try({
		(termVec <- .extractFixedList(nFit, excludeIntercept=FALSE)[[1]])
		(formRes <- .constructLinFormulaString(termVec, "b0"))
		grad <- deriv(parse(text=formRes$form), formRes$coefNames)
	})
	pnFit <- attachVarPrep(nFit, fVarResidual=varResidPower)
	#str(pnFit$varPrep)
	resNlme <- varPredictNlmeGnls(pnFit,nData)	
	checkTrue( all(resNlme > 0))
}

test.categorial <- function(){
	nFit <- update(modExampleStem, fixed=list(b0 ~ author+age, b1~1, b2 ~ 1), random=b1~1|author
		, start=c( rep(as.numeric(fixef(modExampleStem)[1]),length(unique(Wutzler08BeechStem$author))),0, as.numeric(fixef(modExampleStem)[5]), as.numeric(fixef(modExampleStem)[6]) )
	)
	#str(nFit$plist)
	#mtrace(.extractFixedList)
	fixedList <- try( .extractFixedList(nFit) )
	if( inherits(fixedList,"try-error"))
		fixedList <- twNlme:::.extractFixedList(nFit)
	catMapAuthor <- attr(fixedList,"catMap")$author
	checkEquals( as.character(unique(Wutzler08BeechStem$author)[-1]),  catMapAuthor)
	authorVars <- try( .gsubFormTerm(catMapAuthor) )
	if( inherits(authorVars,"try-error"))
		authorVars <- twNlme:::.gsubFormTerm(catMapAuthor)
	checkTrue( all(paste("author",authorVars,sep="") %in% fixedList$b0))

	pnFit <- attachVarPrep(nFit, fVarResidual=varResidPower)
	#str(pnFit$varPrep)
	#tmp.f <- pnFit$varPrep$gradFix; mtrace(tmp.f); pnFit$varPrep$gradFix <- tmp.f
	#mtrace(varPredictNlmeGnls)
	nData$author <- factor("Cienciala", levels=levels(Wutzler08BeechStem$author) )
	resNlme <- varPredictNlmeGnls(pnFit,nData)	
	checkTrue( all(resNlme > 0))
	
}

t_est.formula <- function(){
	# test wheter formula can be retrieved, even if call includes a variable
	# works when executed from workspace
	# fails on RCheck 
	data(Wutzler08BeechStem)
	
	lmStart <- lm(log(stem) ~ log(dbh) + log(height), Wutzler08BeechStem )
	tmpForm <- stem~b0*dbh^b1*height^b2
	nlmeFit <- nlme( tmpForm, data=Wutzler08BeechStem
		,fixed=list(b0 ~ si + age + alt, b1+b2 ~ 1)
		,random=  b0 ~ 1 | author
		,start=c( b0=c(as.numeric(exp(coef(lmStart)[1])),0,0,0), b1=as.numeric(coef(lmStart)[2]), b2=as.numeric(coef(lmStart)[3]) )
		,weights=varPower(form=~fitted(.))
		,method='REML'		# for unbiased error estimates
	)
	formula(nlmeFit)
	summary(nlmeFit)
	
#---- some artificial data for new prediction
	nData <- data.frame( dbh=tmp <- seq(2,80,length.out=40), alt=median(Wutzler08BeechStem$alt), si=median(Wutzler08BeechStem$si) )
	lmHeight <- lm( height ~ dbh, Wutzler08BeechStem)
	nData$height <- predict(lmHeight, nData)
	lmAge <- lm( age ~ dbh, Wutzler08BeechStem)
	nData$age <- predict(lmAge, nData)
	
#---- do the prediction including variance calculation
# automatic derivation with accounting for residual variance model
	#mtrace(attachVarPrep)
	nlmeFit <- attachVarPrep(nlmeFit, fVarResidual=varResidPower)
	resNlme <- varPredictNlmeGnls(nlmeFit,nData)	
	checkTrue( all(resNlme > 0))
	
}







