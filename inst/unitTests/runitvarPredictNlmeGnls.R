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
	resNlme <- varPredictNlmeGnls(pnFit,nData)	
	checkTrue( all(resNlme > 0))
	
}







