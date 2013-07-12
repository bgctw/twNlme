R.methodsS3::setMethodS3("varRanef","lme", function( 
 	### Variance-Covariance of lme random effects (Psi)
	object			##<< the fitted nlme or lme model
	,varName=NULL	##<< for which parameter
	,...
){
	# varRanef.lme
	sigma <- object$sigma
	if( is.null(varName) ){
		reStr <- object$modelStruct$reStruct[[1]]
	}else{
		reStr <- object$modelStruct$reStruct[[varName]]
	}
	#corrm <- nlme::corMatrix(reStr)
	corrm <- corMatrix(reStr)
	stdDev <- sigma * attr(corrm, "stdDev")
	value <- t(corrm*stdDev)*stdDev
	#names are lost in corMatrix, so reconstruct them with array function
	value <- array( value , dim(corrm), attr(reStr,"Dimnames") )
	### named numeric matrix
	##seealso<< \code{\link{twNlme-package}}
})
attr(varRanef.lme,"ex") <- function(){
	fm1 <- lme(distance ~ age, data = Orthodont, random = ~age)
	varRanef(fm1)
	varFixef(fm1)
	fm2 <- gls(distance ~ age, data = Orthodont)
	varFixef(fm2)
	varRanef(fm2)
}

R.methodsS3::setMethodS3("varRanef","gls", function(
		### Empty Variance-Covariance matrix for random effects
		object			##<< the fitted gls or gnls model
		,varName=NULL	##<< for which parameter
		,...
	){
		matrix( numeric(0), nrow=0, ncol=0 )
		### 0 object 0 matrix
		##seealso<< \code{\link{twNlme-package}}
	})


R.methodsS3::setMethodS3("varFixef","lme", function( 
		### Variance-Covariance of lme fixed effects
		object	##<< the fitted model
		,...
	){
		object$varFix
		### named numeric matrix
		##seealso<< \code{\link{twNlme-package}}
	})
R.methodsS3::setMethodS3("varFixef","gls", function( 
		### Variance-Covariance of gls effects
		object	##<< the fitted model
		,...
	){
		object$varBeta
		### named numeric matrix
		##seealso<< \code{\link{twNlme-package}}
	})

R.methodsS3::setMethodS3("fixef","gls", function( 
		### Fixed effects, i.e. coefficients of gls model
		object	##<< the fitted model
		,...
	){
		coef(object)
		### named numeric vector
		##seealso<< \code{\link{twNlme-package}}
	})
R.methodsS3::setMethodS3("ranef","gls", function( 
		### Random effects, i.e. none, of gls model
		object	##<< the fitted model
		,...
	){
		numeric(0)
		### numeric(0)
		##seealso<< \code{\link{twNlme-package}}
	})
R.methodsS3::setMethodS3("seFixef","lme", function( 
		### Standard error of fixed effects 
		object	##<< the fitted model
		,...
	){
		summary(object)$tTable[,"Std.Error"]
		### named numeric vector
		##seealso<< \code{\link{twNlme-package}}
	})
R.methodsS3::setMethodS3("seFixef","gls", function( 
		### Standard error of coefficients 
		object	##<< the fitted model
		,...
	){
		summary(object)$tTable[,"Std.Error"]
		### named numeric vector
		##seealso<< \code{\link{twNlme-package}}
	})





