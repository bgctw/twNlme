nlme:::getParsNlme <- function (plist, fmap, rmapRel, bmap, groups, beta, bvec, b, 
	level, N) 
{
	pars <- array(0, c(N, length(plist)), list(NULL, names(plist)))
	for (nm in names(plist)) {
		if (is.logical(f <- plist[[nm]]$fixed)) {
			if (f) {
				pars[, nm] <- beta[fmap[[nm]]]
			}
		}
		else {
			pars[, nm] <- f %*% beta[fmap[[nm]]]
		}
		if (level > 0) {
			Q <- length(groups)
			for (i in (Q - level + 1):Q) {
				b[[i]][] <- bvec[(bmap[i] + 1):bmap[i + 1]]
				if (is.logical(r <- plist[[nm]]$random[[i]])) {
					if (r) {
						pars[, nm] <- pars[, nm] + b[[i]][rmapRel[[i]][[nm]], 
							groups[[i]]]
					}
				}
				else {
					if (data.class(r) != "list") {
						pars[, nm] <- pars[, nm] + (r * t(b[[i]])[groups[[i]], 
									rmapRel[[i]][[nm]], drop = FALSE]) %*% 
							rep(1, ncol(r))
					}
					else {
						for (j in seq_along(rmapRel[[i]][[nm]])) {
							if (is.logical(rr <- r[[j]])) {
								pars[, nm] <- pars[, nm] + b[[i]][rmapRel[[i]][[nm]][[j]], 
									groups[[i]]]
							}
							else {
								pars[, nm] <- pars[, nm] + (rr * t(b[[i]])[groups[[i]], 
											rmapRel[[i]][[nm]][[j]], drop = FALSE]) %*% 
									rep(1, ncol(rr))
							}
						}
					}
				}
			}
		}
	}
	pars
}


nlme:::nlme.nlsList <- function (model, data = sys.frame(sys.parent()), fixed, random = fixed, 
	groups, start, correlation = NULL, weights = NULL, subset, 
	method = c("ML", "REML"), na.action = na.fail, naPattern, 
	control = list(), verbose = FALSE) 
{
	controlvals <- nlmeControl()
	controlvals[names(control)] <- control
	thisCall <- as.list(match.call())[-1]
	if (any(!is.na(match(names(thisCall), c("fixed", "data", 
					"start"))))) {
		warning(paste("nlme.nlsList will redefine \"fixed\"", 
				"\"data\", and \"start\""))
	}
	method <- match.arg(method)
	REML <- method == "REML"
	last.call <- as.list(attr(model, "call"))[-1]
	last.call$control <- NULL
	last.call$pool <- NULL
	thisCall[names(last.call)] <- last.call
	thisModel <- last.call[["model"]]
	thisCall[["model"]] <- eval(parse(text = paste(deparse(getResponseFormula(thisModel)[[2]]), 
				c_deparse(getCovariateFormula(thisModel)[[2]]), sep = "~")))
	cf <- na.omit(coef(model))
	start <- list(fixed = unlist(lapply(cf, median, na.rm = TRUE)))
	pnames <- names(start$fixed) <- names(cf)
	thisCall[["fixed"]] <- lapply(as.list(pnames), function(el) eval(parse(text = paste(el, 
						1, sep = "~"))))
	if (missing(random)) {
		random <- thisCall[["fixed"]]
	}
	reSt <- reStruct(random, data = NULL)
	if (missing(groups)) {
		thisCall[["groups"]] <- groups <- getGroupsFormula(model)
	}
	if (length(reSt) > 1 || length(groups[[2]]) > 1) {
		stop("Can only fit nlsList objects with single grouping variable")
	}
	ranForm <- formula(reSt)[[1]]
	if (!is.list(ranForm)) {
		ranForm <- list(ranForm)
	}
	mData <- thisCall[["data"]]
	if (is.null(mData)) {
		allV <- unique(unlist(lapply(ranForm, function(el) all.vars(el[[3]]))))
		if (length(allV) > 0) {
			alist <- lapply(as.list(allV), as.name)
			names(alist) <- allV
			alist <- c(as.list(as.name("data.frame")), alist)
			mode(alist) <- "call"
			mData <- eval(alist, sys.parent(1))
		}
	}
	else {
		if (mode(mData) == "name" || mode(mData) == "call") {
			mData <- eval(mData)
		}
	}
	reSt <- reStruct(random, REML = REML, data = mData)
	names(reSt) <- deparse(groups[[2]])
	rnames <- sapply(lapply(ranForm, "[[", 2), deparse)
	if (all(match(rnames, pnames, 0))) {
		madRes <- mad(resid(model), na.rm = TRUE)
		madRan <- unlist(lapply(cf, mad, na.rm = TRUE))
		madRan <- madRan[rnames]
		if (isInitialized(reSt)) {
			warning("Initial value for reStruct overwritten in nlme.nlsList")
		}
		matrix(reSt) <- diag((madRan/madRes)^2, ncol = length(rnames))
	}
	thisCall[["start"]] <- start
	thisCall[["random"]] <- reSt
	val <- do.call("nlme.formula", thisCall)
	val$origCall <- match.call()
	val
}




#nlme:::predict.gnls <-
predict.gnls <-
	function (object, newdata, na.action = na.fail, naPattern = NULL, 
		...) 
{
	if (missing(newdata)) {
		return(fitted(object))
	}
	newdata <- data.frame(newdata, check.names = FALSE)
	mCall <- object$call
	mfArgs <- list(formula = asOneFormula(formula(object), mCall$params, 
			naPattern, omit = c(names(object$plist), "pi", deparse(getResponseFormula(object)[[2]]))), 
		data = newdata, na.action = na.action)
	#mfArgs$drop.unused.levels <- TRUE
	mfArgs$drop.unused.levels <- FALSE
	dataMod <- do.call("model.frame", mfArgs)
	contr <- object$contrasts
	for (i in names(dataMod)) {
		if (inherits(dataMod[, i], "factor") && !is.null(contr[[i]]) && 
			is.matrix(contr[[i]])) {
			levs <- levels(dataMod[, i])
			levsC <- dimnames(contr[[i]])[[1]]
			if (any(wch <- is.na(match(levs, levsC)))) {
				stop(paste("Levels", paste(levs[wch], collapse = ","), 
						"not allowed for", i))
			}
			attr(dataMod[, i], "contrasts") <- contr[[i]][levs, 
				, drop = FALSE]
		}
	}
	N <- nrow(dataMod)
	if (is.null(naPattern)) 
		naPat <- rep(TRUE, N)
	else naPat <- as.logical(eval(asOneSidedFormula(naPattern)[[2]], 
				dataMod))
	plist <- object$plist
	pnames <- names(plist)
	if (is.null(params <- eval(object$call$params))) {
		params <- eval(parse(text = paste(paste(pnames, collapse = "+"), 
					"1", sep = "~")))
	}
	if (!is.list(params)) {
		params <- list(params)
	}
	val <- NULL
	for (i in seq_along(params)) {
		if (is.name(params[[i]][[2]])) {
			val <- c(val, list(params[[i]]))
		}
		else {
			val <- c(val, eval(parse(text = paste("list(", paste(paste(all.vars(params[[i]][[2]]), 
									deparse(params[[i]][[3]]), sep = "~"), collapse = ","), 
							")"))))
		}
	}
	params <- val
	names(params) <- pnames
	prs <- coef(object)
	pn <- names(prs)
	for (nm in pnames) {
		if (!is.logical(plist[[nm]])) {
			plist[[nm]] <- model.matrix(asOneSidedFormula(params[[nm]][[3]]), 
				model.frame(asOneSidedFormula(params[[nm]][[3]]), 
					dataMod))
		}
	}
	modForm <- getCovariateFormula(object)[[2]]
	val <- eval(modForm, data.frame(dataMod, getParsGnls(plist, 
				object$pmap, prs, N)))[naPat]
	names(val) <- row.names(newdata)
	lab <- "Predicted values"
	if (!is.null(aux <- attr(object, "units")$y)) {
		lab <- paste(lab, aux)
	}
	attr(val, "label") <- lab
	val
}


