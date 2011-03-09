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
