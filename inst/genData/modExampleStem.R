.tmp.f <- function(){
	# in zz_archive\09\bef\singleTree
	# .RData
	modExampleStem <- dhcme.beech$stem
	save( modExampleStem, treeMass, file="modExampleStem.RData")
	#move to data directory of this package
}

rm(modExampleStem)
(modExampleStem <- attr(varPredictNlmeGnls,"ex")())

save( modExampleStem, file=file.path("data","modExampleStem.RData"))

