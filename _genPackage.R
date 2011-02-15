#------------ code to help to maintain the package code
# usage: paste single code-lines to the RGui or RTerm
# workspace must be the root of the package, i.e. above R and man

pkg <- "twNlme" 

.tmp.loadPackages <- function(){
	# loading libraries, sourcing the code and loading data
	# usually done on startup library(MyPackage)
	
	# this code uses several packages that need to be installed once
	# install.packages(c("RUnit","inlinedocs","R.methodsS3", "abind","twMisc"), repos = c("http://R-Forge.R-project.org","http://cran.rakanu.com/"), dep = TRUE)
	# in case that you use not the current R-version, you need to download the sources tarball from https://r-forge.r-project.org/R/?group_id=887
	#   upack the sources, and issue from a shell from the folder above extracted folder twMisc "R CMD INSTALL twMisc"
	
	#library(snowfall)
	#sfInit(parallel=TRUE,cpus=4)		# for parallel execution on 4 processors see (library snowfall)
	#library(debug)
	
	library(twMisc) 	# for twStipFileExt
	library(nlme) 	
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
	data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))
}


.tmp.inlinedocs <- function(){
	# generate documentation
	
	# generate RD Files
	library(inlinedocs)
	html_viewer <- function(path) {
		browser <- getOption("browser")
		if(is.null(browser) && .Platform$OS.type == "windows")
			shell.exec(chartr("/", "\\", path))
		else browseURL(paste("file://", URLencode(path), sep=""))
	}
	
	unlink( file.path("man","*.Rd") )
	package.skeleton.dx(".")
	file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man", overwrite = TRUE )	# copy descriptions of data
	file.copy( Sys.glob(file.path("*.Rd")), "man", overwrite = TRUE )	# copy descriptions of data
	
	# generate the HTML  files
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --html ",pkg, sep="") )
	setwd(prevWd)
	
	# show in Browser
	htmlRoot <- file.path( system.file(package = pkg), "html" )
	#html_viewer(file.path(htmlRoot,"00Index.html"))
	html_viewer(file.path(htmlRoot,paste(pkg,"-package.html",sep="") ))
	
	# copy to the generated html into working directory
	#file.copy( htmlRoot, ".", recursive=TRUE)
}


.tmp.UnitTests <- function(){
	# executing Unit Tests before actual installation
	
	library(twMisc) 	# for twUTestF and copy2clip
	twUtestF()											# all unit tests		(displaying only summary output)
	#twUtestF(respTempLlyodAndTaylor94)					# single unit test
	#twUtestF(respTempLlyodAndTaylor94,"test.temp10")		# single test function  (displaying all the output)
	
	# let R check package consistency
	prevWd <- setwd("..")
	system(	cmd <- paste("R CMD check --no-manual ",pkg, sep="") )
	setwd(prevWd)
	copy2clip(cmd)	# for allowing easy pasting into a shell
	
}

# for adding low-level C or Fortran Code to the package ask Thomas (twutz)

.tmp.compile <- function(){
	# compile and test dll on windows
	# only required when providing C or Fortran level code
	# replace howlandInversion with the basename of your dll
	
	dynFilenameLocal <- file.path("src",paste("howlandInversion", .Platform$dynlib.ext, sep = ""))
	sfExport("dynFilenameLocal")
	sfClusterEval( dyn.unload(dynFilenameLocal) )
	dyn.unload(dynFilenameLocal)
	system("R CMD SHLIB -o src/howlandInversion.dll src/icbm1.c src/rc_helpers.c")
	dyn.load(dynFilenameLocal)
	sfExport("dynFilenameLocal")
	sfClusterEval( dyn.load(dynFilenameLocal) )
	
	#needs to be done on host linux machine (pc026):
	system("rm src/*.o")
	system("R CMD SHLIB -o src/howlandInversion.so  src/icbm1.c src/soilmod_fgsd.c src/rc_helpers.c src/mic_c_part.c src/dllinit.c")
}






