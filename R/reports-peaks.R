pack.admin.input.peaks <- function(peaksFile, callsFile, refFile, stutterFile, caseName='dummy',databaseFile=NULL, outputPath=getwd() ) {
	# Packs and verifies administrative information.
	# Documentation in man directory.
    	paths <- c(peaksFile, callsFile, refFile) 
	if(!is.null(databaseFile)) paths <- c(databaseFile, paths, recursive=TRUE)
	if(!is.null(stutterFile)) paths <- c(stutterFile, paths, recursive=TRUE)
	for(path in paths) {
		if(!file.exists(path))
			stop(paste(path, "does not exist."))
		else { 
			info <- file.info(path)
			if(info$isdir) stop(paste(path, "is not a file."))
      		}
    		} # loop over files.
	if(file.exists(outputPath) & !file.info(outputPath)$isdir) 
	stop(paste(outputPath, " exists and is not a directory."))
	admin <- list( caseName=caseName,
                databaseFile=databaseFile,
                peaksFile=peaksFile,
                callsFile=callsFile,
                refFile=refFile,
                stutterFile=stutterFile,
                outputPath=outputPath )
	return(admin)}
