#pack.admin.input
#read.csp.profile
#prosecution.hypotheses
#defence.hypothesis
#optimisation.params
#allele.report
#output.report

# function to subset data
subsetData = function(data,search)
	{
	# subset columns by search parameter
	dataout = data[,grep(search,colnames(data))]
	# add first column (for replicate names etc)
	# can add this back in if necessary
	#dataout = cbind(data[,1],dataout)
	#colnames(dataout)[1] = "info"
	# set rownames to marker names
	rownames(dataout) = data$Marker
	# remove columns that have no data
	index = which(apply(dataout,MARGIN=2,FUN=function(x) all(is.na(x))))
	if(length(index!=0)) dataout = dataout[,-index,drop=FALSE]
	# remove AMEL data
	index = which(rownames(dataout)%in%c("AMEL","amel","Amel"))
	if(length(index)!=0) dataout = dataout[-index,,drop=FALSE]
	return(dataout)
	}

# function to check that the data make sense
checkData1 = function(alleles,heights,sizes)
	{
	# check if dimensions are all the same
	check1 = identical(dim(alleles),dim(heights),dim(sizes))
	# check if all have same loci
	check2 = all(mapply(FUN=function(a,b,c) a==b&b==c,rownames(alleles),rownames(heights),rownames(sizes)))
	return(all(check1,check2))
	}

# function to check data makes sense between replicates
checkData2 = function(alleleReps,heightReps=NULL,sizeReps=NULL)
	{
	# loci names for each replicate
	lociReps = sapply(alleleReps,FUN=function(x) rownames(x),simplify=FALSE)
	foo = function(a,b) all(a==b)
	check1 = all(sapply(lociReps,FUN=function(x) identical(lociReps[[1]],x)))
	return(all(check1))
	}




# function to read in peak height data
read.peaks.profile = function(FILE)
	{
	if(!file.exists(FILE)) stop(paste(FILE, "does not exist."))
	# get data
	data = read.csv(file=FILE,as.is=TRUE,na.strings=c("NA",""))
	# get replicate names
	conditions = names(table(data[,1]))
	# get allele data
	alleles = sapply(conditions, FUN=function(x) likeLTD:::subsetData(data[which(data[,1]==x),],"Allele"),simplify=FALSE)
	# get height data
	heights = sapply(conditions, FUN=function(x) likeLTD:::subsetData(data[which(data[,1]==x),],"Height"),simplify=FALSE)
	# get size data
	# not sure if this is needed
	sizes = sapply(conditions, FUN=function(x) likeLTD:::subsetData(data[which(data[,1]==x),],"Size"),simplify=FALSE)
	# checks
	check1 = !mapply(FUN=checkData1,alleles,heights,sizes)
	if(any(check1)) stop(paste("Error in the CSP provided - replicates ",paste(which(check1),collapse=",")))
	check2 = !checkData2(alleles,heights,sizes)
	if(check2) stop(paste("Error in the CSP provided - replicates do not match"))
	# return data
	return(list(alleles=alleles,heights=heights,sizes=sizes))
	}

# function to read in allele calls
# 0=non-allelic, 1=allelic
read.allelic.calls = function(FILE)
	{
	if(!file.exists(FILE)) stop(paste(FILE, "does not exist."))
	# get data
	data = read.csv(file=FILE,as.is=TRUE,na.strings="")
	# get replicate names
	conditions = names(table(data[,1]))
	# get allele data
	calls = sapply(conditions, FUN=function(x) subsetData(data[which(data[,1]==x),],"Allele"),simplify=FALSE)
	# perform checks
	check = !checkData2(alleleReps=calls)
	if(check) stop(paste("Error in the allele calls provided - replicates do not match"))
	# return data
	return(calls)
	}

# Function to convert from peak heights designation to binary designation
# Requires separate input of designations
# This function is strictly for a single replicate
# Parameters:
# data = CSP alleles (csp$alleles[[x]])
# allelicCalls = alleles designations from file input
# (Coding: 0=nonallelic,1=allelic)
#convert.to.binary = function(data,allelicCalls)
#	{
#	# get rid of spaces
#	data = as.matrix(data)
#	data = matrix(gsub(" ","",data),ncol=ncol(data),dimnames=list(rownames(data),colnames(data)))
#	# replace .0 with nothing
#	data = gsub("[.]0","",data)
#	allelicCalls = as.matrix(allelicCalls)
#	allelicCalls = matrix(gsub(" ","",allelicCalls),ncol=ncol(allelicCalls),dimnames=list(rownames(allelicCalls),colnames(allelicCalls)))
#	# allelic calls
#	allelic = mapply(FUN=function(a,b) a[which(b==1)], split(data,row(data)), split(allelicCalls,row(allelicCalls)))
#	names(allelic) = rownames(data)
#	return(allelic)
#	}

convert.to.binary = function(data)
	{
	# get rid of spaces
	data = as.matrix(data)
	data = matrix(gsub(" ","",data),ncol=ncol(data),dimnames=list(rownames(data),colnames(data)))
	# replace .0 with nothing
	data = gsub("[.]0","",data)
	# allelic calls
	allelic = apply(data,MARGIN=1,FUN=function(x) x[which(!is.na(x))])
	names(allelic) = rownames(data)
	return(allelic)
	}

# combine rare alleles into one joined allele for a single locus
combine.rares.locus.peaks = function(alleleDb,cspProfile,knownProfiles,queriedProfile,rareThreshold=0.05,doDoubleStutter=FALSE)
    {
    # not in any profiles or in overstutter positions of profiles
    inProfile = !rownames(alleleDb)%in%c(unlist(cspProfile),unlist(knownProfiles),unlist(queriedProfile))
    inOverStutter = !as.numeric(rownames(alleleDb))%in%c(as.numeric(unlist(cspProfile))+1,as.numeric(unlist(knownProfiles))+1,as.numeric(unlist(queriedProfile))+1)
    inUnderStutter = !as.numeric(rownames(alleleDb))%in%c(as.numeric(unlist(cspProfile))-1,as.numeric(unlist(knownProfiles))-1,as.numeric(unlist(queriedProfile))-1)
    isRare = alleleDb[,1]<rareThreshold
    # index of alleles not in csp, unc, knowns or queried, and also probability < rareThreshold
    if(doDoubleStutter)
        {
        inDoubleOverStutter = !as.numeric(rownames(alleleDb))%in%c(as.numeric(unlist(cspProfile))+1,as.numeric(unlist(knownProfiles))+1,as.numeric(unlist(queriedProfile))+1)
        inDoubleUnderStutter = !as.numeric(rownames(alleleDb))%in%c(as.numeric(unlist(cspProfile))-1,as.numeric(unlist(knownProfiles))-1,as.numeric(unlist(queriedProfile))-1)
        index = which(inProfile&inOverStutter&inUnderStutter&isRare&inDoubleOverStutter&inDoubleUnderStutter)
        } else {
        index = which(inProfile&inOverStutter&inUnderStutter&isRare)
        }
    if(length(index)>0)
        {
        # remove indexed alleles, new allele has sum of probabilities, and mean of BP
        alleleDb = rbind(alleleDb[-index,],c(sum(alleleDb[index,1]),mean(alleleDb[index,2]))) 
        rownames(alleleDb)[nrow(alleleDb)] = "-1"
        }
    return(alleleDb)   
    }

# combine rare alleles into one joined allele for all loci
combine.rares.peaks = function(alleleDb, cspProfile, knownProfiles, queriedProfile,rareThreshold=0.05,doDoubleStutter=FALSE)
    {
    loci = colnames(cspProfile)
    sapply(loci,FUN=function(x) combine.rares.locus.peaks(alleleDb[[x]],cspProfile[,x],knownProfiles[,x],queriedProfile[,x],rareThreshold=rareThreshold,doDoubleStutter=doDoubleStutter),simplify=FALSE)
    }

add.stutter.index = function(alleleDb)
	{
 	# make i
  	index = which(!as.numeric(rownames(alleleDb))<0)
  	foo=function(x) 2/(max(x)-min(x))*(x-max(x))+1
  	outIndex = vector(length=nrow(alleleDb))
  	outIndex[index] = foo(as.numeric(rownames(alleleDb)[index]))
  	outIndex[-index] = 0
	return(cbind(alleleDb,outIndex))
	}



agnostic.hypothesis.peaks <- function(cspProfile, knownProfiles,
                                queriedProfile, alleleDb, ethnic='EA1',
                                adj=1e0, fst=0.02, combineRare=FALSE, 
                                rareThreshold=0.05,doDoubleStutter=FALSE) {
  # Helper function to figure out the input of most hypothesis.
  #
  # Basically, this groups operations that are done the same by defence and
  # prosection. 

  if(!ethnic%in%colnames(alleleDb)) stop("Chosen race code not included in database")

  # Read database and filter it down to requisite ethnicity and locus. 
  alleleDb = likeLTD:::ethnic.database(ethnic, colnames(cspProfile), alleleDb)

  # Adjust database to contain all requisite alleles
  alleleDb = likeLTD:::missing.alleles.peaks(alleleDb, cspProfile, queriedProfile, knownProfiles)
  alleleDb = likeLTD:::adjust.frequencies( alleleDb, 
                                 queriedProfile[1, colnames(cspProfile), 
                                                drop=FALSE],
                                 adj=adj, fst=fst )
  # combine rare alleles into a single allele
  if(combineRare) alleleDb = likeLTD:::combine.rares.peaks(alleleDb, cspProfile, knownProfiles, queriedProfile[1, colnames(cspProfile), drop=FALSE], rareThreshold,doDoubleStutter)

  # add index for stutter values
  alleleDb = sapply(alleleDb,FUN=add.stutter.index,simplify=FALSE)

  # Construct all profiles as arrays of  
  list(binaryProfile=cspProfile,
       alleleDb=alleleDb,
       queriedProfile=queriedProfile[1, colnames(cspProfile), drop=FALSE],
       knownProfs=knownProfiles) 
}

# function to automatically make allele calls
make.allelic.calls = function(peaksProfile,stutterPercent=0.15)
    {
    out = peaksProfile$alleles
    for(i in 1:length(out))
        {
        tmpAlleles = out[[i]]
        tmpHeights = peaksProfile$heights[[i]]
        for(j in 1:nrow(tmpAlleles))
            {
            for(k in 1:ncol(tmpHeights))
               {
               if(is.na(tmpAlleles[j,k])) break;
               if(k==ncol(tmpHeights)) tmpAlleles[j,k] = 1;
               parentIndex = which(round(as.numeric(tmpAlleles[j,]),1)==round(as.numeric(tmpAlleles[j,k])+1,1))
               if(length(parentIndex)==0)
                    {
                    tmpHeights[j,k] = 1
                    next
                    }
               if(tmpHeights[j,k]<stutterPercent*tmpHeights[j,parentIndex])
                    {
                    tmpHeights[j,k] = 0
                    } else {
                    tmpHeights[j,k] = 1
                    }
               }
            
            }
        out[[i]] = tmpHeights

        }
    return(out)
    }


# Creates prosecution hypothesis.
# Documentation is in man directory.
prosecution.hypothesis.peaks <- function(peaksFile, callsFile=NULL, refFile, ethnic='EA1',
                                   nUnknowns=0, adj=1e0, fst=0.02, linkageFile=NULL,
                                   databaseFile=NULL, relatedness=c(0,0), detectionThresh=30, 
                                   doDropin=FALSE, doDoubleStutter=FALSE,doOverStutter=FALSE, combineRare=TRUE, rareThreshold=0.05, kit=NULL,...) {
  if(is.null(databaseFile)&is.null(kit)) kit = "DNA17"
  alleleDb = load.allele.database(databaseFile,kit)
  peaksProfile = read.peaks.profile(peaksFile)
#  if(is.null(stutterPercent)&is.null(callsFile)) stop("Need either a file with allelic calls, or a stutter percentage to automatically call alleles")
#  if(is.null(callsFile))
#    {
#    allelicCalls = make.allelic.calls(peaksProfile,stutterPercent)
#    } else {
#    allelicCalls = read.allelic.calls(callsFile)
#    }
#  cspProfile = mapply(convert.to.binary,data=peaksProfile$alleles,allelicCalls=allelicCalls,SIMPLIFY=FALSE)
  cspProfile = sapply(peaksProfile$alleles,FUN=likeLTD:::convert.to.binary,simplify=FALSE)
  cspProfile = t(sapply(cspProfile,FUN=function(x) sapply(x,FUN=unlist)))
  if(identical(relatedness,c(0.5,0.25)))
	{
  	linkageInfo = load.linkage.info(linkageFile)
	rownames(linkageInfo) = linkageInfo[,1]
	linkageInfo = linkageInfo[,-1]
	linkIndex = which(colnames(linkageInfo)%in%colnames(cspProfile))
	linkageInfo = linkageInfo[linkIndex,linkIndex]
	}
  knownProfiles = read.known.profiles(refFile)
  if(sum(unlist(knownProfiles[, "queried"])) != 1)
    stop("Expect one queried profile on input.")
  qIndices = which(unlist(knownProfiles[, "queried"]))
  uIndices = which(!unlist(knownProfiles[, "queried"]))
  queriedProfile = knownProfiles[qIndices, , drop=FALSE]  
  # Puts queried profile at the end.
  knownProfiles = knownProfiles[c(uIndices, qIndices), , drop=FALSE] 

  result = likeLTD:::agnostic.hypothesis.peaks(cspProfile, knownProfiles,
                               queriedProfile, alleleDb, ethnic=ethnic,
                               adj=adj, fst=fst, combineRare=combineRare,
			       rareThreshold=rareThreshold, doDoubleStutter=doDoubleStutter)


  	# refIndiv is queried
  	nameRef = rownames(result$queriedProfile)
  	result[["refIndiv"]] = which(rownames(result$knownProfs)==nameRef)


  result = append(result, list(...))
  result[["peaksProfile"]] = peaksProfile$alleles
  result[["heightsProfile"]] = peaksProfile$heights
  result[["sizesProfile"]] = peaksProfile$sizes
  result[["knownProfs"]] = knownProfiles
  result[["nUnknowns"]] = nUnknowns
  result[["hypothesis"]] = "prosecution"
  result[["ethnic"]] = ethnic  
  result[["adj"]] = adj
  result[["fst"]] = fst
  result[["relatedness"]] = relatedness
  result[["doDropin"]] = doDropin
  result[["doDoubleStutter"]] = doDoubleStutter
  result[["doOverStutter"]] = doOverStutter
  result[["peaksFile"]] = peaksFile
  result[["callsFile"]] = callsFile
  result[["refFile"]] = refFile
  result[["databaseFile"]] = databaseFile
  result[["kit"]] = kit
  result[["detectionThresh"]] = detectionThresh
  if(identical(relatedness,c(0.5,0.25)))
	{
    result[["linkageInfo"]] = linkageInfo
    }

  sanity.check.peaks(result) # makes sure hypothesis has right type.
  result
}


# Creates defence hypothesis
# Documentation is in man directory.
defence.hypothesis.peaks <- function(peaksFile, callsFile=NULL, refFile, ethnic='EA1',  nUnknowns=0,
                               adj=1e0, fst=0.02, databaseFile=NULL, stutterPercent=NULL, linkageFile=NULL,
                               relatedness=c(0,0), detectionThresh=30, doDropin=FALSE, doDoubleStutter=FALSE,doOverStutter=FALSE, combineRare=TRUE, 
			       rareThreshold=0.05, kit=NULL,...) {
  if(is.null(databaseFile)&is.null(kit)) kit = "DNA17"
  alleleDb = load.allele.database(databaseFile,kit)
  peaksProfile = read.peaks.profile(peaksFile)
#  if(is.null(stutterPercent)&is.null(callsFile)) stop("Need either a file with allelic calls, or a stutter percentage to automatically call alleles")
#  if(is.null(callsFile))
#    {
#    allelicCalls = make.allelic.calls(peaksProfile,stutterPercent)
#    } else {
#    allelicCalls = read.allelic.calls(callsFile)
#    }
#  cspProfile = mapply(convert.to.binary,data=peaksProfile$alleles,allelicCalls=allelicCalls,SIMPLIFY=FALSE)
  cspProfile = sapply(peaksProfile$alleles,FUN=convert.to.binary,simplify=FALSE)
  cspProfile = t(sapply(cspProfile,FUN=function(x) sapply(x,FUN=unlist)))
  if(identical(relatedness,c(0.5,0.25)))
	{
  	linkageInfo = load.linkage.info(linkageFile)
	rownames(linkageInfo) = linkageInfo[,1]
	linkageInfo = linkageInfo[,-1]
	linkIndex = which(colnames(linkageInfo)%in%colnames(cspProfile))
	linkageInfo = linkageInfo[linkIndex,linkIndex]
	}
  knownProfiles = read.known.profiles(refFile)
  if(sum(unlist(knownProfiles[, "queried"])) != 1)
    stop("Expect one queried profile on input.")
  queriedProfile = knownProfiles[unlist(knownProfiles[, "queried"]), ,
                                 drop=FALSE]
  # Removing queried individual from knownProfiles. 
  knownProfiles = knownProfiles[!unlist(knownProfiles[, "queried"]), ,
                                drop=FALSE]

  result = agnostic.hypothesis.peaks(cspProfile, knownProfiles,
                               queriedProfile, alleleDb, ethnic=ethnic,
                               adj=adj, fst=fst, combineRare=combineRare,
			       rareThreshold=rareThreshold, doDoubleStutter=doDoubleStutter) 
   
  result[["refIndiv"]] = 1

  result = append(result, list(...))
  result[["peaksProfile"]] = peaksProfile$alleles
  result[["heightsProfile"]] = peaksProfile$heights
  result[["sizesProfile"]] = peaksProfile$sizes
  result[["knownProfs"]] = knownProfiles
  result[["nUnknowns"]] = nUnknowns + 1
  result[["hypothesis"]] = "defence"
  result[["ethnic"]] = ethnic 
  result[["adj"]] = adj
  result[["fst"]] = fst
  result[["relatedness"]] = relatedness
  result[["doDropin"]] = doDropin
  result[["doDoubleStutter"]] = doDoubleStutter
  result[["doOverStutter"]] = doOverStutter
  result[["peaksFile"]] = peaksFile
  result[["refFile"]] = refFile
  result[["callsFile"]] = callsFile
  result[["databaseFile"]] = databaseFile
  result[["kit"]] = kit
  result[["detectionThresh"]] = detectionThresh
  result[["stutterPercent"]] = stutterPercent
  if(identical(relatedness,c(0.5,0.25))) result[["linkageInfo"]] = linkageInfo


  sanity.check.peaks(result) # makes sure hypothesis has right type.
  result
}


missing.alleles.peaks = function(alleleDb, cspProfile, queriedProfile, knownProfiles) {
  # Adds missing alleles to database
  #
  # There may be alleles in the crime-scene profile, in the queried individual,
  # or in the individuals subject to dropout, which are not present in the
  # database. Theses alleles are added into the database with count 1 and
  # fragment length 0.
  if(!is.matrix(cspProfile)) stop("input should be a matrix.")
  if(is.null(colnames(cspProfile)))
    stop("input matrix does not have column names.")
  if(!is.matrix(knownProfiles)) stop("input should be a matrix.")
  if(!is.matrix(queriedProfile)) stop("input should be a matrix.")
  for(locus in colnames(cspProfile)) {
    dbAlleles = rownames(alleleDb[[locus]])
    cspAlleles = unique(unlist(cspProfile[, locus]))
    qAlleles = unique(unlist(queriedProfile[, locus]))
    knownAlleles = unique(unlist(knownProfiles[, locus]))
    allAlleles = unique(c(cspAlleles,qAlleles, knownAlleles)) 
    missingAlleles = allAlleles[!allAlleles%in%dbAlleles]
    # At this point, some of the missingAlleles might just be saying move
    # along, nothing to see.
    missingAlleles = setdiff(missingAlleles, c("", NA, "NA"))
    # Now add new rows to database. 
    if(length(missingAlleles)) {
      newrows = matrix(c(1, 0), nrow=length(missingAlleles), ncol=2,
                       byrow=TRUE)
      rownames(newrows) = missingAlleles
      alleleDb[[locus]] = rbind(alleleDb[[locus]], newrows)
    }
  }
  alleleDb
}


# Should check that hypothesis is syntactically correct
# This means matrices as opposed to lists, so on and so forth.
# It does not mean it the input should make any sort of sense.
sanity.check.peaks = function(hypothesis) {
  errors = c()
  if(!is.matrix(hypothesis$binaryProfile)) 
     errors = c(errors, "binaryProfile should be a matrix 'replicates' x 'loci'.")
  if(!is.matrix(hypothesis$queriedProfile))
     errors = c(errors, "queriedProfile should be a matrix 'profiles' x 'loci'.")
  if(ncol(hypothesis$binaryProfile) != ncol(hypothesis$queriedProfile))
      errors = c( errors,
                  "Number of loci of binaryProf and queriedProf do not match." )
  if(!is.null(hypothesis$stutterPercent))
    {
    if(hypothesis$stutterPercent<0|hypothesis$stutterPercent>1)
        errors = c( errors,
            "stutterPercent must be between 0 and 1." )
    if(hypothesis$stutterPercent>0.3)
      print("Warning: stutterPercent is set unusually high.")   
    }
  if(hypothesis$detectionThresh<0)
        errors = c( errors,
            "Detection threshold should not be set lower than 0." )
  if(hypothesis$detectionThresh<0)
    print("Detection threshold set unusually high.")
  if(length(errors) != 0) {
    cat("There seems to be an error with the input hypothesis.\n")
    for(error in errors) cat(paste(c(error, "\n")))
    stop()
  }
}


transform.to.locus.centric.peaks = function(hypothesis) {
  # Transform hypothesis to locus centric modes.
  # 
  # This means we reorganise the data to be in lists of the loci.
  result = list()
  for(locus in colnames(hypothesis$binaryProfile)) {
    # Value of the resulting list for a given locus 
    locusValue = list()
    # Loop over all items in original list.
    for(key in names(hypothesis)) {
      value = hypothesis[[key]]
      # special treatment for peaks profiles (peaks, heights and sizes)
      if(key%in%c("peaksProfile","heightsProfile","sizesProfile"))
        {
        locusValue[[key]] = sapply(value,FUN=function(x) x[locus, , drop=FALSE],simplify=FALSE)
	for(i in 1:length(locusValue[[key]]))
		{
		tmp = data.frame(t(locusValue[[key]][[i]][!is.na(locusValue[[key]][[i]])]),stringsAsFactors=FALSE)
		rownames(tmp) = locus
		if(ncol(tmp)>0) colnames(tmp) = paste("Allele",1:ncol(tmp),sep=".")
		locusValue[[key]][[i]] = tmp
		}
        } else {
      # If a matrix and locus is either in rows or columns, then add only locus
      # part of the matrix. 
      # If a list, then add only the locus specific part of the list.
      if(is.matrix(value)) {
        if(locus %in% colnames(value)) 
          locusValue[[key]] = value[, locus, drop=FALSE]
        else if(locus %in% rownames(value))
          locusValue[[key]] = value[locus, , drop=FALSE]
      } else if(is.list(value) && (locus %in% names(value))) 
        locusValue[[key]] = value[[locus]]
        }
      # If has not been added yet, then not locus dependent, and add it as a
      # whole. 
      if(!key %in% names(locusValue)) locusValue[[key]] = value 
    }
    # Only locus to result if not empty.
    if(length(locusValue) > 0) result[[locus]] = locusValue
  }
  result
}

