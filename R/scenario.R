read.csp.profile = function(path) {
  # Reads crime scene profile from file.
  #
  # Returns a matrix containing the profiles for each locus (columns) and each
  # replicate (rows). Each element is a string or list of string naming the
  # alleles present in the CSP.
  #
  # For simplicity, the return does not contain loci for which replicates are
  # all NA.
  if(!file.exists(path)) stop(paste(path, "does not exist."))
  profile = read.table(path, header=T, colClasses='character', sep=',')
  profile = sapply(profile, function (n) strsplit(n, ','))
  profile = profile[profile[, "X"] != "Uncertain", 5:ncol(profile)]
  allNA = apply(profile, 2, function(n) !all(is.na(n)))
  profile[, allNA]
}
read.unc.profile = function(path) {
  # Reads uncertain alleles in crime scene profile from file.
  #
  # Returns a matrix containing the profiles for each locus (columns) and each
  # replicate (rows). Each element is a string or list of string naming the
  # uncertain alleles present in the CSP.
  if(!file.exists(path)) stop(paste(path, "does not exist."))
  profile = read.table(path, header=T, colClasses='character', sep=',')
  profile = sapply(profile, function (n) strsplit(n, ','))
  profile[profile[, "X"] == "Uncertain", 5:ncol(profile)]
}

read.known.profiles = function(path) {
  # Reads known profiles from file.
  #
  # Returns a matrix where each person is a given row. The column labeled
  # "queried" indicates whether that individual is queried or not. Other
  # columns are loci. Each element in a loci is a vector of two strings naming
  # the allele of the individual.
  if(!file.exists(path)) stop(paste(path, "does not exist."))
  profiles = read.table(path, header=T, colClasses='character', row.names=1,
                        sep=',', quote = "\"")
  result = sapply(profiles, function (n) strsplit(n, ','))
  row.names(result) = row.names(profiles)
  colnames(result)[1] = "queried"
  result[, "queried"] = result[, "queried"] == "queried"
  result
}

determine.dropout = function(knownProfiles, cspProfile) {
  # Adds dropout column to known profiles.
  #
  # An individual is subject to doprout if the individual's allele at one or
  # more locus is not in the crime-scene profile.
  #
  # Parameters:
  #   knownProfiles: A matrix containing the known profiles. 
  #   cspProfile: The crime scene profile.
  # Returns: 
  #   A list with one element per individual. Each is True if that individual
  #   is subject to dropout.
  
  # Alleles present in any replicate of the csp per locus.
  per.indiv = function(indiv) {
    per.locus = function(locus) {
      all(sapply(cspProfile[, locus], function(n) indiv[[locus]] %in% n))
    }
    all(sapply(colnames(cspProfile), per.locus))
  }
  !apply(knownProfiles, 1, per.indiv)
}

masking.and.uncertain.profile = function(uncProfile, maskingProfiles) {
  # Adds masking profiles to each replicate of uncertain alleles.
  # 
  # Masking and uncertain alleles are somehow treated the same.
  if(length(maskingProfiles) == 0) return(uncProfile)
  for(locus in intersect(colnames(uncProfile), colnames(maskingProfiles))) {
    alleles = Reduce(union, maskingProfiles[, locus, drop=FALSE])
    for(i in 1:nrow(uncProfile))
      uncProfile[[i, locus]] = union(uncProfile[[i, locus]], alleles)
  }
  uncProfile
}

missing.alleles = function(alleleDb, cspProfile, noDropoutProfiles) {
  # Adds missing alleles to database
  #
  # There may be alleles in the crime-scene profile, in the queried individual,
  # or in the individuals subject to dropout, which are not present in the
  # database. Theses alleles are added into the database with count 1 and
  # fragment length 0.
  if(!is.matrix(cspProfile)) stop("input should be a matrix.")
  if(is.null(colnames(cspProfile)))
    stop("input matrix does not have column names.")
  if(!is.matrix(noDropoutProfiles)) stop("input should be a matrix.")

  for(locus in colnames(cspProfile)) {
    dbAlleles = rownames(alleleDb[[locus]])
    cspAlleles = unique(unlist(cspProfile[, locus]))
    dropAlleles = unique(unlist(noDropoutProfiles[, locus]))
    missingAlleles = setdiff(cspAlleles, union(dbAlleles, dropAlleles))
    # At this point, some of the missingAlleles might just be saying move
    # along, nothing to see.
    missingAlleles = setdiff(missingAlleles, c("", NA, "NA"))
    # Now add new rows to database. 
    if(length(missingAlleles)) {
      newrows = matrix(c(1, 0), nrow=length(missingAlleles), byrow=T)
      rownames(newrows) = missingAlleles
      alleleDb[[locus]] = rbind(alleleDb[[locus]], newrows)
    }
  }
  alleleDb
}

transform.to.locus.centric = function(scenario) {
  # Transform scenarios to locus centric modes.
  # 
  # This means we reorganise the data to be in lists of the loci.
  result = list()
  for(locus in colnames(scenario$cspProfile)) {
    # Value of the resulting list for a given locus 
    locusValue = list()
    # Loop over all items in original list.
    for(key in names(scenario)) {
      value = scenario[[key]]
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
      # If has not been added yet, then not locus dependent, and add it as a
      # whole. 
      if(!key %in% names(locusValue)) locusValue[[key]] = value 
    }
    # Only locus to result if not empty.
    if(length(locusValue) > 0) result[[locus]] = locusValue
  }
  result
}


agnostic.scenario <- function(cspProfile, uncProfile, knownProfiles,
                              queriedProfile, alleleDb, ethnic='EA1', adj=1e0,
                              fst=0.02) {
  # Helper function to figure out the input of most scenarios.
  #
  # Basically, this groups operations that are done the same by defense and
  # prosection. 

  # Read database and filter it down to requisite ethnicity and locus. 
  alleleDb = ethnic.database(ethnic, colnames(cspProfile), alleleDb)

  # Figure out which profiles show dropout.
  dropout = determine.dropout(knownProfiles, cspProfile)

  # Creates profiles with and without dropout.
  dropoutProfs = knownProfiles[dropout, colnames(cspProfile), drop=FALSE]
  noDropoutProfs = knownProfiles[!dropout, colnames(cspProfile), drop=FALSE]
  # Adjust uncertain profile to also contain masking (no-dropout) profiles.
  uncProf = masking.and.uncertain.profile(uncProfile, noDropoutProfs)

  # Adjust database to contain all requisite alleles.
  alleleDb = missing.alleles(alleleDb, cspProfile, noDropoutProfs)
  alleleDb = adjust.frequencies( alleleDb, 
                                 queriedProfile[1, colnames(cspProfile)],
                                 adj=adj, fst=fst )

  # Construct all profiles as arrays of  
  list(cspProfile=cspProfile, uncProf=uncProf, dropoutProfs=dropoutProfs,
       alleleDb=alleleDb,
       queriedProfile=queriedProfile[1, colnames(cspProfile)])
}

prosecution.scenario <- function(mixedFile, refFile, ethnic='EA1', nUnknowns=0, adj=1e0,
                                 fst=0.02, databaseFile=NULL, ...) {
  # Creates prosecution scenario.
  alleleDb = load.allele.database(databaseFile)
  cspProfile = read.csp.profile(mixedFile)
  uncProfile = read.unc.profile(mixedFile)
  knownProfiles = read.known.profiles(refFile)
  if(sum(unlist(knownProfiles[, "queried"])) != 1)
    stop("Expect one queried profile on input.")
  queriedProfile = knownProfiles[unlist(knownProfiles[, "queried"]), ,
                                  drop=FALSE]
  
  result = agnostic.scenario(cspProfile, uncProfile, knownProfiles,
                             queriedProfile, alleleDb, ethnic=ethnic,
                             adj=adj, fst=fst)

  result = append(result, list(...))
  result[["nUnknowns"]] = nUnknowns
  result[["relatedness"]] = c(0, 0)
  result
}

defense.scenario <- function(mixedFile, refFile, ethnic='EA1',  nUnknowns=0,
                             adj=1e0, fst=0.02, databaseFile=NULL, ...) {
  # Creates defense scenario.
  alleleDb = load.allele.database(databaseFile)
  cspProfile = read.csp.profile(mixedFile)
  uncProfile = read.unc.profile(mixedFile)
  knownProfiles = read.known.profiles(refFile)
  if(sum(unlist(knownProfiles[, "queried"])) != 1)
    stop("Expect one queried profile on input.")
  queriedProfile = knownProfiles[unlist(knownProfiles[, "queried"]), ,
                                  drop=FALSE]
  # Removing queried individual from knownProfiles. 
  knownProfiles = knownProfiles[!unlist(knownProfiles[, "queried"]), ,
                                drop=FALSE]
  
  result = agnostic.scenario(cspProfile, uncProfile, knownProfiles,
                             queriedProfile, alleleDb, ethnic=ethnic,
                             adj=adj, fst=fst)

  result = append(result, list(...))
  result[["nUnknowns"]] = nUnknowns + 1
  result
}
