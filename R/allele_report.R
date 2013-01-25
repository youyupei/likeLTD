unusual.alleles <- function(afreq, profile) {
  # Searches for unusual alleles.
  #
  # Such alleles are likely to be mistakes.
  #
  # Args:
  #   freqencies: Table of allele afreq“
  #   profile: A profile with just the facts. 
  #
  # returns: data frame with the rare alleles. 

  rare =data.frame(name=NULL, locus=NULL, allele=NULL, EA1=NULL, EA3=NULL,
                   EA4=NULL)

  for(name in row.names(profile)){
    for(locus in colnames(profile)){
      # loop over unique allele names.
      uniqueAlleles <- unique(strsplit(profile[name,locus],',')[[1]]) 
      for(allele in uniqueAlleles){

        condition = afreq$Marker==locus & afreq$Allele==allele
        x = afreq[condition, ]
        if(x$EA1<2 | x$EA3<2 | x$EA4<2) {
          frame = data.frame(name=name, locus=locus, allele=allele, EA1=x$EA1,
                             EA3=x$EA3, EA4=x$EA4)
          rare=rbind(rare, frame)
        }
      }  # loop over alleles
    } # loop over loci
  } # loop over people

  return(rare)
}

read.csp.profile <- function(path) {
  # Reads profile from path.
  #
  # Args:
  #   path: Path to file with the profile. 
  # returns: The profile 
  
  raw = read.table(path, header=T, colClasses='character', sep=',')
  return(raw[,5:dim(raw)[2]])
}
read.ref.profile <- function(path) {
  # Reads profile from path.
  #
  # Makes sure everything is formatted as should be, e.g. removes spaces and
  # such.
  #
  # Args:
  #   path: Path to file with the profile. 
  # returns: The profile minus some faff.
  
  raw = read.table(path, header=T, colClasses='character', row.names=1,
                   sep=',', quote = "\"")
  result = raw[,2:ncol(raw)]
  # removes spaces from formatting
  for(locus in 1:ncol(result)){
    for(person in 1:nrow(result)){
      result[person, locus] = gsub(' ', '', result[person, locus])
    }
  }
  return(result)
}
queried.vs.known <- function(path) {
  # Reads profile from path and returns queried vs known column.
  #
  # Args:
  #   path: Path to file with the profile. 
  raw = read.table(path, header=T, colClasses='character', row.names=1,
                   sep=',', quote = "\"")
  return(raw$known.queried == 'queried')
}

internal.representation <- function(profile) {
  # Takes a profile and transforms it to the cprofs format.
  #
  # cprofs is an internal format to do calulations on a crime scene profile.
  #
  # Args:
  #   profile: The crime-scene profile to  process.
  # Returns: The same information in a different format.

  # Half the crime profile is empty for some reason.
  nrep = nrow(profile)/2

  result = list()
  # each column in profile will be a new column in the output
  for(locus in colnames(profile)){	

    # creates a list of reps in each loci
  	result[[locus]]=list()  

    # now creates the items in the list at this locus
	  for(rep in 1:nrep){
		  if (!is.na(profile[(2*rep)-1,locus])) {
        # Figures out csp
			  csp = gsub(' ','', profile[(2*rep)-1, locus]);
        if(csp != '') csp = unlist(strsplit(csp, ',')) 
        # Figures out unc
			  unc = gsub(' ','', profile[(2*rep), locus]);
        if(unc != '') unc = unlist(strsplit(unc, ','))
        # fills it with NULL csp and unc at first
		    result[[locus]][[rep]]=list(csp=csp, unc=unc) 
		  }
      # NAs must be formatted differently
      else result[[locus]][rep]=list(NULL) 
    } # loop over each replicate

  } # loop over loci
  return(result)
}

estimate.csp <- function(ref, cprofs) {
  # Estimate how well each reference profile is represented in the CSP
  #
  # Args: 
  #   ref: reference profile.
  #   cprofs: crime scene profile, internal representation.
  # Returns: 
  #   A data frame where row refer to people and colums to replicates, and the
  #   content is an estimate of how well that person and loci matches the CSP.
  #   There is an additional column for the sum over all replicates.
  nrep = length(cprofs[[1]])
      
  # Constructs the result data frame.
  result = data.frame(array(0, c(nrow(ref), nrep+1)), row.names=rownames(ref))
  colnames(result)[1:nrep] = sapply(1:nrep, function(n) {paste('run', n)})
  colnames(result)[nrep+1] = 'Total' 
  
  # now add data.
  for(person in row.names(ref)){
	  for(rep in 1:nrep){
      # number of alleles in common for given person and CSP replicate, across
      # all loci
			nalleles = 0
      # How many reference loci exists in a given replicate
      nloci = 0
			for(locus in 1:length(cprofs)){
				if(!is.null(cprofs[[locus]][[rep]])) {
          # figure out unique alleles in reference 
			  	alleles  <- unique(unlist(strsplit(ref[person,locus],',')))
          # figure out howmany of these allele are in CSP
          represented = alleles %in% cprofs[[locus]][[rep]]$csp
          nalleles <- nalleles + sum(represented) / length(represented)
          # this locus exists in replicate 
          nloci    <- nloci + 1
		    }
      }
      if(nloci > 0) result[person, rep] = round(100*nalleles/nloci)
	  }
  }
  if(nrep==1)  result[, nrep+1] = result[, nrep]
  else         result[, nrep+1] = round(rowSums(result)/nrep)
   
  # Reorders rows and return
  return(result[order(result[, nrep+1], decreasing=T), ])
}

private_nullify.item <- function(n) {
  # Returns NULL if n is character(0), else returns n.
  #
  # This is a helper function for internal use only.
  if(identical(n, character(0))) return(NULL)
  return(n)
}

flatten.csp <- function(v, subitem='csp') {
  # Flatten and filters elements of cprofs defining one locus.
  # 
  # This function is mostly for internal use. It takes a single lo“flattens the
  # cprofs object, selects its subitems with a given name, and put them all in
  # a single vector.
  #
  # Args:
  #   v: A locus element from cprofs. 
  #   subitem: Name of the subitem
  result <- sapply(v, function(n) n[subitem], simplify=FALSE)
  result <- unlist(result, use.names=FALSE)
  return(Filter(function(n) !(n == ""), result))
}
all.alleles <- function(cprofs, subitem='csp') {
  # Flatten and filters elements of cprofs into a list of vectors.
  # 
  # This function is mostly for internal use. It applies to each element of
  # cprofs the function flatten.csp
  #
  # Args:
  #   v: A locus element from cprofs. 
  # First, get all results across CSP
  result <- sapply(cprofs, function(n) flatten.csp(n, subitem))
  # Then nullifies elements evaluating to character 0. 
  return(sapply(result, private_nullify.item))
}
replicated.alleles <- function(cprofs) {
  # For each locus, lists replicated alleles.
  #  
  # Args:
  #  v: cprofs data object
  # Returns:
  #   A list over loci, where each subitem is a vector of replicated alleles.
  allcsp <- all.alleles(cprofs, subitem='csp')
  allcsp <- sapply(allcsp, function(n) n[duplicated(n)])
  return(sapply(allcsp, private_nullify.item))
}
unreplicated.alleles <- function(cprofs) {
  # For each locus, lists replicated alleles.
  #  
  # Args:
  #  v: cprofs data object
  # Returns: 
  #   A list over loci, where each subitem is a vector of unreplicated alleles.
  allcsp <- all.alleles(cprofs, subitem='csp')
  allcsp <- sapply(allcsp, function(n) setdiff(n, n[duplicated(n)]))
  return(sapply(allcsp, private_nullify.item))
}

allele.table <- function(cprofs) {
  # Allele table of a Crime Scene Profile
  #
  # Creates a better looking table for a CSP. 
  nloc = length(cprofs)
  nrep = length(cprofs[[1]])

  alleles = data.frame(array(0, c(2*nrep,nloc)))
  rowfunc <- function(n) {c(paste('csp', n), paste('unc', n))}
  row.names(alleles) = sapply(1:nrep, rowfunc)
  colnames(alleles)  = names(cprofs)

  for(l in 1:nloc){ 
    for(r in 1:nrep){
      if(!is.null(cprofs[[l]][[r]])){
        alleles[2*r-1,l] = paste(cprofs[[l]][[r]]$csp, collapse=' ')
        alleles[2*r,l]   = paste(cprofs[[l]][[r]]$unc, collapse=' ')
      }
      else {
        alleles[2*r-1,l]='NA'
        alleles[2*r,l] ='NA'
      }
    } # loop over replicates
  } # loop over loci

  return(alleles)
}

summary.generator = function(queried, known, ref, cprofs){
  # Summnarizes the allele table for Q and K
  # 
  # Args:
  #  queried: Names of queried individuals in reference profile (rows).
  #  known: Names of known individuals in reference profile (rows).
  #  ref: the reference profile
  #  cprofs: the crime scene profile, external format
  allNames = c( sapply(queried, function(n) paste(n, '(Q)'), USE.NAMES=FALSE),
                sapply(known, function(n) paste(n, '(K)'), USE.NAMES=FALSE) )
  ncont = length(allNames)
  # index: indexing of ref such that queries come first and known profiles come
  # second.
  index = sapply(c(queried, known), function(n) which(row.names(ref) == n))
  # knownLoci: list over locus where each item contains the alleles with the
  # same indexing as index. 
  knownLoci = list()
  for(locus in names(ref)) {
    string = paste(ref[[locus]][index], collapse=',')
    knownLoci[[locus]] = unlist(strsplit(string, ',')) 
  }

  # Figures out all, replicated, and unreplicated alleles in crime scene.
  cspCombinedAll  <- all.alleles(cprofs)
  replicatedAll   <- replicated.alleles(cprofs) 
  unreplicatedAll <- unreplicated.alleles(cprofs) 

  #  construct summary data frame scaffold.
  summary = data.frame(array(0, c(length(allNames)+1, length(ref))))
  row.names(summary) = c(allNames,'Unattributable')
  colnames(summary)  = names(cprofs)

  # otherRep: per locus, number of unattributable and replicated
  # alleles
  otherRep = list()
  # otherUnrep: per locus, number of unattributable and unreplicated
  # alleles
  otherUnrep = list()
  for(locus in names(ref)){ 
    # loop over people
    for(c in 1:length(allNames)){
      # alleles for this person in reference profile
      cont = knownLoci[[locus]][(2*c-1):(2*c)]
      rep   = paste(cont[ cont %in% replicatedAll[[locus]]],   collapse=' ')
      unrep = paste(cont[ cont %in% unreplicatedAll[[locus]]], collapse=' ')
      none  = paste(cont[!cont %in% cspCombinedAll[[locus]]], collapse=' ')
      summary[c,locus] = paste(rep, '{', unrep,'}','[', none,']',sep='')
    }
    # other alleles not found in the cprofs
    notfound = !cspCombinedAll[[locus]] %in% knownLoci[[locus]]
    unattr = cspCombinedAll[[locus]][notfound]
    unattrRep = unique(unattr[duplicated(unattr)])
    unattrRepString = paste(unattrRep, collapse=' ')
    unattrUnrep = unattr[!unattr %in% unattrRep]
    unattrUnrepString = paste(unattrUnrep, collapse=' ')

    otherRep[locus] = length(unattrRep)
    otherUnrep[locus] = length(unattrUnrep)

    summary[nrow(summary), locus] = paste( unattrRepString, '{',
                                           unattrUnrepString, '}',
                                           sep='' )
  }

  # coerce to numerics
  otherRep = as.numeric(otherRep)
  otherUnrep = as.numeric(otherUnrep)
  return( list(summary=summary, otherRep=otherRep, otherUnrep=otherUnrep) )
}

hypothesis.generator = function(queried, known, otherRep, otherUnrep){
  # function to calculate possible hypotheses
  # 
  # Needs to be done with various permutations of 2 Q's

  otherBoth = otherRep + otherUnrep
  # minU: Minimum number of unknown DNA in CSP 
  minU = ceiling(max(otherRep)/2)
  # maxU: Maximum number of unknown DNA in CSP 
  # TODO: Not sure the next line should have "divide by 2". Might be over
  # replicated only.
  maxU = ceiling(max(otherBoth)/2) 

  hypotheses = c()
  for(unknowns in minU:maxU) {

    extras   = otherBoth - 2*unknowns
    unattributableAlleles = sum(extras[extras>0])

    Q = paste(queried,'(Q)',sep='')
    others = c()
    if(length(known) > 0 ) {
      pasting <- function(n) paste(known[n], '(K', n, ')', sep='')
      others = sapply(1:length(known), pasting)
    }
    if(unknowns > 0) {
      U = sapply(1:unknowns, function(n) paste('U', n, sep=''))
      others = c(others, U, recursive=TRUE)
    }
    # If condition is TRUE, then there might be dropins.
    if(unattributableAlleles !=0) {
      D = paste(unattributableAlleles, 'dropins', sep='')
      others = c(others, D, recursive=TRUE)
    }
         
    others <- paste(others, collapse='+')
    string <- paste(Q, '+', others,' vs X+', others, sep='') 
    hypotheses[unknowns-minU+1] = string
  }
  return(hypotheses[!is.na(hypotheses)])
}


suggested.hypothesis = function(queried, known, ref, cprofs) {
  # Figures out sensible explanation of crime profile.
  #
  # Returns: List of strings with sensible hypothesis. 

  # estimates all possible hypotheses
  summary = summary.generator(queried, known, ref, cprofs)
  otherBoth = summary$otherRep + summary$otherUnrep
  otherRep  = summary$otherRep

  # Each element of hypothesis is a different, hum, hypothesis.
  hypothesis = list(list(queried, known))
  if(length(queried)==2)
    hypothesis = list( hypothesis,
                       list(queried[2], known), 
                       list(queried[1], c(queried[2], known)), 
                       list(queried[2], c(queried[1], known)) )
   
  generate.hypothesis = function(n) {
    #  Generates a hypothesis
    #  Ugly way to call badly designed functions.
    summary = summary.generator(n[[1]], n[[2]], ref, cprofs)
    result = hypothesis.generator(n[[1]], n[[2]], summary$otherRep,
                                  summary$otherUnrep)
    return(result)
  }
  # suggested: The hypothesis generator is applied to each hypothesis. Each
  # element of suggested is a string detailing the hypothesis.
  suggested = sapply(hypothesis, generate.hypothesis, USE.NAMES=FALSE,
                     simplify=FALSE)
  return(suggested)
}

pack.admin.input = function( frequencyFile, mixedFile, refFile,
                             caseName='dummy', outputPath=NULL,
                             checkFiles=TRUE ) {
  # Packs and verifies administrative information.
  #
  # Args:
  #   frequencyFile: File holding the allele afreq
  #   mixedFile: Mixed crime scene profile
  #   refFile: Reference profiles
  #   casename: Name of the current case
  #   outputPath: Path where the output should be stored.
  #   checkFiles: Whether to check file existence and such.
  if(is.null(outputPath)) outputPath = file.path(getwd(), caseName)

  if(checkFiles) {
    for(path in c(frequencyFile, mixedFile, refFile)) {
      if(!file.exists(path))
        stop(paste(path, "does not exist."))
      else { 
        info <- file.info(path)
        if(info$isdir) stop(paste(path, "is not a file."))
      }
    } # loop over files.
    if(file.exists(outputPath) & !file.info(outputPath)$isdir) 
      stop(paste(outputPath, " exists and is not a directory."))
  } # condition whether to check files.
  admin = list( caseName='hammer',
                frequencyFile=frequencyFile,
                mixedFile=mixedFile,
                refFile=refFile,
                outputPath=outputPath )
  return(admin)
}

pack.genetics.input = function(admin, nameK=NULL, nameQ=NULL, dropin=FALSE,
                               unknowns=0, ethnic='EA1', fst=NULL) {
  # Generates and packs genetics input into a list.
  #
  # This list contains all the genetics data needed to perform all subsequent
  # calculations. 
  # Args:
  #   admin: List containing administration information, as packed by
  #          pack.admin.input.
  #   nameK: list of names of known contributors in CSP.
  #   nameQ: list of queried contributors in profile.
  #   dropin: Whether to model drop-ins. 
  #   unknown: Number of unkown contributors in CSP.
  #   ethnic: Ethnicity of contributors.
  #   fst: not sure. 
  afreq        = read.table(admin$frequencyFile, sep="\t", header=T)
  cspData      = read.csp.profile(admin$mixedFile)
  refData      = read.ref.profile(admin$refFile)
  cprofs       = internal.representation(cspData)
  QvK <- queried.vs.known(admin$refFile)
  nameQ        = row.names(refData)[QvK]
  nameK        = row.names(refData)[!QvK]
  estimates    = estimate.csp(refData, cprofs)
  summary      = summary.generator(nameQ, nameK, refData, cprofs)
  if(is.null(fst)) {
    if(ethnic == 'EA1') fst = 0.02
    else                fst = 0.03
  }

  genetics = list( # Allele frequency table.
                   afreq       = afreq,
                   # Crime scene profile data.
                   cspData     = cspData,
                   # Reference contributors data.
                   refData     = refData,
                   # Crime scene profile in internal representation.“
                   cprofs      = cprofs,
                   # Names of queried contributors.
                   nameQ       = nameQ,
                   # Names of known contributors other than queried.
                   nameK       = nameK,
                   # summary
                   summary     = summary,
                   # Estimates on how well each contributor is represented in
                   # CSP.
                   estimates   = estimates, 
                   # Number of unkown contributors.
                   unknowns    = unknowns,
                   # Whether to model dropins.
                   dropin      = dropin,
                   # Contributors' ethnicity.
                   ethnic      = ethnic,
                   # Not sure.
                   fst         = fst,
                   # Number of replicates.
                   nrep        = length(cprofs[[1]]) )
  return(genetics)
}


allele.report <- function(admin, genetics=NULL) {
  # Generates an allele summary from input profiles:
  #
  # Args:
  #   admin: List containing administration data, as packed by
  #          pack.admin.input.
  
  # reads genetics information if not given on input.
  library('gplots')
  if(is.null(genetics)) genetics = pack.genetics.input(admin)


  # Generate 'summary', whether multiple Q or single Q

  #-----------------------------------------------
  # 12. Write the allele report
  # works out sensible font size
  cspTotal = uncTotal = numeric(genetics$nrep)
  for(l in names(genetics$cprofs)){ 
    for(r in 1:length(l)){
      cspTotal[r] = cspTotal[r] + length(genetics$cprofs[[l]][[r]]$csp)
      uncTotal[r] = uncTotal[r] + length(genetics$cprofs[[l]][[r]]$unc)
    }
  }
  tablesize = min(1,25/max(c(cspTotal,uncTotal)))

  outputPath = file.path( admin$outputPath, 
                          paste(admin$caseName, '-allele_report.pdf', sep='') )
  
  # Start plotting proper. 
  pdf(width=8.25,height=11.5,file=outputPath)
  par(mfrow=c(6,1), mai=c(0.7,0.5,0.5,0.5))

  # Allele report
  alleles <- allele.table(genetics$cprofs)
  textplot(alleles, valign='top',cex=tablesize)
  title(paste(admin$caseName,'Allele Report'))

  # Tabular Summary 
  textplot(genetics$summary$summary, valign='top', cex=tablesize)
  title('Summary. {}=unreplicated, []=absent')

  # Unattributable Alleles
  otherBoth = genetics$summary$otherRep + genetics$summary$otherUnrep
  otherRep  = genetics$summary$otherRep
  yaxp=c(0, max(otherBoth), max(otherBoth))
  if(max(otherBoth)==0) yaxp=c(0,1,1)
  barplot(otherBoth, names.arg=names(genetics$cprofs), yaxp=yaxp,
          main='Unattributable alleles', ylab='No. alleles')
  barplot(add=T, col='black', otherRep, yaxt='n')
  if(max(otherBoth) > 3) abline(h=2, lty=2)
  if(max(otherBoth) > 5) abline(h=4, lty=2)
  if(max(otherBoth) > 7) abline(h=6, lty=2) 
  if(max(otherBoth) > 9) abline(h=8, lty=2);
  legend(x=0, y=max(otherBoth), yjust=0, bty='n',
         c('replicated','unreplicated'), xpd=NA, col=c('black','grey'),
         pch=c(15,15))

  # Rare Alleles
  rare <- unusual.alleles(genetics$afreq, genetics$refData)
  if(length(rare)==0) rare = 'No unusual alleles in Q or K profiles'
  textplot(rare, cex=1, valign='top')
  title('Rare alleles')

  # Approximae representation.
  # Estimate how well each reference profile is represented in the CSP
  estimates <- estimate.csp(genetics$refData, genetics$cprofs)
  textplot(estimates, cex=1, valign='top')
  title('Approximate representation %')

  # Suggested hypothesis.
  suggested = suggested.hypothesis(genetics$nameQ, genetics$nameK,
                                   genetics$refData, genetics$cprofs)
  textplot(suggested, cex=1, valign='top')
  title(main='Suggested Hypotheses',
        sub= paste('Note: likeLTD is limited to 2 unknowns. If more than 2 ', 
                   'unknowns are suggested,this may indicate a mixed profile',
                   'with too many alleles to be statisically informative.',
                   sep=''))
  dev.off()
  print(paste('Writing report to', outputPath))
}
