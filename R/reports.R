unusual.alleles <- function(afreq, profile) {
  # Searches for unusual alleles.
  #
  # Such alleles are likely to be mistakes.
  #
  # Args:
  #   freqencies: Table of allele afreqs
  #   profile: The reference profile in the genetics list
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

read.csp.profile.old <- function(path) {
  # Reads profile from path.
  #
  # Args:
  #   path: Path to file with the profile. 
  # returns: The profile 
  
  raw = read.table(path, header=T, colClasses='character', sep=',')
  return(raw[,5:dim(raw)[2]])
}
read.ref.profile.old <- function(path) {
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
  # cprofs is an internal format to do calculations on a crime scene profile.
  #
  # Args:
  #   profile: The crime-scene profile to  process.
  # Returns: The same information in a different format.

  # each rep has two rows: alleleic and uncertain
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
  colnames(result)[1:nrep] = sapply(1:nrep, function(n) {paste('Replicate', n)})
  colnames(result)[nrep+1] = 'Total' 
  
  # now add data.
  for(person in row.names(ref)){
	  for(rep in 1:nrep){
      # number of alleles in common for given person and CSP replicate, across
      # all loci
			represented = c()
			for(locus in 1:length(cprofs)){
				if(!is.null(cprofs[[locus]][[rep]])) {
          # figure out unique alleles in reference 
			  	alleles  <- unique(unlist(strsplit(ref[person,locus],',')))
          # figure out howmany of these allele are in CSP
          represented = c(represented, alleles %in% cprofs[[locus]][[rep]]$csp)
		    }
      }
      if(length(represented) > 0) 
        result[person, rep] = round(100*sum(represented)/length(represented))
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
  return(suggested[[1]])
}

# Packs and verifies administrative information.
# Documentation in man directory.
pack.admin.input = function( mixedFile, refFile, caseName='dummy',
                             databaseFile=NULL, outputPath=getwd()
					 ) {
    paths = c(mixedFile, refFile) 
    if(!is.null(databaseFile)) paths = c(databaseFile, paths, recursive=TRUE)
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
  admin = list( caseName=caseName,
                databaseFile=databaseFile,
                mixedFile=mixedFile,
                refFile=refFile,
                outputPath=outputPath )
  return(admin)
}

# Loads allele database
# Documentation is in man directory.
load.allele.database <- function(path=NULL) {
  if(is.null(path)) { # Load default database
    dummyEnv <- new.env()
    data('lgc-allele-freqs-wbp', package='likeLTD', envir=dummyEnv)
    return(dummyEnv[['lgc-allele-freqs-wbp']])
  }
  if(!file.exists(path)) stop(paste(path, "does not exist."))
  read.table(path, sep="\t", header=TRUE)
}

# Generates and packs genetics input into a list.
# Documentation in man directory.
pack.genetics.input = function(admin, nameK=NULL, nameQ=NULL, dropin=FALSE,
                               unknowns=0, ethnic='EA1', fst=NULL, adj=1) {
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
  #   adj: not sure
  #   fst: not sure.  If NULL and ethnic is 'EA1', then defaults to 0.02,
  #        otherwise to 0.03.
  afreq        = load.allele.database(admin$databaseFile)
  cspData      = read.csp.profile.old(admin$mixedFile)
  refData      = read.ref.profile.old(admin$refFile)
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

  genetics = list( # Allele database table.
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
                   # Not sure.
                   adj         = adj,
                   # Number of replicates.
                   nrep        = length(cprofs[[1]]) )
  return(genetics)
}


# Generates an allele summary from input profiles:
# Documentation in man directory.
allele.report <- function(admin) {
  #
  # Args:
  #   admin: List containing administration data, as packed by
  #          pack.admin.input.
  
  # reads genetics information if not given on input.
  if(!require('gplots')) stop("allele.report reqires the package gplots.")
  genetics = pack.genetics.input(admin)


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




stateHypotheses <- function(prosecutionHypothesis, admin)
	{
	# Table with hypotheses described
	#
	# Parameters:
	# 	args: arguments specified by user
	HP = rownames(read.known.profiles(admin$refFile))[which(read.known.profiles(admin$refFile)[,1]==TRUE)] # nameQ
	HD = c('Unknown (X)')
	nameK = rownames(read.known.profiles(admin$refFile))[which(read.known.profiles(admin$refFile)[,1]!=TRUE)]
	if(length(nameK)>0)for (x in 1:length(nameK))
		{
		HP = paste(HP,nameK[x],sep=' + ')
		HD = paste(HD,nameK[x],sep=' + ')
		}
	if(prosecutionHypothesis$nUnknowns>0)for (x in 1:prosecutionHypothesis$nUnknowns)
		{
		HP = paste(HP,paste('unknown (U',x,')',sep=''),sep=' + ')
		HD = paste(HD,paste('unknown (U',x,')',sep=''),sep=' + ')
		}
	if (prosecutionHypothesis$doDropin ==T)
		{
		HP = paste(HP,' + Dropin')
		HD = paste(HD,' + Dropin')
		}
	hypothesis = t(data.frame(HP,HD))
	colnames(hypothesis) = NA
	row.names(hypothesis) = c('Prosecution(HP)','Defence(HD)')
	return(hypothesis)
	}


locus.likes <- function(hypothesis,results,...) 
	{
	# Generate locus likelihoods from overall likelihood
	#
	# Parameters:
	# 	hypothesis: generated by either defence.hypothesis() or 
        #       prosecution.hypothesis()
	#	results: results from do.call(optim,params)
	model <- create.likelihood.vectors(hypothesis)
	arguments <- relistArguments(results$optim$bestmem, hypothesis, ...)
	likes <- do.call(model,arguments)
	likes <- likes$objectives * likes$penalties
	}


calc.dropout = function(results, hypothesis)
	{
	# Calculates dropout rates for every contributor subject to dropout and for every replicate
	#
	# Parameters:
	# 	hypothesis: generated by either defence.hypothesis() or 
        #       prosecution.hypothesis()
	#	results: results from do.call(optim,params)
	N = nrow(hypothesis$dropoutProfs) + hypothesis$nUnknowns + 1
	# Number of contributors subject to dropout + number of unknowns + 1 (dropin)
	nrep = nrow(hypothesis$cspProfile)
	do = results$optim$bestmem[grep("dropout",names(results$optim$bestmem))]
	rcont = results$optim$bestmem[grep("rcont",names(results$optim$bestmem))]
	rcont = rcontConvert(hypothesis$refIndiv,rcont)
	BB = results$optim$bestmem[grep("power",names(results$optim$bestmem))]
	drout = matrix(0,N-1,nrep)
	if(N>1) for(x in 1:(N-1)) for(z in 1:nrep) drout[x,z] = do[z]/(do[z]+rcont[x]^-BB*(1-do[z]))
	return(drout)
	}
	
ideal = function(hypothesis,rr)
	{
	# Calculates idealised likelihood assuming Q is perfect match
	#
	# Parameters:
	# 	hypothesis: generated by either defence.hypothesis() or 
        #       prosecution.hypothesis()
	#	rr: relatedness arguments from args
	ideal.match = 1
	for(j in 1:ncol(hypothesis$cspProfile))
		{
		af = hypothesis$alleleDb[j][[1]]
		kn = hypothesis$queriedProfile[,j][[1]]
		p1 = af[row.names(af)==kn[1],1]
		p2 = af[row.names(af)==kn[2],1]
		ideal.match = ideal.match/(rr[2] + rr[1]*(p1+p2)/2 + (1-sum(rr))*p1*p2*(1+(kn[1]!=kn[2])))
		}
	return(ideal.match)
	}

rcontConvert <- function(refIndiv,rcont) 
	{
	# Convert rcont to full rcont (including ref individual)
	#
	# Parameters:
	# 	refIndiv: reference individual specified in args
	#	rcont: rcont parameters from do.call(optim,params)
	if(refIndiv == 1) rcont = c(1, rcont)
        else if(refIndiv > length(rcont)) rcont = c(rcont, 1)
        else rcont = c(rcont[1:refIndiv-1], 1,
        rcont[refIndiv:length(rcont)])
	return(rcont)
	}



dropDeg <- function(hypothesis,results,admin) 
	{
	dropout <- calc.dropout(results, hypothesis)
	# Output tables for dropout and degradation
	#
	# Parameters:
	# 	hypothesis: generated by either defence.hypothesis() or 
        #       prosecution.hypothesis()
	#	results: results from do.call(optim,params)
	#	dropout: dropout estimates generated by calc.dropout()
	#	args: arguments specified by user

	# 'known' (in Nkdo, Nknd, knownDropoutsLogical) by definition never includes Q (queried)
	dropoutsLogical = determine.dropout(read.known.profiles(admin$refFile),read.csp.profile(admin$mixedFile))
	Qdrop = dropoutsLogical[names(dropoutsLogical)==rownames(hypothesis$queriedProfile)]
	knownDropoutsLogical = dropoutsLogical[names(dropoutsLogical)!=rownames(hypothesis$queriedProfile)]
	Nkdo = length(which(knownDropoutsLogical))
	Nknd = length(which(!knownDropoutsLogical))

	nrep = nrow(hypothesis$cspProfile)
	nameQ = rownames(read.known.profiles(admin$refFile))[which(read.known.profiles(admin$refFile)[,1]==TRUE)]
	nameK = rownames(read.known.profiles(admin$refFile))[which(read.known.profiles(admin$refFile)[,1]!=TRUE)]
	# names (including 'Q','K','U')in correct order 
	Names=c()
	Names[length(c(nameQ,nameK))] = paste(nameQ,'(Q)')
	if(length(nameK)>0)for (n in 1:length(nameK)) Names[n]=paste(nameK[n],' (K',n,')',sep='')
	if(hypothesis$hypothesis=="defence") Names[length(c(nameQ,nameK))] = 'X'

	condition = (hypothesis$nUnknowns>0)&(hypothesis$hypothesis=="prosecution")
	if(condition)for(n in 1:hypothesis$nUnknowns)Names = c(Names,paste('U',n,sep=''))
	condition = (hypothesis$nUnknowns>1)&(hypothesis$hypothesis=="defence")
	if(condition)for(n in 1:(hypothesis$nUnknowns-1))Names = c(Names,paste('U',n,sep=''))
	runNames = c();for(rName in 1:nrep)runNames[rName]=paste('Replicate',rName)

# Table
	h = h1 = round(dropout,4)	
	Loss = Loss1 = round(10^results$optim$bestmem[grep("degradation",names(results$optim$bestmem))],4)
	if(hypothesis$hypothesis=="prosecution")
		{
		if(Qdrop==F)
			{
			h = rbind(h[0:nrow(hypothesis$dropoutProfs),],0)
			if(hypothesis$nUnknowns>0)
				{
				if(nrep>1) h = rbind(h,h1[(nrow(hypothesis$dropoutProfs)+1):length(Loss1),])
				if(!nrep>1) h = rbind(h,as.matrix(h1[(nrow(hypothesis$dropoutProfs)+1):length(Loss1),]))
				}
			Loss = c(Loss[0:nrow(hypothesis$dropoutProfs)],0)
			if(hypothesis$nUnknowns>0)
				{
				Loss = c(Loss,Loss1[(nrow(hypothesis$dropoutProfs)+1):length(Loss1)])
				}
			}
		}

	if(Nknd>0) for(n in 1:Nknd)
		{
		h = rbind(0,h)
		Loss = c(0,Loss)
		}
	Dropout = format(round(data.frame(h,Loss),3),nsmall=3)
	colnames(Dropout)= c(runNames[1:(nrep)],'Degradation:')
	row.names(Dropout) = Names
	return(Dropout)
	}

output.names <- function(admin,prosecutionHypothesis)
	{
	# Names for output files
	#
	# Parameters:
	# 	admin: output from pack.admin.input()
	#	prosecutionHypothesis: generated by prosecution.hypothesis()
	nameQ = rownames(read.known.profiles(admin$refFile))[which(read.known.profiles(admin$refFile)[,1]==TRUE)]
	nameK = rownames(read.known.profiles(admin$refFile))[which(read.known.profiles(admin$refFile)[,1]!=TRUE)]
	nameCode=c()
	for(n in 1:length(c(nameQ,nameK)))
		{
		nameCode=c(nameCode,strsplit(c(nameQ,nameK)[n],'')[[1]][1:2])
		nameCode = paste(nameCode,collapse='')
		}
	version=1
	FN=paste(admin$caseName,nameCode,prosecutionHypothesis$nUnknowns,prosecutionHypothesis$ethnic,prosecutionHypothesis$doDropin,version,sep='-')
	pdf.name = paste(FN,'pdf',sep='.')
	RData.name = paste(FN,'RData',sep='.')
	while(pdf.name %in% list.files(admin$outputPath))
		{
		version=version+1
		FN=paste(admin$caseName,nameCode,prosecutionHypothesis$nUnknowns,prosecutionHypothesis$ethnic,prosecutionHypothesis$doDropin,version,sep='-')
		pdf.name = paste(FN,'pdf',sep='.')
		RData.name = paste(FN,'RData',sep='.')
		}
	return(list(pdf.name=pdf.name,RData.name=RData.name))
	}



output.report <- function(admin,prosecutionHypothesis,defenceHypothesis, prosecutionResults,defenceResults) {
  #
  # Args:
  #   admin: List containing administration data, as packed by pack.admin.input.
  #   prosecutionHypothesis: generated by prosecution.hypothesis()
  #   defenceHypothesis: generated by defence.hypothesis()
  #   prosecutionResults: results from do.call(optim, prosecutionParams)
  #   defenceResults: results from do.call(optim, defenceParams)
  
  # reads genetics information if not given on input.
  if(!require('gplots')) stop("allele.report requires the package gplots.")
  genetics = pack.genetics.input(admin)

  # Unnatributable alleles
  otherBoth = genetics$summary$otherRep + genetics$summary$otherUnrep


  # Generate 'summary', whether multiple Q or single Q
  
  #-----------------------------------------------
  # Write the output report
  # works out sensible font size
  cspTotal = uncTotal = numeric(genetics$nrep)
  for(l in names(genetics$cprofs)){ 
    for(r in 1:length(l)){
      cspTotal[r] = cspTotal[r] + length(genetics$cprofs[[l]][[r]]$csp)
      uncTotal[r] = uncTotal[r] + length(genetics$cprofs[[l]][[r]]$unc)
    }
  }
  
  # Start plotting pdf
  outputName = output.names(admin,prosecutionHypothesis)
  pdf(paste(admin$outputPath,outputName,sep="/"))

  # report pg1

  par(mai = rep(0.3,times=4))
  heights.pg1 = c(2,2*genetics$nrep,(length(c(genetics$nameQ,genetics$nameK))),3,max(otherBoth)+2)
  heights.pg1 = heights.pg1 +5 # space for headers
  widths.pg1 = c(1)
  layout(matrix(c(1:length(heights.pg1)),nrow=length(heights.pg1)),heights = heights.pg1, widths=widths.pg1) # one extra first row is required

  # 1. Hypotheses used
  hypothesisTable <- stateHypotheses(prosecutionHypothesis, admin)
  textplot(hypothesisTable,valign='top')
  title(paste(admin$caseName,'Output Report'))

  # 2. CSP alleles
  alleles <- allele.table(genetics$cprofs)
  textplot(alleles,valign='top')
  title('Crime Scene Profile')

  # 3. Tabular Summary 
  textplot(genetics$summary$summary,valign='top')
  title(main = 'Reference profiles.{}=unreplicated, []=absent')

  # 4. Locus Likelihood
  prosecuLikes <- locus.likes(prosecutionHypothesis,prosecutionResults)
  defenceLikes <- locus.likes(defenceHypothesis,defenceResults)
  likelihood  = data.frame(signif(prosecuLikes,2),signif(defenceLikes,2),signif(prosecuLikes/defenceLikes,2))
  colnames(likelihood)= c('HP','HD','Ratio')
  row.names(likelihood)= colnames(defenceHypothesis$queriedProfile)
  textplot(t(likelihood),valign='top')
  title('Likelihoods at each locus')

  # 5. Unattributable Alleles
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

  # report pg2
  par(mai = rep(0.3,times=4))

  # Number of people to display
  npeep = genetics$unknowns+length(c(genetics$nameK,genetics$nameQ))
  heightspg2 = c(4,npeep,npeep,npeep-genetics$unknowns,2)
  # add titles space  
  heightspg2 = heightspg2+5
  widthspg2 = c(1)
  layout(matrix(c(1:length(heightspg2)),nrow=length(heightspg2)),heights = heightspg2, widths=widthspg2)


  # Overall likelihood
  overallLikelihood  = data.frame(signif(10^-prosecutionResults$optim$bestval,2),signif(10^-defenceResults$optim$bestval,2),signif(10^(defenceResults$optim$bestval-prosecutionResults$optim$bestval),2),signif(defenceResults$optim$bestval-prosecutionResults$optim$bestval,2))
  colnames(overallLikelihood)= c('HP','HD','Ratio','Log Ratio')
  row.names(overallLikelihood)= ''
  textplot(t(overallLikelihood),valign='top')
  title('Likelihood (overall)')

  # Hp dropout
  pDropout = dropDeg(prosecutionHypothesis,prosecutionResults,admin) 
  textplot(pDropout,valign='top')
  title('HP Drop-out')
  # Hd dropout
  dDropout = dropDeg(defenceHypothesis,defenceResults,admin) 
  textplot(dDropout,valign='top')
  title('HD Drop-out')

  # Approximate representation
  textplot(genetics$estimates,valign='top')
  title('Approximate representation %')

  # Dropin
  if(defenceHypothesis$doDropin){
  	pDropin = prosecutionResults$optim$bestmem[grep("dropin",names(prosecutionResults$optim$bestmem))]
  	dDropin = defenceResults$optim$bestmem[grep("dropin",names(defenceResults$optim$bestmem))]
  	dropin = round(data.frame(pDropin,dDropin),3)
  	colnames(dropin) = c('HP','HD')
  	row.names(dropin) = ''
  	textplot(t(dropin),valign='top',cex=1)
  	title('Drop-in (overall)')
  	}

  if(!defenceHypothesis$doDropin){
  	textplot(' ',valign='top')
  	title('No drop-in modelled')
  	}

  # report pg3
  heights.pg3 = heights = c(13,33,16)
  heights.pg3 = heights.pg3 + 5
  layout(matrix(c(1:3),nrow=3), heights = heights.pg3)
  par(mai = rep(0.3,times=4))
  size = 0.7

  # Operator information
  operator = data.frame(Sys.info())
  colnames(operator) = 'Operator info'
  textplot(operator,valign='top',cex=size)
  title('System info')

  # Parameters
  CSPname = tail(strsplit(admin$mixedFile,split='/')[[1]],1)
  REFname = tail(strsplit(admin$refFile,split='/')[[1]],1)
  ALLELEname = ifelse(is.null(admin$database),"lgc-allele-freqs-wbp.txt",tail(strsplit(admin$databaseFile,split='/')[[1]],1))

  dropoutsLogical = determine.dropout(read.known.profiles(admin$refFile),read.csp.profile(admin$mixedFile))
  Nknd = length(which(!dropoutsLogical))
  Nkdo = length(which(dropoutsLogical))

  # End time
  endTime = as.character(Sys.time())

  # Match probability
  ideal.match=ideal(defenceHypothesis,defenceHypothesis$relatedness)

  administration = data.frame(
  	c(
  	admin$caseName,
	CSPname,
	REFname,
	ALLELEname,
	paste(prosecutionHypothesis$ethnic,defenceHypothesis$ethnic),
	endTime,
	paste(prosecutionHypothesis$adj,defenceHypothesis$adj),
	paste(prosecutionHypothesis$fst,defenceHypothesis$fst),
	paste(prosecutionHypothesis$relatedness[1],defenceHypothesis$relatedness[1]),
	paste(prosecutionHypothesis$relatedness[2],defenceHypothesis$relatedness[2]),
	Nknd,
	Nkdo,
	paste(prosecutionHypothesis$nUnknowns,defenceHypothesis$nUnknowns),
	paste(prosecutionHypothesis$doDropin,defenceHypothesis$doDropin),
	paste(format(round(prosecutionResults$optim$bestmem[grep("dropout",names(prosecutionResults$optim$bestmem))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$optim$bestmem[grep("dropout",names(defenceResults$optim$bestmem))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$optim$bestmem[grep("rcont",names(prosecutionResults$optim$bestmem))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$optim$bestmem[grep("rcont",names(defenceResults$optim$bestmem))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$optim$bestmem[grep("Adjustment",names(prosecutionResults$optim$bestmem))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$optim$bestmem[grep("Adjustment",names(defenceResults$optim$bestmem))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$optim$nfeval, defenceResults$optim$nfeval,3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$optim$iter, defenceResults$optim$iter,3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$upper[grep("locusAdjustment",names(prosecutionResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$lower[grep("locusAdjustment",names(prosecutionResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$upper[grep("power",names(prosecutionResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$lower[grep("power",names(prosecutionResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$upper[grep("dropout",names(prosecutionResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$lower[grep("dropout",names(prosecutionResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$upper[grep("degradation",names(prosecutionResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$lower[grep("degradation",names(prosecutionResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$upper[grep("rcont",names(prosecutionResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(prosecutionResults$member$lower[grep("rcont",names(prosecutionResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$upper[grep("locusAdjustment",names(defenceResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$lower[grep("locusAdjustment",names(defenceResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$upper[grep("power",names(defenceResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$lower[grep("power",names(defenceResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$upper[grep("dropout",names(defenceResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$lower[grep("dropout",names(defenceResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$upper[grep("degradation",names(defenceResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$lower[grep("degradation",names(defenceResults$member$lower))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$upper[grep("rcont",names(defenceResults$member$upper))],3),nsmall=3),collapse=' '),
	paste(format(round(defenceResults$member$lower[grep("rcont",names(defenceResults$member$lower))],3),nsmall=3),collapse=' '),
	signif(as.numeric(ideal.match),3),
	row.names=NULL)
	)
		
  colnames(administration) = "Parameters"
  rownames(administration) = c(
	'Case:',
	'CSP file:',
	'Reference file:',
	'Allele frequency file:',
	'EA ethnic group (HP HD):',
	'Time finished:',
	'adj (HP HD):',
	'fst (HP HD):',
	'rr(1) (HP HD):',
	'rr(2) (HP HD):',
	'Nknd:',
	'Nkdo:',
	'No. unknowns (HP HD):',
	'Drop in (HP HD):',
	'Max HP do:',
	'Max HD do:',
	'Max HP rcont:',
	'Max HD rcont:',
	'Max HP locus adjust:',
	'Max HD locus adjust:',
	'NUmber of function evaluations (HP HD)',
	'NUmber of iterations (HP HD)',
	'Upper bounds of locusAdjutment (HP)',
	'Lower bounds of locusAdjutment (HP)',
	'Upper bounds of power (HP)',
	'Lower bounds of power (HP)',
	'Upper bounds of dropout (HP)',
	'Lower bounds of dropout (HP)',
	'Upper bounds of degradation (HP)',
	'Lower bounds of degradation (HP)',
	'Upper bounds of rcont (HP)',
	'Lower bounds of rcont (HP)',
	'Upper bounds of locusAdjutment (HD)',
	'Lower bounds of locusAdjutment (HD)',
	'Upper bounds of power (HD)',
	'Lower bounds of power (HD)',
	'Upper bounds of dropout (HD)',
	'Lower bounds of dropout (HD)',
	'Upper bounds of degradation (HD)',
	'Lower bounds of degradation (HD)',
	'Upper bounds of rcont (HD)',
	'Lower bounds of rcont (HD)',
	'Ideal match:')
  textplot(administration,valign='top',cex=size)
  title('Parameters')

  # Nomenclature
  nomenclature = '
adj:     sampling adjustment parameter

fst:     coancestry parameter (remote shared ancestry of Q and X).

rr:      probabilities of X and Q sharing 1 and 2 alleles through
         shared inheritance from recent ancestors such as parents 
         and grandparents; for full sibs, rr[1]=0.5 and rr[2]=0.25.

Nknd:    number of profiled potential contributors not subject to dropout.

Nkdo:    number of profiled potential contributors subject to dropout.
'
  textplot(nomenclature,valign='top',cex=size)
  title('Nomenclature')

dev.off()
}
