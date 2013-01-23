UnusualAlleles <- function(frequencies, profile) {
  # Searches for unusual alleles.
  #
  # Such alleles are likely to be mistakes.
  #
  # Args:
  #   freqencies: Table of allele frequencies.“
  #   profile: A profile with just the facts. 
  #
  # returns: data frame with the rare alleles. 

  rare =data.frame(name=NULL, locus=NULL, allele=NULL, EA1=NULL, EA3=NULL,
                   EA4=NULL)

  for(name in row.names(profile)){
    for(locus in colnames(profile)){
      # loop over unique allele names.
      unique.alleles <- unique(strsplit(profile[name,locus],',')[[1]]) 
      for(allele in unique.alleles){

        condition = frequencies$Marker==locus & frequencies$Allele==allele
        x = frequencies[condition, ]
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

ReadCSPProfile <- function(path) {
  # Reads profile from path.
  #
  # Args:
  #   path: Path to file with the profile. 
  # returns: The profile 
  
  raw = read.table(path, header=T, colClasses='character', sep=',')
  return(raw[,5:dim(raw)[2]])
}
ReadRefProfile <- function(path) {
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
QueriedAndKnown <- function(path) {
  # Reads profile from path and returns queried vs known column.
  #
  # Args:
  #   path: Path to file with the profile. 
  raw = read.table(path, header=T, colClasses='character', row.names=1,
                   sep=',', quote = "\"")
  return(raw$known.queried == 'queried')
}

InternalRepresentation <- function(profile) {
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

EstimateCSP <- function(ref, cprofs) {
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

Private_NullifyItem <- function(n) {
  # Returns NULL if n is character(0), else returns n.
  #
  # This is a helper function for internal use only.
  if(identical(n, character(0))) return(NULL)
  return(n)
}

FlattenCSP <- function(v, subitem='csp') {
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
AllCSP <- function(cprofs, subitem='csp') {
  # Flatten and filters elements of cprofs into a list of vectors.
  # 
  # This function is mostly for internal use. It applies to each element of
  # cprofs the function FlattenCSP.
  #
  # Args:
  #   v: A locus element from cprofs. 
  # First, get all results across CSP
  result <- sapply(cprofs, function(n) FlattenCSP(n, subitem))
  # Then nullifies elements evaluating to character 0. 
  return(sapply(result, Private_NullifyItem))
}
ReplicatedAlleles <- function(cprofs) {
  # For each locus, lists replicated alleles.
  #  
  # Args:
  #  v: cprofs data object
  # Returns:
  #   A list over loci, where each subitem is a vector of replicated alleles.
  allcsp <- AllCSP(cprofs, subitem='csp')
  allcsp <- sapply(allcsp, function(n) n[duplicated(n)])
  return(sapply(allcsp, Private_NullifyItem))
}
UnreplicatedAlleles <- function(cprofs) {
  # For each locus, lists replicated alleles.
  #  
  # Args:
  #  v: cprofs data object
  # Returns: 
  #   A list over loci, where each subitem is a vector of unreplicated alleles.
  allcsp <- AllCSP(cprofs, subitem='csp')
  allcsp <- sapply(allcsp, function(n) setdiff(n, n[duplicated(n)]))
  return(sapply(allcsp, Private_NullifyItem))
}

AlleleTable <- function(cprofs) {
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

SummaryGenerator = function(queried, known, ref, cprofs){
  # Summnarizes the allele table for Q and K
  # 
  # Args:
  #  queried: Names of queried individuals in reference profile (rows).
  #  known: Names of known individuals in reference profile (rows).
  #  ref: the reference profile
  #  cprofs: the crime scene profile, external format
  all.names = c( sapply(queried, function(n) paste(n, '(Q)'), USE.NAMES=FALSE),
                 sapply(known, function(n) paste(n, '(K)'), USE.NAMES=FALSE) )
  ncont = length(all.names)
  # index: indexing of ref such that queries come first and known profiles come
  # second.
  index = sapply(c(queried, known), function(n) which(row.names(ref) == n))
  # known.loci: list over locus where each item contains the alleles with the
  # same indexing as index. 
  known.loci = list()
  for(locus in names(ref)) {
    string = paste(ref[[locus]][index], collapse=',')
    known.loci[[locus]] = unlist(strsplit(string, ',')) 
  }

  # Figures out all, replicated, and unreplicated alleles in crime scene.
  csp.combined.all <- AllCSP(cprofs)
  replicated.all   <- ReplicatedAlleles(cprofs) 
  unreplicated.all <- UnreplicatedAlleles(cprofs) 

  #  construct summary data frame scaffold.
  summary = data.frame(array(0, c(length(all.names)+1, length(ref))))
  row.names(summary) = c(all.names,'Unattributable')
  colnames(summary)  = names(cprofs)

  # other.rep: per locus, number of unattributable and replicated
  # alleles
  other.rep = list()
  # other.unrep: per locus, number of unattributable and unreplicated
  # alleles
  other.unrep = list()
  for(locus in names(ref)){ 
    # loop over people
    for(c in 1:length(all.names)){
      # alleles for this person in reference profile
      cont = known.loci[[locus]][(2*c-1):(2*c)]
      rep   = paste(cont[ cont %in% replicated.all[[locus]]],   collapse=' ')
      unrep = paste(cont[ cont %in% unreplicated.all[[locus]]], collapse=' ')
      none  = paste(cont[!cont %in% csp.combined.all[[locus]]], collapse=' ')
      summary[c,locus] = paste(rep, '{', unrep,'}','[', none,']',sep='')
    }
    # other alleles not found in the cprofs
    notfound = !csp.combined.all[[locus]] %in% known.loci[[locus]]
    unattr = csp.combined.all[[locus]][notfound]
    unattr.rep = unique(unattr[duplicated(unattr)])
    unattr.rep.string = paste(unattr.rep, collapse=' ')
    unattr.unrep = unattr[!unattr %in% unattr.rep]
    unattr.unrep.string = paste(unattr.unrep, collapse=' ')

    other.rep[locus] = length(unattr.rep)
    other.unrep[locus] = length(unattr.unrep)

    summary[nrow(summary), locus] = paste( unattr.rep.string, '{',
                                           unattr.unrep.string, '}',
                                           sep='' )
  }

  # coerce to numerics
  other.rep = as.numeric(other.rep)
  other.unrep = as.numeric(other.unrep)
  return( list(summary=summary, other.rep=other.rep, other.unrep=other.unrep) )
}

HypothesisGenerator = function(queried, known, other.rep, other.unrep){
  # function to calculate possible hypotheses
  # 
  # Needs to be done with various permutations of 2 Q's

  other.both = other.rep + other.unrep
  min.u = ceiling(max(other.rep)/2)
  # TODO: Not sure the next line should have "divide by 2". Might be over
  # replicated only.
  max.u = ceiling(max(other.both)/2) 

  hypotheses = c()
  for(unknowns in min.u:max.u) {

    extras   = other.both - 2*unknowns
    unattributable.alleles = sum(extras[extras>0])
    drop.in = unattributable.alleles != 0

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
    if(unattributable.alleles !=0) {
      D = paste(unattributable.alleles, 'dropins', sep='')
      others = c(others, D, recursive=TRUE)
    }
         
    others <- paste(others, collapse='+')
    string <- paste(Q, '+', others,' vs X+', others, sep='') 
    hypotheses[unknowns-min.u+1] = string
  }
  return(hypotheses[!is.na(hypotheses)])
}


SuggestHypothesis = function(queried, known, ref, cprofs) {
  # Figures out sensible explanation of crime profile.
  #
  # Returns: List of strings with sensible hypothesis. 

  # estimates all possible hypotheses
  summary = SummaryGenerator(queried, known, ref, cprofs)
  other.both = summary$other.rep + summary$other.unrep
  other.rep  = summary$other.rep

  # Each element of hypothesis is a different, hum, hypothesis.
  hypothesis = list(list(queried, known))
  if(length(queried)==2)
    hypothesis = list( hypothesis,
                       list(queried[2], known), 
                       list(queried[1], c(queried[2], known)), 
                       list(queried[2], c(queried[1], known)) )
   
  GenerateHypothesis = function(n) {
    #  Generates a hypothesis
    #  Ugly way to call badly designed functions.
    summary = SummaryGenerator(n[[1]], n[[2]], ref, cprofs)
    result = HypothesisGenerator(n[[1]], n[[2]], summary$other.rep,
                                 summary$other.unrep)
    return(result)
  }
  # suggested: The hypothesis generator is applied to each hypothesis. Each
  # element of suggested is a string detailing the hypothesis.
  suggested = sapply(hypothesis, GenerateHypothesis, USE.NAMES=FALSE,
                     simplify=FALSE)
  return(suggested)
}

AlleleReport <- function( frequency.file, mixed.file, ref.file,
                          case.name='dummy', output.path=NA ) {
  # Generates an allele summary from input profiles:
  #
  # Args:
  #   frequency.file: File holding the allele frequencies
  #   mixed.file: Mixed crime scene profile
  #   ref.file: Reference profiles
  #   casename: Name of the current case
  #   output.path: Path where the output should be stored.
  


  # Formats the output path
  if(identical(output.path, NA)) {
    startDir = getwd()
    output.path = file.path(startDir, case.name)
  }

  # Reads raw crime scene profile from  file
  csp.data = ReadCSPProfile(mixed.file)

  # Reads raw reference profile from  file
  ref.data = ReadRefProfile(ref.file)


  # Transform data to an internal represenation
  cprofs <- InternalRepresentation(csp.data)

  # Extracts Queried and Knowns
  queried.vs.known <- QueriedAndKnown(ref.file)
  names.queried = row.names(ref.data)[queried.vs.known]
  names.known   = row.names(ref.data)[!queried.vs.known]



  # Generate 'summary', whether multiple Q or single Q

  #-----------------------------------------------
  # 12. Write the allele report
  # works out sensible font size
  csp.total = unc.total = numeric(nrep)
  for(l in names(cprofs)){ 
    for(r in 1:length(l)){
      csp.total[r] = csp.total[r] + length(cprofs[[l]][[r]]$csp)
      unc.total[r] = unc.total[r] + length(cprofs[[l]][[r]]$unc)
    }
  }
  tablesize = min(1,25/max(c(csp.total,unc.total)))

  case = 'demo/hammer'
  stain = basename(case)
  output.path = file.path(case, paste(stain, '-allele report.pdf', sep=''))
  
  # Start plotting proper. 
  pdf(width=8.25,height=11.5,file=output.path)
  par(mfrow=c(6,1), mai=c(0.7,0.5,0.5,0.5))

  # Allele report
  alleles <- AlleleTable(cprofs)
  textplot(alleles, valign='top',cex=tablesize)
  title(paste(stain,'Allele Report'))

  # Tabular Summary 
  summary = SummaryGenerator(names.queried, names.known, ref.data, cprofs)
  textplot(summary$summary, valign='top', cex=tablesize)
  title('Summary. {}=unreplicated, []=absent')

  # Unattributable Alleles
  other.both = summary$other.rep + summary$other.unrep
  other.rep  = summary$other.rep
  yaxp=c(0, max(other.both), max(other.both))
  if(max(other.both)==0) yaxp=c(0,1,1)
  barplot(other.both, names.arg=names(cprofs), yaxp=yaxp, main='Unattributable
          alleles', ylab='No. alleles')
  barplot(add=T, col='black', other.rep, yaxt='n')
  if(max(other.both) > 3) abline(h=2, lty=2)
  if(max(other.both) > 5) abline(h=4, lty=2)
  if(max(other.both) > 7) abline(h=6, lty=2) 
  if(max(other.both) > 9) abline(h=8, lty=2);
  legend(x=0, y=max(other.both), yjust=0, bty='n',
         c('replicated','unreplicated'), xpd=NA, col=c('black','grey'),
         pch=c(15,15))

  # Rare Alleles
  afreq = read.table(frequency.file,sep="\t",header=T)
  rare <- UnusualAlleles(afreq, ref.data)
  if(length(rare)==0) rare = 'No unusual alleles in Q or K profiles'
  textplot(rare, cex=1, valign='top')
  title('Rare alleles')

  # Approximae representation.
  # Estimate how well each reference profile is represented in the CSP
  estimates <- EstimateCSP(ref.data, cprofs)
  textplot(estimates, cex=1, valign='top')
  title('Approximate representation %')

  # Suggested hypothesis.
  suggested = SuggestHypothesis(names.queried, names.known, ref.data, cprofs)
  textplot(suggested, cex=1, valign='top')
  title(main='Suggested Hypotheses',
        sub= paste('Note: likeLTD is limited to 2 unknowns. If more than 2 ', 
                   'unknowns are suggested,this may indicate a mixed profile',
                   'with too many alleles to be statisically informative.',
                   sep=''))
  dev.off()
  save.image(file = paste(case,'inputs.RData',sep='/'))
}
