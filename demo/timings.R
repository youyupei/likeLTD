source('R/allele_report.R')
source('R/genetics.R')


# Case we are going to be looking at.
caseName = 'hammer'
datapath <- file.path(file.path('inst', 'extdata'), caseName)
# Construct input: frequency file.
databaseFile = NULL #file.path(datapath, 'lgc-allele-freqs-wbp.txt')
# Construct input: crime scene profile
mixedFile = file.path(datapath, 'hammer-CSP.csv')
# Construct input: reference profiles
refFile = file.path(datapath, 'hammer-reference.csv')
# Construct input: output path in the R temp directory for now
outputPath = tempdir()

# Construct list of all administrative input
admin <- pack.admin.input( caseName=caseName,
                           databaseFile=databaseFile,
                           mixedFile=mixedFile,
                           refFile=refFile,
                           outputPath=outputPath )
# Finally call allele.report
genetics <- pack.genetics.input(admin, unknowns=2, dropin=FALSE)

# All known profiles
knownAlleles <- known.alleles( c("Suspect", "Victim 1", "Victim 2"),
                               genetics$refData )
# Add alleles in CSP which are not in database.
alleleDb <- add.missing.alleles( ethnic.database('EA1'), 
                                 genetics$cprofs, knownAlleles[1:2, ] )
# Adjust frequencies depending on case parameters.
alleleDb <- adjust.frequencies(alleleDb, knownAlleles[1:2, ], adj=genetics$adj,
                               fst=genetics$fst)
# Only profiles of people who are not queried and are not subject to dropout.
knownAllelesNoDropouts = known.without.dropouts( c("Victim 1", "Victim 2"),
                                                 genetics$refData,
                                                 genetics$cprofs )
# Only profiles of people who are not queried and are not subject to dropout.
knownAllelesDropouts = known.with.dropouts( c("Victim 1", "Victim 2"),
                                            genetics$refData,
                                            genetics$cprofs )
# Presence matrix of the CSP
cspPresence = presence.matrices(alleleDb, genetics$cprofs)
# Presence matrix of unknown contributors.
if(length(knownAllelesNoDropouts) == 0) {
  uncPresence = presence.matrices(alleleDb, genetics$cprofs, type="unc")
} else {
  uncPresence = presence.matrices(alleleDb, genetics$cprofs, knownAllelesNoDropouts,
                                  type="unc")
}

prof.presence.per.locus <- function(prof, alleleDb) {
  matrix( sapply(prof, function(n) as.integer(rownames(alleleDb) %in% n)), 
          nrow = length(prof), byrow=T )
}
profPresence <- mapply(prof.presence.per.locus, knownAlleles[, names(alleleDb)], alleleDb)


callmeNth <- function(what, times=10, interTimes=10) {
  empty <- function() { return(0) }
  emptyMean <- mean(replicate(times, system.time(replicate(interTimes, empty()))[3]), trim=0.05)
  result <- replicate(times, system.time(replicate(interTimes, what()))[3])
  return(result - emptyMean)
}

timings <- function(times=10, interTimes=10, nunknown=0, dropin=TRUE) {
  loci = names(genetics$cprofs)
  dropins = rep(dropin, length(loci))
  nunknowns = rep(nunknown, length(loci))
  alleles = sapply(loci, function(n) ncol(profPresence[[n]]))
  matsize = sapply(alleles, function(n) (n*(n+1)/2)**nunknown)

  allTimes = list()
  for(locus in loci) {
    missingReps = missing.csp.per.locus(genetics$cprofs[[locus]])
    func <- create.likelihood.per.locus(profPresence[[locus]][1:2, ],
                                        profPresence[[locus]][3:6, ],
                                        cspPresence[[locus]], 
                                        uncPresence[[locus]],
                                        missingReps,
                                        alleleDb[[locus]],
                                        nunknown, dropin)
    newfunc <- function() {
      func(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565,
                           0.108695652173913), 
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )
    }
    allTimes[[locus]] = callmeNth(newfunc, times=times, interTimes=interTimes)
    allTimes[[locus]] = allTimes[[locus]] / interTimes
  }
  allTimes = matrix(unlist(allTimes), nrow=length(allTimes), byrow=T)
  means = apply(allTimes, 1, mean, na.rm=TRUE)
  stdevs = apply(allTimes, 1, sd, na.rm=TRUE)
  return(data.frame(loci=loci, dropin=dropins, unknowns=nunknowns,
                    alleles=alleles, matsize=matsize, mean=means,
                    stdev=stdevs))
}

