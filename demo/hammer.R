library('gplots')
source('R/allele_report.R')
source('R/genetics.R')
source('R/objectives.R')


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
genetics <- pack.genetics.input(admin)

callmeNth <- function(what, times=10, interTimes=10) {
  empty <- function() { return(0) }
  emptyMean <- mean(replicate(times, system.time(replicate(interTimes, empty()))[3]), trim=0.05)
  result <- replicate(times, system.time(replicate(interTimes, what()))[3])
  return(result - emptyMean)
}

timings <- function(times=10, interTimes=10, nUnknowns=0, doDropin=TRUE) {
  loci = names(genetics$cprofs)
  dropins = rep(doDropin, length(loci))
  nunknowns = rep(nUnknowns, length(loci))
  alleles = list()

  allTimes = list()
  objective <- create.likelihood(admin, nUnknowns=0, doDropin=TRUE)
  objectives <- attr(objective, "loci")
  for(locus in loci) {
    print(paste("Working on", locus))
    newfunc <- function() {
      objectives[[locus]](rcont=c(0.923913043478261, 0.565217391304348,
                                  1.000000000000000, 0.543478260869565,
                                  0.108695652173913), 
                          degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                          localAdjustment=1,
                          tvedebrink=-4.35,
                          dropout=c(0.175, 0.105) )
    }
    allTimes[[locus]] = callmeNth(newfunc, times=times, interTimes=interTimes)
    allTimes[[locus]] = allTimes[[locus]] / interTimes
    alleles[[locus]]  = nrow(attr(objectives[[locus]], "alleleDb"))
  }
  alleles = unlist(alleles)
  allTimes = matrix(unlist(allTimes), nrow=length(allTimes), byrow=T)
  means = apply(allTimes, 1, mean, na.rm=TRUE)
  stdevs = apply(allTimes, 1, sd, na.rm=TRUE)
  matsize = sapply(alleles, function(n) (n*(n+1)/2)**nUnknowns)
  return(data.frame(loci=loci, doDropin=dropins, unknowns=nunknowns,
                    alleles=alleles, matsize=matsize, mean=means,
                    stdev=stdevs))
}


trial <- function(nUnknowns=0, doDropin=TRUE, rcont=NULL) {

  if(is.null(rcont)) 
    rcont = c(0.923913043478261, 0.565217391304348, 1.000000000000000,
              0.543478260869565, 0.108695652173913)
  objective <- create.likelihood.vectors(admin, nUnknowns=0, doDropin=TRUE)
  objective(rcont=rcont, degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
            localAdjustment=1, tvedebrink=-4.35,
            dropout=c(0.175, 0.105), beta=-4.35 )
}

