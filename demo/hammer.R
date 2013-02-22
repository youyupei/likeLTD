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
  objective <- create.likelihood.vectors(admin, nUnknowns=nUnknowns, doDropin=doDropin)
  objective(rcont=rcont, degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
            localAdjustment=1, tvedebrink=-4.35,
            dropout=c(0.175, 0.105), beta=-4.35 )
}

args = list(rcont=c(0.923913043478261, 0.565217391304348, 1.000000000000000,
                    0.543478260869565, 0.108695652173913), 
            degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
            localAdjustment=1, tvedebrink=-4.35,
            dropout=c(0.175, 0.105), beta=-4.35 )


plot_rcont <- function(which=1, x=(1:99)/100.0, nUnknowns=0, doDropin=TRUE, ...) {
  require("ggplot2")
  objective = create.likelihood.vectors(admin, nUnknowns=nUnknowns, doDropin=doDropin)
  funcme <- function(i) { args$rcont[which] = x[i]; do.call(objective, args)$objectives}
  y = sapply(1:length(x), funcme) 
  total = array(1, ncol(y))
  for(i in 1:nrow(y)) total = total * y[i, ]
  y = rbind(y, total=total)
  data = data.frame(prob=c(y), locus=rep(rownames(y), ncol(y)),
                    rcont=rep(x, rep(nrow(y), length(x))))
  return(ggplot(data, aes(x=rcont, y=prob, group=locus, colour=locus)) + geom_line())
}

plot_deg <- function(which=1, range=c(1e-4, 1e-2), n=100, nUnknowns=0, doDropin=TRUE, ...) {
  objective = create.likelihood(admin, nUnknowns=nUnknowns, doDropin=doDropin)
  funcme <- function(i) { args$degradation[which] = x[i]; do.call(objective, args) }
  x = 1:n/n * (max(range) - min(range)) + min(range)
  y = sapply(1:length(x), funcme)
  plot(x, y, ...)
}

modify.args = function(index, newvalue, args) { 
  for(n in names(args)) {
    if(index == 1 & length(args[[n]]) == 1) { args[[n]] = newvalue; break }
    else if( index < length(args[[n]]) ) { args[[n]][index] = newvalue; break }
    else index = index - length(args[[n]])
  }
  return(args)
}

plot_all.log <- function(which=1, x=(1:99)/100.0, nUnknowns=0, doDropin=TRUE) {
  require("ggplot2")
  objective = create.likelihood.vectors(admin, nUnknowns=nUnknowns, doDropin=doDropin)
  funcme <- function(i) { 
    do.call(objective, modify.args(which, x[i], args))$objectives
  }
  y = log(sapply(1:length(x), funcme) )
  total = colSums(y)
  total = total - min(total) 
  y = rbind(y, total=total)
  data = data.frame(prob=c(y), locus=rep(rownames(y), ncol(y)),
                    rcont=rep(x, rep(nrow(y), length(x))))
  return(ggplot(data, aes(x=rcont, y=prob, group=locus, colour=locus)) + geom_line())
}

plot_all2d.log <- function(which=c(1, 2), x=(1:99)/100.0, y=(1:99)/100.0, nUnknowns=0, doDropin=TRUE) {
  require("ggplot2")
  objective = create.likelihood.log(admin, nUnknowns=nUnknowns, doDropin=doDropin)
  
  result = matrix(1, nrow=length(x), ncol=length(y))
  for(i in 1:length(x)) 
    for(j in 1:length(y))
      result[i, j] = do.call(objective, 
                             modify.args(which[2], y[j],
                                         modify.args(which[1], x[i], args)))
  amap = data.frame(x=rep(x, length(y)), y=rep(y, rep(length(x), length(y))), z=c(result))
  return(ggplot(amap, aes(x, y, z=z)))
}
# plot_all2d.log(which=c(1,3), nUnknowns=1, x=(1:20/10.0), y=(1:20/10.0)) + 
#          geom_tile(aes(fill=z))                                         + 
#          stat_contour()
