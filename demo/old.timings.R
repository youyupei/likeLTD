# Demo for running allele report.

# Load library!
source('R/allele_report.R')
source('R/genetics.R')

# Case we are going to be looking at.
caseName = 'hammer'
datapath <- file.path(file.path('inst', 'extdata'), caseName)
# Construct input: frequency file.
databaseFile = NULL #file.path(datapath, 'lgc-allele-freqs-wbp.txt')
# Construct input: crime scene profile
cspFile = file.path(datapath, 'hammer-CSP.csv')
# Construct input: reference profiles
refFile = file.path(datapath, 'hammer-reference.csv')
# Construct input: output path in the R temp directory for now
outputPath = tempdir()

# Construct list of all administrative input
admin <- pack.admin.input( caseName=caseName,
                           databaseFile=databaseFile,
                           cspFile=mixedFile,
                           refFile=refFile,
                           outputPath=outputPath )
# Finally call allele.report
genetics <- pack.genetics.input(admin, unknowns=2, dropin=FALSE)

create.original <- function(unknowns=1, dropin=FALSE) {
  # Recreates original env
  originalEnv = new.env()
  evalq(source('R/control-functions.R', TRUE), envir=originalEnv)
  originalEnv$afreq <- load.allele.database()
  originalEnv$onlyCSP <- genetics$cspData
  originalEnv$REF <- genetics$refData
  originalEnv$cprofs <- genetics$cprofs
  originalEnv$estimates <- genetics$estimates
  originalEnv$nameK <- genetics$nameK
  originalEnv$nameQ <- genetics$nameQ
  originalEnv$nameQK <- c(genetics$nameQ, genetics$nameK)
  originalEnv$remainingNames = originalEnv$nameQK[originalEnv$nameQK != originalEnv$nameQ]
  originalEnv$nrep <- 2
  originalEnv$nloc <- length(originalEnv$cprofs)
  originalEnv$NU = unknowns
  originalEnv$Drin = dropin
  originalEnv$fst = 0.02
  originalEnv$ethnic = 'EA1'
  originalEnv$adj = 1
  originalEnv$rr = c(0, 0)/4
  originalEnv$known <- FALSE
  originalEnv$Nknd <- FALSE
  originalEnv$Nkdo <- FALSE
  originalEnv$nN <- 0
  originalEnv$Qdrop <- TRUE
  originalEnv$frq <- FALSE
  originalEnv$initiated <- FALSE
  originalEnv$BB <- FALSE
  originalEnv$nRef <- FALSE
  originalEnv$dN <- FALSE
  originalEnv$dRef <- FALSE
  originalEnv$mfunpr <- FALSE
  originalEnv$mfunpr <- FALSE
  originalEnv$otherBoth <- genetics$summary$otherUnrep + 
                           genetics$summary$otherRep
  evalq(library(gtools), envir=originalEnv)
  evalq(global.objects(), envir=originalEnv)
  evalq({
    nu = de = list(); for(j in c(1:10)) {
      print(paste('Currently processing fixed objects for',names(cprofs)[j]))

      # sampling (and other) adjustments, creates a list of objects for each locus (fixed across iterations)
  #   nu[[j]] = prepro(Qdrop,known[[j]][c(1:(2*Nknd+2))],cprofs[[j]],frq[[j]],adj,fst) 
      de[[j]] = prepro(T,known[[j]][c(1:(2*Nknd+2))],cprofs[[j]],frq[[j]],adj,fst) 

      # joins on some more objects that are calculated from the previous ones
  #   nu[[j]] = c(nu[[j]], calc.fixed(Nkdo,NU,nu[[j]]$af,known[[j]][mfunpr],nu[[j]]$csp,Drin,1+Qdrop))
      de[[j]] = c(de[[j]], calc.fixed(Nkdo, NU, de[[j]]$af, known[[j]][mfunpr], de[[j]]$csp, Drin, 0))
      }
      names(de) <- names(cprofs)[1:length(de)]
      }, envir=originalEnv)

  evalq( { start = start.values()
           start$depa$beta=-4.35 
           maxd = depa = start$depa
           depa$deg <- rep(2e-2, 4)
           depa$deg <- c(0.00723217060006922, 0.00569441925951047,
                         0.00216652022387600, 0.00131485405088635)

           names(depa$l) = names(cprofs)[1:length(depa$l)]
           names(depa$locadj) = names(cprofs)[1:length(depa$locadj)]
          }, originalEnv )

  print("Done computing evaluation functions.")
  callmeOnce <- function() {
    evalq({
           de.tmp=Calclik.1(de[[j]]$kpdo, depa$do, Drin, de[[j]]$af, de[[j]]$csp,
                            de[[j]]$unc, nrep, c(depa$locadj[[j]], depa$beta),
                            de[[j]]$pUall, NU, depa$rcont, 1+depa$deg,
                            de[[j]]$index, de[[j]]$fragments, de[[j]]$v.index)
#          depa$l[j] <- Adjust.Like(de.tmp, de[[j]]$pUall, de[[j]]$af, NU, rr)
          }, envir=originalEnv) 
  }
  return(list(env=originalEnv, func=callmeOnce))
}

callmeNth <- function(what, times=10, interTimes=10) {
  empty <- function() { return(0) }
  emptyMean <- mean(replicate(times, system.time(replicate(interTimes, empty()))[3]), trim=0.05)
  result <- replicate(times, system.time(replicate(interTimes, what()))[3])
  return(result - emptyMean)
}

timings <- function(times=10, interTimes=10, nunknown=0, dropin=TRUE) {
  a <- create.original(nunknown, dropin)
  loci = names(a$env$cprofs)
  dropins = rep(a$env$Drin, length(loci))
  nunknowns = rep(a$env$NU, length(loci))
  alleles = sapply(loci, function(n) nrow(a$env$de[[n]]$af))
  matsize = sapply(loci, function(n) nrow(a$env$de[[n]]$pUall))
  print(c("matsize", matsize))
  allTimes = sapply(loci, function(n) {
                      a$env$j <- n
                      callmeNth(a$func, times=times, interTimes=interTimes) 
                     }) 
  allTimes = allTimes / interTimes
  means = apply(allTimes, 2, mean, na.rm=TRUE)
  stdevs = apply(allTimes, 2, sd, na.rm=TRUE)
  return(data.frame(loci=loci, dropin=dropins, unknowns=nunknowns,
                    alleles=alleles, matsize=matsize, mean=means,
                    stdev=stdevs))
}
