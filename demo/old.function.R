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

create.original <- function(nUnknowns=1, dropin=FALSE) {
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
  originalEnv$NU = nUnknowns
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
           depa$deg <- c(0.00723217060006922, 0.00569441925951047,
                         0.00216652022387600, 0.00131485405088635)
           depa$deg <- rep(3e-3, 4)
           #depa$deg <- rep(2e-2, 4)

           names(depa$l) = names(cprofs)[1:length(depa$l)]
           names(depa$locadj) = names(cprofs)[1:length(depa$locadj)]
          }, originalEnv )

  print("Done computing evaluation functions.")
  callmeOnces = list()
  for(j in names(originalEnv$cprofs)) {
    args1 = list(originalEnv$de[[j]]$kpdo, originalEnv$depa$do,
                 originalEnv$Drin, originalEnv$de[[j]]$af,
                 originalEnv$de[[j]]$csp, originalEnv$de[[j]]$unc,
                 originalEnv$nrep,
                 c(originalEnv$depa$locadj[[j]], originalEnv$depa$beta),
                 originalEnv$de[[j]]$pUall, originalEnv$NU,
                 originalEnv$depa$rcont, 1+originalEnv$depa$deg,
                 originalEnv$de[[j]]$index, originalEnv$de[[j]]$fragments,
                 originalEnv$de[[j]]$v.index)
    args2 = list(originalEnv$de[[j]]$pUall,
                 originalEnv$de[[j]]$af, originalEnv$NU, originalEnv$rr)
    fuckmelocally <- function(a, b) {
      a = eval(a)
      b = eval(b)
      function() {
        de.tmp= do.call(originalEnv$Calclik.1, a, envir=originalEnv)
        do.call(originalEnv$Adjust.Like, append(list(de.tmp), b), envir=originalEnv)
      }
    }
    callmeOnces[[j]] <- fuckmelocally(args1, args2)
  }
  return(callmeOnces)
}
