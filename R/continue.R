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
# All known profiles
knownAlleles <- known.alleles( c("Suspect", "Victim 1", "Victim 2"),
                               genetics$refData )
# Only profiles of people who are not queried and are not subject to dropout.
knownAllelesNoDropouts = known.without.dropouts( c("Victim 1", "Victim 2"),
                                                 genetics$refData,
                                                 genetics$cprofs )
# Only profiles of people who are not queried and are not subject to dropout.
knownAllelesDropouts = known.with.dropouts( c("Victim 1", "Victim 2"),
                                            genetics$refData,
                                            genetics$cprofs )
# Whether queried profile has dropout.
hasDropouts <- has.dropouts(c("Suspect"), genetics$refData, genetics$cprofs)
# Number of not-queried profiles which are not subject to dropouts.
nbNoDrop <- length(knownAllelesNoDropouts)

# Add alleles in CSP which are not in database.
alleleDb <- add.missing.alleles( ethnic.database('EA1'), 
                                 genetics$cprofs, knownAlleles[1:2, ] )
# Adjust frequencies depending on case parameters.
alleleDb <- adjust.frequencies(alleleDb, knownAlleles[1:2, ], adj=genetics$adj,
                          fst=genetics$fst)
# Presence matrix of the CSP
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

# allProfiles = mapply( function(csp, prof, freq, alleles, n, d) {
#                                 allProfiles(csp, prof, rownames(freq), d, n)
#                               }, 
#                               cspPresence, profPresence, alleleDb,
#                               MoreArgs=list(n=genetics$dropin, d=genetics$unknowns) )

# fragments = mapply( function(freq, perms) t(matrix(freq[perms,2],nrow=nrow(perms))), 
#                     alleleDb, allProfiles )

indices <- rep( 1: max(1, genetics$unknowns), 
                rep(2, max(1, genetics$unknowns)) ) + nrow(knownAlleles) / 2 

degradation <- rep(1.003, nrow(knownAlleles) / 2 + genetics$unknowns)
rcont <- rep(0.1, nrow(knownAlleles) / 2 + genetics$unknowns)


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
originalEnv$NU = 1
originalEnv$Drin = TRUE
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
  nu = de = list(); for(j in 7:7) { #c(3:7, 9:10)) {
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

         j = 7
         de.tmp=Calclik.1(de[[j]]$kpdo, depa$do, Drin, de[[j]]$af, de[[j]]$csp,
                          de[[j]]$unc, nrep, c(depa$locadj[[j]], depa$beta),
                          de[[j]]$pUall, NU, depa$rcont, 1+depa$deg,
                          de[[j]]$index, de[[j]]$fragments, de[[j]]$v.index)
         depa$l[j] <- Adjust.Like(de.tmp, de[[j]]$pUall, de[[j]]$af, NU, rr)
       }, envir=originalEnv)

get_permorder <- function(oldPerms, newPerms) {
  # Figure out mapping from newPerms to oldPerms
  # 
  # The return should be such that newPerms[return, ] == oldPerms. 
  check <- function(i) which(apply(newPerms, 1, function(n)  all(i == n))) 
  return(apply(oldPerms, 1, check))
} 

missingReps      = missing.csp.per.locus(genetics$cprofs$D18)
cspPresence      = cspPresence$D18
uncPresence      = uncPresence$D18
queriedPresence  = profPresence$D18[1:2, ]
profPresence     = profPresence$D18[3:6, ]
alleleDb         = alleleDb$D18
genetics$dropin = TRUE
genetics$unknowns = 0
objective.function <- create.likelihood.per.locus(queriedPresence,
                                                  profPresence, cspPresence,
                                                  uncPresence, missingReps,
                                                  alleleDb,
                                                  genetics$unknowns,
                                                  genetics$dropin)
arguments = list(rcont           = originalEnv$depa$rcont,
                 degradation     = originalEnv$depa$deg,
                 localAdjustment = originalEnv$depa$locadj[7],
                 tvedebrink      = originalEnv$BB,
                 dropout         = originalEnv$depa$do )
newRes <- do.call(objective.function, arguments)
