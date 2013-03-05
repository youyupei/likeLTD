#-------------------------------------------------------------------------
# Main control script
# Calls all required functions and scripts
# Runs main loop through all iterations
#-------------------------------------------------------------------------
# 1. Runs  scripts containing all required functions
# source('R/control-functions.R') 
# source('output-report.R')

afreq <- load.allele.database()
onlyCSP <- genetics$cspData
REF <- genetics$refData
cprofs <- genetics$cprofs
estimates <- genetics$estimates
nameK <- genetics$nameK
nameQ <- genetics$nameQ
nameQK <- c(genetics$nameQ, genetics$nameK)
remainingNames = nameQK[nameQK != nameQ]
nrep <- 2
nloc <- length(cprofs)
NU = 1
Drin = TRUE
fst = 0.02
ethnic = 'EA1'
adj = 1
rr = c(0.5, 0.75)/4
known <- FALSE
Nknd <- FALSE
Nkdo <- FALSE
nN <- 0
Qdrop <- TRUE
frq <- FALSE
initiated <- FALSE
BB <- FALSE
nRef <- FALSE
dN <- FALSE
dRef <- FALSE
mfunpr <- FALSE
mfunpr <- FALSE
if(Drin) DI <- 1 else DI = 0
otherBoth <- genetics$summary$otherUnrep + genetics$summary$otherRep
#-------------------------------------------------------------------------
# 2. Global objects created in workspace. Quick and dirty solution, ideally all objects that are constant across iterations and loci should be values of a single object.
# Runs GUI script for user defined case specific inputs, which creates most global objects. NOTE: requires 'allele.report.R' to have previously been run. Work is required to incorporate correct file path for allele.report.R if it hasn't yet been run, probably easier to add an 'optimised' version.
# Calculates other required global objects. 
# source('GUI-optimised.R',echo=T)
global.objects()

#-------------------------------------------------------------------------
# 3. Objects that are constant across iterations, but vary across loci, therefore one object per locus is generated.
nu = de = list(); for(j in 1:nloc){ 
  print(paste('Currently processing fixed objects for',names(cprofs)[j]))

  # sampling (and other) adjustments, creates a list of objects for each locus (fixed across iterations)
  nu[[j]] = prepro(Qdrop,known[[j]][c(1:(2*Nknd+2))],cprofs[[j]],frq[[j]],adj,fst) 
  de[[j]] = prepro(T,known[[j]][c(1:(2*Nknd+2))],cprofs[[j]],frq[[j]],adj,fst) 

  # joins on some more objects that are calculated from the previous ones
  nu[[j]] = c(nu[[j]], calc.fixed(Nkdo,NU,nu[[j]]$af,known[[j]][mfunpr],nu[[j]]$csp,rr,1+Qdrop,Drin))
  de[[j]] = c(de[[j]], calc.fixed(Nkdo,NU+1,de[[j]]$af,known[[j]][mfunpr],de[[j]]$csp,rr,0,Drin))
  }

#-------------------------------------------------------------------------
# 4. Initialises start values for nupa and depa
start = start.values()
start$nupa$beta=-4.35 # = append(start$nupa, list(beta=-4.35))
start$depa$beta=-4.35 # = append(start$depa, list(beta=-4.35))
maxn = nupa = start$nupa
maxd = depa = start$depa


#-------------------------------------------------------------------------
# 5. Main loop.
shrink <<- 1

niter = 1
for(itr in 1:niter){ # start simulated annealing loop
	
	# copy of inital objects before updates are proposed, in case updates are rejected
	nupa.old = nupa
	depa.old = depa


	new = propose.new(nupa,depa,itr) # convenient here for a simple loop, but note this slightly changes the initial values before they are used by Calclik.1()
	nupa = new$nupa
	depa = new$depa

# callmenow <- function() {
depa$l <- array(0, nloc)

# CHANGE: this loop
for(j in 1:nloc){ # for each locus
  nu.tmp=Calclik.1(nu[[j]]$kpdo,nupa$do,rev(nupa$rcont)[1],nu[[j]]$af,nu[[j]]$csp,nu[[j]]$unc,nrep,c(nupa$locadj[[j]],nupa$beta),nu[[j]]$pUall,NU,nupa$rcont,1+nupa$deg,nu[[j]]$index,nu[[j]]$fragments,nu[[j]]$v.index)
  if(NU!=0) {
  	nupa$l[j] <- Adjust.Like(nu.tmp,nu[[j]]$pUall,nu[[j]]$af,NU,rr)
  	} 
  else {
  	if(Drin==0) nupa$l[j] <- nu.tmp
  	else nupa$l[j] <- Adjust.Like(nu.tmp,nu[[j]]$pUall,nu[[j]]$af,NU,rr)
  	}
  de.tmp=Calclik.1(de[[j]]$kpdo,depa$do,rev(depa$rcont)[1],de[[j]]$af,de[[j]]$csp,de[[j]]$unc,nrep,c(depa$locadj[[j]],depa$beta),de[[j]]$pUall,NU+1,depa$rcont,1+depa$deg,de[[j]]$index,de[[j]]$fragments,de[[j]]$v.index)
  depa$l[j] <- Adjust.Like(de.tmp,de[[j]]$pUall,de[[j]]$af,NU+1,rr)
  
  
#   Adjustment before relatedness (stays the same if rr=c(0,0))
  depa$l[j] = depa$l[j] * (1-sum(rr))
# Relatedness between Q and X only taken into account in defense
# N.B. de[[j]]$pUall.rel[[i]],   de[[j]]$fragments.rel[[i]],    de[[j]]$v.index.rel[[i]]
# Also flag at end is 1
	if(rr[1]>0) {
		for (i in 1:2) {
			de.tmp = Calclik.1(de[[j]]$kpdo,depa$do,rev(depa$rcont)[1],de[[j]]$af,de[[j]]$csp,de[[j]]$unc,nrep,c(depa$locadj[[j]],depa$beta),de[[j]]$pUall.rel[[i]],NU+1,depa$rcont,1+depa$deg,de[[j]]$index,de[[j]]$fragments.rel[[i]],de[[j]]$v.index.rel[[i]])
			depa$l[j] = depa$l[j] + Adjust.Like(de.tmp,de[[j]]$pUall.rel[[i]],de[[j]]$af,NU+1,rr,1)
			}
		}
# N.B. 	de[[j]]$pUall.rel[[3]],   de[[j]]$fragments.rel[[3]],    de[[j]]$v.index.rel[[3]]
# Flag at end is 2
  if(rr[2]>0) {
  	de.tmp = Calclik.1(de[[j]]$kpdo,depa$do,rev(depa$rcont)[1],de[[j]]$af,de[[j]]$csp,de[[j]]$unc,nrep,c(depa$locadj[[j]],depa$beta),de[[j]]$pUall.rel[[3]],NU+1,depa$rcont,1+depa$deg,de[[j]]$index,de[[j]]$fragments.rel[[3]],de[[j]]$v.index.rel[[3]])
  	depa$l[j] = depa$l[j] + Adjust.Like(de.tmp,de[[j]]$pUall.rel[[3]],de[[j]]$af,NU+1,rr,2)
  	}		
}
# CHANGE: Moved these two out of loop
	nupa$L = calc.L(nupa,nN)
	depa$L = calc.L(depa,dN)
# } # callmenow
# callmenow()
  names(depa$l) <- names(cprofs)

# # decide whether to accept proposed updates (separately for pros and def)
# decision = decide(depa.old,nupa.old,depa,nupa,itr)
# nupa = decision$nupa
# depa = decision$depa



	# Generate reports (interim or final)
#	if((itr)%in%interims | itr==niter) outputs() 
	}

#----------------------------------------------------------------------
