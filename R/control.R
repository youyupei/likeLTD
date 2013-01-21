#-------------------------------------------------------------------------
# Main control script
# Calls all required functions and scripts
# Runs main loop through all iterations
#-------------------------------------------------------------------------
# 1. Runs  scripts containing all required functions
source('control-functions.R') 
source('output-report.R')

#-------------------------------------------------------------------------
# 2. Global objects created in workspace. Quick and dirty solution, ideally all objects that are constant across iterations and loci should be values of a single object.
# Runs GUI script for user defined case specific inputs, which creates most global objects. NOTE: requires 'allele.report.R' to have previously been run. Work is required to incorporate correct file path for allele.report.R if it hasn't yet been run, probably easier to add an 'optimised' version.
# Calculates other required global objects. 
source('GUI-optimised.R',echo=T)
global.objects()

#-------------------------------------------------------------------------
# 3. Objects that are constant across iterations, but vary across loci, therefore one object per locus is generated.
nu = de = list(); for(j in 1:nloc){ 
	print(paste('Currently processing fixed objects for',names(cprofs)[j]))

	# sampling (and other) adjustments, creates a list of objects for each locus (fixed across iterations)
	nu[[j]] = prepro(Qdrop,known[[j]][c(1:(2*Nknd+2))],cprofs[[j]],frq[[j]],adj,fst) 
	de[[j]] = prepro(T,known[[j]][c(1:(2*Nknd+2))],cprofs[[j]],frq[[j]],adj,fst) 

	# joins on some more objects that are calculated from the previous ones
	nu[[j]] = c(nu[[j]], calc.fixed(Nkdo,NU,nu[[j]]$af,known[[j]][mfunpr],nu[[j]]$csp,Drin,1+Qdrop))
	de[[j]] = c(de[[j]], calc.fixed(Nkdo,NU+1,de[[j]]$af,known[[j]][mfunpr],de[[j]]$csp,Drin,0))
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
for(itr in 1:niter){ # start simulated annealing loop
	
	if(itr==1) time.begin <- proc.time()[3]


	# copy of inital objects before updates are proposed, in case updates are rejected
	nupa.old = nupa
	depa.old = depa


	new = propose.new(nupa,depa,itr) # convenient here for a simple loop, but note this slightly changes the initial values before they are used by Calclik.1()
	nupa = new$nupa
	depa = new$depa


# CHANGE: this loop
for(j in 1:nloc){ # for each locus
	#if(DI*max(1-DO)*max(afbp[setdiff(rownames(afbp)[apply(unc,2,prod)==0],kpdo),1])>1) return(0) # don't allow dropin rate in any replicate to exceed 1/max(af)
if(NU!=0) {
	if(Drin*max(1-nupa$do)*max(nu[[j]]$af[setdiff(rownames(nu[[j]]$af)[apply(nu[[j]]$unc,2,prod)==0],known[[j]][mfunpr]),1])>1) {nupa$l[j]=0} else {
		nu.tmp=Calclik.1(known[[j]][mfunpr],nupa$do,Drin,nu[[j]]$af,nu[[j]]$csp,nu[[j]]$unc,nrep,c(nupa$locadj[[j]],nupa$beta),nu[[j]]$pUall,NU,nupa$rcont,1+nupa$deg,nu[[j]]$tprof,nu[[j]]$index,nu[[j]]$fragments,nu[[j]]$v.index) 
		nupa$l[j] = Adjust.Like(nu.tmp,nu[[j]]$pUall,nu[[j]]$af,NU,rr)
		}
	} else {
	if(Drin*max(1-nupa$do)*max(nu[[j]]$af[setdiff(rownames(nu[[j]]$af)[apply(nu[[j]]$unc,2,prod)==0],known[[j]][mfunpr]),1])>1) {nupa$l[j]=0} else {
		nupa$l[j] = zero.cont(Drin,nupa$do,nu[[j]]$csp,nu[[j]]$hyptadinit*nupa$beta,nu[[j]]$af[,1],nu[[j]]$unc,nrep,c(nupa$locadj[j],nupa$beta))
		}
	}

	if(Drin*max(1-depa$do)*max(de[[j]]$af[setdiff(rownames(de[[j]]$af)[apply(de[[j]]$unc,2,prod)==0],known[[j]][mfunpr]),1])>1) {depa$l[j]=0} else {
		de.tmp=Calclik.1(known[[j]][mfunpr],depa$do,Drin,de[[j]]$af,de[[j]]$csp,de[[j]]$unc,nrep,c(depa$locadj[[j]],depa$beta),de[[j]]$pUall,NU+1,depa$rcont,1+depa$deg,de[[j]]$tprof,de[[j]]$index,de[[j]]$fragments,de[[j]]$v.index)
		depa$l[j] = Adjust.Like(de.tmp,de[[j]]$pUall,nu[[j]]$af,NU+1,rr)
	}
}
# CHANGE: Moved these two out of loop
	nupa$L = calc.L(nupa,nN)
	depa$L = calc.L(depa,dN)


	# decide whether to accept proposed updates (separately for pros and def)
	decision = decide(depa.old,nupa.old,depa,nupa,itr)
	nupa = decision$nupa
	depa = decision$depa


if(itr==1) {time.end <- proc.time()[3]; print(paste("Approx run time =", round(((time.end-time.begin)*niter)/60), "minutes"))}


	# Generate reports (interim or final)
	if((itr)%in%interims | itr==niter) outputs() 
	}

#----------------------------------------------------------------------
