#------------------------------------------
# script in development for new docx output format
# should contain all revised functions to completely replace the 'reports' script
# IMPORTANT:
# this new script (specifically the allele.report() and outpout.report() functions are dependent on R2DOC and R2DOCX
# however, these packages are not available on the CRAN (yet?), so instead we have removed them from the 
# 'Depends' list in 'DESCRIPTION', and instead coded their installation here:
if(!require(R2DOCX)){
	install.packages(c("devtools","rJava"),repos="http://cran.ma.imperial.ac.uk/")
	devtools::install_github('R2DOC', 'davidgohel') # installs from the github
	devtools::install_github('R2DOCX', 'davidgohel')  # installs from the github
	Sys.setenv(NOAWT=1) #prevents usage of awt - required on Mac
	}

#------------------------------------------
latex.cleaner <- function(text){
# cleans up the text produced by any of the ...to.latex() functions below
	text <- gsub('{\\bf}','',text,fixed=T)
	text <- gsub(',\\fx{}','',text,fixed=T)
	text <- gsub('{\\em}','',text,fixed=T)
	text <- gsub(',,',',',text,fixed=T)
	text <- gsub(',&','&',text,fixed=T)
	text <- gsub('&,','&',text,fixed=T)
	text <- gsub('}}','}',text,fixed=T)
	text <- gsub(',\\\\\\hline','\\\\\\hline',text,fixed=T)
return(text)}

latex.table.header <- function(genetics){
	N <- ncol(genetics$summary$table)+1
	text <- c('
	\\begin{sidewaystable}[p]\n
	%\\setcounter{table}{1}
	',
	paste('\\begin{center}\\begin{tabular}{|',paste(rep('c|',N),collapse=''),'}\\hline',sep=''),
	paste('&',paste(names(genetics$summary$table),collapse='&'),'\\\\\\hline',sep=''),
	paste('\\multicolumn{',N,'}{|l|}{Crime scene profiles}\\\\\\hline',sep=''))
return(text)}

csp.table.to.latex <- function(genetics){
	# formats the CSP table of alleles for latex
	table.csp <- table.collapser(genetics$cspData)
	table.unc <- table.collapser(genetics$uncData)
	table <- NULL; 
	for(n in 1:genetics$nrep){
		table <- rbind(table,table.csp[n,])
		table <- rbind(table,table.unc[n,])
		}
	colnames(table) <- colnames(genetics$cspData)
	row.names(table) <- rep(c('csp','unc'),genetics$nrep)
	text = c()
	N.col <- ncol(table)
	N.row <- nrow(table)
	for(row in 1:N.row){
		text.line <- character(N.col)
		for(col in 1:N.col){
			if(as.character(table[row,col])=='')text.line[col] <- '--'
			if(as.character(table[row,col])!='')text.line[col] <- gsub(' ',',',table[row,col])
			}
		if(row%%2==1)text=c(text,paste(paste(row.names(table)[row],paste(text.line,collapse='&'),sep='&'),'\\\\',sep=''))
		if(row%%2==0)text=c(text,paste(paste(row.names(table)[row],paste(text.line,collapse='&'),sep='&'),'\\\\\\hline',sep=''))
		}
return(text)}

ref.table.to.latex <- function(genetics){
# table: Reference table, from genetics$summary$table
	table <- genetics$summary$table
	text = paste('\\multicolumn{',ncol(table)+1,'}{|l|}{SUMMARY:}\\\\\\hline',sep='')
	N.col <- ncol(table)
	N.row <- nrow(table)
	for(row in 1:N.row){
		text.line = character(N.col)
		for(col in 1:N.col){
			tmp = sub('{','},{\\em',table[row,col],fixed=T)
			tmp = sub('[',',\\fx{',tmp,fixed=T)
			tmp = sub(']','}',tmp,fixed=T)

			text.line[col] = paste('{\\bf',tmp,'}',sep='')
			}
		text=c(text,paste(paste(row.names(table)[row],paste(text.line,collapse='&'),sep='&'),'\\\\\\hline',sep=''))
		}
return(text)}

latex.maker <- function(genetics,filename){
# table1: CSP table produced by allele.table() 
# table2: Reference table, from genetics$summary$table
	text.1 <- latex.table.header(genetics)
	text.2 <- csp.table.to.latex(genetics)
	text.3 <- ref.table.to.latex(genetics)
	text <- latex.cleaner(c(text.1,text.2,text.3))

	file <- file(filename)
	writeLines(text,file)
	close(file)
	}


pack.admin.input <- function(cspFile, refFile, caseName='dummy',databaseFile=NULL, outputPath=getwd() ) {
	# Packs and verifies administrative information.
	# Documentation in man directory.
    	paths <- c(cspFile, refFile) 
	if(!is.null(databaseFile)) paths <- c(databaseFile, paths, recursive=TRUE)
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
	admin <- list( caseName=caseName,
                databaseFile=databaseFile,
                cspFile=cspFile,
                refFile=refFile,
                outputPath=outputPath )
	return(admin)}


load.allele.database <- function(path=NULL) {
	# Loads allele database
	# Documentation is in man directory.
	if(is.null(path)) { # Load default database
	dummyEnv <- new.env()
	data('lgc-allele-freqs-wbp', package='likeLTD', envir=dummyEnv)
	return(dummyEnv[['lgc-allele-freqs-wbp']])
	}
	if(!file.exists(path)) stop(paste(path, "does not exist."))
	read.table(path, sep="\t", header=TRUE)
	}


unattributable.plot.maker <- function(genetics){
	plot <- ggplot(genetics$summary$counts, aes(x=loci,y=counts,fill=status))+
		geom_bar(stat='identity')+
		scale_fill_grey()
return(plot)}

table.collapser <- function(table){
	# collapses a table so that fields stored as lists become comma separated strings
	result <- array(,dim(table))
	for(row in 1:nrow(table)){
		for(locus in 1:ncol(table)){
			result[row,locus] <- paste(unlist(table[row,locus]),collapse=',')
			}}
return(result)}

unusual.alleles.per.table <- function(table,afreq){
	# finds unusual alleles in any specific table, (provided in original listed format)
	rare <- data.frame(locus=NULL,allele=NULL,frequency=NULL)
	loci <- colnames(table); loci <- loci[loci!='queried']# ref table includes 'queried'
	for(row in 1:nrow(table)){
		for(locus in loci){
			alleles <- unique(unlist(table[row,locus]))
			for(allele in alleles){
				condition <- afreq$Marker==locus & afreq$Allele==allele
				x <- afreq[condition,]
				if(nrow(x)==1){ # if the allele is present once in the database (should be!)
					if(x$EA1<2 | x$EA3<2 | x$EA4<2) {
            				frame <- data.frame(locus=locus, allele=allele, frequency=x[,4:6])
            				rare <- rbind(rare, frame)}
          					}
				if(nrow(x)==0){ # if the allele is absent from database it is probably a typo	
						frame <- data.frame(locus=locus, allele=allele, frequency=0)
            				rare <- rbind(rare, frame)
          					}
				if(nrow(x)>1){ # if the allele is more than once there is a problem with the database!	
						frame <- data.frame(locus=locus,allele=allele,frequency='present multiple times in database')
            				rare <- rbind(rare, frame)
          					}
        			}  # loop over alleles
     			} # loop over loci
		} # loop over profiles				
return(unique(rare))}

unusual.alleles <- function(genetics){
	#  creates and formats a combined table of unusual alleles

	t1.tmp <- unusual.alleles.per.table(genetics$refData,genetics$afreq)
	t1 <- cbind(data.frame(source=rep('Reference profiles',nrow(t1.tmp))),t1.tmp)
	t2.tmp <- unusual.alleles.per.table(genetics$cspData,genetics$afreq)
	t2 <- cbind(data.frame(source=rep('Crime scene certain',nrow(t2.tmp))),t2.tmp)
	t3.tmp <- unusual.alleles.per.table(genetics$uncData,genetics$afreq)
	t3 <- cbind(data.frame(source=rep('Crime scene uncertain',nrow(t3.tmp))),t3.tmp)
	table <- rbind(t1,t2,t3)
	if(nrow(table)==0)table <- data.frame(status='No unusual alleles present')
return(table)}

csp.table.reformatter <- function(genetics){
	table.csp <- table.collapser(genetics$cspData)
	table.unc <- table.collapser(genetics$uncData)
	table <- NULL; 
	for(n in 1:genetics$nrep){
		table <- rbind(table,table.csp[n,])
		table <- rbind(table,table.unc[n,])
		}
	colnames(table) <- colnames(genetics$cspData)
	extra <- data.frame(rep=rep(1:genetics$nrep,each=2),status=rep(c("certain","uncertain"),genetics$nrep))
	combined <- cbind(extra,table)
return(combined)}

reference.table.reformatter <- function(genetics){
	table <- genetics$summary$table
	extra <- data.frame(profile=row.names(table))
	combined <- cbind(extra,table)
return(combined)}

local.likelihood.table.reformatter <- function(prosecutionHypothesis,defenceHypothesis,prosecutionResults,defenceResults){
	P <- log10(locus.likes(prosecutionHypothesis,prosecutionResults))
	D <- log10(locus.likes(defenceHypothesis,defenceResults))
	table  <- t(data.frame(Prosecution.log10=P,Defence.log10=D,Ratio.log10=(P-D),Ratio=10^(P-D)))
	extra <- data.frame(Likelihood=row.names(table))
	combined <- cbind(extra,table)
return(combined)}

hypothesis.generator <- function(genetics){
	# replicated alleles may not be explained by dropin, therefore this dictates minimum number of Us (minU)
	# maximum Us determined by whichever locus has the most unattributables (rep and unrep combined)
	# reports all U + dropin combos between minU and maxU
	table <- genetics$summary$counts
	rep <- subset(table,table$status=='replicated')$counts
	unrep <- subset(table,table$status=='unreplicated')$counts
	minU <- ceiling(max(rep)/2)
	maxU <- ceiling(max(rep+unrep)/2) 
	unknowns <- minU:maxU
	N <- length(unknowns)
	dropins <- numeric(N)
	recommendation <- character(N)
	for(n in 1:N){
		dropins[n] <- sum(pmax(0,rep+unrep-2*unknowns[n]))
		if(dropins[n]<=2)recommendation[n]<-"strongly recommended"
		if(dropins[n]==3)recommendation[n]<- "worth considering"
		if(unknowns[n]==3)recommendation[n]<-"Can only be evaluated by removing the additional U from defence"
		if(unknowns[n]>3)recommendation[n]<-"Too many U's to evaluate"
		}
	result <- data.frame(nUnknowns=unknowns, doDropin=dropins, Recommendation=recommendation)
return(result)}

system.info <- function(){
	date <- date()
	package.info <- sessionInfo()$otherPkgs$likeLTD
	sys.info <- Sys.info()
	Details <- t(data.frame(Details=c(date,package.info,sys.info) ))
	Details  <- gsub("\n","",Details ,fixed=T)
	Details  <- gsub("   ","",Details ,fixed=T)
	Type <- data.frame(Type=c('Date report generated:',names(package.info),names(sys.info)))
	all <- cbind(Type,Details)
return(all)}

filename.maker <- function(outputPath,file,type=NULL){
	# type: report type, one of 'allele' or 'results'
	if(is.null(file)){
		if(is.null(type))report.type <- 'report'
		if(type=='allele')report.type <- 'allele report'
		if(type=='results')report.type <- 'results report'
		
		n <- 1
		file <- file.path(outputPath,paste(report.type,n,'docx',sep='.'))
		while(file.exists(file)){
			n <- n + 1
			file <- file.path(outputPath,paste(report.type,n,'docx',sep='.'))
			}
		}
return(file)}

hyp.P <- function(genetics){
	Q <- paste(genetics$nameQ,'(Q)') 
	HP <- paste('Prosecution Hypothesis:',paste(c(Q,genetics$nameK),collapse=' + ')  )
return(HP)}

hyp.D <- function(genetics){
	X <- 'Unknown (X)'
	HD <- paste('Defence Hypothesis:',paste(c(X,genetics$nameK),sep=' + ') )
return(HD)}

locus.likes <- function(hypothesis,results,...) 
	{
	# Generate locus likelihoods from overall likelihood
	#
	# Parameters:
	# 	hypothesis: generated by either defence.hypothesis() or 
      #     prosecution.hypothesis()
	# results: results from do.call(optim,params)
	model <- create.likelihood.vectors(hypothesis)
	arguments <- relistArguments(results$optim$bestmem, hypothesis, ...)
	likes <- do.call(model,arguments)
	likes <- likes$objectives * likes$penalties
	}

pack.genetics.for.allele.report <- function(admin){
	# packs together all the genetic information required for the allele.report()

	# primary objects
	cspData <- read.csp.profile(admin$cspFile)
	uncData <- read.unc.profile(admin$cspFile)
	refData <- read.known.profiles(admin$refFile)
	afreq <- load.allele.database(admin$databaseFile)

	# secondary objects
	summary <- summary.generator(refData, cspData)
	estimates <- estimate.csp(refData, cspData)
	nrep	<- nrow(cspData)

	allele.report.genetics <- list( 
	cspData = cspData, 
	uncData = uncData,
	refData = refData,
	afreq = afreq,
	summary = summary,
	estimates = estimates,
	nrep = nrep)
return(allele.report.genetics)}
	
pack.genetics.for.output.report <- function(P.hyp,D.hyp){
	# packs together all the genetic information required for the output.report()

	# Comparison checks between P.hyp and D.hyp 
	compare.hypothesis.inputs(P.hyp,D.hyp)

	# primary objects
	cspData <- read.csp.profile(P.hyp$cspFile)
	uncData <- read.unc.profile(P.hyp$cspFile)
	refData <- read.known.profiles(P.hyp$refFile)
	afreq <- load.allele.database(P.hyp$databaseFile)

	# secondary objects
	summary <- summary.generator(refData, cspData)
	estimates <- estimate.csp(refData, cspData)
	nrep	<- nrow(cspData)

	QvK <- queried.vs.known(P.hyp$refFile)
	nameQ <- row.names(refData)[QvK]
	nameK <- row.names(refData)[!QvK]

	output.report.genetics <- list( 
	cspData = cspData, 
	uncData = uncData,
	refData = refData,
	afreq = afreq,
	summary = summary,
	estimates = estimates,
	nrep = nrep,
	nameQ = nameQ,
	nameK = nameK, 
	P.hyp = P.hyp,
	D.hyp = D.hyp)  
return(output.report.genetics)}

compare.hypothesis.inputs <- function(P.hyp,D.hyp){
	# Checks the input data and parameters were the same for P and D
	if(!identical(P.hyp$cspFile,D.hyp$cspFile))warning("P and D hypotheses were constructed using two different crime scene files!!")
	if(!identical(P.hyp$refFile,D.hyp$refFile))warning("P and D hypotheses were constructed using two different reference files!!")
	if(!identical(P.hyp$databaseFile,D.hyp$databaseFile))warning("P and D hypotheses were constructed using two different allele database files!!")
	if(!identical(P.hyp$outputPath,D.hyp$outputPath))warning("P and D hypotheses were given two different output paths!!")
	if(!identical(P.hyp$ethnic,D.hyp$ethnic))warning("P and D hypotheses were constructed using two different ethnic codes!!")
	if(!identical(P.hyp$ethnic,D.hyp$ethnic))warning("P and D hypotheses were constructed using two different ethnic codes!!")
	if(!identical(P.hyp$adj,D.hyp$adj))warning("P and D hypotheses were constructed using two different locus adjustment parameter vectors!!")
	if(!identical(P.hyp$fst,D.hyp$fst))warning("P and D hypotheses were constructed using two different fsts!!")
	}


queried.vs.known <- function(path) {
	# Reads profile from path and returns queried vs known column.
	# Args:
	#	path: Path to file with the profile. 
	raw <- read.table(path, header=T, colClasses='character', row.names=1, sep=',', quote = "\"")
	if(is.null(raw$known.queried))stop("Reference csv must contain a column 'known/queried'. This column must contain one field 'queried'. ")
	if(sum(raw$known.queried=='queried')!=1)stop("The reference csv column 'known/queried' must contain one field 'queried'. ")
	return(raw$known.queried == 'queried')
	}

estimate.csp <- function(refData, cspData) {
  # Estimate how well each reference profile is represented in the CSP

	nrep <- nrow(cspData)
	# Constructs the result data frame.
	result <- data.frame(array(0, c(nrow(refData), nrep+1)), row.names=rownames(refData))
	colnames(result)[1:nrep] = sapply(1:nrep, function(n) {paste('Rep', n)})
	colnames(result)[nrep+1] = 'Total' 
  
	# now add data.
	for(person in row.names(refData)){
		for(rep in 1:nrep){
			# number of alleles in common for given person and CSP replicate, across all loci
			represented <- c()
			for(locus in 1:ncol(cspData)){
				if(!is.null(cspData[rep,locus])) {
					# figure out unique alleles in reference 
					ref.alleles  <- unique(unlist(refData[person,locus+1]))# +1 ignores first column(queried)
					csp.alleles <- unique(unlist(cspData[rep,locus]))
					# figure out how many of these allele are in CSP
					represented <- c(represented, ref.alleles %in% csp.alleles)
					}
				}
			if(length(represented) > 0) 
			result[person, rep] <- round(100*sum(represented)/length(represented))
			}
  		}
	if(nrep==1)  result[, nrep+1] <- result[, nrep]
	else         result[, nrep+1] <- round(rowSums(result)/nrep)
   
	# Reorders rows and return
return(result[order(result[, nrep+1], decreasing=T), ])}

estimates.reformatter <- function(genetics){
	table <- genetics$estimates
	extra <- data.frame(Contributor=row.names(table))
	combined <- cbind(extra,table)
return(combined)}

summary.generator <- function(refData, cspData){

	# summary table for Q and K, showing which of their alleles are replicated, unreplicated, or absent
	table <- as.data.frame(array(,c(nrow(refData)+1,ncol(cspData))))
	colnames(table) <- colnames(cspData)
	row.names(table) <- c(row.names(refData),'Unattributable')

	# count of unattributable alleles, to establish if they are replicated or not
	rep.counts <- unrep.counts <- c()

	# generate the results 
	for(locus in colnames(cspData)){
		# for unattributable alleles are rep, unrep, absent
		cspAlleles <- unlist(cspData[,locus])
		refAlleles <- refData[,locus][[1]]
		unattributableAlleles <- cspAlleles[!cspAlleles%in%refAlleles]
		table['Unattributable',locus] <- summary.helper(unattributableAlleles,cspAlleles)$string

		# for unattributable alleles calculate how many are replicated /unreplicated
		rep.counts <- c(rep.counts, summary.helper(unattributableAlleles,cspAlleles)$rep )
		unrep.counts <- c(unrep.counts, summary.helper(unattributableAlleles,cspAlleles)$unrep )

		# for each contributor calculate which alleles are rep, unrep, absent
		for(name in row.names(refData)){
			refAlleles <- refData[name,locus][[1]]
			table[name,locus] <- summary.helper(refAlleles,cspAlleles)$string
			}}	
		
		# reformat the counts table (c.)
		c.rep  <- data.frame(loci=colnames(cspData), counts=rep.counts, status='replicated')
		c.unrep <- data.frame(loci=colnames(cspData), counts=unrep.counts, status='unreplicated')
		counts <- rbind(c.rep,c.unrep)		

		summary <- list(table=table,counts=counts)
return(summary)}

summary.helper <- function(refAlleles,cspAlleles){
	# intricate and irritating operations required to check each allele 
	rep <- unrep <- absent <- c()
	for(allele in unique(refAlleles)){
		condition <- sum(cspAlleles%in%allele)
		if(condition>1)rep <- c(rep,allele)
		if(condition==1)unrep <- c(unrep,allele)
		if(condition==0)absent <- c(absent,allele)
		}
	# then add correct parentheses etc
	rep1 <- paste(rep,collapse=',')
	unrep1 <- paste(unrep,collapse=',')
	absent1 <- paste(absent,collapse=',')
	string <- paste(rep1, '{', unrep1,'}','[', absent1,']',sep='')
	# remove parentheses if they have nothing in them
	string <- gsub('{}','',string,fixed=T)
	string <- gsub('[]','',string,fixed=T)
return(list(string=string,rep=length(rep),unrep=length(unrep)))}

overall.likelihood.table.reformatter <- function(prosecutionResults,defenceResults){
	P <- -prosecutionResults$optim$bestval
	D <- -defenceResults$optim$bestval
	table  <- data.frame(estimate=t(data.frame(Prosecution.log10=P,Defence.log10=D,Ratio.log10=P-D,Ratio=10^(P-D))))
	extra <- data.frame(calculation=row.names(table))
	combined <- cbind(extra,table)
return(combined)}

rcontConvert <- function(refIndiv,rcont){
	# Convert rcont to full rcont (including ref individual)
	# Parameters:
	# 	refIndiv: reference individual specified in args
	#	rcont: rcont parameters from do.call(optim,params)
	if(refIndiv == 1) rcont = c(1, rcont)
      else if(refIndiv > length(rcont)) rcont = c(rcont, 1)
      else rcont = c(rcont[1:refIndiv-1], 1, rcont[refIndiv:length(rcont)])
return(rcont)}

calc.dropout = function(results, hypothesis){ 
	# Calculates dropout rates for every contributor subject to dropout and for every replicate
	# Parameters:
	# 	hypothesis: generated by either defence.hypothesis() or  prosecution.hypothesis()
	#	results: results from do.call(optim,params)

	# Number of contributors subject to dropout + number of unknowns + 1 (dropin)
	N <- nrow(hypothesis$dropoutProfs) + hypothesis$nUnknowns + 1
	nrep <- nrow(hypothesis$cspProfile)
	do <- results$optim$bestmem[grep("dropout",names(results$optim$bestmem))]
	rcont <- results$optim$bestmem[grep("rcont",names(results$optim$bestmem))]
	rcont <- rcontConvert(hypothesis$refIndiv,rcont)
	BB <- results$optim$bestmem[grep("power",names(results$optim$bestmem))]
	drout <- matrix(0,N-1,nrep)
	if(N>1) for(x in 1:(N-1)) for(z in 1:nrep) drout[x,z] <- do[z]/(do[z]+rcont[x]^-BB*(1-do[z]))
	return(drout)
	}

dropDeg <- function(hypothesis,results,genetics){
	# Output tables for dropout and degradation
	# Parameters:
	# 	hypothesis: generated by either defence.hypothesis() or prosecution.hypothesis() 
	#	results: results from do.call(optim,params)

	dropoutsLogical <- determine.dropout(genetics$refData,genetics$cspData)
	Qdrop <- dropoutsLogical[names(dropoutsLogical)==genetics$nameQ]
	knownDropoutsLogical <- dropoutsLogical[names(dropoutsLogical)!=rownames(hypothesis$queriedProfile)]
	# Number of Known contributors with No Dropout
	Nknd <- length(which(!knownDropoutsLogical)) 

	# create the row names, are arranged into the correct order: K, Q/X, U
	names.Dropout <- names(which(knownDropoutsLogical))
	names.NoDropout <- names(which(!knownDropoutsLogical))
	nameK <- c(names.NoDropout,names.Dropout)	
	Names=c()
	if(length(nameK)>0)for (n in 1:length(nameK)) Names[n]=paste(nameK[n],' (K',n,')',sep='')
	Names[length(c(genetics$nameQ,nameK))] = paste(genetics$nameQ,'(Q)')
	nU <- hypothesis$nUnknowns
	if(hypothesis$hypothesis=="defence"){
		Names[length(c(genetics$nameQ,nameK))] = 'X' 
		nU <- nU -1 # it already contains an extra one for X
		}
	if(nU>0)for(n in 1:nU)Names <- c(Names,paste('U',n,sep=''))

	# column names for the table
	runNames = c();for(rName in 1:genetics$nrep)runNames[rName]=paste('Replicate',rName)

	# dropout values
	h <- h1 <- calc.dropout(results, hypothesis)	
	# degradation values
	d <- d1 <- 10^results$optim$bestmem[grep("degradation",names(results$optim$bestmem))]
	# under pros, if Q is not subject to dropout, an extra 0 needs to be added to both dropout and degradation for Q
	if(hypothesis$hypothesis=="prosecution"){
		if(Qdrop==F){
			h = rbind(h[0:nrow(hypothesis$dropoutProfs),,drop=F],0)
			d = c(d[0:nrow(hypothesis$dropoutProfs)],0)
			if(hypothesis$nUnknowns>0){
				h = rbind(h,h1[(nrow(hypothesis$dropoutProfs)+1):length(d1),,drop=F])
				d = c(d,d1[(nrow(hypothesis$dropoutProfs)+1):length(d1)])
				}
			}
		}

	# Dropout is assumed zero for fully represented profiles, therefore suffixed on the dropout table
	if(Nknd>0) for(n in 1:Nknd){
		h = rbind(0,h)
		d = c(0,d)
		}
	Dropout = data.frame(h,d)
	colnames(Dropout)= c(runNames[1:genetics$nrep],'degradation')
	row.names(Dropout) = Names
	return(Dropout)
	}

overall.dropout.table.reformatter <- function(prosecutionHypothesis,defenceHypothesis,prosecutionResults,defenceResults,genetics){
	P.dropDeg <- dropDeg(prosecutionHypothesis,prosecutionResults,genetics)
	P.extra <- data.frame(hypothesis=rep('Prosecution',nrow(P.dropDeg)),contributor=rownames(P.dropDeg))
	P.combined <- cbind(P.extra,P.dropDeg)
	D.dropDeg <- dropDeg(defenceHypothesis,defenceResults,genetics)
	D.extra <- data.frame(hypothesis=rep('Defence',nrow(D.dropDeg)),contributor=rownames(D.dropDeg))
	D.combined <- cbind(D.extra,D.dropDeg)
	combined <- rbind(P.combined,D.combined)
return(combined)}

overall.dropin.table.reformatter <- function(prosecutionResults,defenceResults){
	P.dropin <- prosecutionResults$optim$bestmem["dropin"]
	D.dropin <- defenceResults$optim$bestmem["dropin"]
	table <- data.frame(hypothesis=c('Prosecution','Defence'),dropin=c(P.dropin,D.dropin))
return(table)}

optimised.parameter.table.reformatter <- function(hypothesis,result){
	table <- t(rbind(result$optim$bestmem,result$member$upper,result$member$lower))
	extra <- data.frame(parameter=rownames(table))
	combined <- cbind(extra,table)
	colnames(combined) <- c('parameter','estimate','upper bound','lower bound')
return(combined)}

chosen.parameter.table.reformatter <- function(prosecutionHypothesis){
	keep <- c(7,8,9,11,12,13,14,15)
	table <- as.data.frame(unlist(prosecutionHypothesis[keep]))
	extra <- data.frame(Parameter= rownames(table))
	combined <- cbind(extra,table)
	colnames(combined) <- c('Parameter','User input')
return(combined)}

ideal <- function(hypothesis,rr){
	# Calculates idealised likelihood assuming Q is perfect match
	# Parameters:
	# 	hypothesis: generated by either defence.hypothesis() or prosecution.hypothesis(), although rr only makes sense under defence
	#	rr: relatedness arguments from args
	ideal.match <- 1
	for(j in 1:ncol(hypothesis$cspProfile)){
		af = hypothesis$alleleDb[j][[1]]
		kn = hypothesis$queriedProfile[,j][[1]]
		p1 = af[row.names(af)==kn[1],1]
		p2 = af[row.names(af)==kn[2],1]
		ideal.match = ideal.match/(rr[2] + rr[1]*(p1+p2)/2 + (1-sum(rr))*p1*p2*(1+(kn[1]!=kn[2])))
		}
	result <- data.frame(calculation =c('likelihood ratio','Log10 likelihood ratio'),estimate=c(ideal.match,log10(ideal.match)))
	return(result)
	}

allele.report <- function(admin=NULL,file=NULL){

 	# admin: List containing administration data, as packed by pack.admin.input()
	# file: defaults to creating its own sequential file names (to avoid overwriting)

	if(is.null(admin))stop('missing argument: admin')

	# create genetics information 
	genetics <- pack.genetics.for.allele.report(admin)

  	# Latex output
	latex.maker(genetics,(paste(admin$outputPath,"/",file," table.tex",sep="")))

	# Create a new Docx. We need landscape which must be adapted from an existing template (seems to inherit funny properties)
	template.file <- file.path( find.package("likeLTD"), "extdata/landscape_template.docx" )
	myformat = tableProperties(integer.digit=1,fraction.double.digit=0)

	doc = new("Docx", basefile = template.file)

	# page 1
	doc = addParagraph( doc, value = c("Allele Report",admin$caseName), stylename = "TitleDoc" )
	doc = addPageBreak( doc )

	# page 2 contents
	doc = addHeader(doc, "Table of content", 1);
	doc = addTOC(doc)
	doc = addPageBreak( doc )

	doc = addHeader(doc, "Data provided by forensic scientist", 1)
	doc = addHeader(doc, "Crime scene profiles (CSP)", 2)
	doc = addTable(doc, data = csp.table.reformatter(genetics), span.columns = "rep",formats = myformat)
	doc = addHeader(doc, "Reference profiles", 2)
	doc = addParagraph( doc, value = "{unreplicated alleles}, [absent alleles]", stylename = "Caption" )
	doc = addTable(doc, data = reference.table.reformatter(genetics) ,formats = myformat)
	doc = addParagraph( doc, value = "Table 2.2 above uses the 'certain' allelic designations only.", stylename="Normal")
	doc = addPageBreak( doc )

	doc = addHeader(doc, "Summary", 1)
	doc = addHeader(doc, "Unattributable alleles", 2)
	doc = addPlot( doc, fun = print, x = unattributable.plot.maker(genetics) , width = 12, height = 4)
	doc = addParagraph( doc, value = "Table 3.1 The number of 'certain' alleles that cannot be attributed to a known profile. This is intended to assist in choosing how many unknowns are required, and whether dropin is required. Suggested values are given in Table 3.4 below.", stylename="Normal")
	doc = addPageBreak( doc )

	doc = addHeader(doc, "Unusual alleles", 2)
	doc = addTable(doc, unusual.alleles(genetics),format=myformat)
	doc = addHeader(doc, "Approximate representation", 2)
	doc = addTable(doc, estimates.reformatter(genetics), format=tableProperties(integer.digit=1,fraction.double.digit=1))
	doc = addParagraph( doc, value = "Table 3.3 The fraction of an individual's alleles (as a percentage) that have been designated as 'certain' alleles in each replicate. This estimate is not used by likeLTD, and is intended to assist informal assessments of possible known contributors to the CSP (other than the queried contributor). A more formal approach is to do a likeLTD run to compute the likelihood ratio (LR) for that individual contributor.",stylename="Normal")
	doc = addHeader(doc, "Suggested parameter values", 2)
	doc = addTable(doc, hypothesis.generator(genetics),format=myformat)
	doc = addParagraph( doc, value = "Table 3.4 Recommended values for 'nUnknowns' (0,1 or 2 under the prosecution hypothesis) and 'doDropin' (TRUE or FALSE). All the attributable alleles (see Table 3.1 above) must either come from an unknown or dropin. likeLTD automatically adds and additional unknown (X) to the defence hypothesis in place of the queried profile (Q).",stylename="Normal")
	doc = addPageBreak( doc )

	doc = addHeader(doc, "System information", 1)
	doc = addTable(doc,  system.info(),format=myformat)

writeDoc( doc, filename.maker(admin$outputPath,file,type='allele') )}

output.report <- function(prosecutionHypothesis=NULL,defenceHypothesis=NULL,prosecutionResults=NULL,defenceResults=NULL,file=NULL){

 	# prosecutionHypothesis: generated by prosecution.hypothesis()
 	# defenceHypothesis: generated by defence.hypothesis()
	# prosecutionResults: results from do.call(optim, prosecutionParams)
  	# defenceResults: results from do.call(optim, defenceParams)
	# file: defaults to creating its own sequential file names (to avoid overwriting)

	if(is.null(prosecutionResults))stop('missing argument: prosecutionResults')
	if(is.null(defenceResults))stop('missing argument: defenceResults')
	if(is.null(prosecutionHypothesis))stop('missing argument: prosecutionHypothesis')
	if(is.null(defenceHypothesis))stop('missing argument: defenceHypothesis')

	# create genetics information 
	genetics <- pack.genetics.for.output.report(prosecutionHypothesis,defenceHypothesis)

	# Create a new Docx. We need landscape which must be adapted from an existing template (seems to inherit funny properties)
	template.file <- file.path( find.package("likeLTD"), "extdata/landscape_template.docx" )
	myformat = tableProperties(integer.digit=1,fraction.double.digit=0)
	doc = new("Docx", basefile = template.file)

	# page 1
	doc = addParagraph( doc, value = c("Statistical Evaluation report",prosecutionHypothesis$caseName), stylename = "TitleDoc" )
	doc = addLineBreak( doc )
	doc = addParagraph( doc, value = hyp.P(genetics), stylename = "BulletList"  )
	doc = addParagraph( doc, value = hyp.D(genetics), stylename = "BulletList"  )
	doc = addParagraph( doc, value = "", stylename = "TitleDoc" )
	doc = addPageBreak( doc )

	# page 2 contents
	doc = addHeader(doc, "Table of content", 1);
	doc = addTOC(doc)
	doc = addPageBreak( doc )

	doc = addHeader(doc, "Data provided by forensic scientist", 1)
	doc = addHeader(doc, "Crime scene profiles (CSP)", 2)
	doc = addTable(doc, data = csp.table.reformatter(genetics), span.columns = "rep",formats = myformat)
	doc = addHeader(doc, "Reference profiles", 2)
	doc = addParagraph( doc, value = "{unreplicated alleles}, [absent alleles]", stylename = "Caption" )
	doc = addTable(doc, data = reference.table.reformatter(genetics) ,formats = myformat)
	doc = addParagraph( doc, value = "Table 2.2 above uses the 'certain' allelic designations only.", stylename="Normal")
	doc = addPageBreak( doc )

	doc = addHeader(doc, "Evaluation results", 1)
	doc = addHeader(doc, "Unattributable alleles", 2)
	doc = addPlot( doc, fun = print, x = unattributable.plot.maker(genetics) , width = 12, height = 4)
	doc = addParagraph( doc, value = "Table 3.1 above uses the 'certain' allelic designations to estimate how many alleles cannot be attributed to a known profile. See the Allele report for further details.", stylename="Normal")
	doc = addPageBreak( doc )

	doc = addHeader(doc, "Unusual alleles", 2)
	doc = addTable(doc, unusual.alleles(genetics),format=myformat)
	doc = addHeader(doc, "Approximate representation", 2)
	doc = addTable(doc, estimates.reformatter(genetics), format=tableProperties(integer.digit=1,fraction.double.digit=1))
	doc = addHeader(doc, "Likelihoods at each locus", 2)
	doc = addTable(doc, data = local.likelihood.table.reformatter(prosecutionHypothesis,defenceHypothesis,prosecutionResults,defenceResults) ,
		formats = tableProperties(integer.digit=1,fraction.double.digit=2))
	doc = addHeader(doc, "Overall Likelihood", 2)
	doc = addTable(doc, data = overall.likelihood.table.reformatter(prosecutionResults,defenceResults) ,
		formats = tableProperties(integer.digit=1,fraction.double.digit=2))
	doc = addHeader(doc, "Theoretical maximum LR", 2)
	doc = addTable(doc, data = ideal(defenceHypothesis,defenceHypothesis$relatedness) ,
		formats = tableProperties(integer.digit=1,fraction.double.digit=2))

	doc = addHeader(doc, "Dropout and degradation parameter estimates", 2)
	doc = addTable(doc, data = overall.dropout.table.reformatter(prosecutionHypothesis,defenceHypothesis,prosecutionResults,defenceResults,genetics) ,
		span.columns = "hypothesis", formats = tableProperties(integer.digit=1,fraction.double.digit=3))
	doc = addHeader(doc, "Dropin parameter estimates", 2)
	doc = addTable(doc, data = overall.dropin.table.reformatter(prosecutionResults,defenceResults) ,
		formats = tableProperties(integer.digit=1,fraction.double.digit=4))
	doc = addParagraph( doc, value = "NA indicates no dropout was modelled.", stylename="Normal")
	doc = addPageBreak( doc )

	doc = addHeader(doc, "User defined parameters", 2)
	doc = addTable(doc, data = chosen.parameter.table.reformatter(prosecutionHypothesis), formats = myformat)
	
	doc = addHeader(doc, "Optimised parameters", 2)
	doc = addHeader(doc, "Prosecution parameters", 3)
	doc = addTable(doc, data = optimised.parameter.table.reformatter(prosecutionHypothesis,prosecutionResults),
		formats = tableProperties(integer.digit=1,fraction.double.digit=3))
	doc = addHeader(doc, "Defence parameters", 3)
	doc = addTable(doc, data = optimised.parameter.table.reformatter(defenceHypothesis,defenceResults),
		formats = tableProperties(integer.digit=1,fraction.double.digit=3))
	doc = addPageBreak( doc )

	doc = addHeader(doc, "System information", 1)
	doc = addTable(doc,  system.info(),format=myformat)

writeDoc( doc, filename.maker(prosecutionHypothesis$outputPath,file,type='results') )}


