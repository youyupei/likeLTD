#------------------------------------------------
# Generates an allele summary, suggested hypotheses, and sundry info
# from the mixed and reference profile .csv,
# using various interactive user inputs
#------------------------------------------------
# NOT NECESSARY if packages are already installed!!!!
#install.packages('tcltk')
#install.packages('gplots')

#----------------------------------------------
require(tcltk)
tclRequire("BWidget")
library('gplots')
#-----------------------------------------------
startDir = getwd()

#-----------------------------------------------
# 1. choose allele frequency file
alleleFreq = tk_choose.files(caption='Please choose the ALLELE FREQUENCY file...')
afreq = read.table(alleleFreq,sep="\t",header=T)

#-----------------------------------------------
# 1. choose stain folder
case = tk_choose.dir(startDir, 'Please choose the STAIN folder...') 
caseName = rev(strsplit(case,'/')[[1]])[1]
setwd(case)

#-----------------------------------------------
# 2. choose mixed CSP file
CSPfile = tk_choose.files(caption='Please choose the MIXED CSP file...')

#-----------------------------------------------
# 3. choose reference file
REFfile = tk_choose.files(caption='Please choose the REFERENCE file...')

#-----------------------------------------------
# 4. Read mixed profile, and reference profiles
rawCSP = read.table(CSPfile,header=T,colClasses='character',sep=',')
onlyCSP = rawCSP[,5:dim(rawCSP)[2]]

REF.all = read.table(REFfile,header=T,colClasses='character',row.names=1,sep=',',quote = "\"")
REF = REF.all[,2:dim(REF.all)[[2]]]

# removes spaces from formatting
for(locus in 1:dim(REF)[2]){
	for(person in 1:dim(REF)[1]){
		REF[person,locus] = gsub(' ','',REF[person,locus])}}

caseName = strsplit((strsplit(case,'cases/')[[1]])[2],'/')[[1]][1]

#-----------------------------------------------
# 4. Search for unusual alleles that may be mistakes
rare =data.frame(name=NULL,locus=NULL,allele=NULL,EA1=NULL,EA3=NULL,EA4=NULL)
for(name in 1:dim(REF)[1]){
	for(locus in 1:dim(REF)[2]){
		for(allele in unique(strsplit(REF[name,locus],',')[[1]])){
			x = afreq[afreq$Marker==names(REF)[locus]&afreq$Allele==allele,]
			if(x$EA1<2|x$EA3<2|x$EA4<2)rare=rbind(rare,data.frame(name=row.names(REF)[name],locus=names(REF)[locus],allele=allele,EA1=x$EA1,EA3=x$EA3,EA4=x$EA4))
}}}
if(length(rare)==0)rare='No unusual alleles in Q or K profiles'
#-----------------------------------------------
# 5. Transform CSP into 'cprofs' format
nrep = dim(onlyCSP)[1]/2
cprofs = as.list(onlyCSP[0,]) # creates a list of all loci
for(locus in 1:length(cprofs)){	
	cprofs[[locus]]=list()  # creates a list of reps in each loci
	for(rep in 1:nrep){
		cprofs[[locus]][[rep]]=list(csp=NULL,unc=NULL) # fills it with NULL csp and unc
		if (is.na(onlyCSP[(2*rep)-1,locus])) cprofs[[locus]][rep]=list(NULL) # NAs must be formatted differently
		else {
			csp.value = gsub(' ','', onlyCSP[(2*rep)-1,locus]); if(csp.value!='')cprofs[[locus]][[rep]]$csp = strsplit(csp.value,',')[[1]] # values into csp 
			unc.value = gsub(' ','', onlyCSP[(2*rep),locus]); if(unc.value!='')cprofs[[locus]][[rep]]$unc = strsplit(unc.value,',')[[1]] # values into unc 
		}}}
#-----------------------------------------------
# 6. Estimate how well each reference profile is represented in the CSP
allRepCSP = array(0,c(dim(REF)[1],nrep+1))
for(person in 1:dim(REF)[1]){
	for(rep in 1:nrep){
		representedCSP=c()
		for(locus in 1:length(cprofs)){
			alleles = unique(strsplit(REF[person,locus],',')[[1]])
			if(!is.null(cprofs[[locus]][[rep]]))representedCSP = c(representedCSP,alleles%in%cprofs[[locus]][[rep]]$csp)
			}
		allRepCSP[person,rep] = round(100*sum(representedCSP)/length(representedCSP))
		}}
if(nrep==1)allRepCSP[,nrep+1]= allRepCSP[,nrep]
if(nrep>1)allRepCSP[,nrep+1]= round(rowSums(allRepCSP)/nrep)
runNames = c(); for(n in 1:nrep)runNames[n]=paste('run',n)
estimates = data.frame(allRepCSP,row.names=rownames(REF));colnames(estimates)=c(runNames,'total')
estimates = estimates[order(estimates[,nrep+1],decreasing=T),]
names = rownames(estimates)

#-----------------------------------------------
# 7. allele table for CSP
nloc = length(cprofs); nrep = length(cprofs[[1]])
alleles = array(0,c(2*nrep,nloc))
cspCombined.all = replicated.all = unreplicated.all = vector('list',nloc)
for(l in 1:nloc){ 
	cspCombined = c()
	for(r in 1:nrep){
		if(is.null(cprofs[[l]][[r]])){
			alleles[2*r-1,l]='NA'
			alleles[2*r,l] ='NA'
			}
		if(!is.null(cprofs[[l]][[r]])){
			alleles[2*r-1,l] = paste(cprofs[[l]][[r]]$csp,collapse=' ')
			alleles[2*r,l] = paste(cprofs[[l]][[r]]$unc,collapse=' ')
			}
		cspCombined=c(cspCombined,cprofs[[l]][[r]]$csp)
		}
	if(!is.null(cspCombined)){
		cspCombined.all[[l]] = cspCombined 
		replicated.all[[l]] = unique(cspCombined[duplicated(cspCombined)])
		unreplicated.all[[l]] = cspCombined[!cspCombined%in%replicated.all[[l]]]}}

alleleRows = c();for(r in 1:nrep)alleleRows=c(alleleRows,paste('csp',r),paste('unc',r))
alleles = data.frame(alleles);row.names(alleles)=alleleRows;colnames(alleles)=names(cprofs)

#-----------------------------------------------
# 8. Extracts Queried and Knowns
nameQ = row.names(REF)[REF.all$known.queried=='queried']
nameK = row.names(REF)[REF.all$known.queried=='known']

#-----------------------------------------------
# 9. Function to calculate 'summary', the allele table for Q and K
# This is re-used in the wrapper once Q and K is specified, because 'summary' varies accordingly.
# It also gets used to calculate possible hypotheses when there are 2 Q's, using otherRep and other Unrep
summary.generator = function(nameQ,nameK){
	knownNames = c()
	for (name in nameQ)knownNames =c(knownNames,paste(name,'(Q)'))
	for (name in nameK)knownNames =c(knownNames,paste(name,'(K)'))
	ncont = length(knownNames)

	index = c(); for (x in c(nameQ,nameK))index = c(index,which(row.names(REF)==x))
	known = as.list(onlyCSP[0,])
	for(locus in 1:length(cprofs))known[[locus]] = strsplit(paste(REF[[locus]][index],collapse=','),',')[[1]] # converts into a vector

	summary = array(0,c(ncont+1,nloc))
	otherRep <<- numeric(nloc); otherUnrep = numeric(nloc)
	for(l in 1:nloc){ 
		contCombined = c()
		for(c in 1:ncont){
			cont = known[[l]][(2*c-1):(2*c)]
			contCombined = c(contCombined,cont)
			summary[c,l] = paste(paste(cont[cont%in%replicated.all[[l]]],collapse=' '),'{',paste(cont[cont%in%unreplicated.all[[l]]],collapse=' '),'}','[',paste(cont[!cont%in%cspCombined.all[[l]]],collapse=' '),']',sep='')
			}
		# other alleles not found in the cprofs
		unattributable = cspCombined.all[[l]][!cspCombined.all[[l]]%in%contCombined]
		unattributableRep = unique(unattributable[duplicated(unattributable)]); otherRep[l] = length(unattributableRep)
		unattributableUnrep = unattributable[!unattributable%in%unattributableRep]; otherUnrep[l] = length(unattributableUnrep)
		summary[dim(summary)[1],l] = paste(paste(unattributableRep,collapse=' '),'{',paste(unattributableUnrep,collapse=' '),'}',sep='')
		}
	summary = data.frame(summary);row.names(summary)=c(knownNames,'Unattributable');colnames(summary)=names(cprofs)
return(list(summary=summary,otherRep=otherRep,otherUnrep=otherUnrep))}



#--------------------------------------------------------------
# 10. function to calculate possible hypotheses
# needs to be done with various permutations of 2 Q's
hypothesis.generator = function(nameQ,nameK,otherRep,otherUnrep){
	
	otherBoth = otherRep+otherUnrep
	minU = ceiling(max(otherRep)/2)
	maxU = ceiling(max(otherBoth)/2) 

	hypotheses = c()
	for(n in 1:(maxU-minU+1)){
		Unknowns=minU+n-1
		extras = otherBoth-(2*(Unknowns))
		Unattributable.alleles = sum(extras[extras>0])
		if(Unattributable.alleles==0)Drop.in=F;if(Unattributable.alleles!=0)Drop.in=T

		Q = paste(nameQ,'(Q)',sep='')
		K = c(); if(length(nameK)>0)for(n in 1:length(nameK))K=paste(K,paste(nameK[n],'(K',n,')',sep=''),sep='+')
		U = c(); if(Unknowns>0)for(n in 1:Unknowns)U=paste(U,paste('U',n,sep=''),sep='+')
		D = c(); if(Drop.in)D = paste('+',Unattributable.alleles,'dropins',sep='')
		hypotheses[n] = paste(Q,K,U,D,' vs X',K,U,D,sep='')
		}
return(hypotheses[!is.na(hypotheses)])}

#--------------------------------------------------------------
# 10.1 estimates all possible hypotheses
if(length(nameQ)==1){
	s = summary.generator(nameQ,nameK)
	suggested = hypothesis.generator(nameQ,nameK,s$otherRep,s$otherUnrep)}
if(length(nameQ)==2){
	suggested=c()
	s1 = summary.generator(nameQ[1],nameK)
	suggested=c(suggested,hypothesis.generator(nameQ[1],nameK,s1$otherRep,s1$otherUnrep))

	s2 = summary.generator(nameQ[2],nameK)
	suggested=c(suggested,hypothesis.generator(nameQ[2],nameK,s2$otherRep,s2$otherUnrep))

	nameK1 = c(nameQ[2],nameK)
	s3 = summary.generator(nameQ[1],nameK1)
	suggested=c(suggested,hypothesis.generator(nameQ[1],nameK1,s3$otherRep,s3$otherUnrep))

	nameK2 = c(nameQ[1],nameK)
	s4 = summary.generator(nameQ[2],nameK2)
	suggested=c(suggested,hypothesis.generator(nameQ[2],nameK2,s4$otherRep,s4$otherUnrep))
	}	
text = 
'Note: likeLTD is limited to 2 unknowns. If more than 2 
unknowns are suggested,this may indicate a mixed profile
with too many alleles to be statisically informative.'

#-----------------------------------------------
# 11. generate 'summary', whether multiple Q or single Q
summary = summary.generator(nameQ,nameK)
otherBoth = summary$otherRep+summary$otherUnrep
otherRep = summary$otherRep

#-----------------------------------------------
# 12. Write the allele report
# works out sensible font size
CSPtotal = UNCtotal = numeric(nrep)
for(l in 1:nloc){ 
	for(r in 1:nrep){
		CSPtotal[r] = CSPtotal[r]+length(cprofs[[l]][[r]]$csp)
		UNCtotal[r] = UNCtotal[r]+length(cprofs[[l]][[r]]$unc)
		}}
tablesize = min(1,25/max(c(CSPtotal,UNCtotal)))

stain = rev(strsplit(case,'/')[[1]])[1]
pdf(width=8.25,height=11.5,file=paste(case,'/',stain,'-allele report.pdf',sep=''))
par(mfrow= c(6,1),mai = c(0.7,0.5,0.5,0.5))
textplot(alleles,valign='top',cex=tablesize);title(paste(stain,'Allele Report'))
textplot(summary$summary,valign='top',cex=tablesize);title('Summary. {}=unreplicated, []=absent')
yaxp=c(0,max(otherBoth),max(otherBoth));if(max(otherBoth)==0)yaxp=c(0,1,1)
barplot(otherBoth,names.arg=names(cprofs),yaxp=yaxp,main='Unattributable alleles',ylab='No. alleles')
barplot(add=T,col='black',otherRep, yaxt='n')
if(max(otherBoth)>3)abline(h=2,lty=2);if(max(otherBoth)>5)abline(h=4,lty=2);if(max(otherBoth)>7)abline(h=6,lty=2);if(max(otherBoth)>9)abline(h=8,lty=2);
legend(x=0,y=max(otherBoth),yjust=0,bty='n',c('replicated','unreplicated'),xpd=NA,col=c('black','grey'),pch=c(15,15))
textplot(rare,cex=1,valign='top');title('Rare alleles')
textplot(estimates,cex=1,valign='top');title('Approximate representation %')
textplot(suggested,cex=1,valign='top');title(main='Suggested Hypotheses',sub=text)
dev.off()
save.image(file = paste(case,'inputs.RData',sep='/'))

