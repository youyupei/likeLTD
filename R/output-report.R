#-----------------------------------------------------------------------------------
# output() 
# function called to generate both a final pdf report and (if requested) interim pdf reports
#-----------------------------------------------------------------------------------

outputs = function(){ 
#-----------------------------------------------------------------------------------	
# 1.1. couple of other required objects, rushed over from control. Functions found in control functions
ideal.match=ideal() # Calculates idealised likelihood assuming Q is perfect match
all.dropouts = calc.all.dropouts(); hpdrout=all.dropouts$hpdrout; hddrout=all.dropouts$hddrout # Calculates dropout rates for every contributor subject to dropout and for every replicate

#-----------------------------------------------------------------------------------	
# 1. Status selection for interim or final report
if((itr)%in%interims)status='Interim'
if(itr==niter)status='Complete'
endTime = as.character(Sys.time())

#-----------------------------------------------------------------------------------
# 2. Table for reports (likelihood for each locus)
likelihood  = data.frame(signif(maxn$l,2),signif(maxd$l,2),signif(maxn$l/maxd$l,2))
colnames(likelihood)= c('HP','HD','Ratio')
row.names(likelihood)= names(alleles)

#-----------------------------------------------------------------------------------
# 3. Table for reports (likelihood overall)
overallLikelihood  = data.frame(signif(maxn$L,2),signif(maxd$L,2),signif(maxn$L/maxd$L,2),signif(log10(maxn$L/maxd$L),2))
colnames(overallLikelihood)= c('HP','HD','Ratio','Log Ratio')
row.names(overallLikelihood)= ''

#-----------------------------------------------------------------------------------
# 4. Table for reports (DROPOUT and DEGRADATION)

# names (including 'Q','K','U')in correct order 
prosNames=c()
prosNames[length(c(nameQ,nameK))] = paste(nameQ,'(Q)')
if(length(nameK)>0)for (n in 1:length(nameK)) prosNames[n]=paste(nameK[n],' (K',n,')',sep='')

defNames = prosNames
defNames[length(c(nameQ,nameK))] = 'X'

if(NU>0)for(n in 1:NU)prosNames = c(prosNames,paste('U',n,sep=''))
if(NU+1>1)for(n in 1:(NU+1-1))defNames = c(defNames,paste('U',n,sep=''))
runNames = c();for(rName in 1:nrep)runNames[rName]=paste('Run',rName)

# pros
hp = hp1 = round(hpdrout,4);	pLoss = pLoss1 = round(maxn$deg,4)
if(Qdrop==F){
	hp = rbind(hp[0:Nkdo,],0);if(NU>0){
		if(nrep>1) hp = rbind(hp,hp1[(Nkdo+1):length(pLoss1),])
		if(!nrep>1) hp = rbind(hp,as.matrix(hp1[(Nkdo+1):length(pLoss1),]))
		}
	pLoss = c(pLoss[0:Nkdo],0);if(NU>0)pLoss = c(pLoss,pLoss1[(Nkdo+1):length(pLoss1)])
	}
if(Nknd>0) for(n in 1:Nknd){
				hp = rbind(0,hp)
				pLoss = c(0,pLoss)
				}
pDropout = format(round(data.frame(hp,pLoss),3),nsmall=3)
colnames(pDropout)= c(runNames[1:(nrep)],'Degradation:')
row.names(pDropout) = prosNames

#def
hd = round(hddrout,4); dLoss = round(maxd$deg,4)
if(Nknd>0) for(n in 1:Nknd){
				hd = rbind(0,hd)
				dLoss = c(0,dLoss)
				}
dDropout = format(round(data.frame(hd,dLoss),3),nsmall=3)
colnames(dDropout)= c(runNames[1:(nrep)],'Degradation:')
row.names(dDropout) = defNames

#-----------------------------------------------------------------------------------
# 5. Table for reports (DROP-IN)
dropin = round(data.frame(maxn$rcont[length(maxn$rcont)],maxd$rcont[length(maxd$rcont)]),3)
colnames(dropin) = c('HP','HD')
row.names(dropin) = ''

#-----------------------------------------------------------------------------------
# 6. Table for reports (SYSTEM INFO)
operator = data.frame(Sys.info())
colnames(operator) = 'Operator info'

#-----------------------------------------------------------------------------------
# 7. Table for reports (PARAMETERS)
CSPname = tail(strsplit(CSPfile,split='/')[[1]],1)
REFname = tail(strsplit(REFfile,split='/')[[1]],1)
ALLELEname = tail(strsplit(alleleFreq,split='/')[[1]],1)

administration = data.frame(status,caseName,stain,CSPname,REFname,ALLELEname,ethnic,initiated,endTime,
adj,fst,rr[1],rr[2],Nknd,Nkdo,NU+1,NU,Drin,
maxn$best,maxn$better,maxd$best,maxd$better,
paste(format(round(start$nupa$do,2),nsmall=2),collapse=' '),paste(format(round(start$depa$do,2),nsmall=2),collapse=' '),
paste(format(round(maxn$do,2),nsmall=2),collapse=' '),paste(format(round(maxd$do,2),nsmall=2),collapse=' '),
paste(format(round(maxn$deg,2),nsmall=2),collapse=' '),
paste(format(round(start$nupa$rcont,2),nsmall=2),collapse=' '),paste(format(round(start$depa$rcont,2),nsmall=2),collapse=' '),
paste(format(round(maxn$rcont,2),nsmall=2),collapse=' '),paste(format(round(maxd$rcont,2),nsmall=2),collapse=' '),
paste(format(round(maxn$locadj,2),nsmall=2),collapse=' '),paste(format(round(maxd$locadj,2),nsmall=2),collapse=' '),
niter,signif(ideal.match,3),row.names=NULL)

colnames(administration) = c('Report status:','Case:','Stain:','CSP file:','Reference file:','Allele frequency file:','EA ethnic group:','Time started:','Time finished:',
'adj:','fst:','rr(1):','rr(2):','Nknd:','Nkdo:','NU+1:','NU:','Drin:',
'Best iteration(HP):','No. improvements to ML(HP):','Best iteration(HD):','No. improvements to ML(HD):',
'Initial HP do:','Initial HD do:','Max HP do:','Max HD do:','Initial deg:','Initial HP rcont:','Initial HD rcont:','Max HP rcont:','Max HD rcont:',
'Max HP local adjust:','Max HD local adjust:',
'End iteration:','Ideal match:')
row.names(administration) = 'Parameters'

#-----------------------------------------------------------------------------------
# 8. Creates unique file name (FN) that is unique to this particular combination of inputs
nameCode=c();for(n in 1:length(c(nameQ,nameK)))nameCode=c(nameCode,strsplit(c(nameQ,nameK)[n],'')[[1]][1:2]);nameCode = paste(nameCode,collapse='')
version=1; FN=paste(stain,nameCode,NU,ethnic,Drin,itr,version,sep='-'); pdf.name = paste(FN,'pdf',sep='.'); RData.name = paste(FN,'RData',sep='.')
while(pdf.name %in% list.files(case)){
	version=version+1
	FN=paste(stain,nameCode,NU,ethnic,Drin,itr,version,sep='-')
	pdf.name = paste(FN,'pdf',sep='.')
	RData.name = paste(FN,'RData',sep='.')
}

#-----------------------------------------------------------------------------------
# 9. State overall general hypothesis
HP = nameQ; HD = c('Unknown (X)')
if(length(nameK)>0)for (x in 1:length(nameK)){
	HP = paste(HP,nameK[x],sep=' + ')
	HD = paste(HD,nameK[x],sep=' + ')
	}
if(NU>0)for (x in 1:NU){
	HP = paste(HP,paste('unknown (U',x,')',sep=''),sep=' + ')
	HD = paste(HD,paste('unknown (U',x,')',sep=''),sep=' + ')
	}
if (Drin ==T){
	HP = paste(HP,' + Dropin')
	HD = paste(HD,' + Dropin')
	}
hypothesis = t(data.frame(HP,HD))
colnames(hypothesis) = NA
row.names(hypothesis) = c('Prosecution(HP)','Defence(HD)')

#-----------------------------------------------------------------------------------
# 10. Nomenclature
nomenclature = '
adj:		sampling adjustment parameter
fst:		coancestry parameter (remote shared ancestry of Q and X).
rr:		probabilities of X and Q sharing 1 and 2 alleles through
		shared inheritance from recent ancestors such as parents 
		and grandparents; for full sibs, rr[1]=0.5 and rr[2]=0.25.
Nknd:		number of profiled potential contributors not subject to dropout.
Nkdo:		number of profiled potential contributors subject to dropout.
NU+1:		number of unknown contributors under Hd.
NU:		number of unknown contributors under Hp.
Drin:   TRUE if dropin is being modelled.'

#-----------------------------------------------------------------------------------
# 11. Prints report
# report pg1
pdf(width=8.25,height=11.5,file=pdf.name);size = 1;hsize = 1.3
layout(mat = c(1,2,3,4,5,6), heights = c(1,1,3,1,1,3))
par(mai = c(0.5,0.5,0.5,0.5))
textplot(paste(status,'Report for ',stain,' in ',caseName),cex=1.5)
textplot(hypothesis,valign='top',cex=size);title('Hypothesis',cex=hsize)
textplot(alleles,valign='top',cex=tablesize);title('Crime scene profiles',cex=hsize)
textplot(summary$summary,valign='top',cex=tablesize);title('Summary. {}=unreplicated, []=absent',cex=hsize)
textplot(t(likelihood),valign='top',cex=size);title('likelihood',cex=hsize)
yaxp=c(0,max(otherBoth),max(otherBoth));if(max(otherBoth)==0)yaxp=c(0,1,1)
barplot(otherBoth,names.arg=names(cprofs),yaxp=yaxp,main='Unattributable alleles',ylab='No. alleles')
barplot(add=T,col='black',otherRep, yaxt='n')
if(max(otherBoth)>3)abline(h=2,lty=2);if(max(otherBoth)>5)abline(h=4,lty=2);if(max(otherBoth)>7)abline(h=6,lty=2);if(max(otherBoth)>9)abline(h=8,lty=2);
legend(x=0,y=max(otherBoth),yjust=0,bty='n',c('replicated','unreplicated'),xpd=NA,col=c('black','grey'),pch=c(15,15))

# report pg2
par(mfrow= c(6,1),mai = c(0.5,0.5,0.5,0.5))
textplot(t(overallLikelihood),valign='top',cex=size);title('Likelihood (overall)',cex=hsize)
textplot(pDropout,valign='top',cex=size);title('HP Drop-out',cex=hsize)
textplot(dDropout,valign='top',cex=size);title('HD Drop-out',cex=hsize)
textplot(estimates,cex=1,valign='top');title('Approximate representation %',cex=hsize)
textplot(t(dropin),valign='top',cex=size);title('Drop-in (overall)',cex=hsize)
textplot(operator,valign='top',cex=size);title('System info',cex=hsize)

# report pg3
layout(mat = c(1,2), heights = c(2,1))
par(mai = c(0.5,0.5,0.5,0.5));size = 0.7
textplot(t(administration),valign='top',cex=size);title('Parameters')
textplot(nomenclature,valign='top',cex=size);title('Nomenclature')

dev.off()

#-----------------------------------------------------------------------------------
# 11. Saves workspace
save.image(file=RData.name)

}
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------


