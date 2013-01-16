#------------------------------------------------
# Uses various interactive user inputs,
# both generated from this script and from 'allele table.R',
# and then calls the wrapper to run the main likeLTD program,
# and produce final reports

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
# 3. choose stain folder
case = tk_choose.dir(startDir, 'Please choose the STAIN folder...') 

#-----------------------------------------------
# 4. load data from previously generating allele report
setwd(case)
if(length(list.files(pattern='inputs.RData'))==1)load(file='inputs.RData')
if(length(list.files(pattern='inputs.RData'))==0){setwd(startDir);source(paste(startDir,'allele report.R',sep='/'))}

#-----------------------------------------------
# 5. choose Queried
nameQK =c(nameQ,nameK)
tt = tktoplevel()
tl = tklistbox(tt,height=length(names)+1)
tkgrid(tklabel(tt,text='Please choose the queried profile (Q)'))
tkgrid(tl)
for (i in 1:length(nameQ)) tkinsert(tl,'end',nameQ[i])
tkselection.set(tl,0)
OnOK = function(){
	nameQ <<- nameQ[as.numeric(tkcurselection(tl))+1]
	tkdestroy(tt)
	}
OK.but = tkbutton(tt,text='   OK   ',command=OnOK)
tkgrid(OK.but)
tkfocus(tt)
tkwait.window(tt)

#-----------------------------------------------
# 6. choose Known(s)

remainingNames = nameQK[nameQK!=nameQ]
if(length(remainingNames)>0){
	tt = tktoplevel()
	tl = tklistbox(tt,height=length(remainingNames)+1,selectmode='multiple')
	tkgrid(tklabel(tt,text='Please choose the known profiles (K), can be none, one or several'))
	tkgrid(tl)
	for (i in 1:length(remainingNames)) tkinsert(tl,'end',remainingNames[i])
	OnOK = function(){
		nameK <<- remainingNames[as.numeric(tkcurselection(tl))+1]
		tkdestroy(tt)
		}
	OK.but = tkbutton(tt,text='   OK   ',command=OnOK)
	tkgrid(OK.but)
	tkfocus(tt)
	tkwait.window(tt)
	}

#-----------------------------------------------
# 7. choose number of Unknowns
tt = tktoplevel()
tl = tklistbox(tt,height=3)
tkgrid(tklabel(tt,text='Please choose the number of UNKNOWN CONTRIBUTORS. Please refer to the Allele report for recommendations'))
tkgrid(tl)
for (i in 0:2) tkinsert(tl,'end',i)
OnOK = function(){
	NU <<- as.numeric(tkcurselection(tl))
	tkdestroy(tt)
	}
OK.but = tkbutton(tt,text='   OK   ',command=OnOK)
tkgrid(OK.but)
tkfocus(tt)
tkwait.window(tt)

#-----------------------------------------------
# 8. choose dropin
tt = tktoplevel()
tl = tklistbox(tt,height=2)
tkgrid(tklabel(tt,text='Please choose DROP-IN model'))
tkgrid(tl)
choiceDropin = c('TRUE','FALSE')
for (i in 1:2) tkinsert(tl,'end',choiceDropin[i])
tkselection.set(tl,0)
OnOK = function(){
	Drin <<- c(T,F)[as.numeric(tkcurselection(tl))+1]
	tkdestroy(tt)
	}
OK.but = tkbutton(tt,text='   OK   ',command=OnOK)
tkgrid(OK.but)
tkfocus(tt)
tkwait.window(tt)

#-----------------------------------------------
# 9. choose ethnicity
tt = tktoplevel()
tl = tklistbox(tt,height=3)
tkgrid(tklabel(tt,text='Please choose ETHNICITY: '))
tkgrid(tl)
choiceEthnicity = c('EA1 (Caucasian)','EA3 (Afro-Caribbean)','EA4 (South Asian)')
for (i in 1:3) tkinsert(tl,'end',choiceEthnicity[i])
tkselection.set(tl,0)
OnOK = function(){
	ethnic <<- c('EA1','EA3','EA4')[as.numeric(tkcurselection(tl))+1]
	tkdestroy(tt)
	}
OK.but = tkbutton(tt,text='   OK   ',command=OnOK)
tkgrid(OK.but)
tkfocus(tt)
tkwait.window(tt)
if(ethnic=='EA1')fst=0.02; if(ethnic!='EA1')fst=0.03

#-----------------------------------------------
# 10. select number of iterations
tt = tktoplevel()
tl = tklistbox(tt,height=6)
tkgrid(tklabel(tt,text='Please choose the number of iterations:'))
tkgrid(tl)
choiceIterations = c('1 (Test run)','100','5000 (Full run)','10,000')
for (i in 1:4) tkinsert(tl,'end',choiceIterations[i])
tkselection.set(tl,0)
OnOK = function(){
	niter <<- c(1,100,5000,10000)[as.numeric(tkcurselection(tl))+1]
	tkdestroy(tt)
	}
OK.but = tkbutton(tt,text='   OK   ',command=OnOK)
tkgrid(OK.but)
tkfocus(tt)
tkwait.window(tt)

#-----------------------------------------------
# 11. choose number of interim reports
tt = tktoplevel()
tl = tklistbox(tt,height=4)
tkgrid(tklabel(tt,text='Please choose the number of interim reports:'))
tkgrid(tl)
choiceIterations = c(0,4)
for (i in 1:2) tkinsert(tl,'end',choiceIterations[i])
tkselection.set(tl,0)
OnOK = function(){
	Ninterim <<- c(0,4)[as.numeric(tkcurselection(tl))+1]
	tkdestroy(tt)
	}
OK.but = tkbutton(tt,text='   OK   ',command=OnOK)
tkgrid(OK.but)
tkfocus(tt)
tkwait.window(tt)
interims = round(seq(0,niter,length.out=(Ninterim+2)))[2:(Ninterim+1)]
if(niter<100)interims=NA

#-----------------------------------------------
# 13. needs coding...sampling, relatedness adjustments
adj=1; rr=c(0,0)/4




