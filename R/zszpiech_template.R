#!/usr/bin/env Rscript --vanilla

# can probably incorporate these into autoload
library(doParallel)
library(qqman)
library(coda)


# notes:
# this script will eventually become several files with R code
# when coding these scripts, DO NOT use explicit covariate names 
# you MUST rename them in-house, e.g. "age => cov1", and use "cov1" in the scripts 

# functions
linear.lsa=function(in.lsa,in.pheno) {
    #n.nonlsa=ncol(in.pheno)
    #pheno.lsa=merge(in.pheno,in.lsa,by.x=0,by.y=0)
    tempresult=lm(in.pheno$delta1~in.lsa+in.pheno$age.months+in.pheno$Male+in.pheno$bmi)
    coef=summary(tempresult)$coef[2,]
    return(coef)
}

# variables
inNamePheno  = "PR.delta.brd.pheno.txt"
roh.file.pfx = "GALA2_mergedLAT-LATP_noParents_030816_PR."
rho.file.sfx = ".ROH.R.out"

# call this from doParallel for parallelizing calculations over chromosomes
registerDoParallel(cores=22)



pheno=read.table(inNamePheno,header=T,check.names=F,stringsAsFactors=F)
pheno$Male[pheno$Male == "Male"]=0
pheno$Male[pheno$Male == "Female"]=1

for (i in 1:nrow(pheno[,c(3,4)])){
	pheno$delta1[i]=max(pheno[i,c(3,4)],na.rm=T)
}

#foreach (i=1:22) %do%{
foreach (i=1:22) %dopar%{

	inNameLSA=paste('GALA2_mergedLAT-LATP_noParents_030816_PR.',i,'.ROH.R.out',sep="")
	outName=paste('GALA2_mergedLAT-LATP_noParents_030816_PR.',i,'.ROH.R.out.results',sep="")
    in.name.lsa = paste(roh.file.pfx, i, roh.file.sfx, sep = "")
    out.name    = paste(in.name.lsa, ".results", sep = "")

    ### replace this with a call to fread() and make a data.table
	lsa=read.table(in.name.lsa, header=T, check.names=F, stringsAsFactors = FALSE)

	new.lsa=lsa[rowSums(lsa[,3:dim(lsa)[2]]) > 0,3:dim(lsa)[2]]
	new.lsa=new.lsa[,pheno$SubjectID]
	name.order=order(colnames(new.lsa))
	new.lsa=new.lsa[,name.order]

	pos=lsa[rowSums(lsa[,3:dim(lsa)[2]]) > 0,2]
	result=t(apply(new.lsa,1,linear.lsa,in.pheno=pheno))

	result.final=cbind(pos,result)
	colnames(result.final)=c('Probe','beta','se','t','p')

    ### replace this with fwrite()
	write.table(result.final,out.name, quote=F, sep="\t", row.names=F)
}

# merge all data.frames by row
# start with chromosome 1
# attach chromosome 2 results to bottom
# do same with 3
# proceed for all 22 chromosomes
chr=1
file=paste("GALA2_mergedLAT-LATP_noParents_030816_PR.",chr,".ROH.R.out.results",sep="")
data=read.table(file,header = T)
gwas=data.frame(chr=rep(chr,dim(data)[1]),data)

for(chr in 2:22){
	file=paste("GALA2_mergedLAT-LATP_noParents_030816_PR.",chr,".ROH.R.out.results",sep="")
	data=read.table(file,header = T)
	df=data.frame(chr=rep(chr,dim(data)[1]),data)
	gwas=rbind(gwas,df)
	rm(df)
	rm(data)
}

# better plotting tools?
png(file="PR.roh.bdr.manhattan.png")
manhattan(gwas,chr="chr",p="p",bp="Probe",main="PR")



### the code below repeats what is done above
### goal: make completely reusable code
### then apply it to MX data, PR data
### the data will be previously split by pop before running

#registerDoParallel(cores=22)
#
#linear.lsa=function(in.lsa,in.pheno) {
#				     #n.nonlsa=ncol(in.pheno)
#				     #pheno.lsa=merge(in.pheno,in.lsa,by.x=0,by.y=0)
#				     tempresult=lm(in.pheno$delta1~in.lsa+in.pheno$age.months+in.pheno$Male+in.pheno$bmi)
#				     coef=summary(tempresult)$coef[2,]
#				     return(coef)
#}
#
#inNamePheno='MX.delta.brd.pheno.txt'
#pheno=read.table(inNamePheno,header=T,check.names=F,stringsAsFactors=F)
#pheno$Male[pheno$Male == "Male"]=0
#pheno$Male[pheno$Male == "Female"]=1
#
#for (i in 1:nrow(pheno[,c(3,4)])){
#	pheno$delta1[i]=max(pheno[i,c(3,4)],na.rm=T)
#}
#
##foreach (i=1:22) %do%{
#foreach (i=1:22) %dopar%{
#
#	inNameLSA=paste('GALA2_mergedLAT-LATP_noParents_030816_MX.',i,'.ROH.R.out',sep="")
#	outName=paste('GALA2_mergedLAT-LATP_noParents_030816_MX.',i,'.ROH.R.out.results',sep="")
#
#	lsa=read.table(inNameLSA,header=T,check.names=F,stringsAsFactors=F)
#
#	new.lsa=lsa[rowSums(lsa[,3:dim(lsa)[2]]) > 0,3:dim(lsa)[2]]
#	new.lsa=new.lsa[,pheno$SubjectID]
#	name.order=order(colnames(new.lsa))
#	new.lsa=new.lsa[,name.order]
#
#	pos=lsa[rowSums(lsa[,3:dim(lsa)[2]]) > 0,2]
#	result=t(apply(new.lsa,1,linear.lsa,in.pheno=pheno))
#
#	result.final=cbind(pos,result)
#	colnames(result.final)=c('Probe','beta','se','t','p')
#
#	write.table(result.final,outName,quote=F,sep='\t',row.names=F)
#}
#
#chr=1
#file=paste("GALA2_mergedLAT-LATP_noParents_030816_MX.",chr,".ROH.R.out.results",sep="")
#data=read.table(file,header = T)
#gwas=data.frame(chr=rep(chr,dim(data)[1]),data)
#
#for(chr in 2:22){
#	file=paste("GALA2_mergedLAT-LATP_noParents_030816_MX.",chr,".ROH.R.out.results",sep="")
#	data=read.table(file,header = T)
#	df=data.frame(chr=rep(chr,dim(data)[1]),data)
#	gwas=rbind(gwas,df)
#	rm(df)
#	rm(data)
#}
#
#png(file="MX.roh.bdr.manhattan.png")
#	manhattan(gwas,chr="chr",p="p",bp="Probe",main="MX")
#dev.off()
