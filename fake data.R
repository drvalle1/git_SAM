rm(list=ls(all=TRUE))
library('mvtnorm')
set.seed(1)

nloc=2679
nspp=443
nparam=6
ngroup1=32

xmat=matrix(rnorm(nparam*nloc),nloc,nparam)
alpha.true=alpha=rnorm(nspp,mean=0,sd=0.4)
betas.true=betas=t(rmvnorm(ngroup1,mean=rep(0,nparam),sigma=diag(1,nparam)))
cs.true=cs=sample(1:ngroup1,size=nspp,replace=T)
    
omega=matrix(NA,nloc,nspp)
for (i in 1:nspp){
  media=alpha[i]+xmat%*%betas[,cs[i]]
  omega[,i]=rnorm(nloc,mean=media,sd=1)
}
y=omega.true=omega
y[omega>0]=1
y[omega<0]=0
    
setwd('U:\\GIT_models\\git_cluster_rcurve_gamma')
nome=paste0('fake ',c('data ','xmat ','cs ','betas ','alpha '),'sim ',ngroup1,'ng','.csv')
write.csv(y,nome[1],row.names=F)
write.csv(xmat,nome[2],row.names=F)
write.csv(cs,nome[3],row.names=F)
write.csv(betas,nome[4],row.names=F)
write.csv(alpha,nome[5],row.names=F)
