theta.estim=res$theta[nrow(res$theta),]
plot(theta.estim,type='h')
sum(theta.estim>0.005)

cs.estim=res$cs[nrow(res$cs),]
betas.estim=matrix(res$betas[nrow(res$betas),],6,50)
alpha.estim=res$alpha[nrow(res$alpha),]

setwd('U:\\GIT_models\\git_cluster_rcurve_gamma')
alpha.true=read.csv('fake alpha sim 32ng.csv',as.is=T)$x
betas.true=read.csv('fake betas sim 32ng.csv',as.is=T)
cs.true=read.csv('fake cs sim 32ng.csv',as.is=T)$x

plot(alpha.true,alpha.estim)

combo=data.frame(true=cs.true,estim=cs.estim)
z=table(combo)
ind=numeric()
for (i in 1:nrow(z)){
  tmp=which(z[i,]==max(z[i,]))
  ind=c(ind,colnames(z)[tmp])
}

plot(betas.estim[,as.numeric(ind)],data.matrix(betas.true))

plot(theta,type='h')
