tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  #qnorm can give some imprecise results
  cond=z<lo;    z[cond] = lo[cond]
  cond=z==-Inf; z[cond] = lo[cond]
  cond=z>hi;    z[cond] = hi[cond]
  cond=z==Inf;  z[cond] = hi[cond]
  z
}
#----------------------------------------------------------------------------------------------
sample.omega=function(y,nspp,nloc,xmat,betas,cs,alpha){
  omega=matrix(NA,nloc,nspp)
  
  media=xmat%*%betas
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  media1=alpha1+media[,cs]
  cond=y==1
  omega[ cond]=tnorm(sum(cond) ,lo=0,hi=Inf ,mu=media1[cond] ,sig=1)
  omega[!cond]=tnorm(sum(!cond),lo=-Inf,hi=0,mu=media1[!cond],sig=1)
  omega
}
#-----------------------------------------------
sample.alpha=function(nloc,xmat,betas,omega,cs,nspp){
  prec=nloc+(1/10)
  var1=1/prec
  media=xmat%*%betas
  media1=media[,cs]
  pmedia=colSums(omega-media1)
  rnorm(nspp,mean=var1*pmedia,sd=sqrt(var1))
}
#-----------------------------------------------
sample.betas=function(ngroups,cs,nparam,xtx,t.xmat,alpha,nloc,nspp,omega){
  i1=diag(1,nparam)
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  err=omega-alpha1
  betas=matrix(NA,nparam,ngroups)
  for (i in 1:ngroups){
    cond=cs==i
    nc=sum(cond)
    if (nc==0) betas[,i]=rmvnorm(1,mean=rep(0,nparam),sigma=i1)
    if (nc> 0) {
      prec=nc*xtx+i1
      var1=solve(prec)
      if (nc==1) soma.err=err[,cond]
      if (nc!=1) soma.err=rowSums(err[,cond])
      pmedia=t.xmat%*%soma.err
      betas[,i]=rmvnorm(1,mean=var1%*%pmedia,sigma=var1)
    }      
  }
  betas
}
#-----------------------------------------------
sample.cs=function(ngroups,omega,xmat,alpha,betas,theta,nspp,nloc){
  #calculate probabilities for each species for each group
  prob=matrix(NA,ngroups,nspp)
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  for (i in 1:ngroups){
    media=xmat%*%betas[,i]
    media1=alpha1+matrix(media,nloc,nspp)
    prob[i,]=colSums(dnorm(omega,mean=media1,sd=1,log=T))+log(theta[i])
  }
  
  #make max=0
  max1=apply(prob,2,max)
  max2=matrix(max1,ngroups,nspp,byrow=T)
  tmp=prob-max2
  
  #normalize probabilities
  tmp1=exp(tmp)
  soma=matrix(colSums(tmp1),ngroups,nspp,byrow=T)
  prob=tmp1/soma
  
  #sample cs
  rmultinom1(prob=t(prob),randu=runif(nspp))+1
}
#-----------------------------------------------
sample.theta=function(cs,ngroups,gamma,burnin,gibbs.step,betas,theta){
  #re-order thetas. Based on that, re-order everything else
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(theta,decreasing=T)
    theta=theta[ind]
    betas=betas[,ind]
    
    #get cs.new
    cs.new=cs; cs.new[]=NA
    for (i in 1:ngroups){
      cond=cs==ind[i]
      cs.new[cond]=i
    }
    cs=cs.new
  }
  
  #regular function based on FCD
  n=rep(0,ngroups)
  tmp=table(cs)
  n[as.numeric(names(tmp))]=tmp
  
  v=theta=rep(NA,ngroups)
  prod=1
  for (i in 1:(ngroups-1)){
    n.greater.k=n[-(1:i)]
    v[i]=rbeta(1,n[i]+1,sum(n.greater.k)+gamma)
    theta[i]=v[i]*prod
    prod=prod*(1-v[i])
  }
  theta[ngroups]=prod

  #to avoid numerical issues
  cond=v>0.99999
  v[cond]=0.99999
  
  list(theta=theta,betas=betas,cs=cs,v=v)
}
#----------------------------
get.logl=function(y,omega,nspp,nloc,xmat,betas,cs,alpha){
  media=xmat%*%betas
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  media1=alpha1+media[,cs]
  p=y*pnorm(media1,log=T)+(1-y)*pnorm(media1,lower.tail=F,log.p=T) #log(1-pnorm(media1))
  sum(p)
}
#----------------------------
sample.gamma=function(v,ngroups,gamma.possib){
  ngamma=length(gamma.possib)
  res=rep(NA,ngamma)
  for (i in 1:ngamma){
    res[i]=sum(dbeta(v,1,gamma.possib[i],log=T))
  }
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}