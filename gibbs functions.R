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
sample.cs=function(ngroups,omega,xmat,alpha,betas,theta,nspp,nloc,
                   InvSigma.precalc,Sigma.precalc,lds,cs,nparam){
  #pre-calculate stuff when beta is integrated out
  p1=nloc*log(2*pi)
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  err=omega-alpha1
  p2s=colSums(err^2)
  q1s=t(xmat)%*%err  
  mu=Sigma.precalc%*%q1s
  p3s=-t(mu)%*%InvSigma.precalc%*%mu; 
  p3s=as.numeric(diag(p3s))
  p4=-lds
  ltheta=log(theta)
  p5k=-2*ltheta
  
  #which groups exist and which don't?
  tab=rep(0,ngroups)
  tmp=table(cs)
  tab[as.numeric(names(tmp))]=tmp

  #calculate probabilities if group already exists or if group is brand new
  lprob.exist=matrix(NA,nspp,ngroups)
  lprob.not.exist=matrix(NA,nspp,ngroups)
  k=log(2*pi)
  resto.media=xmat%*%betas  
  for (i in 1:ngroups){
    #lprob if group exists
    resto.media1=matrix(resto.media[,i],nloc,nspp)
    err2=(omega-(alpha1+resto.media1))^2
    dnorm1=-0.5*(k+err2)
    lprob.exist[,i]=colSums(dnorm1)+ltheta[i]
    # z=colSums(dnorm(omega,mean=alpha1+resto.media1,sd=1,log=T))
    
    #lprob if group does not exist
    lprob.not.exist[,i]=-0.5*(p1+p2s+p3s+p4+p5k[i])
  }

  #for each species, sample the cs's
  for (i in 1:nspp){
    tab[cs[i]]=tab[cs[i]]-1 #all but this species
    
    #get probabilities for each group
    lprob=rep(NA,ngroups)
    cond=tab==0
    lprob[ cond]=lprob.not.exist[i,cond]
    lprob[!cond]=lprob.exist[i,!cond]

    #normalize probabilities
    max1=max(lprob)
    tmp=lprob-max1
    tmp1=exp(tmp)
    prob=tmp1/sum(tmp1)
    
    #sample cs
    ind=rmultinom(1,size=1,prob)
    ind1=which(ind==1)
    cs[i]=ind1
    tab[ind1]=tab[ind1]+1
  }
  cs
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
  soma=sum(log(1-v[-ngroups]))
  k=(ngroups-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  # sum(dbeta(v[-ngroups],1,gamma.possib[5],log=T))
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}