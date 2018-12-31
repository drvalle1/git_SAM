gibbs.clust.response.curve=function(y,xmat,ngroups,ngibbs,burnin){
t.xmat=t(xmat)
xtx=t.xmat%*%xmat

#useful stuff
nloc=nrow(y)
nspp=ncol(y)
nparam=ncol(xmat)
gamma.possib=seq(from=0.1,to=1,by=0.05)

#get invSigma, Sigma, lds for different values of n_k
InvSigma.precalc=xtx+diag(1,nparam)
Sigma.precalc=solve(InvSigma.precalc)
lds=determinant(Sigma.precalc)$modulus[1] #log of determinant of covariance matrix

#initial parameter values
betas=matrix(0,nparam,ngroups)
alpha=rep(0,nspp)
omega=matrix(-1,nloc,nspp)
omega[y==1]=1
cs=sample(1:ngroups,size=nspp,replace=T)
theta=rep(1/ngroups,ngroups)
gamma=0.1

#gibbs stuff
vec.betas=matrix(NA,ngibbs,nparam*ngroups)
vec.alpha=matrix(NA,ngibbs,nspp)
vec.theta=matrix(NA,ngibbs,ngroups)
vec.logl=matrix(NA,ngibbs,1)
vec.gamma=rep(NA,ngibbs,1)
vec.cs=matrix(NA,ngibbs,nspp)

for (i in 1:ngibbs){
  print(i)
  
  cs=sample.cs(ngroups=ngroups,omega=omega,xmat=xmat,alpha=alpha,betas=betas,theta=theta,
               nspp=nspp,nloc=nloc,InvSigma.precalc=InvSigma.precalc,
               Sigma.precalc=Sigma.precalc,lds=lds,cs=cs,nparam=nparam)
  # cs=cs.true
  
  betas=sample.betas(ngroups=ngroups,cs=cs,nparam=nparam,xtx=xtx,t.xmat=t.xmat,alpha=alpha,
                     nloc=nloc,nspp=nspp,omega=omega)
  # betas=cbind(betas.true,matrix(0,nparam,ngroups-ncol(betas.true)))
  
  omega=sample.omega(y=y,nspp=nspp,nloc=nloc,xmat=xmat,betas=betas,cs=cs,alpha=alpha)
  # omega=omega.true
    
  alpha=sample.alpha(nloc=nloc,xmat=xmat,betas=betas,omega=omega,cs=cs,nspp=nspp)
  # alpha=alpha.true
  
  tmp=sample.theta(cs=cs,ngroups=ngroups,gamma=gamma,burnin=burnin,gibbs.step=i,
                   betas=betas,theta=theta)
  betas=tmp$betas
  cs=tmp$cs
  theta=tmp$theta
  v=tmp$v[-ngroups]
  # theta=rep(1/ngroups,ngroups)
  
  gamma=sample.gamma(v=v,ngroups=ngroups,gamma.possib=gamma.possib)
  
  logl=get.logl(y=y,omega=omega,nspp=nspp,nloc=nloc,xmat=xmat,betas=betas,cs=cs,alpha=alpha)
  
  #store results
  vec.betas[i,]=betas
  vec.alpha[i,]=alpha
  vec.theta[i,]=theta
  vec.logl[i]=logl
  vec.cs[i,]=cs
  vec.gamma[i]=gamma
}
  seq1=burnin:ngibbs
  list(theta=vec.theta[seq1,],
       logl=vec.logl[seq1],
       betas=vec.betas[seq1,],
       alpha=vec.alpha[seq1,],
       cs=vec.cs[seq1,],
       gamma=vec.gamma[seq1])
}