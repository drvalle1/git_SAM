compare=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,xlim=rango,ylim=rango)
  lines(rango,rango,col='red')
}

compare(omega,omega.true)

compare(alpha,alpha.true)

k=data.frame(estim.cs=cs,true.cs=cs.true)
table(k)

ind=c(2,8,9,7,1)

compare(betas[,ind],betas.true)

plot(theta,type='h')
