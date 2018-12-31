rm(list=ls(all=TRUE))
set.seed(2)
library('mvtnorm')
library('Rcpp')

setwd('U:\\GIT_models\\git_cluster_rcurve_gamma')
source('gibbs functions.R')
sourceCpp('aux1.cpp')
source('gibbs clust response curve gamma.R')

ngroups=32
gamma=0.1
ngibbs=1000

nome=paste0('fake ',c('data ','xmat '),'sim ',ngroups,'ng','.csv')
y=data.matrix(read.csv(nome[1],as.is=T))
xmat=data.matrix(read.csv(nome[2],as.is=T))

res=gibbs.clust.response.curve(y=y,xmat=xmat,ngroups=50,
                               ngibbs=ngibbs,burnin=ngibbs/2)

    
