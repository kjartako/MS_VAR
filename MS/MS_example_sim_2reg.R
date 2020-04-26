######################################################################################
######################################################################################
######################################################################################
# Simulated data

library(MASS)
mu = c(0.2,0.3);mu2 = c(-0.3,0.45);rho = 0.2;sd=c(0.15,0.15)
sigma = matrix(c(sd[1]^2,rep(rho*sd[1]*sd[2],2),sd[2]^2),nrow=2,byrow=TRUE)
Q = matrix(c(0.97,0.03,0.1,0.9),nrow=2,byrow=TRUE)
delta = solve(t(diag(2)-Q +1),rep(1,2))

n=750
S = numeric(n)
Y = matrix(NA,nrow=n,ncol=2)

set.seed(123)
S[1] = sample(c(0,1),size=1,prob=delta)
Y[1,] = mvrnorm(n=1,mu=(1-S[1])*mu+S[1]*mu2,Sigma=sigma)

for (i in 2:n)
{
  S[i]=sample(c(0,1),size=1,prob=Q[S[i-1]+1,])
  Y[i,] =mvrnorm(n=1,mu=(1-S[i])*mu+S[i]*mu2,Sigma=sigma)
}

######################################################################################
######################################################################################
######################################################################################

library(coda)
library(rstan)
options(mc.cores = parallel::detectCores())

chains=4
it_burnin = 500
it_sample = 1000

sm = stan_model('Stan/MS_lkj.stan',auto_write=TRUE)
#sm = stan_model('Stan/MS_wishart.stan',auto_write=TRUE)

######################################################################################

n_dim=ncol(Y)
n_reg=2
reg_m=TRUE
reg_s=FALSE

stanfit = sampling(sm,data=list(T=n,dim=n_dim,y=Y,nreg=n_reg,mean_reg=reg_m,sigma_reg=reg_s,Q_alpha=as.array(rep(1,n_reg)),
                   mu_mean=0,mu_sd=0.2,eta=1,gamma_alpha=2,gamma_beta=0.1),
                   iter=(it_burnin+it_sample), warmup = it_burnin, chains=chains,
                   pars=c("Q","mu","sigma","S"))

# stanfit = sampling(sm,data=list(T=n,dim=n_dim,y=Y,nreg=n_reg,mean_reg=reg_m,sigma_reg=reg_s,Q_alpha=as.array(rep(1,n_reg)),
#                                 mu_mean=0,mu_sd=1,eye=diag(n_dim),L=chol(30*diag(n_dim)),nu=n_dim+1),
#                    iter=(it_burnin+it_sample), warmup = it_burnin, chains=chains,
#                    pars=c("Q","mu","sigma","S"))

######################################################################################
# Label-switching

pars =mcmc.list(lapply(1:ncol(stanfit), function(x) mcmc(as.array(stanfit)[,x,])))

source("labsw_2reg_mu.R")
pars_ls = labsw_2reg(pars,chains,n_dim,n_reg,reg_m,reg_s)

######################################################################################
# Results

options(scipen=999)
summary(pars_ls$pars)$statistics[1:pars_ls$index[1],]
plot(pars_ls$pars[,c(1,2,pars_ls$index[2],pars_ls$index[2]+1)])

states=as.matrix(pars_ls$pars)[,(pars_ls$index[1]+1):pars_ls$index[3]]
plot(colMeans(states),ylab="Mean state")
plot(round(colMeans(states),0),ylab="Most likely state",main="True states in red, estimated states in black")
points(S+1.01,col='red')
