######################################################################################
######################################################################################
######################################################################################
# Simulated data

library(MASS)
mu = c(0,0);mu2 = c(0.01,0);rho = 0.2;rho2 = 0.2;sd=c(0.15,0.15);sd2=c(0.3,0.25)
sigma = matrix(c(sd[1]^2,rep(rho*sd[1]*sd[2],2),sd[2]^2),nrow=2,byrow=TRUE)
sigma2 = matrix(c(sd2[1]^2,rep(rho2*sd2[1]*sd2[2],2),sd2[2]^2),nrow=2,byrow=TRUE)
phi = matrix(c(0.6,0.3,0.01,0.98),nrow=2,byrow=TRUE)
phi2 = matrix(c(0.98,0.01,0.01,0.99),nrow=2,byrow=TRUE)
Q = matrix(c(0.97,0.03,0.1,0.9),nrow=2,byrow=TRUE)
delta = solve(t(diag(2)-Q +1),rep(1,2))

n=750
S = numeric(n)
Y = matrix(rep(NA,2*n),nrow=n,byrow=TRUE)

set.seed(123)
S[1] = sample(c(0,1),size=1,prob=delta)
Y[1,] = mvrnorm(n=1,mu=(1-S[1])*mu+S[1]*mu2,Sigma=(1-S[1])*sigma+S[1]*sigma2)

for (i in 2:n)
{
  S[i]=sample(c(0,1),size=1,prob=Q[S[i-1]+1,])
  Y[i,] =mvrnorm(n=1,mu=(1-S[i])*(mu+phi%*%Y[i-1,])+S[i]*(mu2+phi2%*%Y[i-1,]),Sigma=(sigma*(1-S[i])+sigma2*S[i]))
}

######################################################################################
######################################################################################
######################################################################################

library(coda)
library(rstan)
options(mc.cores = parallel::detectCores())

chains=8
it_burnin = 500
it_sample = 1000

sm = stan_model('Stan/MS_VAR_wishart.stan',auto_write=TRUE)
#sm = stan_model('Stan/MS_VAR_lkj.stan',auto_write=TRUE)

######################################################################################

n_dim=ncol(Y)
n_AR=1
n_reg=2
reg_m=TRUE
reg_s=TRUE

stanfit = sampling(sm,data=list(T=n,dim=n_dim,y=Y,ARdim=n_AR,nreg=n_reg,mean_reg=reg_m,sigma_reg=reg_s,Q_alpha=as.array(rep(1,n_reg)),
                   mu_mean=0,mu_sd=0.2,phi_mean=0,phi_sd=1,eye=diag(n_dim),L=chol(diag(n_dim)),nu=n_dim+1),
                   iter=(it_burnin+it_sample), warmup = it_burnin, chains=chains,
                   pars=c("Q","mu","phi","sigma","S"))

# stanfit = sampling(sm,data=list(T=n,dim=n_dim,y=Y,ARdim=n_AR,nreg=n_reg,mean_reg=reg_m,sigma_reg=reg_s,Q_alpha=as.array(rep(1,n_reg)),
#                    mu_mean=0,mu_sd=0.2,phi_mean=0,phi_sd=1,eta=1,gamma_alpha=2,gamma_beta=0.1),
#                    iter=(it_burnin+it_sample), warmup = it_burnin, chains=chains,
#                    pars=c("Q","mu","phi","sigma","S"))

######################################################################################
# Label-switching

pars =mcmc.list(lapply(1:ncol(stanfit), function(x) mcmc(as.array(stanfit)[,x,])))

source("labsw_2reg.R")
pars_ls = labsw_2reg(pars,chains,n_dim,n_AR,n_reg,reg_m,reg_s)

######################################################################################
# Results

options(scipen=999)
summary(pars_ls$pars)$statistics[1:pars_ls$index[1],]
plot(pars_ls$pars[,c(1,2,pars_ls$index[2],pars_ls$index[2]+1)])

states=as.matrix(pars_ls$pars)[,(pars_ls$index[1]+1):pars_ls$index[3]]
plot(colMeans(states),ylab="Mean state")
plot(round(colMeans(states),0),ylab="Most likely state",main="True states in red, estimated states in black")
points(S+1.01,col='red')
