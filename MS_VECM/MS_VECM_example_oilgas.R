library(coda)
library(rstan)
options(mc.cores = parallel::detectCores())

chains=8
it_burnin = 500
it_sample = 1000

sm = stan_model('Stan/MS_VECM_wishart.stan',auto_write=TRUE)

######################################################################################
Y_0 = log(unname(as.matrix(read.table("Data/Oilgas.txt",sep="\t"))))
Y = diff(Y_0)
n=nrow(Y)
z=(Y_0[2:(n+1),1]-Y_0[2:(n+1),2])
######################################################################################
n_AR=1
n_reg=2
reg_m=TRUE
reg_s=TRUE
reg_a=TRUE

stanfit = sampling(sm,data=list(T=n,y=Y,z=z,ARdim=n_AR,nreg=n_reg,mean_reg=reg_m,sigma_reg=reg_s,alpha_reg=reg_a,Q_alpha=as.array(rep(1,n_reg)),
                   mu_mean=0,mu_sd=0.2,phi_mean=0,phi_sd=1,alpha_mean=0,alpha_sd=0.2,eye=diag(2),L=chol(30*diag(2)),nu=3),
                   iter=(it_burnin+it_sample), warmup = it_burnin, chains=chains,
                   pars=c("Q","mu","phi","alpha","sigma","S"))

######################################################################################
# Label-switching

pars =mcmc.list(lapply(1:ncol(stanfit), function(x) mcmc(as.array(stanfit)[,x,])))

source("labsw.R")
pars_ls = labsw_VECM(pars,chains,n_AR,n_reg,reg_m,reg_s)

######################################################################################
# Results

options(scipen=999)
summary(pars_ls$pars)$statistics[1:pars_ls$index[1],]
plot(pars_ls$pars[,c(1,2,pars_ls$index[2],pars_ls$index[2]+1)])

states=as.matrix(pars_ls$pars)[,(pars_ls$index[1]+1):pars_ls$index[3]]
plot(colMeans(states),ylab="Mean state")
plot(round(colMeans(states),0),ylab="Most likely state")
