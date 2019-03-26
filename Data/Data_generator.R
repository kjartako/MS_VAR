library(MASS)
# Simulated data
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

write.table(Y,file=paste("Y_",n,".txt",sep=""),sep="\t", col.names = F, row.names = F)
write.table(S,file=paste("S_",n,".txt",sep=""),sep="\t", col.names = F, row.names = F)

