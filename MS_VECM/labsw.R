labsw_VECM = function(pars,chains,n_AR,n_reg,reg_m,reg_s,beta=FALSE)
{
  # Reorder Q[2,1] and Q[2,2], to make permutation easier
  for(i in 1:chains)
  {
    pars[[i]][,c(2,4)]=pars[[i]][,c(4,2)]
    dimnames(pars[[i]])$parameters[c(2,4)]=dimnames(pars[[i]])$parameters[c(4,2)]
  }
  # Identify parameters indexes
  n_m= max(1,reg_m*n_reg)
  n_s= max(1,reg_s*n_reg)
  
  q_range = c(1:n_reg^2)
  mean_range = c(1:(2*(n_m*2)+(n_m*n_AR*2^2)))+max(q_range)
  sigma_range= c(1:(n_s*2^2))+max(mean_range) 
  
  # Identify which parameters to permutate
  permu_range = q_range
  if(reg_m){permu_range=c(permu_range,mean_range)}
  if(reg_s){permu_range=c(permu_range,sigma_range)}
  
  index_lastpar=max(sigma_range)+beta
  index_laststate=index_lastpar+n-n_AR
  
  # Impose identifiability constraint, sigma[1,1,1]<sigma[2,1,1]
  index_s=index_lastpar-beta-(n_s*2^2)+1
  for(i in 1:chains)
  {
    chainvals=pars[[i]]
    
    if(chainvals[1,index_s]>chainvals[1,index_s+1])
    {
      # permute parameters
      for(j in permu_range[seq(from=1,to=length(permu_range),by=2)])
      {
        pars[[i]][,c(j,j+1)]=chainvals[,c(j+1,j)]
      }
      # permute states
      tmp=pars[[i]][,(index_lastpar+1):index_laststate]
      tmp[(pars[[i]][,(index_lastpar+1):index_laststate])==1]=2
      tmp[(pars[[i]][,(index_lastpar+1):index_laststate])==2]=1
      
      pars[[i]][,(index_lastpar+1):index_laststate]=tmp
    }
  }
  
  return(list(pars=pars,index=c(index_lastpar,index_s,index_laststate)))
}