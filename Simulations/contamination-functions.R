contaminacaoA=function(eta.mu,eta.phi,ep,seed){
  linkobj=robustbetareg:::set.link()
  mu=linkobj$linkfun.mu$inv.link(eta.mu)
  phi=linkobj$linkfun.phi$inv.link(eta.phi)
  n=length(mu) 
  set.seed(seed)
  result=list()
  ind.c=0

  y=y2=rbeta(n,mu*phi,(1-mu)*phi)
  if(ep>0){
    ind.c=which(mu %in% sort(mu,decreasing = FALSE)[1:ep])
    mu_l=(1+mu[ind.c])/2 
    phi_l=phi[ind.c] 
    y[ind.c]=rbeta(ep,mu_l*phi_l,(1-mu_l)*phi_l)
  }
  result$sample=pmax(pmin(y,1-.Machine$double.eps),.Machine$double.eps)
  result$ind=ind.c
  result$y=y2 
  return(result)
}

contaminacaoB=function(eta.mu,eta.phi,ep,seed,x){
  linkobj=robustbetareg:::set.link()
  mu=linkobj$linkfun.mu$inv.link(eta.mu)
  phi=linkobj$linkfun.phi$inv.link(eta.phi)
  n=length(mu) 
  set.seed(seed)
  result=list()
  ind.c=0

  y=y2=rbeta(n,mu*phi,(1-mu)*phi)
  if(ep>0){
    ind.c=which(x[,2] %in% sort(x[,2],decreasing = F)[1:ep])
    mu_l=0.001
    phi_l=phi[ind.c] 
    y[ind.c]=rbeta(ep,mu_l*phi_l,(1-mu_l)*phi_l)
  }
  result$sample=pmax(pmin(y,1-.Machine$double.eps),.Machine$double.eps)
  result$ind=c(ind.c) 
  result$y=y2 
  return(result)
}

contaminacaoC=function(eta.mu,eta.phi,ep,seed){
  linkobj=robustbetareg:::set.link()
  mu=linkobj$linkfun.mu$inv.link(eta.mu)
  phi=linkobj$linkfun.phi$inv.link(eta.phi)
  n=length(mu) 
  set.seed(seed)
  result=list()
  y=y2=rbeta(n,mu*phi,(1-mu)*phi)
  if(ep>0){
    ind.c=which(mu %in% sort(mu,decreasing = T)[1:ep])
    mu_l=0.05
    phi_l=20 
    y[ind.c]=rbeta(ep,mu_l*phi_l,(1-mu_l)*phi_l)
  }
  result$sample=pmax(pmin(y,1-.Machine$double.eps),.Machine$double.eps)
  result$ind=c(ind.c) 
  result$y=y2 
  return(result)
}
