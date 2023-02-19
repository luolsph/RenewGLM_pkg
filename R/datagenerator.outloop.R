datagenerator_out <-
function(beta, b, n, family, construct, rho, intercept = TRUE, categorical=FALSE, seed=NA){
  
  if(!is.na(seed)){ set.seed(b+seed) }
  set.seed(b+123456)
  
  p <- length(beta)
  if(construct=='ind'){Sigma<- diag(rep(1, p))}
  if(construct=='cs'){Sigma<-matrix(rho, p, p); diag(Sigma)<-1}
  if(construct=='ar1'){Sigma<-rho^(abs(outer(1 : p, 1 : p, "-")))}
  X <- mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma)
  
  #change half of the X to be categorical variables
  if(categorical==TRUE){some<-sample(1 : p,round(p/2),replace=FALSE); X[,some] <- (X[,some]>0)}
  
  if(intercept == TRUE){X[, 1] <- 1}
  
  if(family=="gaussian"){ y <- X %*% beta + rnorm(n) }
  if(family=="binomial"){ y <- rbinom(n, 1, prob=exp(X%*%beta)/(exp(X%*%beta)+1) ) }
  if(family=="poisson"){y <- rpois(n, exp(X%*%beta)) }
  
  simdata <- cbind(y, X)
  return(simdata)
}
