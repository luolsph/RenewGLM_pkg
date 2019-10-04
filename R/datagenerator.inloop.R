datagenerator_in <-
function(beta, n, p, B, family, construct, rho, tempdatadir, categorical=FALSE, seed=NA){
  
  if(length(n)!=1 & length(n)!=B){ stop("n must be a number or a vector of length B.") }
  if(!is.na(seed)){ set.seed(seed) }
  if(length(n)==1){ n <- rep(n,B) } #same sample size
  
  dir.create(tempdatadir)
  for(b in 1:B){
    set.seed(b+123456)
    if(construct=='ind'){Sigma<- diag(rep(1,p-1))}
    if(construct=='cs'){Sigma<-matrix(rho,p-1,p-1); diag(Sigma)<-1}
    if(construct=='ar1'){Sigma<-rho^(abs(outer(1:(p-1), 1:(p-1), "-")))}
    X <- mvrnorm(n=n[b], mu=rep(0,p-1), Sigma=Sigma)
    
    #change half of the X to be categorical variables
    if(categorical==TRUE){some<-sample(1:(p-1),round((p-1)/2),replace=FALSE); X[,some]<-(X[,some]>0)}
    
    X<-cbind(1,X)
    
    if(family=="gaussian"){ y <- X%*%beta + rnorm(n[b]) }
    if(family=="binomial"){ y <- rbinom(n = n[b], size = 1, prob=exp(X%*%beta)/(exp(X%*%beta)+1)) }
    if(family=="poisson") { y <- rpois(n = n[b], lambda = exp( X%*%beta )) }
    X <-X[,-1]
    
    save(y, X, file = paste(tempdatadir, "/Simdata", b, ".RData", sep=""))
  }
}
