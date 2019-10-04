RenewGLM.inloop <-
function(B, tempdatadir, type, init=NA, p, intercept=TRUE){
  sum2<-diag(0,p,p);
  s<-0
  phi<-0
  tol=1e-6;
  max_iter=100;
  
  if(is.na(init)){init<-rep(0,p)}
  betahat<-init;
  
  for(b in 1:B){
    X<-y<-NULL;
    load(paste(tempdatadir,"/Simdata",b,".RData",sep=""))
    if(intercept==TRUE){X<-cbind(1,X)}
    
    betahat_old=betahat;
    
    #record W in a list form rather than a matrix to simplify calculation
    W<-invlinkdiv(X,betahat_old,type=type)
    H<-cp(X,y,W)
    
    U=chol(sum2+H)
    L=t(U)
    
    for (r in 1:max_iter){
      
      g_0=t(crossprod((y-invlink(X, betahat, type)),X));
      g_1=-t(crossprod((betahat-betahat_old),sum2));
      g=g_0+g_1;
      
      d_beta=backsolve(U,forwardsolve(L,g))
      
      df_beta=crossprod(g,d_beta);
      if (abs(df_beta)<tol){
        break
      }else {
        betahat=betahat+d_beta;
      }
    }
    if(type=="gaussian"){
      s<-s+length(y)
      phi<-(s-length(y)-p)/(s-p)*phi+1/(s-p)*t(betahat_old)%*%sum2%*%(betahat_old-betahat)+t(y)%*%(y-X%*%betahat)/(s-p);
    }
    
    W<-invlinkdiv(X,betahat,type=type)
    H_new<-cp(X,y,W)
    sum2<-sum2+H_new
  }
  
  if (type=="gaussian"){
    sd<-sqrt(diag(solve(sum2))*phi);
  }else{
    sd<-sqrt(diag(solve(sum2)));
  }
  pvalue<-2*pnorm(-abs(betahat)/sd)
  result<-cbind(betahat=betahat,sd=sd,pvalue=pvalue)
  colnames(result)<-c("Estimates","Std.Errors","p-values")
  return(result)
  
}

cp <-
function(X, y, w){
  if (length(y)==1){
    H<-w*tcrossprod(X,X)
  }else{
    H<-crossprod(sqrt(w)*X)
  }
}

invlink <-
function(X, beta, type){
  if (length(beta)==1){
    eta<-X*beta
  }else{
    eta<-drop(X%*%beta)
  }
  if(type=="gaussian"){ out <- eta }
  if(type=="binomial"){ out <- exp(eta)/(1 + exp(eta)) }
  if(type=="poisson") { out <- exp(eta) }
  out
}

invlinkdiv <-
function(X, beta, type){
  if (length(beta)==1){
    eta<-X*beta
  }else{
    eta<-drop(X%*%beta)
  }
  if(type=="gaussian"){ out <- 1 }
  if(type=="binomial"){ out <- exp(eta)/(1 + exp(eta))^2 }
  if(type=="poisson") { out <- exp(eta) }
  out
}


