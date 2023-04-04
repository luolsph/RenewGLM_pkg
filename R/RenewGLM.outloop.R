RenewGLM_out <-
function(X, y, type, betahat, infomats, intercept, s, phi){
  tol <- 1e-6;
  max_iter <- 100;
  
#  if(intercept==TRUE){X<-cbind(1,X)};
  
  p <- dim(as.matrix(X))[2];
  
  betahat_old <- betahat;
  
  #record W in a list form rather than a matrix to simplify calculation
  W <- invlinkdiv(X, betahat_old, type=type)
  H <- cp(X,y,W)
  
  U=chol(infomats+H)
  L=t(U)
  
  for (r in 1:max_iter){
    
    g_0=t(crossprod((y-invlink(X, betahat, type)),X));
    g_1=-t(crossprod((betahat-betahat_old),infomats));
    g <- g_0 + g_1;
    
    d_beta <- backsolve(U,forwardsolve(L,g))
    
    df_beta <- crossprod(g,d_beta);
    if (abs(df_beta)<tol){
      break
    }else {
      betahat <- betahat + drop(d_beta);
    }
  }
  if(type=="gaussian"){
    s <- s+length(y)
    phi<-(s-length(y)-p)/(s-p)*phi+1/(s-p)*t(betahat_old)%*%infomats%*%(betahat_old-betahat)+t(y)%*%(y-X%*%betahat)/(s-p);
  }
  
  W <- invlinkdiv(X,betahat,type=type)
  H_new <- cp(X,y,W)
  infomats<-infomats+H_new
  
  if(type=="gaussian"){
    return(list(betahat,infomats,s,phi))
  }else{
    return(list(betahat,infomats))
  }
}

cp <-
function(X, y, w){
  if (length(y)==1){
    H <- w*tcrossprod(X,X)
  }else{
    H <- crossprod(sqrt(w)*X)
  }
}

invlink <-
function(X, beta, type){
#  if (length(beta)==1){
#    eta <- X * beta
#  }else{
    eta <- drop(X %*% beta)
#  }
  if(type=="gaussian"){ out <- eta }
  if(type=="binomial"){ out <- exp(eta)/(1 + exp(eta)) }
  if(type=="poisson") { out <- exp(eta) }
  out
}

invlinkdiv <-
function(X, beta, type){
#  if (length(beta)==1){
#    eta <- X * beta
#  }else{
    eta <- drop(X %*% beta)
#  }
  if(type=="gaussian"){ out <- 1 }
  if(type=="binomial"){ out <- exp(eta)/(1 + exp(eta))^2 }
  if(type=="poisson") { out <- exp(eta) }
  out
}

