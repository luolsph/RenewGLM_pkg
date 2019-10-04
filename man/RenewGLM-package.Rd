\name{RenewGLM-package}
\alias{RenewGLM-package}
\alias{RenewGLM}
\docType{package}
\title{
\packageTitle{RenewGLM}
}
\description{
\packageDescription{RenewGLM}
}
\details{
This package aims to update the regression coefficients as data batches arrive sequentially. There are two main functions in package. \code{RenewGLM.inloop} is used to processing a sequence of datasets that are stored in a given directory, and the major input is the name of the directory. \code{RenewGLM.outloop} is applied in the case where data batches are imported externally, and the form of the input is \code{X} and \code{y} of the current data batch.
}

\author{
\packageAuthor{RenewGLM} \cr
Maintainer: \packageMaintainer{RenewGLM}
}
%\references{
%~~ Literature or other references for background information ~~
%}

\keyword{Renewable estimator}
\keyword{Streaming datasets}
\keyword{Generalized linear models}

\examples{
#Processing data batches internally
N=1000
B=10
p=5
n=N/B
beta<-c(0.2,-0.2,0.2,-0.2,0.2)

tempdatadir<-"~/Desktop/tempdata"
datagenerator.inloop(beta=beta,n=n, p=p, B=B, family="binomial", construct="cs",
    rho=0.5, tempdatadir=tempdatadir)
RenewGLM.inloop(B, tempdatadir=tempdatadir, "binomial", p=p, intercept=TRUE)
unlink(tempdatadir)

#Processing data batches externally
N=1000
B=10
p=5
n=N/B
beta<-c(0.2,-0.2,0.2,-0.2,0.2)
infomats<-diag(0,p,p);
betahat<-rep(0,p)
for(b in 1:10){
 data<-datagenerator.outloop(beta,b,n,"binomial","cs",0.5)
 y<-data[,1]
 X<-data[,-1]
 summary<-RenewGLM.outloop(X,y,"binomial",betahat,infomats,intercept=TRUE,s,phi)
 betahat<-summary[[1]]
 infomats<-summary[[2]]
 rm(data)
}
sd<-sqrt(diag(solve(infomats)));
pvalue<-2*pnorm(-abs(betahat)/sd)
result<-cbind(betahat=betahat,sd=sd,pvalue=pvalue)
colnames(result)<-c("Estimates","Std.Errors","p-values")

}
