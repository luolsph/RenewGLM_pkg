% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datagenerator_out.R
\name{datagenerator.outloop}
\alias{datagenerator.outloop}
\title{This function is used to generate one data batch each time to illustrate the usage of RenewGLM function}
\usage{
datagenerator_out(beta, b, n, family, construct, rho, categorical = FALSE,
  seed = NA)
}
\arguments{
\item{beta}{designed coefficients for simulated data batches, including the intercept}

\item{b}{index for the current data batch}

\item{n}{sample size of the bth data batch}

\item{family}{exponential family of the reponses, including \code{c("gaussian","binomial","poisson")}}

\item{construct}{structure of covariance matrix for generating the covariate matrix \code{X}}

\item{rho}{the correlation coefficient in the covariance matrix, choices inlcude \code{c("ind", "cs", "ar1")}}

\item{categorical}{set to TRUE to let some \code{X} be dichotomized, the default is \code{FALSE}}

\item{seed}{set random seed and the default is \code{NA}}
}
\value{
a data matrix containing response vector \code{y} and covariate matrix \code{X}
}
\description{
This function is used to generate one data batch each time to illustrate the usage of \code{RenewGLM.outloop} function
}
