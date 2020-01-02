#' @title BCa confidence interval
#' @useDynLib SC19009
#' @description A function calculating BCa confidence interval for population mean using R
#' @import Rcpp
#' @importFrom stats pnorm qnorm quantile
#' @param x a vector composed of random sample points
#' @param alpha confidence level
#' @param B number of bootstrap replicates
#' @return a named vector representing the confidence interval whose first and second elements are respectively lower bound and upper bound
#' @examples
#' \dontrun{
#' x<-rt(200,3)
#' BCa(x,0.95,500)
#' }
#' @export
BCa<-function(x,alpha,B){
  n<-length(x)
  bar<-mean(x)
  
  bar.boot<-numeric(B)
  for (i in 1:B) {
    j<-sample(1:n,n,replace=TRUE)
    bar.boot[i]<-mean(x[j])
  }
  
  bar.jack<-(sum(x)-x)/(n-1)
  bar0.jack<-mean(bar.jack)
  
  z<-qnorm((1-alpha)/2)
  z0<-qnorm(mean(bar.boot<bar))
  alpha0<-sum((bar0.jack-bar.jack)^3)/(6*sum(abs(bar0.jack-bar.jack)^3))
  
  alpha1<-pnorm(z0+(z0+z)/(1-alpha0*(z0+z)))
  alpha2<-pnorm(z0+(z0-z)/(1-alpha0*(z0-z)))
  
  lb.BCa<-quantile(bar.boot,alpha1,names=FALSE)
  ub.BCa<-quantile(bar.boot,alpha2,names=FALSE)
  
  out<-c(lb.BCa,ub.BCa)
  names(out)<-c("lb","ub")
  return(out)
}