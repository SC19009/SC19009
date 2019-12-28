#' @title Kernel density estimator
#' @useDynLib SC19009
#' @description A kernel density estimator using R
#' @import Rcpp
#' @importFrom stats density
#' @importFrom graphics hist lines
#' @param x a vector composed of random sample points
#' @return a histogram for sample frequency whith the estimated density curve
#' @examples
#' \dontrun{
#' x<-rt(2000,3)
#' kde(x)
#' }
#' @export
kde<-function(x){
  hist(x,xlab="x",ylab="f(x)",main="Kernel density estimation",freq=F)
  f<-density(x)
  lines(f,xlab="x",ylab="f(x)",col="blue",sub="")
}