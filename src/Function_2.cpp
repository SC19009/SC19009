#include <Rcpp.h>
using namespace Rcpp;

//' @title  Jackknife estimator for mean and its variance
//' @description Jackknife estimator for mean and its variance using Rcpp
//' @param x a vector composed of sample points
//' @return a 2-dimension vector represent the jackknife estimate of mean and its variance
//' @examples
//' \dontrun{
//' x<-rnorm(1000,5,5)
//' jack(x)
//' }
//' @export
// [[Rcpp::export]]
NumericVector jack(NumericVector x){
  int n=sizeof(x); double s=sum(x); double bar0=s/n;
  double bar[n]; double y[n]; double sum1=0; double sum2=0;
  
  for (int i=0;i<=n-1;i++){
    bar[i]=(s-x[i])/(n-1);
    sum1=sum1+bar[i];
  }
  
  double bar1=sum1/n;
  double jack=n*bar0-(n-1)*bar1;
  
  for (int i=0;i<=n-1;i++){
    double z=bar[i]-bar1;
    y[i]=z*z;
    sum2=sum2+y[i];
  }
  
  double var=(n-1)*sum2/n;
  
  NumericVector out= NumericVector::create(Named("jack",jack),Named("var",var));
  return (out);
}