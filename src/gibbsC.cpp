#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param x0 initialization
//' @param y0 initialization
//' @param n default is 10
//' @param a default is 1
//' @param b default is 2
//' @return a random sample of size n
//' @examples
//' \dontrun{
//' N <- 5000
//' x0 <- 0.5*n
//' y0 <- (x0+a) / (n+a+b)
//' Chain_Cpp <- gibbsC(N, x0, y0)
//' }
//' @import microbenchmark
//' @importFrom Rcpp evalCpp
//' @importFrom stats rbinom rbeta
//' @useDynLib WAP
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, double x0, double y0, int n=10, int a=1, int b=2){
  NumericMatrix mat(N, 2);
  double x = x0, y = y0;
  mat(0, 0) = x;
  mat(0, 1) = y;
  
  for(int i = 1; i < N; i++){
    y = mat(i-1, 1);
    mat(i, 0) = rbinom(1, n, y)[0];
    x = mat(i, 0);
    mat(i, 1) = rbeta(1, x+a, n-x+b)[0];
  }
  return(mat);
}
