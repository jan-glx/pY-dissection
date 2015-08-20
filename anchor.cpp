#include <Rcpp.h>
#include <string>
using namespace Rcpp;

extern "C" {
#include "anchor.c"
}

// [[Rcpp::export]]
std::string anchor(std::string seq) {
  const char* seq_name = "BLA";
  
  anchor_(seq_name, seq.c_str(), seq.length(), 1, 0);
  return out;
}