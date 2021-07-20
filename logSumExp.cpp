#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

// from: https://stackoverflow.com/questions/45943750/calculating-log-sum-exp-function-in-c

template <typename Iter>
std::iterator_traits<Iter>::value_type
  logSumExp(Iter begin, Iter end)
  {
    using VT = std::iterator_traits<Iter>::value_type{};
    if (begin==end) return VT{};
    using std::exp;
    using std::log;
    auto max_elem = *std::max_element(begin, end);
    auto sum = std::accumulate(begin, end, VT{}, 
                               [max_elem](VT a, VT b) { return a + exp(b - max_elem); });
    return max_elem + log(sum);
  }