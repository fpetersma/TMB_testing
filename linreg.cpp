// Simple linear regression.
#include <TMB.hpp>

template<class Type>

Type objective_function<Type>::operator() () 
{
  // Add data macros
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  
  // Add parameter macros
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  
  // whatever this does
  ADREPORT(exp(2*logSigma));
  
  // Version 1
  // Define the negative log likelihood
  // Type nll = -sum(dnorm(Y, a+b*x, exp(logSigma), true));
  
  // Version 2
  int Y_size = Y.size();

  vector<Type> out(Y_size);

  for(int i = 0; i < Y_size; i++) {
    out(i) = dnorm(Y(i), a + b * x(i), exp(logSigma), true);
  }

  Type nll = -sum(out);
  
  // Return the nll
  return nll;
}