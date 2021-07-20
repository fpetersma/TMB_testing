// The first line includes the source code for the whole TMB package (and all its 
// dependencies). The objective function is a templated class where <Type> is 
// the data type of both the input values and the return value of the objective 
// function. This allows us to evaluate both the objective function and its 
// derivatives using the same chunk of C++ code (via the AD package CppAD). The 
// technical aspects of this are hidden from the user. There is however one 
// aspect that surprises the new TMB user. When a constant like “1.2” is used
// in a calculation that affects the return value it must be “cast” to Type.

// R code             C++/TMB code
// 
// Comments      #                  //                          // Comment symbol
// Constants     3.4                Type(3.4);                  // Explicit casting recommended in TMB
// Scalar        x = 5.2            Type x = Type(5.2);         // Variables must have type
// Arrays        x = numeric(10)    vector<Type> x(10);         // C++ code here does NOT initialize to 0
// Indexing      x[1]+x[10]         x(0)+x(9);                  // C++ indexing is zero-based
// Loops         for(i in 1:10)     for(int i=1;i<=10;i++)      // Integer i must be declared in C++
// Increments    x[1] = x[1] + 3    x(0) += 3.0;                // += -= *= /= incremental operators in C++


// pi = M_PI
// inf = INFINITY
// -inf = -INFINITY


// Simple ASCR with constant density.
#include <TMB.hpp>
// #include <fstream> // doesnt work
#include <cmath>
#include <math.h>
// #include <tiny_ad/bessel/bessel.h>
// #include <iostream>


template<class Type>
Type objective_function<Type>::operator() () 
{
  // Add data macros
  DATA_MATRIX(W);               // detection histories
  DATA_MATRIX(R);               // received levels
  DATA_MATRIX(Y_rec);           // recorded bearings
  DATA_MATRIX(Y_grid);          // grid bearings
  DATA_MATRIX(X);               // distances
  DATA_VECTOR(A);               // areas of grid cells of mesh
  DATA_SCALAR(trunc_level);     // truncation level for received levels
  // DATA_MATRIX(design_matrix);   // the design matrix from the GAM object
  
  // Add parameter macros
  PARAMETER(log_D);
  PARAMETER(logit_g0); 
  PARAMETER(log_beta_r);
  PARAMETER(log_sd_r);
  PARAMETER(log_mu_s);
  
  // PARAMETER(log_kappa);
  
  PARAMETER(log_kappa_low);
  PARAMETER(log_kappa_high);
  PARAMETER(logit_mix_bear);
  
  // Convert to real scale
  Type D = exp(log_D);
  Type g0 = exp(logit_g0) / (exp(logit_g0) + 1);
  Type beta_r = exp(log_beta_r);
  Type sd_r = exp(log_sd_r);
  Type mu_s = exp(log_mu_s);
  
  // Type kappa = exp(log_kappa);
  
  Type kappa_low = exp(log_kappa_low);
  Type kappa_high = exp(log_kappa_high) + kappa_low;
  Type mix_bear = exp(logit_mix_bear) / (exp(logit_mix_bear) + 1);
  
  
  // Set constants
  // Type pi = Type(3.14159265358979323846);
  
  // whatever this does
  // ADREPORT(exp(log_D));
  // ADREPORT(exp(2 * log_sd_r));
  // ADREPORT(exp(logit_g0) / (exp(logit_g0) + 1));
  // ADREPORT(exp(log_kappa));
  // ADREPORT(exp(log_beta_r);
  // ADREPORT(exp(log_mu_s));
  
  // Extract dimensions
  int n_call = W.rows();
  int n_det = W.cols();
  int n_grid = X.rows();
  
  std::cout << "Number of observations: " << n_call << std::endl;
  std::cout << "Number of detectors: " << n_det << std::endl; 
  std::cout << "Grid size of the mesh: " << n_grid << std::endl;
  std::cout << "kappa_high parameter: " << kappa_high << std::endl;
  
  // std::cout<< "1 + INFINITY: " << 1 + INFINITY << std::endl; // Is still inf!
  
  // std::cout << "test: " << - pow(10, 16) << std::endl;
  
  std::cout << "X(0, 0): " << X(0, 0) << std::endl;
  
  // ===========================================================================
  // The non-conditional part of the log likelihood
  // ===========================================================================
  
  // Define the likelihood of n_calls with constant density
  vector<Type> lambda_vector(n_grid);

  for (int m = 0; m < n_grid; m++) { // loop through the mesh
    vector<Type> probs_m(n_det);
    vector<Type> non_probs_m(n_det);
    for (int j = 0; j < n_det; j++) { // loop through the detectors
      // std::cout << j << std::endl; // works
      
      Type expected_rl = mu_s - beta_r * log10(X(m, j));
      // std::cout << expected_rl << std::endl; // works
      
      Type xx = (trunc_level - expected_rl) / sd_r;
      Type p = g0 * (Type(1) - pnorm(xx, Type(0), Type(1)));
      probs_m(j) = p;
      // std::cout << p << std::endl; // works
      
      non_probs_m(j) = Type(1) - probs_m(j);
    }
    
    // std::cout << "probs_m(0): " << probs_m(0) << std::endl; //  works!
  
    // Get probability of at least two detections
    vector<Type> v1(n_det + 1);
    for (int a = 0; a < n_det; a++) { // derive all probabilities of one detection
      // Vector for configurations of a single detection
      vector<Type> config_vector(n_det);
      // Vector to store the probabilities associated with v2 configurations
      vector<Type> probability_vector(n_det); 
      for (int iii = 0; iii < n_det; iii++) {
        config_vector(iii) = 0; // fill the vector v2 with zeros
      }
      config_vector(a) = 1; // set one position equal to
      // std::cout << "config_vector: " << config_vector << std::endl;
      for (int ii = 0; ii < n_det; ii++) {
        if (config_vector(ii) == 1) probability_vector(ii) = probs_m(ii);
        if (config_vector(ii) == 0) probability_vector(ii) = non_probs_m(ii);
      }
      // std::cout << "probability_vector: " << probability_vector << std::endl;
      // std::cout << "probs_m: " << probs_m << std::endl;
      // std::cout << "non_probs_m: " << non_probs_m << std::endl;
      v1(a) = probability_vector.prod();
    }
    v1(n_det) = non_probs_m.prod(); // probability of no detections
    Type p_twice = 1 - v1.sum(); // p_twice is complement of p_none and p_once
    // std::cout << "p_twice: " << p_twice << " and v1: " << v1 << std::endl;

    // store lambda given x in lambda vector
    lambda_vector(m) = D * p_twice * A(m);
  }

  Type lambda = lambda_vector.sum();
  
  std::cout << "The value for lambda is: " << lambda << std::endl;
  
  // Type llk_n = dpois(Type(n_call), lambda, true); 
  Type llk_n = -lambda; // Only need part due to logarithm in llk
  
  if (llk_n  == -INFINITY) llk_n = -pow(10, 16);
  if (llk_n  == INFINITY) llk_n = pow(10, 16);
  
  // ===========================================================================
  // The conditional part of the log likelihood 
  // ===========================================================================
  
  // Define the conditional negative log likelihood
  Type llk_cond = Type(0.0); // define vector to store likelihoods per individual call

  for (int i = 0; i < n_call; i++) { // loop over the calls
    vector<Type> log_probs_i(n_grid);
    for (int m = 0; m < n_grid; m++) { // loop through the mesh
      // Likelihood of capture histories
      Type log_prob_w_im = Type(0.0);

      for (int j = 0; j < n_det; j++) {
        Type expected_rl = mu_s - beta_r * log10(X(m, j));
        Type xx = (trunc_level - expected_rl) / sd_r;
        Type prob = g0 * (Type(1) - pnorm(xx, Type(0), Type(1)));
        
        Type test = Type(1) * W(i, j); // Multiply W(i, j)  by one so to not use DATA directly
        
        if (test == Type(1)) log_prob_w_im += log(prob); 
        else log_prob_w_im += log(Type(1) - prob);
      }
      if (log_prob_w_im == -INFINITY) log_prob_w_im = -pow(10, 16);
      if (log_prob_w_im == INFINITY) log_prob_w_im = pow(10, 16);

      // Likelihood of bearings
      Type log_prob_y_im = Type(0.0);
      
      for (int j = 0; j < n_det; j++) {
        if (W(i, j) == 1) {
          Type obs_minus_exp = Y_rec(i, j) - Y_grid(m, j);
          if (false) { // REPLACE LATER FOR MIXTURE CHECK
            // log_prob_y_im += kappa * cos(obs_minus_exp) - 
            //   log(2 * M_PI * besselI(kappa, Type(0)));// log version
          } else {
            log_prob_y_im += log(mix_bear * exp(kappa_low * cos(obs_minus_exp)) / 
              (2 * M_PI * besselI(kappa_low, Type(0))) +
              (1 - mix_bear) * exp(kappa_high * cos(obs_minus_exp)) / 
              (2 * M_PI * besselI(kappa_high, Type(0))));
          }
        }
      }
      if (log_prob_y_im == -INFINITY) log_prob_y_im = -pow(10, 16);
      if (log_prob_y_im == INFINITY) log_prob_y_im = pow(10, 16);
      
      // Likelihood of received levels
      Type log_prob_r_im = Type(0.0); // sum the log probabilities
      
      // Loop through all detectors and sum the log probabilities
      for (int j = 0; j < n_det; j++) {
        if (W(i, j) == 1) {
          Type expected_rl = mu_s - beta_r * log10(X(m, j));
          Type recorded_rl = Type(1) * R(i, j);
  
          // Probability of exp_rl - rec_rl
          Type log_prob = dnorm(recorded_rl, expected_rl, sd_r, 1); // log probabilities
          // std::cout << "log_prob1: " << log_prob << std::endl;
          
          log_prob -= log(sd_r * (Type(1.0) - pnorm(trunc_level, expected_rl, sd_r))); // DIFFERENT FROM R CODE
          
          // std::cout << "log_prob2: " << log_prob << std::endl;
          
          log_prob_r_im += log_prob;
        }
      }
      // Set log_prob_r_im to -pow(10, 16) if it is -inf
      if (log_prob_r_im == -INFINITY) log_prob_r_im = -pow(10, 16);
      if (log_prob_r_im == INFINITY) log_prob_r_im = pow(10, 15);
      
      
      // std::cout << "log(D): " << log(D) << std::endl; // works
      // std::cout << "prob_w_im: " << prob_w_im << std::endl; // works for now
      // std::cout << "prob_y_im: " << prob_y_im << std::endl; // works
      // std::cout << "prob_r_im: " << prob_r_im << std::endl; // works
      // log likelihood of a single call given a single location
      Type cond_llk_im = log(A(m)) + log_D + log_prob_w_im + log_prob_y_im + log_prob_r_im;
      // IMPORTANT: In the above line, it crashes when D is inf, as the log(D) is then nan. 
      //            thus it is better to just log_D. However, this will need to change
      //            when creating a linear predictor. 
      
      
      bool test = std::isnan(asDouble(cond_llk_im)); // isnan from cmath requires a double
      if (test)  {
        std::cout << "log_prob_w_im: " << log_prob_w_im <<
          "\nlog_prob_y_im: " << log_prob_y_im <<
            "\nlog_prob_r_im: " << log_prob_r_im << std::endl;
      }
      
      log_probs_i(m) = cond_llk_im;
    }
    
    // sum the values in log_probs_i for every grid point, in a similar way to logsumexp()
    // Can we use the function logspace_add() ?
    Type aa = log_probs_i(0);
    for (int ii = 1; ii < n_grid; ii++) {
      aa = logspace_add(aa, log_probs_i(ii));
    }
    
    llk_cond += aa;
    // std::cout << "llk_cond after " << i << " additions: " << llk_cond << std::endl;
  }
  // Combine the two log likelihood elements and make it negative
  Type nll = -(llk_n + llk_cond);
  
  std::cout << "The total n part of the log likelihood: " << llk_n << std::endl; // nan
  std::cout << "The conditional part of the log likelihood: " << llk_cond << std::endl; // nan
  std::cout << "The negative log likelihood: " << nll << std::endl;

  REPORT(D);
  // ADREPORTs
  // ADREPORT(D);
  ADREPORT(g0);
  ADREPORT(beta_r);
  ADREPORT(pow(sd_r, 2));
  // ADREPORT(sd_r ^ 2);
  // ADREPORT(mu_s);
  
  // ADREPORT(mix_bear);
  // ADREPORT(kappa_low);
  // ADREPORT(kappa_high);
  
  // Return the nll
  return nll;
}