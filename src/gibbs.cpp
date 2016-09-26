#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;


/***
 * Convenience functions to compute log normal density and log sum of exponentials
 */

// [[Rcpp::export]]
double log_d_norm(double x, double mu, double sigma) {
  double ldn = -0.5 * log(2 * PI);
  ldn -= log(sigma) + 1 / (2 * sigma * sigma) * (x - mu) * (x - mu);
  return ldn;
}

// [[Rcpp::export]]
double log_sum_exp(NumericVector x) {
  double y;
  y = log(sum(exp(x - max(x)))) + max(x);
  return(y);
}

// [[Rcpp::export]]
IntegerVector r_bernoulli_mat(NumericMatrix pi) {
  int N = pi.nrow();
  int S = pi.ncol(); // n objects
  
  IntegerVector gamma(N); // store results
  
  NumericVector rands = runif(N); // n random numbers
  
  for(int i = 0; i < N; i++) {
    NumericVector pi_i = pi(i,_);

    double cum_prob = 0;    
    for(int s = 0; s < S; s++) {
      cum_prob += pi_i[s];
      if(rands[i] < cum_prob) {
        gamma[i] = s;
        break;
      } 
    }
  }
  return gamma;
}

/*******
 * K sampling here
 ******/

//[[Rcpp::export]]
NumericVector calculate_nuk(NumericMatrix y, NumericVector pst, NumericVector c,
                       NumericVector tau, NumericVector theta, NumericVector tau_k,
                       LogicalVector which_l) {
  int N = y.nrow();
  int G = y.ncol();
  
  NumericVector nu_k(G);

  for(int g = 0; g < G; g++) {
    double tmp = 0;
    for(int i = 0; i < N; i++) {
      if(which_l[i])
        tmp += pst[i] * ( y(i,g) - c[g] );
    }
    nu_k[g] = tau_k[g] * theta[g] + tau[g] * tmp;
  }  
  
  return nu_k;
}

//[[Rcpp::export]]
NumericVector calculate_lamk(NumericVector tau_k, NumericVector tau, 
                             NumericVector pst, LogicalVector which_l) {
  int G = tau.size();
  int N = pst.size();
  
  NumericVector lamk(G);

  double pst_sum = 0;  
  for(int i = 0; i < N; i++) {
    if(which_l[i])
      pst_sum += pst[i] * pst[i];
  }
  
  lamk = tau * pst_sum + tau_k;

  return lamk;
}

// [[Rcpp::export]]
NumericVector sample_k(NumericMatrix y, NumericVector pst, NumericVector c,
                       NumericVector tau, NumericVector theta, NumericVector tau_k,
                       LogicalVector which_l) {
  int G = y.ncol();
  
  NumericVector nuk = calculate_nuk(y, pst, c, tau, theta, tau_k, which_l);
  NumericVector lamk = calculate_lamk(tau_k, tau, pst, which_l);
  
  for(int g = 0; g < G; g++) 
    nuk[g] /= lamk[g];
  
  NumericVector k_new(G);
  
  for(int g = 0; g < G; g++)
    k_new[g] = as<double>(rnorm(1, nuk[g], 1 / sqrt(lamk[g])));
  
  return k_new;
}


/*******
 * C sampling here
 ******/

//[[Rcpp::export]]
NumericVector calculate_nuc(NumericMatrix y, NumericVector pst, NumericVector k,
                            NumericVector tau, double eta, double tau_c,
                            LogicalVector which_l) {
  int G = k.size();
  int N = pst.size();
  
  NumericVector nuc(G);
  
  for(int g = 0; g < G; g++) {
    double inner_sum = 0;
    for(int i = 0; i < N; i++) {
      if(which_l[i])
        inner_sum += y(i,g) - k[g] * pst[i];
    }
    nuc[g] = tau_c * eta + tau[g] * inner_sum;
  }
  
  
  return nuc;
}

// [[Rcpp::export]]
NumericVector calculate_lamc(NumericVector tau, double tau_c, int N) {
  int G = tau.size();
  NumericVector lamc(G);
  for(int g = 0 ; g < G; g++)
    lamc[g] =  tau_c  + N * tau[g];
  return lamc;
}

// [[Rcpp::export]]
NumericVector sample_c(NumericMatrix y, NumericVector pst, NumericVector k,
                       NumericVector tau, double eta, double tau_c,
                       LogicalVector which_l, int N) {
  int G = k.size();
  NumericVector nuc = calculate_nuc(y, pst, k, tau, eta, tau_c, which_l);
  NumericVector lamc = calculate_lamc(tau, tau_c, N);
  
  for(int g = 0; g < G; g++) 
    nuc[g] /= lamc[g];
  
  NumericVector c_new(G);
  
  for(int g = 0; g < G; g++)
    c_new[g] = as<double>(rnorm(1, nuc[g], 1 / sqrt(lamc[g])));
  
  return c_new;
}

/** 
 * Pseudotime updates
 * */

/* 
 *First column is lambda, second is mean. We explicitly calculate this to check 
 *consistency with previous results.
 */

NumericMatrix pst_update_par(NumericMatrix y, NumericMatrix c, NumericMatrix k, 
                          double r, NumericVector gamma, NumericVector tau) {
  int N = y.nrow();
  int G = y.ncol();

  NumericMatrix pst_parameters(N, 2); // what we return
  
  // Placehold vectors for current branch
  NumericVector k_(G);
  NumericVector c_(G);
  
  double lam_ti, nu_ti;
  
  for(int i = 0; i < N; i++) {
    nu_ti = 0;
    
    k_ = k(_, gamma[i] - 1); // C++ as 0 indexing
    c_ = c(_, gamma[i] - 1);

    lam_ti = pow(r, 2) + sum(tau * pow(k_, 2));
    for(int g = 0; g < G; g++) 
      nu_ti += tau[g] * k_[g] * (y(i,g) - c_[g]);
    
    nu_ti /= lam_ti;
    pst_parameters(i, 0) = nu_ti;
    pst_parameters(i, 1) = lam_ti;
  }
  
  return pst_parameters;
}

// [[Rcpp::export]]
NumericVector sample_pst(NumericMatrix y, NumericMatrix c, NumericMatrix k, 
                         double r, NumericVector gamma, NumericVector tau) {
  int N = y.nrow();

  NumericMatrix pst_pars(N, 2);
  
  pst_pars = pst_update_par(y, c, k, r, gamma, tau);
  
  NumericVector pst_new(N);
  for(int i = 0; i < N; i++) {
    pst_new[i] = as<double>(rnorm(1, pst_pars(i, 0), 1 / sqrt(pst_pars(i, 1))));
  }
  
  return pst_new;
}


NumericMatrix tau_params(NumericMatrix y, NumericMatrix c, NumericMatrix k,
                         NumericVector gamma, NumericVector pst, double alpha, double beta) {
  int N = y.nrow();
  int G = y.ncol();
  
  NumericMatrix alpha_beta(G, 2); // new alpha and beta: first column alpha, second beta
  NumericMatrix mu(N, G);
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      mu(i,g) = c(g, gamma[i] - 1) + k(g, gamma[i] - 1) * pst[i];
    }
  }
  
  for(int g = 0; g < G; g++) {
    alpha_beta(g, 0) = alpha + N / 2;
    double beta_new = beta;
    for(int i = 0; i < N; i++)
      beta_new += pow(y(i,g) - mu(i,g), 2) / 2;
    
    alpha_beta(g, 1) = beta_new;
  }
  
  return alpha_beta;
  
  
}

// [[Rcpp::export]]
NumericVector sample_tau(NumericMatrix y, NumericMatrix c, NumericMatrix k,
                         NumericVector gamma, NumericVector pst, double alpha, double beta) {
  // N = y.nrow();
  int G = y.ncol();
  
  NumericVector tau(G);
  NumericMatrix alpha_beta(G, 2);
  alpha_beta = tau_params(y, c, k, gamma, pst, alpha, beta);
  
  for(int g = 0; g < G; g++) {
    tau[g] = as<double>(rgamma(1, alpha_beta(g, 0), 1 / alpha_beta(g, 1))); // !!! RCPP gamma parametrised by shape - scale
  }
  
  return tau;
}



// [[Rcpp::export]]
NumericMatrix calculate_pi(NumericMatrix y, NumericMatrix c, NumericMatrix k, NumericVector pst, NumericVector tau, 
                           NumericVector eta, double tau_c, bool collapse, NumericVector log_w) {
  int N = y.nrow();
  int G = y.ncol();
  int b = c.ncol(); // number of branches

  NumericMatrix pi(N,b); // probability of class membership for cell (row) and branch (column)


  if(collapse == 0) {
    for(int i = 0; i < N; i++) {
      NumericVector comp_ll(b, 0.0); // log-likelihood vector for each branch, default to 0.0

      for(int branch = 0; branch < b; branch++) {

        for(int g = 0; g < G; g++) {
          double y_ = y(i,g);
          double sd = 1 / sqrt(tau[g]);
        
          double comp_mean = c(g, branch) + k(g, branch) * pst[i];
          // std::cout << comp_mean << " ";
          comp_ll[branch] += log_d_norm(y_, comp_mean, sd);
        }
        // std::cout << std::endl;
        comp_ll[branch] += log_w[branch];
      }

      for(int branch = 0; branch < b; branch++) {
        // std::cout << comp_ll[branch] << " ";
        pi(i, branch) = exp(comp_ll[branch] - log_sum_exp(comp_ll));
      }
      // std::cout << std::endl;
    }
  } else {
    for(int i = 0; i < N; i++) {
      NumericVector comp_ll(b, 0.0); // log-likelihood vector for each branch, default to 0.0
      
      for(int branch = 0; branch < b; branch++) {
        for(int g = 0; g < G; g++) {
          double y_ = y(i,g);
          double sd = 1/sqrt( (1/tau_c) + (1/tau[g]) );
        
          double comp_mean = eta[branch] + k(g, branch) * pst[i];
          comp_ll[branch] += log_d_norm(y_, comp_mean, sd);
        }
        comp_ll[branch] += log_w[branch]; 
      }
      for(int branch = 0; branch < b; branch++)
        pi(i, branch) = exp( comp_ll[branch] - log_sum_exp(comp_ll));
    }
  }
  return pi;
}


/***
 * x sampling
 ***/

// [[Rcpp::export]]
NumericMatrix sample_x(NumericMatrix x, LogicalMatrix is_dropout,
                       NumericVector c0, NumericVector c1, NumericVector k0, NumericVector k1,
                           NumericVector gamma, NumericVector pst, NumericVector tau, double lambda) {
  int N = x.nrow();
  int G = x.ncol();
  
  NumericVector k(G);
  NumericVector c(G);
  
  NumericMatrix x_new(N, G);
  

  for(int i = 0; i < N; i++) {
    if(gamma[i] == 0) {
      k = k0;
      c = c0;
    } else {
      k = k1;
      c = c1;
    }
    for(int g = 0; g < G; g++) {
      if(is_dropout(i,g) == true) {
        double mu_ig = c[g] + k[g] * pst[i]; 
        x_new(i,g) = as<double>(rnorm(1, mu_ig + lambda / (N * tau[g]), 1 / sqrt(tau[g])));
      } else {
        x_new(i,g) = x(i,g);
      }
    }
  }
  
  
  return x_new;
}
