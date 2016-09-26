## Gibbs sampling for mixture of factor analyzers
## kieran.campbell@sjc.ox.ac.uk

rbernoulli <- function(pi) sapply(pi, function(p) sample(c(0,1), 1, 
                                                         prob = c(1-p,p)))

#' Turn a matrix's columns into informative names
mcmcify <- function(m, name) {
  colnames(m) <- paste0(name, "[", seq_len(ncol(m)),"]")
  return(m)
}


log_sum_exp <- function(x) log(sum(exp(x - max(x)))) + max(x)

map_branch <- function(g) {
  gamma <- g$gamma_trace
  df <- apply(gamma, 2, function(gam) {
    tab <- table(gam)
    max <- which.max(tab)
    max_n <- as.integer(names(max))
    prop <- mean(gam == max_n)
    return(c(max_n, prop))
  })
  df <- t(df) %>% as_data_frame()
  names(df) <- c("max", "prop")
  return(df)
}


posterior <- function(y, c, k, pst, tau, gamma, theta, eta, chi, tau_c, r, alpha, beta,
                      theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_k, beta_k) {
  G <- ncol(y)
  N <- nrow(y)
  b <- ncol(k)
  
  branch_likelihoods <- sapply(seq_len(b), function(branch) {
    ll_branch <- sapply(seq_len(N), function(i) {
        sum(dnorm(y[i,], c[,branch] + k[branch] * pst[i], 1 / sqrt(tau), log = TRUE))
    })
    sum(ll_branch[gamma == branch])
  })
    
  ll <- sum(branch_likelihoods)
  
  prior <- 
    sum(dnorm(theta, theta_tilde, 1 / sqrt(tau_theta), log = TRUE)) +
    sum(dnorm(eta, eta_tilde, 1 / sqrt(tau_eta), log = TRUE)) +
    sum(dgamma(tau, alpha, beta, log = TRUE)) +
    sum(dnorm(pst, 0, 1 / r, log = TRUE)) +
    sum(dgamma(chi, alpha_k, beta_k, log = TRUE)) 
  
  k_prior <- sum( apply(k, 2, function(k_b) sum(dnorm(k_b, theta, 1 / sqrt(chi), log = TRUE))) )
  c_prior <- sum( sapply(seq_len(b), function(branch) sum(dnorm(c[,branch], eta[branch], 1 / sqrt(tau_c), log = TRUE)))) 
  
  prior <- prior + k_prior + c_prior
  
  return( ll + prior )
}



to_ggmcmc <- function(g) {
  x <- do.call(cbind, g)
  mcmcl <- mcmc.list(list(mcmc(x)))
  return(ggs(mcmcl))
}


#' Fit a MFA object
#' 
#' @export
#' @return Something horrific
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib WMUtils
mfa <- function(y, iter = 2000, thin = 1, burn = iter / 2, b = 2,
                pc_initialise = 1, collapse = FALSE, seed = 123L,
                eta_tilde = mean(y)) {
  
  # set.seed(seed)
  N <- nrow(y)
  G <- ncol(y)
  message(paste("Sampling for", N, "cells and", G, "genes"))
  
  ## branching hierarchy
  theta_tilde <- 0
  # eta_tilde <- 5
  tau_eta <- tau_theta <- 1e-2
  alpha_k <- beta_k <- 1e-2
  
  
  ## precision parameters
  alpha <- 2
  beta <- 1
  tau <- rep(1, G)
  
  ## pseudotime parameters
  r <- 1
  pst <-  prcomp(y)$x[,pc_initialise] # rep(0.5, N) # rnorm(N, 0, 1 / r^2)
  pst <- pst / sd(pst) # make it a little more consistent with prior assumptions
  
  
  ## c & k parameters
  lms <- apply(y, 2, function(gex) coef(lm(gex ~ pst)))
  theta <- theta0 <- lms[2, ]
  k <- matrix(NA, nrow = G, ncol = b)
  for(i in seq_len(b)) k[,i] <- lms[2,]
  
  c <- matrix(NA, nrow = G, ncol = b)
  for(i in seq_len(b)) c[,i] <- lms[1,]
  
  eta <- colMeans(c)
  
  chi <- rep(1, G) # rgamma(G, alpha_k, beta_k)
  tau_c <- 1 # 0.1 # rgamma(G, alpha_c, beta_c)
  
  
  ## assignments for each cell
  gamma <- sample(seq_len(b), N, replace = TRUE) # as.numeric( pst < mean(pst)  ) #

  nsamples <- floor((iter - burn) / thin)
  G_dim <- c(nsamples, G)
  N_dim <- c(nsamples, N)
  
  eta_trace <- mcmcify(matrix(NA, nrow = nsamples, ncol = b), "eta")
  theta_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "theta")
  chi_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "chi")
  
  k_trace <- array(dim = c(nsamples, b, G))
  c_trace <- array(dim = c(nsamples, b, G))
  
  tau_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "tau")
  gamma_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "gamma")
  pst_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "pst")
  lp_trace <- matrix(NA, nrow = nsamples, ncol = 1)
  colnames(lp_trace) <- "lp__"
  
  rownames(y) <- colnames(y) <- NULL
  
  for(it in 1:iter) {
    

    
    k_new <- sapply(seq_len(b), function(branch) sample_k(y, pst, c[, branch], tau, theta, chi, gamma == branch))
    
    c_new <- sapply(seq_len(b), function(branch) sample_c(y, pst, k_new[, branch], tau, eta[branch], tau_c, gamma == branch, sum(gamma == branch)))
    
    
    # ## remove
    # calculate_nuc(y, pst, k_new[,branch], tau, eta, tau_c, which())
    # 
    ####

    ## update for pseudotimes
    pst_new <- sample_pst(y, c_new, k_new, r, gamma, tau);
  
    tau_new <- sample_tau(y, c_new, k_new, gamma, pst_new, alpha, beta)
    
    ## updates for theta (k)
    lambda_theta <- 2 * chi + tau_theta
    nu_theta <- tau_theta * theta_tilde + chi * rowSums(k_new)
    nu_theta <- nu_theta / lambda_theta
    
    theta_new <- rnorm(G, nu_theta, 1 / sqrt(lambda_theta))

    ## updates for eta (c)
    lambda_eta <- tau_eta + G * tau_c
    nu_eta <- tau_eta * eta_tilde + tau_c * colSums(c_new)
    nu_eta <- nu_eta / lambda_eta
    
    eta_new <- rnorm(length(nu_eta), nu_eta, 1 / sqrt(lambda_eta))
    
    ## update for tau_k
    alpha_new <- alpha_k + 1
    beta_new <- beta_k + 0.5 * rowSums( (k_new - theta_new)^2 )
    chi_new <- rgamma(G, alpha_new, beta_new)

    pi <-  calculate_pi(y, c_new, k_new, pst_new, tau_new, eta_new, tau_c, collapse)
    gamma <- r_bernoulli_mat(pi) + 1 # need +1 to convert from C++ to R
    # print(table(gamma))

    k <- k_new
    c <- c_new
    pst <- pst_new
    tau <- tau_new
    eta <- eta_new
    theta <- theta_new
    chi <- chi_new
    # tau_c <- tau_c_new
    
    if((it > burn) && (it %% thin == 0)) {
      sample_pos <- (it - burn) / thin
      # k0_trace[sample_pos,] <- k0
      # k1_trace[sample_pos,] <- k1
      # c0_trace[sample_pos,] <- c0
      # c1_trace[sample_pos,] <- c1
      tau_trace[sample_pos,] <- tau
      gamma_trace[sample_pos,] <- gamma
      pst_trace[sample_pos,] <- pst
      theta_trace[sample_pos,] <- theta
      
      eta_trace[sample_pos,] <- eta
      # tau_k_trace[sample_pos,] <- tau_k
      # tau_c_trace[sample_pos,] <- tau_c
      
      post <- posterior(y, c, k, pst,
                        tau, gamma, theta, eta, chi, tau_c, r,
                        alpha, beta, theta_tilde, 
                        eta_tilde, tau_theta, tau_eta,
                        alpha_k, beta_k)
      lp_trace[sample_pos,] <- post 
    }
  }
  return(list(#k0_trace = k0_trace, k1_trace = k1_trace,
              #c0_trace = c0_trace, c1_trace = c1_trace,
              tau_trace = tau_trace, gamma_trace = gamma_trace,
              pst_trace = pst_trace, theta_trace = theta_trace,
              eta_trace = eta_trace, lp_trace = lp_trace))
}
