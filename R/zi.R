## Gibbs sampling for mixture of factor analyzers

zi_rbernoulli <- function(pi) sapply(pi, function(p) sample(c(0,1), 1, 
                                                         prob = c(1-p,p)))

#' Turn a matrix's columns into informative names
zi_mcmcify <- function(m, name) {
  colnames(m) <- paste0(name, "[", seq_len(ncol(m)),"]")
  return(m)
}

zi_posterior <- function(y, c0, c1, k0, k1, pst, tau, gamma, theta, eta, tau_k, tau_c, r, alpha, beta,
                      theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_k, beta_k) {
  G <- ncol(y)
  N <- nrow(y)
  
  zero_ll <- sapply(seq_len(N), function(i) {
    sum(dnorm(y[i,], c0 + k0 + pst[i], 1 / sqrt(tau), log = TRUE))
  })
  
  one_ll <- sapply(seq_len(N), function(i) {
    sum(dnorm(y[i,], c1 + k1 + pst[i], 1 / sqrt(tau), log = TRUE))
  })
        
  ll <- sum(zero_ll[gamma == 0]) + sum(one_ll[gamma == 1])

  prior <- sum(dnorm(k0, theta, 1 / sqrt(tau_k), log = TRUE)) +
    sum(dnorm(k1, theta, 1 / sqrt(tau_k), log = TRUE)) +
    sum(dnorm(c0, eta[1], 1 / sqrt(tau_c), log = TRUE)) +
    sum(dnorm(c1, eta[2], 1 / sqrt(tau_c), log = TRUE)) +
    sum(dnorm(theta, theta_tilde, 1 / sqrt(tau_theta), log = TRUE)) +
    sum(dnorm(eta, eta_tilde, 1 / sqrt(tau_eta), log = TRUE)) +
    sum(dgamma(tau, alpha, beta, log = TRUE)) +
    sum(dnorm(pst, 0, 1 / r, log = TRUE)) +
    sum(dgamma(tau_k, alpha_k, beta_k, log = TRUE)) 
  
  
  return( ll + prior )
}


zi_to_ggmcmc <- function(g) {
  g['x_mean_trace'] <- NULL
  x <- do.call(cbind, g)
  mcmcl <- mcmc.list(list(mcmc(x)))
  return(ggs(mcmcl))
}


mfa_zi <- function(y, iter = 2000, thin = 1, burn = iter / 2, lambda = 1,
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
  
  ## REMOVE ---------------
  # pst <- true_t
  
  ## c & k parameters
  lms <- apply(y, 2, function(gex) coef(lm(gex ~ pst)))
  theta <- theta0 <- k0 <- k1 <- lms[2, ]
  c0 <- c1 <- lms[1, ]
  eta <- c(mean(c0), mean(c1))
  
  tau_k <- rep(1, G) # rgamma(G, alpha_k, beta_k)
  tau_c <- 1 # rgamma(G, alpha_c, beta_c)
  
  
  ## assignments for each cell
  gamma <- sample(0:1, N, replace = TRUE) # as.numeric( pst < mean(pst)  ) #
  # gamma <- branch # REMOVE-------
  
  nsamples <- floor((iter - burn) / thin)
  G_dim <- c(nsamples, G)
  N_dim <- c(nsamples, N)
  
  eta_trace <- mcmcify(matrix(NA, nrow = nsamples, ncol = 2), "eta")
  theta_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "theta")
  tau_k_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "tau_k")
  # tau_c_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "tau_c")
  k0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "k0")
  k1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "k1")
  c0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "c0")
  c1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "c1")
  tau_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "tau")
  gamma_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "gamma")
  pst_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "pst")
  lp_trace <- matrix(NA, nrow = nsamples, ncol = 1)
  colnames(lp_trace) <- "lp__"
  

  rownames(y) <- colnames(y) <- NULL
  
  x <- y
  is_dropout <- x == 0
  
  ## Just store the posterior mean of x to avoid huge overheads
  x_mean_trace <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  
  for(it in 1:iter) {
    
    ## set which and N parameters
    which_0 <- which(gamma == 0)
    which_1 <- which(gamma == 1)
    
    which_0_l <- gamma == 0
    which_1_l <- !which_0_l
    
    N_0 <- length(which_0)
    N_1 <- length(which_1)
    
    k0_new <- sample_k(x, pst, c0, tau, theta, tau_k, which_0_l)
    k1_new <- sample_k(x, pst, c1, tau, theta, tau_k, which_1_l)
    
    
    ## update for c0
    c0_new <- sample_c(x, pst, k0_new, tau, eta[1], tau_c, which_0_l, N_0);
    c1_new <- sample_c(x, pst, k1_new, tau, eta[2], tau_c, which_1_l, N_1);
    
    
    ## update for pseudotimes
    pst_new <- sample_pst(x, k0_new, k1_new, c0_new, c1_new, r, gamma, tau);
  
    tau_new <- sample_tau(x, c0_new, c1_new, k0_new, k1_new, gamma, pst_new, alpha, beta)
    
    ## updates for theta (k)
    lambda_theta <- 2 * tau_k + tau_theta
    nu_theta <- tau_theta * theta_tilde + tau_k * (k0_new + k1_new)
    nu_theta <- nu_theta / lambda_theta
    
    theta_new <- rnorm(G, nu_theta, 1 / sqrt(lambda_theta))
    # theta_new <- theta0 # REMOVE ------
    
    ## updates for eta (c)
    lambda_eta <- tau_eta + G * tau_c
    nu_eta <- tau_eta * eta_tilde + tau_c * c(sum(c0_new), sum(c1_new))
    nu_eta <- nu_eta / lambda_eta
    
    eta_new <- rnorm(2, nu_eta, 1 / sqrt(lambda_eta))
    
    ## update for tau_k
    alpha_new <- alpha_k + 1
    beta_new <- beta_k + 0.5 * ( (k0_new - theta_new)^2 + (k1_new - theta_new)^2 )
    tau_k_new <- rgamma(G, alpha_new, beta_new)
    
    pi <-  calculate_pi(x, c0_new, c1_new, k0_new, k1_new, gamma, 
                       pst_new, tau_new, eta_new, tau_c, collapse)
    gamma_new <- r_bernoulli_vec(1 - pi)
    
    ## update for x
    x_new <- sample_x(x, is_dropout, c0_new, c1_new, k0_new, k1_new,
                      gamma, pst_new, tau_new, lambda)
    
    
    ## accept new parameters
    # gamma <- gamma_new
    k0 <- k0_new
    k1 <- k1_new
    c0 <- c0_new
    c1 <- c1_new
    pst <- pst_new
    tau <- tau_new
    eta <- eta_new
    theta <- theta_new
    tau_k <- tau_k_new
    gamma <- gamma_new
    x <- x_new
    
    if((it > burn) && (it %% thin == 0)) {
      sample_pos <- (it - burn) / thin
      k0_trace[sample_pos,] <- k0
      k1_trace[sample_pos,] <- k1
      c0_trace[sample_pos,] <- c0
      c1_trace[sample_pos,] <- c1
      tau_trace[sample_pos,] <- tau
      gamma_trace[sample_pos,] <- gamma
      pst_trace[sample_pos,] <- pst
      theta_trace[sample_pos,] <- theta
      eta_trace[sample_pos,] <- eta
      tau_k_trace[sample_pos,] <- tau_k
      # tau_c_trace[sample_pos,] <- tau_c
      
      post <- posterior(y, c0, c1, k0, k1, pst,
                        tau, gamma, theta, eta, tau_k, tau_c, r,
                        alpha, beta, theta_tilde, 
                        eta_tilde, tau_theta, tau_eta,
                        alpha_k, beta_k)
      lp_trace[sample_pos,] <- post 
      
      x_mean_trace <- x_mean_trace + x
    }
  }
  return(list(k0_trace = k0_trace, k1_trace = k1_trace,
              c0_trace = c0_trace, c1_trace = c1_trace,
              tau_trace = tau_trace, gamma_trace = gamma_trace,
              pst_trace = pst_trace, theta_trace = theta_trace,
              eta_trace = eta_trace, lp_trace = lp_trace,
              tau_k_trace = tau_k_trace, 
              x_mean_trace = x_mean_trace / nsamples))
}

zi_log_sum_exp <- function(x) log(sum(exp(x - max(x)))) + max(x)


empirical_lambda <- function(Y) {
  means <- colMeans(Y)
  pdrop <- colMeans(Y == 0)
  fit <- nls(pdrop ~ exp(-lambda * means), start = list(lambda = 1))
  coef(fit)['lambda']
}

plot_dropout_relationship <- function(Y, lambda = empirical_lambda(Y)) {
  ff <- function(x, lambda) exp(- lambda * x)
  d <- data_frame(mean = colMeans(Y), pdrop = colMeans(Y == 0)) %>% 
    mutate(fit = ff(mean, lambda))
  
  ggplot(arrange(d, desc(fit)), aes(x = mean)) + 
    geom_point(aes(y = pdrop), size = 1, alpha = 0.5) +
    geom_line(aes(y = fit), color = 'red') +
    xlab("Mean expression") + ylab("Proportion dropout")
  
}
