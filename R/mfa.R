## Gibbs sampling for mixture of factor analyzers
## kieran.campbell@sjc.ox.ac.uk

rbernoulli <- function(pi) sapply(pi, function(p) sample(c(0,1), 1, 
                                                         prob = c(1-p,p)))

#' Turn a matrix's columns into informative names
#' 
#' @param m A fit returned from \code{mfa}
#' @param name The name of the parameter
mcmcify <- function(m, name) {
  colnames(m) <- paste0(name, "[", seq_len(ncol(m)),"]")
  return(m)
}


log_sum_exp <- function(x) log(sum(exp(x - max(x)))) + max(x)


#' Calculate the log-posterior during inference
#' 
#' @importFrom stats dgamma dnorm
posterior <- function(y, c, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
                      theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi,
                      zero_inflation = FALSE, lambda = NULL) {
  G <- ncol(y)
  N <- nrow(y)
  b <- ncol(k)
  
  ll_cell <- sapply(seq_len(N), function(i) {
    sum(dnorm(y[i,], c[,gamma[i]] + k[,gamma[i]] * pst[i], 1 / sqrt(tau), log = TRUE))
  })
  
  ll <- sum(ll_cell)
  
  prior <- 
    sum(dnorm(theta, theta_tilde, 1 / sqrt(tau_theta), log = TRUE)) +
    sum(dnorm(eta, eta_tilde, 1 / sqrt(tau_eta), log = TRUE)) +
    sum(dgamma(tau, alpha, beta, log = TRUE)) +
    sum(dnorm(pst, 0, 1 / r, log = TRUE)) +
    sum(dgamma(chi, alpha_chi, beta_chi, log = TRUE)) +
    log(MCMCpack::ddirichlet(w, rep(1/b, b)))
  
  k_prior <- sum( apply(k, 2, function(k_b) sum(dnorm(k_b, theta, 1 / sqrt(chi), log = TRUE))) )
  c_prior <- sum( sapply(seq_len(b), function(branch) sum(dnorm(c[,branch], eta, 1 / sqrt(tau_c), log = TRUE)))) 
  
  prior <- prior + k_prior + c_prior
  
  if(zero_inflation) {
    means <- colMeans(y)
    p_drop <- exp(-lambda * means)
  }
  

  return( ll + prior )
}


#' Turn a trace list to a \code{ggmcmc} object
#' 
#' @param g A list of trace matrices
#' 
to_ggmcmc <- function(g) {
  x <- do.call(cbind, g)
  mcmcl <- coda::mcmc.list(list(coda::mcmc(x)))
  return(ggmcmc::ggs(mcmcl))
}


#' Fit a MFA object
#' 
#' Perform Gibbs sampling inference for a hierarchical Bayesian mixture of factor analysers
#' to identify bifurcations in single-cell expression data.
#' 
#' @param y A cell-by-gene single-cell expression matrix
#' @param iter Number of MCMC iterations
#' @param thin MCMC samples to thin
#' @param burn Number of MCMC samples to throw away
#' @param b Number of branches to model
#' @param zero_inflated Logical, should zero inflation be enabled?
#' @param lambda The dropout parameter - by default estimated using the function \code{empirical_lambda}
#' @param pc_initialise Which principal component to initialise pseudotimes to
#' @param collapse Collapsed Gibbs sampling of branch assignments
#' @param seed Random seed to set
#' @param scale_input Logical. If true, input is scaled to have mean 0 variance 1
#' @param eta_tilde Hyperparameter
#' @param alpha Hyperparameter
#' @param beta Hyperparameter
#' @param theta_tilde Hyperparameter
#' @param tau_eta Hyperparameter
#' @param tau_theta Hyperparameter
#' @param alpha_chi Hyperparameter
#' @param beta_chi Hyperparameter
#' @param w_alpha Hyperparameter
#' @param clamp_pseudotimes This clamps the pseudotimes to their initial values and
#' doesn't perform sampling. Should be FALSE except for diagnostics.
#' 
#' @details 
#' The column names of Y are used as feature (gene/transcript) names while the row names
#' are used as cell names. If either of these is undefined then the corresponding names
#' are set to cell_x or feature_y.
#' 
#' It is recommended the form of Y is analogous to log-expression to mitigate the impact of 
#' outliers.
#' 
#' In the absence of prior information, three valid local maxima in the posterior likelihood
#' exist (see manuscript). Setting the initial values to a principal component typically
#' fixes sampling to one of them, analogous to specifying a root cell in similar methods.
#' 
#' The hyper-parameter \code{eta_tilde} represents the expected expression in the absence of
#' any actual expression measurements. While a Bayesian purist might reason this based on 
#' knowledge of the measurement technology, simply taking the mean of the input matrix in
#' an Empirical Bayes style seems reasonable.
#' 
#' The degree of shrinkage of the factor loading matrices to a common value is given by the
#' gamma prior on \code{chi}. The mean of this is \code{alpha_chi / beta_chi} while the variance 
#' \code{alpha_chi / beta_chi^2}. Therefore, to obtain higher levels of shrinkage increase
#' \code{alpha_chi} with respect to \code{beta_chi}.
#' 
#' The collapsed Gibbs sampling option given by \code{collapse} involves marginalising out
#' \code{c} (the factor loading intercepts) when updating the branch assignment parameters
#' \code{gamma} which tends to soften the branch assignments.
#' 
#' If zero inflation is enabled using the \code{zero_inflation} parameter then scaling should 
#' *not* be enabled.
#' 
#' @export
#' @return 
#' An S3 structure with the following entries:
#' \itemize{
#' \item \code{traces} A list of iteration-by-dim trace matrices for several important variables
#' \item \code{iter} Number of iterations
#' \item \code{thin} Thinning applied
#' \item \code{burn} Burn period at the start of MCMC
#' \item \code{b} Number of branches modelled
#' \item \code{prop_collapse} Proportion of updates for gamma that are collapsed
#' \item \code{N} Number of cells
#' \item \code{G} Number of features (genes/transcripts)
#' \item \code{feature_names} Names of features
#' \item \code{cell_names} Names of cells
#' }
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom stats prcomp nls coef sd lm rnorm rgamma
#' @useDynLib mfa
mfa <- function(y, iter = 2000, thin = 1, burn = iter / 2, b = 2,
                zero_inflation = FALSE,
                pc_initialise = 1, prop_collapse = 0,
                scale_input = !zero_inflation,
                lambda = empirical_lambda(y),
                eta_tilde = mean(y), alpha = 1e-1, beta = 1e-1,
                theta_tilde = 0, tau_eta = 1, tau_theta = 1, tau_c = 1,
                alpha_chi = 1e-2, beta_chi = 1e-2, w_alpha = 1 / b,
                clamp_pseudotimes = FALSE) {
  
  N <- nrow(y)
  G <- ncol(y)
  message(paste("Sampling for", N, "cells and", G, "genes"))
  
  feature_names <- colnames(y)
  cell_names <- rownames(y)
  
  if(is.null(feature_names)) feature_names <- paste0("feature_", seq_len(G))
  if(is.null(cell_names)) cell_names <- paste0("cell_", seq_len(N))
  
  if(scale_input) y <- scale(y)
  
  ## precision parameters
  tau <- rep(1, G)
  
  ## pseudotime parameters
  r <- 1
  pst <-  prcomp(y, scale = TRUE)$x[,pc_initialise] # rep(0.5, N) # rnorm(N, 0, 1 / r^2)
  pst <- pst / sd(pst) # make it a little more consistent with prior assumptions
  
  
  ## c & k parameters
  lms <- apply(y, 2, function(gex) coef(lm(gex ~ pst)))
  theta <- theta0 <- lms[2, ]
  k <- matrix(NA, nrow = G, ncol = b)
  for(i in seq_len(b)) k[,i] <- lms[2,]
  
  c <- matrix(NA, nrow = G, ncol = b)
  for(i in seq_len(b)) c[,i] <- lms[1,]
  
  eta <- mean(c)
  
  chi <- rep(1, G) # rgamma(G, alpha_chi, beta_chi)

  ## assignments for each cell
  w <- rep(1/b, b) # prior probability of each branch
  gamma <- sample(seq_len(b), N, replace = TRUE, prob=w) # as.numeric( pst < mean(pst)  ) 

  nsamples <- floor((iter - burn) / thin)
  G_dim <- c(nsamples, G)
  N_dim <- c(nsamples, N)
  
  eta_trace <- mcmcify(matrix(NA, nrow = nsamples, ncol = 1), "eta")
  theta_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "theta")
  lambda_theta_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "lambda_theta")
  
  chi_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "chi")
  
  k_trace <- array(dim = c(nsamples, G, b))
  c_trace <- array(dim = c(nsamples, G, b))
  
  tau_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "tau")
  
  gamma_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "gamma")
  pst_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "pst")
  lp_trace <- matrix(NA, nrow = nsamples, ncol = 1)
  colnames(lp_trace) <- "lp__"
  
  rownames(y) <- colnames(y) <- NULL
  
  ## Sort out zero-inflation
  is_dropout <- x_mean_trace <- NULL
  x <- y
  
  if(zero_inflation) {
    is_dropout <- x == 0
    
    ## Just store the posterior mean of x to avoid huge overheads
    x_mean_trace <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  }
  
  # We need to be somewhat careful here:
  # If zero_inflated == FALSE then x == y always and we delete y
  # If zero_inflated == TRUE then we gibbs sample x and keep y
  if(!zero_inflation) rm(y)
  
  for(it in 1:iter) {
    
    # Factor loading matrix sampling
    k_new <- sapply(seq_len(b), function(branch) sample_k(x, pst, c[, branch], tau, theta, chi, gamma == branch))
    c_new <- sapply(seq_len(b), function(branch) sample_c(x, pst, k_new[, branch], tau, eta, tau_c, gamma == branch, sum(gamma == branch)))
    
    # Pseudotime sampling
    if(!clamp_pseudotimes) {
      pst_new <- sample_pst(x, c_new, k_new, r, gamma, tau);
    } else {
      pst_new <- pst
    }
    
    # Precision sampling
    tau_new <- sample_tau(x, c_new, k_new, gamma, pst_new, alpha, beta)
    
    # Theta sampling
    lambda_theta <- 2 * chi + tau_theta
    nu_theta <- tau_theta * theta_tilde + chi * rowSums(k_new)
    nu_theta <- nu_theta / lambda_theta
    theta_new <- rnorm(G, nu_theta, 1 / sqrt(lambda_theta))

    # Eta sampling
    lambda_eta <- tau_eta + G * b * tau_c
    nu_eta <- tau_eta * eta_tilde + tau_c * sum(c_new)
    nu_eta <- nu_eta / lambda_eta
    eta_new <- rnorm(1, nu_eta, 1 / sqrt(lambda_eta))
    
    # Chi sampling
    alpha_new <- alpha_chi + b / 2
    beta_new <- beta_chi + 0.5 * rowSums( (k_new - theta_new)^2 )
    chi_new <- rgamma(G, alpha_new, beta_new)

    # Gamma sampling
    collapse <- runif(1) < prop_collapse
    pi <-  calculate_pi(x, c_new, k_new, pst_new, tau_new, eta_new, tau_c, collapse, log(w),
                        log_result = FALSE)
    gamma <- r_bernoulli_mat(pi) + 1 # need +1 to convert from C++ to R
    
    # update prior probabilities of each branch
    n_gamma <- tabulate(gamma, nbins = b)
    w <- MCMCpack::rdirichlet(1, n_gamma + w_alpha) 
    
    # Gibbs sample x
    if(zero_inflation) {
      x <- sample_x(y, is_dropout, c, k, gamma, pst, tau, lambda)
    }
    
    # Gibbs sampling - accept all parameters
    k <- k_new
    c <- c_new
    pst <- pst_new
    tau <- tau_new
    eta <- eta_new
    theta <- theta_new
    chi <- chi_new

    # Add some relevant variables to trace    
    if((it > burn) && (it %% thin == 0)) {
      sample_pos <- (it - burn) / thin
      tau_trace[sample_pos,] <- tau
      gamma_trace[sample_pos,] <- gamma
      pst_trace[sample_pos,] <- pst
      theta_trace[sample_pos,] <- theta
      chi_trace[sample_pos,] <- chi
      lambda_theta_trace[sample_pos,] <- lambda_theta
      eta_trace[sample_pos,] <- eta
      k_trace[sample_pos,,] <- k
      c_trace[sample_pos,,] <- c

      post <- posterior(x, c, k, pst,
                        tau, gamma, theta, eta, chi, w, tau_c, r,
                        alpha, beta, theta_tilde, 
                        eta_tilde, tau_theta, tau_eta,
                        alpha_chi, beta_chi)
      lp_trace[sample_pos,] <- post 
      
      if(zero_inflation) x_mean_trace <- x_mean_trace + x
    }
  }
  traces <- list(tau_trace = tau_trace, gamma_trace = gamma_trace,
                      pst_trace = pst_trace, theta_trace = theta_trace, lambda_theta_trace = lambda_theta_trace, chi_trace = chi_trace,
                      eta_trace = eta_trace, k_trace = k_trace, c_trace = c_trace,
                 lp_trace = lp_trace)
  if(zero_inflation) traces$x_mean <- x_mean_trace / nsamples
  
  mfa_res <- structure(list(traces = traces, iter = iter, thin = thin, burn = burn,
                            b = b, collapse = collapse, N = N, G = G,
                            feature_names = feature_names, cell_names = cell_names), 
                       class = "mfa")
}

#' Print an mfa fit
#' 
#' @param x An MFA fit returned by \code{mfa}
#' @param ... Additional arguments
#' 
#' @export
print.mfa <- function(x, ...) {
  msg <- paste("MFA fit with\n",
               x$N, "cells and", x$G, "genes\n",
               "(", x$iter, "iterations )")
  cat(msg)
}

#' Plot MFA trace
#' 
#' @param m A fit returned from \code{mfa}
#' 
#' @export
#' 
#' @importFrom methods is
plot_mfa_trace <- function(m) {
  stopifnot(is(m, "mfa"))
  lp <- m$traces$lp_trace[,1]
  qplot(seq_along(lp), lp, geom = 'line') + stat_smooth(se = F) +
    xlab("Iteration") + ylab("log-probability")
}


#' Plot MFA autocorrelation
#' 
#' @param m A fit returned from \code{mfa}
#' 
#' @export
#' @importFrom dplyr data_frame
#' @importFrom methods is
plot_mfa_autocorr <- function(m) {
  stopifnot(is(m, "mfa"))
  lp <- m$traces$lp_trace[,1]
  lp_df <- ggmcmc::ggs(coda::mcmc.list(list(coda::mcmc(data.frame(lp))))) # I long for a simple life
  
  ggmcmc::ggs_autocorrelation(lp_df)
}


map_branch <- function(g) {
  gamma <- g$gamma_trace
  df <- apply(gamma, 2, function(gam) {
    tab <- table(gam)
    max <- which.max(tab)
    max_n <- as.integer(names(max))
    prop <- mean(gam == max_n)
    return(c(max_n, prop))
  })
  df <- dplyr::as_data_frame(t(df))
  names(df) <- c("max", "prop")
  df$max <- factor(df$max)
  return(df)
}


#' Summarise mfa fit
#' 
#' @param object An MFA fit returned by a call to \code{mfa}
#' @param ... Additional arguments
#' @export
#' 
#' @importFrom MCMCglmm posterior.mode
summary.mfa <- function(object, ...) {
  
  ## Please someone fix R:
  branch <- branch_certainty <- pseudotime_lower <- pseudotime_upper <- NULL
  
  # map branching
  df <- map_branch(object$traces)
  
  # pseudotimes
  tmap <- posterior.mode(coda::mcmc(object$traces$pst_trace))
  hpd_credint <- coda::HPDinterval(coda::mcmc(object$traces$pst_trace))
  
  df$pseudotime <- tmap
  df$pseudotime_lower <- hpd_credint[,1]
  df$pseudotime_upper <- hpd_credint[,2]
  
  df <- dplyr::rename(df, branch = max, branch_certainty = prop)
  df <- dplyr::select(df, pseudotime, branch, branch_certainty, pseudotime_lower, pseudotime_upper)
  return( df )
}

#' Calculate posterior precision parameters
#' 
#' @param m A fit returned from \code{mfa}
#' 
#' @export
#' @importFrom MCMCglmm posterior.mode
calculate_chi <- function(m) {
  chi_map <- posterior.mode(coda::mcmc(m$traces$chi_trace))
  return(
    dplyr::data_frame(feature = m$feature_names, chi_map)
  )
}

#' Plot posterior precision parameters
#' 
#' @param m A fit returned from \code{mfa}
#' @param nfeatures Top number of 
#' 
#' @export
#' @import ggplot2
plot_chi <- function(m, nfeatures = m$G) {
  chi <- calculate_chi(m)
  chi <- dplyr::mutate(chi, chi_map_inverse = 1 / chi_map)
  chi <- dplyr::arrange(chi, chi_map_inverse)
  chi <- chi[seq_len(nfeatures), ]
  
  chi$feature <- factor(chi$feature, levels = chi$feature)
  
  ggplot(chi, aes_string(x = "feature", y = "chi_map_inverse")) +
    geom_bar(stat = 'identity') + coord_flip() +
    ylab(expression(paste("[MAP ", chi[g] ,"]" ^ "-1"))) 
}


#' Estimate the dropout parameter
#' 
#' @param y A cell-by-gene expression matrix
#' @param lower_limit The limit below which expression counts as 'dropout'
#' 
#' @export
#' @return The estimated lambda
empirical_lambda <- function(y, lower_limit = 0) {
  means <- colMeans(y)
  pdrop <- colMeans(y == 0)
  fit <- nls(pdrop ~ exp(-lambda * means), start = list(lambda = 1))
  coef(fit)['lambda']
}

