context("Conditional distributions through Gwecke tests")



## This implements a series of Gwecke tests 
## Where we check the conditional distributions for
## Gibbs sampling are correct (assuming that our calculation
## of the log-posterior is correct and vice-versa)

log_dnorm <- function(x, mean, precision) dnorm(x, mean, 1 / sqrt(precision), log = TRUE)


# Make some synthetic data ------------------------------------------------


N <- 80
G <- 8
b <- B <- 3

alpha = 1 
beta = 1
theta_tilde = 0
eta_tilde = 0
tau_eta = 1
tau_theta = 1
tau_c = 1
alpha_chi = 1
beta_chi = 1
w_alpha = 1 / b
r <- 1

w <- c(0.2, 0.3, 0.5)

gamma <- sample(seq_len(B), N, prob = w, replace = TRUE)

eta <- rnorm(1)
theta <- rnorm(G, 0, 1)

## precision parameters
tau <- rep(1, G)

chi <- rgamma(G, alpha_chi, beta_chi)

c <- matrix(rnorm(G*B), ncol = B)
k <- matrix(rnorm(G*B), ncol = B)

## pseudotime parameters
pst <-  rnorm(N)

## parameter initialisation
y <- sapply(1:N, function(i) {
  m <- c[,gamma[i]] + k[,gamma[i]] * pst[i]
  rnorm(G, m, 1 / sqrt(tau))
})
y <- t(y)

test_that("Conditional distribution of k is correct", {
  nu_k <- calculate_nuk(y, pst, c[,1], tau, theta, chi, gamma == 1)
  lam_k <- calculate_lamk(chi, tau, pst, gamma == 1)
  
  kp <- k; kp[1] <- 0
  
  rhs <- posterior(y, c, kp, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
                   theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) -
    posterior(y, c, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
              theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi)
  
  lhs <- sum(log_dnorm(kp[,1], nu_k / lam_k, lam_k)) -
    sum(log_dnorm(k[,1], nu_k / lam_k, lam_k))
  
  expect_equal(lhs, rhs)
})

test_that("Conditional distirbution of c is correct", {
  nu_c <- calculate_nuc(y, pst, k[,1], tau, eta, tau_c, gamma == 1)
  lam_c <- calculate_lamc(tau, tau_c, sum(gamma == 1))
  
  cp <- c; cp[1] <- 0
  
  rhs <- posterior(y, cp, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
                   theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) -
    posterior(y, c, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
              theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi)
  
  lhs <- sum(log_dnorm(cp[,1], nu_c / lam_c, lam_c)) -
    sum(log_dnorm(c[,1], nu_c / lam_c, lam_c))
  
  expect_equal(lhs, rhs)
})

test_that("Conditional distribution of t is correct", {
  pst_par <- pst_update_par(y, c, k, r, gamma, tau)
  
  pstp <- pst; pstp[1] <- 0
  
  
  rhs <- posterior(y, c, k, pstp, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
                   theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) -
    posterior(y, c, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
              theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi)
  
  lhs <- log_dnorm(pstp[1], pst_par[1,1], pst_par[1,2]) -
    log_dnorm(pst[1], pst_par[1,1], pst_par[1,2])
  
  expect_equal(lhs, rhs)
})

test_that("Conditional distribution of tau is correct", {
  tau_par <- tau_params(y, c, k, gamma, pst, alpha, beta)
  
  taup <- tau
  taup[1] <- 2
  
  rhs <- posterior(y, c, k, pst, taup, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
                   theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) -
    posterior(y, c, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
              theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi)
  lhs <- dgamma(taup[1], tau_par[1,1], tau_par[1,2], log = TRUE) -
    dgamma(tau[1], tau_par[1,1], tau_par[1,2], log = TRUE)
  
  expect_equal(lhs, rhs)
})

test_that("Conditional distribution of theta is correct", {
  lambda_theta <- b * chi + tau_theta
  nu_theta <- tau_theta * theta_tilde + chi * rowSums(k)
  nu_theta <- nu_theta / lambda_theta
  
  thetap <- theta; thetap[1] <- 0
  
  
  rhs <- posterior(y, c, k, pst, tau, gamma, thetap, eta, chi, w, tau_c, r, alpha, beta,
                   theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) -
    posterior(y, c, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
              theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi)
  
  lhs <- log_dnorm(thetap[1], nu_theta[1], lambda_theta[1]) -
    log_dnorm(theta[1], nu_theta[1], lambda_theta[1])
  
  expect_equal(lhs, rhs)
})

test_that("Conditional distribution of eta is correct", {
  lambda_eta <- tau_eta + G * b * tau_c
  nu_eta <- tau_eta * eta_tilde + tau_c * sum(c)
  nu_eta <- nu_eta / lambda_eta
  
  etap <- eta
  etap[1] <- 0
  
  rhs <- posterior(y, c, k, pst, tau, gamma, theta, etap, chi, w, tau_c, r, alpha, beta,
                   theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) -
    posterior(y, c, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
              theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi)
  
  lhs <- log_dnorm(etap[1], nu_eta[1], lambda_eta[1]) -
    log_dnorm(eta[1], nu_eta[1], lambda_eta[1])
  
  expect_equal(lhs, rhs)
})

test_that("Conditional distribution of chi is correct", {
  alpha_new <- alpha_chi + b / 2
  beta_new <- beta_chi + 0.5 * rowSums( (k - theta)^2 )
  
  chip <- chi
  chip[1] <- 2
  
  rhs <- posterior(y, c, k, pst, tau, gamma, theta, eta, chip, w, tau_c, r, alpha, beta,
                   theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) -
    posterior(y, c, k, pst, tau, gamma, theta, eta, chi, w, tau_c, r, alpha, beta,
              theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi)
  
  lhs <- dgamma(chip[1], alpha_new, beta_new[1], log = TRUE) -
    dgamma(chi[1], alpha_new, beta_new[1], log = TRUE)
  
  expect_equal(lhs, rhs)
})

test_that("Conditional distribution of pi is correct", {
  pi <- calculate_pi(y, c, k, pst, tau, eta, tau_c, FALSE, log(w), log_result = TRUE)
  
  gammap <- gamma
  gammap[1] <- 2
  
  lhs <- pi[1, gammap[1]] - pi[1, gamma[1]]
  
  rhs <- posterior(y, c, k, pst, tau, gammap, 
                   theta, eta, chi, w, tau_c, r, alpha, beta,
                   theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) -
    posterior(y, c, k, pst, tau, gamma, 
              theta, eta, chi, w, tau_c, r, alpha, beta,
              theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) 
  
  expect_equal(lhs, rhs)
})




