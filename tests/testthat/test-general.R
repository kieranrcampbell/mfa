context("General package functionality")

test_that("A call to MFA returns a valid object", {
  iter <- 200
  burn <- 100
  thin <- 2
  B <- 2
  C <- 40
  G <- 10
  
  synth <- create_synthetic(C = C, G = G)
  
  m <- mfa(synth$X, iter = iter, burn = burn, thin = thin)
  
  expect_is(m, "mfa")
  expect_equal(iter, m$iter)
  expect_equal(burn, m$burn)
  expect_equal(thin, m$thin)
  expect_equal(C, m$N)
  expect_equal(G, m$G)
})

test_that("A call to create_synthetic returns a valid object", {
  C <- 40
  G <- 10
  
  synth <- create_synthetic(C = C, G = G)
  
  expect_is(synth, 'list')
  expect_is(synth$X, 'matrix')  

  expect_equal(dim(synth$X), c(C, G))
  
  expect_is(synth$branch, 'integer')  
  expect_is(synth$pst, 'numeric')  

  expect_equal(length(synth$branch), C)
  expect_equal(length(synth$pst), C)
})


