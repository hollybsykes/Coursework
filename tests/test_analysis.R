source("../R/likelihood_function.R")

#test that function sums to one
#install.packages("testthat")
library("testthat")
test_that("y1 probs sum to one", {
  n<-39
  theta<-0.5
  s <- sum(prob_y1(0:n, n, theta))
  expect_equal(s, 1)
})



