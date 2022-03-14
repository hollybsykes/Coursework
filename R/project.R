# set the seed to make sure code is repeatable
set.seed(2305679)

# install.packages("styler") - use "style selection

# time how long optimising the parameters lambda and gamma takes
ptm <- proc.time()
# initial guess of the sample size
initial <- function(sigma, delta, alpha, beta) {
  n <- sigma * (qnorm(1 - alpha) - qnorm(beta))^2 / delta^2
  return(n)
}
# using parameters given
initial(0.25, 0.2, 0.05, 0.2) # n1=39, n2=78


# function to find type I and type II errors
# a function which allows a grid of values for lambda and gamma
# change the function to allow for a grid search
evaluate_design_gridsearch <- function(x, n1, n2, theta) {

  # extract values of gamma and lambda from data frame
  gamma <- x[1]
  lambda <- x[2]


  # stage 1

  # set number of simulations
  M1 <- 10^4

  # simulation observations using prior distribution
  y1 <- rbinom(M1, n1, theta)

  # set new posterior distribution parameters
  a1 <- 0.5 + y1
  b1 <- 0.5 + n1 - y1

  # find probability of futility
  fut1 <- pbeta(0.5, a1, b1)

  # decision value
  c1 <- 1 - lambda * (n1 / n2)^gamma

  # stage 2

  # simulation y2
  y2 <- rbinom(M1, n2 - n1, theta)

  # set new posterior parameters
  a2 <- 0.5 + y1 + y2
  b2 <- 0.5 + n2 - y2 - y1

  # find probability of futility
  fut2 <- pbeta(0.5, a2, b2)

  # find decision value
  c2 <- 1 - lambda * (n2 / n2)^gamma

  # find the type I and type II error using monte carlo estimates
  return(c(typeI = sum(fut2 < c2 & fut1 < c1) / M1))
}

# make a grid of possible combinations of decision parameters
gamma <- seq(0.5, 5, 0.5)
lambda <- seq(0.05, 0.95, 0.05)
df <- expand.grid(gamma = gamma, lambda = lambda)

# apply function to all pairs of decision parameters

# type I error under null hypothesis value
type_I <- apply(df, 1, evaluate_design_gridsearch, n1 = 39, n2 = 78, theta = 0.5)

# type II error under alternative hypothesis value
typeII <- apply(df, 1, evaluate_design_gridsearch, n1 = 39, n2 = 78, theta = 0.7)
type_II <- 1 - typeII

# make matrix of type I and type II errors
type_I_II <- cbind(type_I, type_II)

# find which values of lambda and gamma satisfy the constraints on errors
v <- which(type_I_II[, 1] <= 0.05 & type_I_II[, 2] <= 0.2)

# make a new data frame of corresponding appropriate gamma and lambdas
df2 <- df[v, ]


# use function of probability of observing responses in n trials
source("./R/likelihood_function.R")

# use a unit test to check function sums to one for all values of y
# library("testthat")
test_dir("./tests")

# function for the sample size
sample_size <- function(x, n1, n2, theta) {

  # extract values of gamma and lambda from data frame
  gamma <- x[1]
  lambda <- x[2]

  # find the decision rule for stage 1
  C1 <- 1 - lambda * (n1 / n2)^gamma

  # Vector of possible stage 1 outcomes.
  y_1s <- 0:n1

  # Vector of corresponding progression decisions - a vector of true or false
  stops <- pbeta(0.5, y_1s + 0.5, n1 - y_1s + 0.5) < C1

  # For each outcome, calculate its probability
  y_1_probs <- prob_y1(y_1s, n1, theta)

  # if stops=TRUE, gives n1
  # if stops=FALSE, gives n2
  sum(n1 * stops * y_1_probs + n2 * (!stops) * y_1_probs)
}

# apply function to all possible pairs of decision parameters under null
samplesize <- apply(df2, 1, sample_size, n1 = 20, n2 = 40, theta = 0.5)

# round up sample sizes to integer value
round_up <- ceiling(samplesize)

# make a matrix of decision parameter values and corresponding sample sizes
final <- cbind(df2, round_up)

# stop timing
proc.time() - ptm

# evaluation

# change sample size function to now take different values of theta
sample_size_theta <- function(gamma, lambda, n1, n2, x) {

  # extract value of theta from data frame
  theta <- x[1]

  # find the decision rule for stage 1
  C1 <- 1 - lambda * (n1 / n2)^gamma

  # Vector of possible stage 1 outcomes.
  y_1s <- 0:n1

  # Vector of corresponding progression decisions - a vector of true or false
  stops <- pbeta(0.5, y_1s + 0.5, n1 - y_1s + 0.5) < C1

  # For each outcome, calculate its probability
  y_1_probs <- prob_y1(y_1s, n1, theta)

  # if stops=TRUE, gives n1
  # if stops=FLASE, gives n2
  sum(n1 * stops * y_1_probs + n2 * (!stops) * y_1_probs)
}


# create a grid of different values of theta
theta_grid <- seq(0.1, 0.99, 0.01)
theta_expand <- expand.grid(theta = theta_grid)

# apply values of theta to the function of sample size
# with optimized decision parameters
fun1 <- apply(theta_expand, 1, sample_size_theta, gamma = 5, lambda = 0.95, n1 = 39, n2 = 78)

# plot theta against the expected sample sizes
plot(theta_grid, fun1, xlab = "theta", ylab = "expected sample size")


# sensitivity of the error rates to the null/alternative hypothesis values
evaluate_design_gridtheta <- function(gamma, lambda, n1, n2, x) {

  # extract values of gamma and lambda form data frame
  theta <- x[1]

  # stage 1

  # set number of simulations
  M1 <- 10^4

  # simulation observations using prior distribution
  y1 <- rbinom(M1, n1, theta)

  # find posterior distribution parameters
  a1 <- 0.5 + y1
  b1 <- 0.5 + n1 - y1

  # find probability of futility (posterior distribution)
  fut1 <- pbeta(0.5, a1, b1)

  # decision value
  c1 <- 1 - lambda * (n1 / n2)^gamma

  # stage 2

  # simulation y2
  y2 <- rbinom(M1, n2 - n1, theta)

  # set new posterior parameters
  a2 <- 0.5 + y1 + y2
  b2 <- 0.5 + n2 - y2 - y1

  # find probability of futility
  fut2 <- pbeta(0.5, a2, b2)
  # find decision value
  c2 <- 1 - lambda * (n2 / n2)^gamma

  # find the type I and type II error using monte carlo estimates
  return(c(typeI = sum(fut2 < c2 & fut1 < c1) / M1))
}

# find type I error for different theta values
type_I <- apply(theta_expand, 1, evaluate_design_gridtheta, n1 = 39, n2 = 78, gamma = 5, lambda = 0.95)

# find type II error for different theta values
type_II <- 1 - type_I

# make a matrix of theta and errors
cbind(type_I, type_II, theta_grid)

# plot the errors against theta values
par(mfrow = c(1, 2))
plot(theta_grid, type_I, xlab = "theta", ylab = "Type I error")
plot(theta_grid, type_II, xlab = "theta", ylab = "Type II error")

# trial now uses theta_0 = 0.4
# find an initial guess of sample size
initial(0.24, 0.3, 0.05, 0.2)
# n1=17, n2=34

# apply gridsearch function with different sample sizes and null value
# find type I error under null hypothesis value
type_I_eval <- apply(df, 1, evaluate_design_gridsearch, n1 = 17, n2 = 34, theta = 0.4)

# find type II error under alternative hypothesis value
typeII_eval <- apply(df, 1, evaluate_design_gridsearch, n1 = 17, n2 = 34, theta = 0.7)
type_II_eval <- 1 - typeII_eval

# make a matrix of type I and type II values
type_I_II_eval <- cbind(type_I_eval, type_II_eval)

# find which values of lambda and gamma satisfy the constraints on errors
v_eval <- which(type_I_II_eval[, 1] <= 0.05 & type_I_II_eval[, 2] <= 0.2)

# find the errors which correspond to these values of the decision parameters
type_eval <- type_I_II_eval[v_eval, ]

# make a new data frame of appropiate gamma and lambdas
df2_eval <- df[v_eval, ]

# use sample size function to find the estimates of sample size of these pairs of decision parameters
samplesize_eval <- apply(df2_eval, 1, sample_size, n1 = 17, n2 = 34, theta = 0.4)

# round up to integer value
round_up_eval <- ceiling(samplesize_eval)

# make a matrix of decision parameters, errors and sample size
final_eval <- cbind(df2_eval, round_up_eval, type_eval)

# find corresponding parameters which give the minimum sample size
yes <- which(round_up_eval == 20)
table <- final_eval[yes, ]

# repeat for theta_0= 0.45
initial(0.2475, 0.25, 0.05, 0.2)
# n1=25, n2=50

# grid search with different sample size and null value
# type I error under the null hypothesis
type_I_eval2 <- apply(df, 1, evaluate_design_gridsearch, n1 = 25, n2 = 50, theta = 0.45)

# type II error under the alternative hypothesis
typeII_eval2 <- apply(df, 1, evaluate_design_gridsearch, n1 = 25, n2 = 50, theta = 0.7)
type_II_eval2 <- 1 - typeII_eval2

# make a matrix of the errors
type_I_II_eval2 <- cbind(type_I_eval2, type_II_eval2)


# find which values of lambda and gamma satisfy the constraints on errors
v_eval2 <- which(type_I_II_eval2[, 1] <= 0.05 & type_I_II_eval2[, 2] <= 0.2)

# find the errors which correspond to these values of the decision parameters
type_eval2 <- type_I_II_eval2[v_eval2, ]

# make a new data frame of appropriate gamma and lambdas
df2_eval2 <- df[v_eval2, ]

# find the sample size function to find the sample sizes
samplesize_eval2 <- apply(df2_eval2, 1, sample_size, n1 = 25, n2 = 50, theta = 0.45)

# round up to an integer value
round_up_eval2 <- ceiling(samplesize_eval2)

# find th minimal sample size
min(round_up_eval2)

# make a matrix with the errors, sample size and decision parameters
final_eval2 <- cbind(df2_eval2, round_up_eval2, type_eval2)

# find corresponding parameters for the minimal sample size
yes2 <- which(round_up_eval2 == 27)
final_eval2[yes2, ]

