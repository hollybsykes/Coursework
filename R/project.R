set.seed(2305679)
#install.packages("styler") - use "style selection

# initial guess of the sample size
# we know that alpha>=0.05 and beta>=0.2
# we know delta = 0.7- 0.5 = 0.2
initial <- function(sigma, delta, alpha, beta) {
  n <- sigma * (qnorm(1 - alpha) - qnorm(beta))^2 / delta^2
  return(n)
}
# using the parameters of the trial
initial(0.25, 0.2, 0.05, 0.2)
# 38.64098 , round up no 39=n1
# choose n2=78 = 39*2


# now we want to find an appropriate gamma and alpha
# finding the expected sample size for any hypothesis
#could it be 2 values of theta
#evaluate_design <- function(gamma, lambda, n1, n2, theta) {

  # stage 1

  # set number of simulations
  #M1 <- 10^4

  # simulation observations using prior distribution
  #y1 <- rbinom(M1, n1, theta)

  # find posterior distribution parameters
  #a1 <- 0.5 + y1
  #b1 <- 0.5 + n1 - y1

  # find probability of futility (posterior distribution)
  #fut1 <- pbeta(0.5, a1, b1)

  # decision value
 # c1 <- 1 - lambda * (n1 / n2)^gamma

  # number of successes, is used in stage 2
  #M2 <- sum(fut1 < c1)

  # simulation y2
  #y2 <- rbinom(M2, n2, theta)

  # set new posterior parameters
  #a2 <- 0.5 + y2
 # b2 <- 0.5 + n2 - y2

  # find propbability of futility
  #fut2 <- pbeta(0.5, a2, b2)

  # find decision value
  #c2 <- 1 - lambda * (n2 / n2)^gamma

  # find the type I and type II error using monte carlo estimates
 # return(c(typeI = sum(fut2 <c2) / M2))
#}

#evaluate_design(0.4, 0.75, 50, 70, 0.5)


# we want to do a grid search using this function fro type I and type II errors
# create a data frame of all different pairs of gamma and lambda
ptm<-proc.time()
gamma <- seq(0.5, 5, 0.5)
lambda <- seq(0.05, 0.95, 0.05)

df <- expand.grid(gamma = gamma, lambda = lambda)
df

# change the function to allow for a grid search
evaluate_design_gridsearch <- function(x, n1, n2, theta) {

  # extract values of gamma and lambda form data frame
  gamma <- x[1]
  lambda <- x[2]

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

  # number of successes, is used in stage 2
  #M2 <- sum(fut1 < c1)

  # simulation y2
  y2 <- rbinom(M1, n2 - n1, theta)

  # set new posterior parameters
  a2 <- 0.5 + y1 + y2 
  b2 <- 0.5 + n2- y2 - y1

  # find probability of futility
  fut2 <- pbeta(0.5, a2, b2)
  # find decision value
  c2 <- 1 - lambda * (n2 / n2)^gamma

  # find the type I and type II error using monte carlo estimates
  return(c(typeI = sum(fut2 < c2 & fut1 < c1) / M1))
}

type_I <- apply(df, 1, evaluate_design_gridsearch, n1 = 39, n2 = 78, theta = 0.5)
type_I
typeII <- apply(df, 1, evaluate_design_gridsearch, n1=39, n2 = 78, theta = 0.7)
type_II <- 1-typeII

type_I_II <- cbind(type_I, type_II)
#head(type_I_II)
#tail(type_I_II)
# find which values of lambda and gamma satisfy the constraints on errors
v <- which(type_I_II[,1] <= 0.05 & type_I_II[,2] <= 0.2)
v

type_I_II


# make a new data frame of appropiate gamma and lambdas
df2 <- df[v, ]
df2

source("./R/likelihood_function.R")
#library("testthat")
test_dir("./tests")

# function for the sample size - using likelihood function
sample_size <- function(x, n1, n2, theta) {
  gamma <- x[1]
  lambda <- x[2]
  #find the decision rule for stage 1
  C1 <- 1 - lambda * (n1 / n2)^gamma

  # Vector of possible stage 1 outcomes.
  y_1s <- 0:n1

  # Vector of corresponding progression decisions - a vector of true or false
  stops <- pbeta(0.5, y_1s + 0.5, n1 - y_1s + 0.5) < C1

  # For each outcome, calculate its probability
  y_1_probs <- prob_y1(y_1s, n1, theta)

  #if stops=TRUE, gives n1
  #if stops=FLASE, gives n2
  sum(n1 * stops * y_1_probs + n2 * (!stops) *y_1_probs)
}

samplesize <- apply(df2, 1, sample_size, n1 = 20, n2 = 40, theta = 0.5)

round_up=ceiling(samplesize)

final<- cbind(df2, round_up)
final
proc.time()-ptm

# question3 - how does the changing the hypothesis value change the error values
# create a grid of different values of theta
theta_grid <- seq(0.1, 0.99, 0.01)
theta_expand <- expand.grid(theta = theta_grid)

# change function to except different values of theta
sample_size_theta <- function(gamma, lambda, n1, n2, x) {
  theta<-x[1]
  #find the decision rule for stage 1
  C1 <- 1 - lambda * (n1 / n2)^gamma
  
  # Vector of possible stage 1 outcomes.
  y_1s <- 0:n1
  
  # Vector of corresponding progression decisions - a vector of true or false
  stops <- pbeta(0.5, y_1s + 0.5, n1 - y_1s + 0.5) < C1
  
  # For each outcome, calculate its probability
  y_1_probs <- prob_y1(y_1s, n1, theta)
  
  #if stops=TRUE, gives n1
  #if stops=FLASE, gives n2
  sum(n1 * stops * y_1_probs + n2 * (!stops) *y_1_probs)
}
# apply values of theta to the function
fun1 <- apply(theta_expand, 1, sample_size_theta, gamma = 5, lambda = 0.95, n1=39, n2=78)
fun1


# plot theta against the the errors setting gamma and lambda to be constant
par(mfrow=c(1,2))
plot(theta_grid, fun1, xlab="theta", ylab="expected sample size")



#sensitivity of the error rates to the null/alternative hypothesis values
evaluate_design_gridtheta <- function(gamma, lambda, n1, n2, x) {
  
  # extract values of gamma and lambda form data frame
  theta<-x[1]
  
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
  
  # number of successes, is used in stage 2
  #M2 <- sum(fut1 < c1)
  
  # simulation y2
  y2 <- rbinom(M1, n2 - n1, theta)
  
  # set new posterior parameters
  a2 <- 0.5 + y1 + y2 
  b2 <- 0.5 + n2- y2 - y1
  
  # find probability of futility
  fut2 <- pbeta(0.5, a2, b2)
  # find decision value
  c2 <- 1 - lambda * (n2 / n2)^gamma
  
  # find the type I and type II error using monte carlo estimates
  return(c(typeI = sum(fut2 < c2 & fut1 < c1) / M1))
}

type_I <- apply(theta_expand, 1, evaluate_design_gridtheta, n1 = 39, n2 = 78, gamma=5, lambda=0.95)
length(type_I)
type_II <- 1-type_I

par(mfrow=c(1,2))
plot(theta_grid, type_I, xlab="theta", ylab = "Type I error")
plot(theta_grid, type_II, xlab="theta", ylab = "Type II error")




