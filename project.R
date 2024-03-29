set.seed(2305679)
#initial guess of the sample size
#we know that alpha>=0.05 and beta>=0.2
#we know delta = 0.7- 0.5 = 0.2
  initial<-function(sigma, delta, alpha, beta){
    n<-sigma*(qnorm(1-alpha)-qnorm(beta))^2/delta^2
    return(n)
  }
#using the parameters of the trial
initial(0.25,0.2,0.05,0.2)
#38.64098 , round up no 39=n1
#choose n2=78 = 39*2

#now we want to find an appropriate gamma and alpha
#finding the expected sample size for any hypothesis
evaluate_design <- function(gamma, lambda, n1, n2, theta) {
  
  #stage 1
  
  #set number of simulations 
  M1 <- 10^4
  
  #simulation observations using prior distribution
  y1 <- rbinom(M1, n1, theta)
  
  #find posterior distribution parameters
  a1 <- 0.5 + y1
  b1 <- 0.5 + n1 - y1
  
  #find probability of futility (posterior distribution)
  fut1 <- pbeta(theta, a1, b1)
  
  #decision value
  c1 <- 1 - lambda * (n1 / n2)^gamma
  
  #number of successes, is used in stage 2
  M2<- sum(fut1 > c1)
  
  #simulation y2 
  y2 <- rbinom(M2, n2, theta)
  
  #set new posterior parameters
  a2 <- 0.5 + y2
  b2 <- 0.5 + n2 - y2
  
  #find propbability of futility
  fut2 <- pbeta(theta, a2, b2)
  
  #find decision value
  c2<- 1 - lambda * (n2 / n2)^gamma
  
  #find the type I and type II error using monte carlo estimates
  return(c(typeI = sum(fut1>c1)/M1, typeII=sum(fut2>c2)/M2))
}
#example - test it works
evaluate_design(0.4,0.75,50,70,0.5)

#we want to do a grid search using this function fro type I and type II errors
#create a data frame of all different pairs of gamma and lambda
gamma<- seq(0.05,1,0.05)
lambda<- seq(0.05,1,0.05)
theta=0.5
df<- expand.grid(gamma=gamma,lambda=lambda)
df

#change the function to allow for a grid search
evaluate_design_gridsearch <- function(x, n1, n2, theta) {
  
  #extract values of gamma and lambda form data frame
  gamma<- x[1]
  lambda<- x[2]
  
  #stage 1
  
  #set number of simulations 
  M1 <- 10^4
  
  #simulation observations using prior distribution
  y1 <- rbinom(M1, n1, theta)
  
  #find posterior distribution parameters
  a1 <- 0.5 + y1
  b1 <- 0.5 + n1 - y1
 
  #find probability of futility (posterior distribution)
  fut1 <- pbeta(theta, a1, b1)
  
  #decision value
  c1 <- 1 - lambda * (n1 / n2)^gamma
  
  #number of successes, is used in stage 2
  M2<- sum(fut1>c1)
  
  #simulation y2 
  y2 <- rbinom(M2, n2, theta)
  
  #set new posterior parameters
  a2 <- 0.5 + y2
  b2<- 0.5 + n2 - y2
  
  #find propbability of futility
  fut2 <- pbeta(theta, a2, b2)
  #find decision value
  c2<- 1 - lambda * (n2 / n2)^gamma
  
  #find the type I and type II error using monte carlo estimates
  return(c(typeI = sum(fut1>c1)/M1, typeII=sum(fut2>c2)/M2))
}

fun<- apply(df, 1, evaluate_design_gridsearch, n1=38, n2=78, theta=0.5)
fun 
rownames(fun)
#find which values of lambda and gamma satisfy the constraints on errors
v<-which(fun[1,]<=0.05 & fun[2,]<=0.2)

#make a new data frame of appropiate gamma and lambdas
df2<-df[v,]

#create a function for the likelihood
prob_y1 <- function(y,n,theta){
  choose(n, y)*theta^(y)*(1-theta)^(n-y)
}
prob_y(10,30,0.5)
#check this function works - 
n<-30
theta<-0.5
sum(prob_y(0:n, n, 0.5)) # this sums to one -so okay!

#function for the sample size
sample_size <- function(x, n1, n2, theta) {
  gamma<-x[1]
  lambda<-x[2]
  C1 <- 1 - lambda * (n1 / n2)^gamma
  
  # Vector of possible stage 1 outcomes.
  y_1s <- 0:n1
  
  # Vector of corresponding progression decisions.
  stops <- pbeta(0.5, y_1s + 0.5, n1 - y_1s + 0.5) < C1
  
  # For each outcome, calculate its probability.
  y_1_probs <- prob_y1(y_1s, n1, theta)
  
  sum(n1 * stops * y_1_probs + n2 * (!stops) * y_1_probs)
}

samplesize = apply(df2, 1, sample_size, n1=39, n2=78, theta=0.5)
min(samplesize)

t<-which(round(samplesize,5)==40.03841)
t
# how do you decide which one
df2
df3<-df2[t,]
df3
#question3 - how does the changing the hypothesis value change the error values

#create a grid of different values of theta
theta_grid<- seq(0.1,0.99,0.01)
theta_expand<-expand.grid(theta=theta_grid)

#change function to except different values of theta
evaluate_design_theta <- function(gamma, lambda, n1, n2, x) {
  
  #stage 1
  theta<-x[1]
  #set number of simulations 
  M1 <- 10^4
  
  #simulation observations using prior distribution
  y1 <- rbinom(M1, n1, theta)
  
  #find posterior distribution parameters
  a1 <- 0.5 + y1
  b1 <- 0.5 + n1 - y1
  
  #find probability of futility (posterior distribution)
  fut1 <- pbeta(theta, a1, b1)
  
  #decision value
  c1 <- 1 - lambda * (n1 / n2)^gamma
  
  #number of successes, is used in stage 2
  M2<- sum(fut1 > c1)
  
  #simulation y2 
  y2 <- rbinom(M2, n2, theta)
  
  #set new posterior parameters
  a2 <- 0.5 + y2
  b2 <- 0.5 + n2 - y2
  
  #find propbability of futility
  fut2 <- pbeta(theta, a2, b2)
  
  #find decision value
  c2<- 1 - lambda * (n2 / n2)^gamma
  
  #find the type I and type II error using monte carlo estimates
  return(c(typeI = sum(fut1>c1)/M1, typeII=sum(fut2>c2)/M2))
}

#apply values of theta to the function
fun1<- apply(theta_expand, 1, evaluate_design_theta, gamma =0.4, lambda=0.05, n1=38, n2=78)
fun1

#plot theta against the the errors setting other parameters to be constant
par(mfrow=c(1,2))
plot(theta_grid,fun1[1,], xlab="theta", ylab="type I" )
plot(theta_grid, fun1[2,], xlab="theta", ylab="typeII")
