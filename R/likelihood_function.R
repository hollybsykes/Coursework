#create a function for the likelihood
prob_y1 <- function(y,n,theta){
  choose(n, y)*theta^(y)*(1-theta)^(n-y)
}

prob_y1(0,39,0.5)
