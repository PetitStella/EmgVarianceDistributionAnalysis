# EMG variance distribution estimation

library(stats)

# Calculate log-marginal likelihood function
log.marginal.likelihood <- function(x, nu, lambda)
{
  alpha <- nu * 0.5
  beta <- lambda * alpha
  
  len <- length(x)
  
  return(0.5*len*log(2.0*pi) - len*alpha*log(beta) + len*log(gamma(alpha)) - 
           len*log(gamma(alpha+0.5)) + (alpha+0.5)*sum(log(beta + x^2*0.5)))
}

# Calculate likelihood function for update of nu
nu.likelihood <- function(x, lambda, omega)
{
  len <- length(x)
  return(function(nu){
    return(-digamma(nu*0.5) + log(nu*0.5) + 1.0 + mean(log(omega)-omega) + 
             mean(digamma((nu+1.0)*0.5) - log((nu+1.0)*0.5) ))
  })
}

# Estimate variance distribution
estimate.var.dist <- function(x, maxIter = 1000, tol = 0.0000001)
{
  xLength <- length(x)
  
  # Initialize
  nu <- runif(1, min = 1.0, max = 20.0)
  lambda <- var(x)
  
  preJ <- log.marginal.likelihood(x, nu, lambda)
  
  for(i in 1:maxIter){
    
    # E-step
    delta <- x^2/lambda
    omega <- (nu+1.0)/(nu+delta)
    
    # M-step
    lambda <- mean(omega*x^2)
    nuOptim <- uniroot(nu.likelihood(x, lambda, omega), interval=c(0.1,200), tol = 0.001)
    nu <- nuOptim$root
    
    # Calculate log-marginal likelihood
    J <- log.marginal.likelihood(x, nu, lambda)
    cat(sprintf("%f\n", J))

    if(abs(preJ - J) < tol){
      alpha <- nu * 0.5
      beta <- lambda * alpha
      return(list(alpha,beta))
    }
    preJ <- J
  }
  warning("reach maximum iteration")
  return(list(alpha,beta))
}

