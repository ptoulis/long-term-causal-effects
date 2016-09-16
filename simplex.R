## Distribution/Diffusion on the simplex.
## Has the likelihood for temporal evolution of continuous proportions.
##
##
library(mvtnorm)

cos.angle <- function(a, b) {
  # angle between two vectors.
  a2 = sum(a^2)
  b2 = sum(b^2)
  if(a2 == 0 | b2 == 0) return(1)
  sum(a * b) / (sqrt(a2) * sqrt(b2))
}

simplex.fix.proportion <- function(v, tol=1e-3) {
  i = which(v < tol)
  if(length(i) > 0) {
    s = sum(v[i])
    v[i] <- tol
  }
  return(v / sum(v))
}

simplex.logit <- function(proportion) {
  # Returns (p-1) vector
  proportion = as.vector(proportion)
   X1 = 0
  Z1 = proportion[1]
  Xi = log(tail(proportion, length(proportion)-1) / Z1)
  Xi[Xi < -25] <- -25
  Xi[Xi > 25] <- 25
  return(c(Xi))
}

simplex.logistic <- function(v) {
  # assumes it was produced by simplex.logit
  v = c(0, v) 
  v[v > 20] <- 20
  v[v < -20] <- -20
  u = exp(v) / sum(exp(v))
  i = which(u < 1e-7)
  if(length(i) > 0) {
    s = sum(u[i])
    u[i] <- 0
    u[-i] <- u[-i] / sum(u[-i])
  }
  return(u)
}

check.psi <- function(psi) {
  stopifnot(psi$sigma >= 0, psi$beta >= 0, psi$beta < 1)
  stopifnot(length(psi$mu) > 0, all(psi$mu >= 0), abs(1- sum(psi$mu))< 1e-3)
}

diffusion.logLikelihood <- function(B_0T, psi) {
  # Beta 0:T is nxp matrix of proportions.
  #
  # Assume Xt+1 = fi * Xt + (1- fi) * (mu + et)
  #
  check.psi(psi)
  
  n = nrow(B_0T) # number of steps.
  p = ncol(B_0T) # dimension of the simplex

  # unload params.
  mu.X = simplex.logit(psi$mu)
  sigma = psi$sigma
  beta = psi$beta
  
  ll <- 0
  
  for(i in 2:n) {
    bt.next = B_0T[i, ]  # next proportion
    bt.prev = B_0T[i-1, ]  # previous proportion.
    stopifnot(abs(1-sum(bt.next)) < 1e-3, abs(1-sum(bt.prev)) < 1e-3)
    # transform
    Xt.next = simplex.logit(bt.next)  # NOTE: (p-1) dimension
    Xt.prev = simplex.logit(bt.prev)
    
    pivot = (Xt.next - beta * Xt.prev) - mu.X  # standard normal.
    # get log-likelihood.
    ll <- ll + dmvnorm(pivot, mean=rep(0, p-1), sigma=sigma * diag(p-1), log=T)
    if(is.na(ll)) {
      stop("LL in simplex is NA")
    }
  }
  return(ll)
}

sample.diffusion <- function(p=3, n=20, beta=0.5, sigma=1.1,
                             mu = seq(1,p)/sum(seq(1,p))) {

  psi = list(sigma=sigma, beta=beta, mu=mu)
  check.psi(psi)
  
  print("True mu")
  print(round(mu, 3))
  print(sprintf("Beta = %.3f sigma= %.3f", beta, sigma))
  # Compute params
  S = sigma * diag(p-1)
  X = matrix(NA, nrow=n, ncol=p)
  B = matrix(NA, nrow=n, ncol=p)
  B[1,] <- random.proportion(p)
  X[1,] <- simplex.logit(B[1,]) ## because first is 0
  
  mu.X = simplex.logit(mu)
  # create dataset
  for(t in 2:n) {
    epsilon = c(rmvnorm(1, mean=rep(0, p-1), sigma=S))
    X[t,] = beta * X[t-1,] + (mu.X + epsilon) 
    stopifnot(X[t, 1] == 0)
    B[t,] = simplex.logistic(X[t,])
  }
  return(list(B=B, true.psi=list(beta=beta, sigma=sigma, mu=mu)))
}

diffusion.mle <- function(B0t) {
  f = function(par) {
    p = ncol(B0t)
    mu = simplex.logistic(tail(par, p-1))
    psi = list(beta=par[1], sigma=par[2], mu=mu)
    diffusion.logLikelihood(B0t, psi)
  }
  
  translate.par <- function(par) {
    print(sprintf("Beta=%.3f sigma=%.3f", par[1], par[2]))
    mu = simplex.logistic(tail(par, p-1))
    print(sprintf("mu = %s", paste(round(mu, 3), collapse=", ")))
    return(list(beta=par[1], sigma=par[2], mu=mu))
  }
  
  p = ncol(B0t)
  fit = optim(par=c(0.5, 1, runif(p-1)), fn=f, method="L-BFGS-B",
              lower=c(1e-5, 1e-4, rep(-10, p-1)), upper=c(1-1e-4, rep(9, p+1)),
              control=list(fnscale=-1))
  x = translate.par(fit$par)  
  return(x)
}

simplex.impute <- function(nsamples, B0t, psi) {
  ## Impute behaviors.
  n = nrow(B0t)
  p = ncol(B0t)
  out = matrix(NA, nrow=nsamples, ncol=p)
  last = B0t[n, ]
  Xt = simplex.logit(last)
  mu.X = simplex.logit(psi$mu)
  for(i in 1:nsamples) {
    X.next = psi$beta * Xt + mu.X + rnorm(p-1, mean=0, sd=psi$sigma)
    out[i, ] = simplex.logistic(X.next)
  }
  return(out)
}