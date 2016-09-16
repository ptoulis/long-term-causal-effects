## New experiments for long-term causal effects.
#
#  Glossary
# -------------------
# b = behavior 
# beta = distribution over behaviors, 
# beta_0t = matrix of betas (column is time.)
#
# beta_free
#
# a = action
# alpha = distribution over actions.
# alpha_0t = matrix of alphas.

rm(list=ls())
dataset = data.frame(Game = c(1, 1, 1, 1, 2, 2, 2, 2),
                     Session = c(1, 1, 2, 2, 1, 1, 2, 2),
                     Period = c("early", "late", "early", "late", "early", "late", "early", "late"),
                     Obs = c(1, 1, 1, 0, 1, 1, 1, 0),
                     A1 = c(0.308, 0.293, 0.273, 0.295, 0.258, 0.290, 0.355, 0.323),
                     A2 = c(0.307, 0.272, 0.350, 0.292, 0.367, 0.347, 0.313, 0.270),
                     A3 = c(0.113, 0.162, 0.103, 0.113, 0.105, 0.118, 0.082, 0.093),
                     A4 = c(0.120, 0.100, 0.123, 0.135, 0.143, 0.110, 0.100, 0.105),
                     A5 = c(0.152, 0.173, 0.151, 0.165, 0.127, 0.135, 0.150, 0.209),
                     B1 = c(0.350, 0.333, 0.353, 0.372, 0.332, 0.355, 0.355, 0.343),
                     B2 = c(0.218, 0.177, 0.133, 0.192, 0.115, 0.198, 0.215, 0.243),
                     B3 = c(0.202, 0.190, 0.258, 0.222, 0.245, 0.208, 0.187, 0.168),
                     B4 = c(0.092, 0.140, 0.102, 0.063, 0.140, 0.108, 0.110, 0.107),
                     B5 = c(0.138, 0.16, 0.154, 0.151, 0.168, 0.131, 0.133, 0.139))

# Parameters
# lambda = QLk model,  phi = Dirichlet prior, psi = VAR model.
source("games.R")
library(mvtnorm)
# library(gtools) # for Dirichlet

logdirichlet <- function(x, alpha) {
  stopifnot(all(alpha > 0))
  stopifnot(all(x > 0) && all(x < 1))
  A1 = sum((alpha - 1) * log(x))
  A2 = sum(log(gamma(alpha))) - log(gamma(sum(alpha)))
  A1 - A2
}

Behaviors = c(1, 2, 3) # QL3 model.
numBehaviors = length(Behaviors) # number of behaviors.
Actions = c(1, 2, 3, 4, 5) # actions.
numActions = length(Actions) # number of actions.
a_unif = rep(1, numActions)/numActions

CHECK_simplex <- function(x) {
  stopifnot(abs(sum(x)-1) < 1e-8 && all(x >= 0))
}

## Helper functions
logistic <- function(v) {
  v[v > 25] <- 25
  v[v < -25] <- -25
  exp(v) / sum(exp(v))
}


alpha_b <- function(lambda, b, Pj) {
  # Get action frequency given behavior and QLk model.
  #
  alpha = rep(0, numActions)
  
  if(b==1) {
    alpha = a_unif
  } else if(b==2){
    # actions of other player
    alpha_other = a_unif
    util = Pj %*% alpha_other
    lam = lambda[1]
    alpha = logistic(util * lam)
  } else {
    lam12 = lambda[2]
    # get the adversarial actions
    alpha_other = alpha_b(c(lam12), b=2, -t(Pj))
    util = Pj %*% alpha_other
    lam = lambda[2]
    alpha = logistic(util * lam)
  }
  CHECK_simplex(alpha)
  return(alpha)
}

alpha_beta <- function(lambda, beta, Pj) {
  # Given a distribution over behaviors, return the distribution over actions.
  #
  # checks.
  beta = as.numeric(beta)
  CHECK_simplex(beta)
  stopifnot(length(beta) == numBehaviors)
  alpha = rep(0, numActions)
  
  for(b in 1:numBehaviors) {
    b_freq = beta[b]
    alpha = alpha + b_freq * alpha_b(lambda, b, Pj)
  }  
  CHECK_simplex(alpha)
  return(alpha)
}

alpha0t_beta0t <- function(lambda, beta0t, Pj) {
  to = ncol(beta0t)
  alpha0t = matrix(0, nrow=numActions, ncol=to)
  
  for(t in 1:to) {
    alpha0t[, t] <- alpha_beta(lambda, beta0t[,t], Pj)
  }
  return(alpha0t)
}

logdensity.alpha0t_beta0t <- function(lambda, alpha0t, beta0t, Pj) {
  expected_alpha0t = alpha0t_beta0t(lambda, beta0t, Pj)
  
  logd = 0
  to = ncol(alpha0t)
  stopifnot(to > 1 || ncol(beta0t) != ncol(alpha0t))
  N = 1000
  
  for(t in 1:to) {
    counts = alpha0t[,t] * N
    prob_t = expected_alpha0t[, t]
    logd <- logd + dmultinom(counts, size = N, prob = prob_t, log = T)
  }
  return(logd)
}

sample.alpha_beta0t <- function(lambda, psi, beta0t, Pj) {
  
  to = ncol(beta0t)
  stopifnot(to > 1 || ncol(beta0t) != ncol(alpha0t))
  beta_t = beta0t[, to]
  
  betaFree.old = betaFree_beta(beta_t)
  mu = matrix(head(psi, 2), ncol=1)

  p = length(betaFree.old)
  I = 0.5 * diag(p)
  
  x0 = (0.2 * betaFree.old + 0.8 * mu) 
 # print(dim(x0))
  et = t(rmvnorm(1, mean=rep(0, p), sigma=I))
#  print(dim(et))
  betaFree.new = x0 + et

  beta_next = beta_betaFree(betaFree.new)
  expected_alpha = alpha_beta(lambda, beta_next, Pj)
  
  N = 1000
  alpha_next = rmultinom(1, size = N, prob = expected_alpha)/N
  return(alpha_next)
}


###      Time-series model.
betaFree_beta <- function(beta) {
  stopifnot(length(beta) == numBehaviors)
  CHECK_simplex(beta)
  
  if(beta[1]==0) {
    beta[1] = 1e-5
    # beta[2:numBehaviors] = (1-beta[1]) * beta[2:numBehaviors]
    beta = beta / sum(beta)
  }
  CHECK_simplex(beta)
  
  # Transform to free parameters.
  betaFree = log(beta[2:numBehaviors] / beta[1])
  return(betaFree)
}

beta_betaFree <- function(betaFree) {
  x = c(0, betaFree)
  z = exp(x) / sum(exp(x))
  # sanitize
  f = 1e-5
  index = which(z < f)
  if(length(index) > 0) {
    z[index] <- f 
    z[-index] <- z[-index] * (1-length(index) * f) / sum(z[-index])
  }
  return(z)
  
}

logdensity.VAR <- function(psi, beta.new, beta.old) {
  #
  # beta.new = psi[1, 2] + psi[3] * beta.old + psi[4] * e
  #
  stopifnot(length(psi)==4)
  #stopifnot(psi[3] != 0)
  if(psi[4]==0) return(-1e10)
  
  betaFree.new = betaFree_beta(beta.new)
  betaFree.old = betaFree_beta(beta.old)
  mu = matrix(head(psi, 2), ncol=1)
  # residual
  # not going to matter.
  # TODO(ptoulis): Is this ok?
  # psi[4] = .1
  # et = as.numeric((betaFree.new - psi[3] * betaFree.old - mu) / psi[4])
  et = as.numeric(betaFree.new - (0.2 * betaFree.old + 0.8 * mu))
  p = length(betaFree.new)
  I = 0.5 * diag(p)
  # print(et)
  lt = dmvnorm(et, mean=rep(0, p), sigma=I, log=T)
  # print(lt)
  return(lt)
}

logdensity.beta0t <- function(psi, phi, beta0t) {
  to = ncol(beta0t)
  stopifnot(to > 1)
  beta_0 = beta0t[, 1]
  
  # Prior
  logd = logdirichlet(beta_0, phi)
  
  # Iterate from t=2 to t=to
  if(to > 1) {
    for(t in 2:to) {
      beta.new = beta0t[, t]
      beta.old = beta0t[, t-1]
      # print(sprintf("t=%d logD(%d)=%.2f", t, t-1, logd))
      logd <- logd + logdensity.VAR(psi, beta.new, beta.old)
    } 
  }
  return(logd)
}

viterbi <- function(alpha0t, Payoff) {
  # Finds marginal-likelihood MLE
  #
  #
  # par = (2 * to, lambda>0, psi, phi)
  to = ncol(alpha0t)
  numFlatBeta = (numBehaviors-1) * to
  
  flat_beta0t <- function(beta0t) {
    # transform beta0t to free parameters.
    #
    par = c()
    for(t in 1:ncol(beta0t)) {
      betaFree_t = betaFree_beta(beta = beta0t[,t])
      par = c(par, as.numeric(betaFree_t))
    }
    return(par)
  }
  
  beta0t_flat <- function(beta0t_flat) {
    # transform from free parameter
    #
    stopifnot(to == length(beta0t_flat)/2)
    beta0t = matrix(0, ncol=to, nrow=numBehaviors)
    for(t in 1:to) {
      index = c(2*t-1, 2*t)
      beta0t[, t] = beta_betaFree(betaFree = beta0t_flat[index])
    }
    return(beta0t)
  } 
  
  translate.par <- function(x) {
    beta0t = beta0t_flat(head(x, numFlatBeta))
    lambda = x[seq(numFlatBeta+1, numFlatBeta+3)]
    psi = x[seq(numFlatBeta+4, numFlatBeta+7)]
    phi = x[seq(numFlatBeta+8, numFlatBeta+10)]
    return(list(beta0t=beta0t, lambda=lambda, psi=psi, phi=phi))
  }
  
  objective <- function(par) {
    params = translate.par(par)
    # log-likelihoods.
    loglik_a_b = logdensity.alpha0t_beta0t(lambda = params$lambda, 
                                           alpha0t = alpha0t, 
                                           beta0t = params$beta0t, 
                                           Pj=Payoff)
    
    loglik_b = logdensity.beta0t(params$psi, params$phi, params$beta0t)
    # Shrink to equilbrium
    alpha.Eq = c(0.375, 0.25, rep(0.125, 3))  # By 
    beta.Eq = beta_betaFree(betaFree = head(params$psi, 2))
    alpha.Eq_pred = alpha_beta(lambda = params$lambda, beta = beta.Eq, Pj = Payoff)
    penalty = 100 * sum((alpha.Eq - alpha.Eq_pred)**2)

    
    print(penalty)
    print(alpha.Eq_pred)
    print(alpha.Eq)
    print("loglik a_b")
    print(loglik_a_b)
    print("loglik b")
    print(loglik_b)
    #
    stopifnot(is.numeric(loglik_a_b) || is.numeric(loglik_b))
    return(-loglik_a_b - loglik_b + penalty)
  }
  
  ## Optimization
  print(numFlatBeta)
  par0 = c(rep(0, numFlatBeta), rep(0, 3), rep(1, 4), rep(1, 3))
  x = optim(par = par0, fn = objective, method = "L-BFGS-B",
            lower=c(rep(-100, numFlatBeta), rep(1e-3, 3), rep(-10, 4), rep(1e-3, 3)),
             upper=c(rep(100, numFlatBeta), rep(20, 3), rep(10, 4), rep(1, 3)))$par
  print(-objective(x))
  print(-objective(par0))
  
  
  return(translate.par(x))
}

viterbi.both <- function(alpha0t_A, Payoff_A,
                         alpha0t_B, Payoff_B) {
  # Finds marginal-likelihood MLE
  #
  #
  # par = (2 * to, lambda>0, psi, phi)
  to = ncol(alpha0t_A)
  numFlatBeta = (numBehaviors-1) * to
  
  flat_beta0t <- function(beta0t) {
    # transform beta0t to free parameters.
    #
    par = c()
    for(t in 1:ncol(beta0t)) {
      betaFree_t = betaFree_beta(beta = beta0t[,t])
      par = c(par, as.numeric(betaFree_t))
    }
    return(par)
  }
  
  beta0t_flat <- function(beta0t_flat) {
    # transform from free parameter
    #
    stopifnot(to == length(beta0t_flat)/2)
    beta0t = matrix(0, ncol=to, nrow=numBehaviors)
    for(t in 1:to) {
      index = c(2*t-1, 2*t)
      beta0t[, t] = beta_betaFree(betaFree = beta0t_flat[index])
    }
    return(beta0t)
  } 
  
  translate.par <- function(x) {
    beta0t_A = beta0t_flat(x[1:numFlatBeta])
    beta0t_B = beta0t_flat(x[seq(numFlatBeta + 1, 2 * numFlatBeta)])
    
    offset = 2 * numFlatBeta
    lambda = x[seq(offset + 1, offset+3)]
    psi = x[seq(offset + 4, offset + 7)]
    phi = x[seq(offset + 8, offset + 10)]
    return(list(beta0t_A=beta0t_A, beta0t_B=beta0t_B, 
                lambda=lambda, psi=psi, phi=phi))
  }
  
  objective <- function(par, player) {
    params = translate.par(par)
    
    alpha0t = NA
    if(player=="A") {
      alpha0t = alpha0t_A
      beta0t = params$beta0t_A
      P = Payoff_A
    } else {
      alpha0t = alpha0t_B
      beta0t = params$beta0t_B
      P = Payoff_B
    }
  
    # log-likelihoods.
    loglik_a_b = logdensity.alpha0t_beta0t(lambda = params$lambda, 
                                           alpha0t = alpha0t, 
                                           beta0t = beta0t, 
                                           Pj=P)
    
    loglik_b = logdensity.beta0t(params$psi, params$phi, beta0t)
    # Shrink to equilbrium
    alpha.Eq = c(0.375, 0.25, rep(0.125, 3))  # By 
    beta.Eq = beta_betaFree(betaFree = head(params$psi, 2))
    alpha.Eq_pred = alpha_beta(lambda = params$lambda, beta = beta.Eq, Pj = P)
    
    shrinkage_one = 1e4 * sum((alpha.Eq - alpha.Eq_pred)**2) 
    shrinkage_two = 5  * sum((params$lambda - c(0.5, 0.2, 1.9))**2)
    shrinkage_three = -1 * sum(diff(params$beta0t[1, ]))
    penalty = shrinkage_one + shrinkage_two + shrinkage_three
    
        
#       print(shrinkage_one)  
#     print("shrinkage from lambda")
#     print(params$lambda)
#     print(shrinkage_two)
#     
#       print(shrinkage_three)
#         print("penalty")
#         print(penalty)
#         # print(alpha.Eq_pred)
#         # print(alpha.Eq)
#         print("loglik a_b")
#         print(loglik_a_b)
#         print("loglik b")
#         print(loglik_b)
    
    stopifnot(is.numeric(loglik_a_b) || is.numeric(loglik_b))
    return(-loglik_a_b - loglik_b + penalty)
  }
  
  ## Optimization
  par0 = c(rep(0, numFlatBeta), rep(0, numFlatBeta),
           rep(0, 3), rep(1, 4), rep(1, 3))
  
  total.objective <- function(par) {
    objective(par, "A") + objective(par, "B")
  }
  
  x = optim(par = par0, fn = total.objective, method = "L-BFGS-B",
            lower=c(rep(-20, 2 * numFlatBeta), rep(1e-3, 3), rep(-10, 4), rep(1e-3, 3)),
            upper=c(rep(20, 2 * numFlatBeta), rep(20, 3), rep(10, 4), rep(1, 3)))$par
  print(-total.objective(x))
  print(-total.objective(par0))
  
  
  return(translate.par(x))
}


fit.rapoport <- function() {
  alpha0t_rowPlayer = t(dataset[1:3, 5:9])
  alpha_A_true =  t(dataset[4, 5:9])
  # alpha_A_2_true =  t(dataset[4, 5:9])
  
  alpha0t_colPlayer =  t(dataset[1:3, 10:14])
  alpha_B_true =  t(dataset[4, 10:14])
  D = cbind(alpha0t_rowPlayer, alpha0t_colPlayer)
  colnames(D) <- c()
  
  Payoff_A = game.rapoport(1, -1)$A
  Payoff_B = game.rapoport(1, -1)$B
  
  params <- NA
  if(file.exists("model.rda")) {
    print("Loading model")
    print(params)
    load("model.rda")
  } else {
    print("Fitting model.")
    params = viterbi.both(alpha0t_rowPlayer, Payoff_A, alpha0t_colPlayer, Payoff_B)
    save(params, file="model.rda")
    print("Model saved.")
  }
  
  ##  Estimation process.
  taus = c()
  dids = c()
  naives = c()
  longs = c()
  
  for(i in 1:25) {
    # coefficient
    cR = runif(10, min=0, max=1)
    
    get.estimand <- function(alpha_A, alpha_B) {
      sum(cR * c(alpha_A, alpha_B))
    }
    alpha.Eq = c(0.375, 0.25, rep(0.125, 3))
    tau = get.estimand(alpha.Eq, alpha.Eq)
    did = get.estimand(alpha0t_rowPlayer[,3], alpha0t_colPlayer[,3])
    naive = get.estimand(alpha0t_rowPlayer[,1], alpha0t_rowPlayer[,1])
    taus = c(taus, tau)
    dids = c(dids, did)
    naives = c(naives, naive)
    # our method
    longterm = c()
    for(j in 1:100) {
      alpha_A_sample = sample.alpha_beta0t(params$lambda, params$psi, beta0t = params$beta0t_A, Pj = Payoff_A)
      alpha_B_sample = sample.alpha_beta0t(params$lambda, params$psi, beta0t = params$beta0t_B, Pj = Payoff_B)  
      # print(alpha_A_sample)
      longterm = c(longterm, get.estimand(alpha_A_sample, alpha_B_sample))
    }
    longs = c(longs, median(longterm))
  }
  plot(taus, taus, type="l", xlab =expression(paste("ground truth  (", tau, ") ")), 
       ylab=expression(paste("ground truth  (", tau, ") ")), col="white", xlim=c(0.5, 1.5), ylim=c(0.5, 1.5))
  sub = seq(1, length(taus), length.out=10)
  z = seq(0.5, 1.5, length.out=100)
  lines(z, z, lty=3, col="black")
  points(taus, jitter(dids, amount=.0), pch=15, col="red")
  print(summary(longs))
  points(taus, longs, pch=17, col="blue")
  points(taus, jitter(naives, amount=.0), pch=16, col="orange")
  
  print(sqrt(sum(dids-taus)**2))
  print(sqrt(sum(longs-taus)**2))
  print(sqrt(sum(naives-taus)**2))
  legend(1.0, 0.8, legend=c("naive", "DID", "LACE     "),
         pch=c(16, 15, 17), col=c("orange", "red", "blue"))
}

## SOME TESTS.
test.viterbi <- function() {
  alpha0t = matrix(1/numActions, nrow=numActions, ncol=5)
  Payoff = game.rapoport(1, -1)$A
  print(alpha0t)
  x = viterbi(alpha0t, Payoff)
  print(x)
}


plot.behaviors <- function() {
  load("model.rda")
  b0 = params$beta0t_A[1, ]
  print(b0)
}

