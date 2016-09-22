## Code for the paper "Long-term causal effects via behavioral game theory" (NIPS, 2016)
## 
## Panagiotis (Panos) Toulis (U Chicago, Booth) and David Parkes (Harvard, SEAS)
## All rights reserved.

rm(list=ls())
library(compositions)
source("definitions.R")

sample_params_prior <- function(W=1, L=-1) {
  # Samples QLk parameters from their prior distribution.
  # Returns an object of type params (see definitions.R)
  # TODO: Could think about p(ψ, φ| λ, obs_actions) i.e, define the prior for temporal/behavioral model c
  #                                                  conditional on the QLk parameters and the observed actions. 
  Payoff_row <- 1 * matrix(c(W, L, L, L, L,
                             L, L, W, W, W,
                             L, W, L, L, W, 
                             L, W, L, W, L, 
                             L, W, W, L, L), nrow=5, byrow=T)
  Payoff_col <- -t(Payoff_row)
  #
  lambda = runif(kNumBehaviorLevels, min=-5, max=5)
  
  # Could do an Empirical Bayes here.
  phi = as.numeric(rDirichlet.acomp(1, alpha = rep(1, kNumBehaviorLevels)))
  psi = c(runif(2, min=-3, max=3), runif(2))
  
  
  return(list(psi=psi, phi=phi, lambda=lambda, 
              Payoff_ColPlayer=Payoff_col, Payoff_RowPlayer=Payoff_row))
}

population_behavior <- function(behaviorLevel_prevalence, params, is_rowplayer) {
  # Computes population behavior as a mixture of behavior levels given their prevalence.
  #   behaviorLevel_prevalence = 3 x 1 vector of prevalence for the 3 behavior-levels of QLk.
  #   params = QLk parameters.
  #
  # Returns behavior (px1 distribution vector over actions).
  actions_per_behaviorLevels = param_behaviorLevels(params, is_rowplayer = is_rowplayer)
  CHECK_EQ(ncol(actions_per_behaviorLevels), length(behaviorLevel_prevalence))

  population_behavior = actions_per_behaviorLevels %*% matrix(behaviorLevel_prevalence, ncol=1)
  CHECK_NEAR(sum(population_behavior), 1, tol=1e-5)
  population_behavior = population_behavior / sum(population_behavior)
  return(population_behavior)
}

likelihood_obs_actions <- function(obs_action_counts, behaviorLevel_prevalence, params,
                                   is_rowplayer, is_log=F) {
  # Calculates p(a | b, θ) = density of aggregate action α given aggregate behavior β, params θ.
  #   obs_actions = p x 1 aggregate action profile (counts); p=#actions.
  #   behaviorLevel_prevalence = how prevalent a behavior level is (on simplex)
  #   params = parameters of behavioral model (see definitions.R)
  #   is_rowplayer = observed actions from row players or not?
  #
  # Returns a positive scalar for density.
  #
  CHECK_TRUE(all(obs_action_counts >= 0))
  
  population_behavior = population_behavior(behaviorLevel_prevalence, params, is_rowplayer=is_rowplayer)
  # obs_alpha = obs_actions / sum(obs_actions)
  dmultinom(obs_action_counts, prob=population_behavior, log=is_log)
}

sample_actions <- function(num_samples, behaviorLevel_prevalence, params, is_rowplayer, log=F) {
  # Samples from p(a | b, θ), as defined above.
  #
  # Returns a vector of positive numbers that sums up to num_samples.
  #
  population_behavior = population_behavior(behaviorLevel_prevalence, params, is_rowplayer = is_rowplayer)

  # obs_alpha = obs_actions / sum(obs_actions)
  rmultinom(1, size=num_samples, prob=population_behavior)
}

## Code to sample behavior pervalence over time.
sample_B0t <- function(to_t, params, verbose=F) {
  # Samples prevalence of behavior levels in the population given parameters.
  #   to_t = positive integer denoting the long-term horizon we want to extrapolate to (quantity "T" in the paper)
  #   params = QLK parameters (see definitions.R)
  #
  # Returns B = 3 x (T+1) matrix of behavior prevalence such that B_ij = is prevalence of behavior i at time j.
  # Reqs: colSums(B)=1, B>0.
  CHECK_params(params)
  B0t = matrix(NA, nrow=3, ncol=to_t+1)
  B0t[, 1] = as.numeric(rDirichlet.acomp(1, alpha=param_beta0(params)))
  
  psi = param_temporal(params)
  b_infty = c(psi[1], psi[2])
  #
  if(to_t > 0) {
    for(i in seq(2, to_t+1)) {
      # previous prevalence
      b = logitSimplex(B0t[, i-1])
      # noise
      z_t = rnorm(2)
      # next prevalence
      b_next = b_infty + psi[3] * b + psi[4] * rnorm(2)
      #
      if(verbose) {
        print("b = ")
        print(b)
        print("b_next = ")
        print(b_next)
      }
      B0t[, i] = invLogitSimplex(b_next)
    }
  }
  if(verbose) {
    print("b_infty = ")
    print(b_infty)
    print("beta_infty = ")
    print(invLogitSimplex(b_infty / (1-psi[3])))
  }
  return(B0t)
}

expected_longTerm_actions <- function(obs_actions_history, to_t, is_rowplayer, nsamples=1000) {
  W <- rep(0, nrow(obs_actions_history))
  num_agents <- 20
  
  for(k in 1:nsamples) {
    params = sample_params_prior()
    B0t = sample_B0t(to_t, params, verbose=F)
    loglik <- 0
    for(i in 1:ncol(obs_actions_history)) {
      actions_i = as.integer(num_agents * obs_actions_history[, i])
      beta_i = B0t[, i]
      loglik <- loglik + likelihood_obs_actions(actions_i, beta_i, params, is_rowplayer = is_rowplayer, is_log=T)
     # print(sprintf("k=%d, i=%d", k, i))
    #  print(obs_actions_history[, i])
     # print(actions_i)
    #  print(beta_i)
     # print(loglik)
    }
    beta_infty = B0t[, to_t]
    alpha_infty = sample_actions(20, beta_infty, params, is_rowplayer=is_rowplayer)
    
    W <- W + (alpha_infty / sum(alpha_infty)) * exp(loglik)
  }
  M = W/nsamples
  return(M / sum(M))
}

run_experiment <- function() {
  At = t(subset(Rapoport_Data, Game==1, select=c(A1, A2, A3, A4, A5)))
  rs = seq(0, 1, length.out=5)
  # equilibrium <- c(0.375, 0.250, 0.125, 0.125, 0.125) # nash
  equilibrium <- c(0.286, 0.302, 0.138, 0.138, 0.138) # QRE
  med = rowMeans(At)
  
  performance_lte = sapply(rs, function(r) {
    y = r * med + (1-r) * equilibrium
    y_pred = expected_longTerm_actions(At[, 1:4], to_t = 8, is_rowplayer = T, nsamples = 3000)
    # y_pred_naive = At[, 4]
    sum((y_pred-y)**2)
  })
  performance_naive = sapply(rs, function(r) {
    y = r * med + (1-r) * equilibrium
    y_pred = At[, 4]
    # y_pred_naive = At[, 4]
    sum((y_pred-y)**2)
  })
  
  plot(rs, performance_lte, type="l", cex=1.2, ylim=c(0, max(performance_lte)))
  lines(rs, performance_naive, col="red", lty=3)
}




