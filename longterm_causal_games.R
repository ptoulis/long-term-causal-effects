## Code for the paper "Long-term causal effects via behavioral game theory" (NIPS, 2016)
## 
## Panagiotis (Panos) Toulis (U Chicago, Booth) and David Parkes (Harvard, SEAS)
## All rights reserved.

rm(list=ls())
library(compositions)
source("definitions.R")

sample_params <- function(W=1, L=-1) {
  # Samples QLk parameters.
  # Returns an object of type params (see definitions.R)
  # 
  Payoff_row <- 1 * matrix(c(W, L, L, L, L,
                             L, L, W, W, W,
                             L, W, L, L, W, 
                             L, W, L, W, L, 
                             L, W, W, L, L), nrow=5, byrow=T)
  Payoff_col <- -t(Payoff_row)
  #
  lambda = runif(3, min=-5, max=5)
  phi = as.numeric(rDirichlet.acomp(1, alpha = rep(.1, nrow(Payoff_row))))
  psi = runif(3, min=-10, max=10)
  
  return(list(psi=psi, phi=phi, lambda=lambda, 
              Payoff_ColPlayer=Payoff_col, Payoff_RowPlayer=Payoff_row))
}

population_behavior <- function(behaviorLevel_prevalence, params) {
  # Computes population behavior as a mixture of behavior levels given their prevalence.
  #   behaviorLevel_prevalence = 3 x 1 vector of prevalence for the 3 behavior-levels of QLk.
  #   params = QLk parameters.
  #
  # Returns behavior (px1 distribution vector over actions).
  behaviorLevels = param_behaviorLevels(params, is_rowplayer = is_rowplayer)
  CHECK_EQ(ncol(behaviorLevels), length(behaviorLevel_prevalence))
  population_behavior = behaviorLevels %*% matrix(behaviorLevel_prevalence, ncol=1)
  CHECK_EQ(sum(population_behavior), 1)
  return(population_behavior)
}

likelihood_obs_actions <- function(obs_actions, behaviorLevel_prevalence, params,
                                   is_rowplayer, log=F) {
  # Calculates p(a | b, θ) = density of aggregate action α given aggregate behavior β, params θ.
  #   obs_actions = p x 1 aggregate action profile (counts); p=#actions.
  #   behaviorLevel_prevalence = how prevalent a behavior level is (on simplex)
  #   params = parameters of behavioral model (see definitions.R)
  #   is_rowplayer = observed actions from row players or not?
  #
  # Returns a positive scalar for density.
  #
  CHECK_TRUE(all(obs_actions >= 0))
  CHECK_EQ(nrow(behaviorLevels), length(obs_actions))
  
  population_behavior = population_behavior(behaviorLevel_prevalence, params)
  # obs_alpha = obs_actions / sum(obs_actions)
  dmultinom(obs_actions, prob=population_behavior, log=log)
}

sample_obs_actions <- function(num_samples, behaviorLevel_prevalence, params, is_rowplayer, log=F) {
  # Samples from p(a | b, θ), as defined above.
  #
  # Returns a vector of positive numbers that sums up to n.
  #
  population_behavior = population_behavior(behaviorLevel_prevalence, params)

  # obs_alpha = obs_actions / sum(obs_actions)
  rmultinom(1, size=num_samples, prob=population_behavior)
}