## Code for the paper "Long-term causal effects via behavioral game theory" (NIPS, 2016)
## 
## Panagiotis (Panos) Toulis (U Chicago, Booth) and David Parkes (Harvard, SEAS)
## All rights reserved.

# THIS FILE: definitions. R
#
# Here, we define the objects used in the code and also provide sanity checks.

source("../r-toolkit/checks.R")
CHECK_simplex <- function(x) {
  stopifnot(abs(sum(x)-1) < 1e-8 && all(x >= 0))
}

CHECK_params <- function(params) {
  CHECK_SETEQ(names(params), c("lambda", "psi", "phi", "Payoff_RowPlayer", "Payoff_ColPlayer"))
  CHECK_EQ(length(params$lambda), 3, msg = "3 params for behavioral model")
}

expit <- function(x) {
  S = sum(exp(x))
  exp(x) / S
}

param_rowplayer_Payoff <- function(params) {
  return(params$Payoff_RowPlayer)
}
param_colplayer_Payoff <- function(params) {
  return(params$Payoff_ColPlayer)
}
param_numactions <- function(params) {
  nrow(param_rowplayer_Payoff(params))
}
param_lam1 <- function(params) {
  params$lambda[1]  
}
param_lam1_2 <- function(params) {
  params$lambda[2]
}
param_lam2 <- function(params) {
  params$lambda[3]
}

param_behaviorLevels <- function(params, is_rowplayer) {
  # Calculate the behaviors for this particular QLk parameter.
  #
  # Returns p x 3 matrix of the form  (p=#actions, 3=total behavior levels in QLk)
  #     [b0, b1, b2] where b0 = random behavior, 
  #                        b1 = level-1 playing against b0 with precision lam1
  #                   and  b2 = level-2 playing with precision lam2 against level-1 players with precision lam1_2.
  #
  CHECK_params(params)
  nActions = param_numactions(params)
  # level-0 players
  b0 = rep(1, nActions) / nActions
  
  Payoff = param_rowplayer_Payoff(params)
  if(!is_rowplayer) {
    Payoff = param_colplayer_Payoff(params)
  }
  # level-1 players
  b1 = expit(param_lam1(params) * Payoff %*% b0)
  # level-2 players
  b1_2 = expit(param_lam1_2(params) * Payoff %*% b0)
  b2 = expit(param_lam2(params) * Payoff %*% b1_2)
  
  return(matrix(c(b0, b1, b2), ncol=3))
}





