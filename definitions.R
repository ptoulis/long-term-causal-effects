## Code for the paper "Long-term causal effects via behavioral game theory" (NIPS, 2016)
## 
## Panagiotis (Panos) Toulis (U Chicago, Booth) and David Parkes (Harvard, SEAS)
## All rights reserved.

# THIS FILE: definitions. R
#
# Here, we define the objects used in the code and also provide sanity checks.
#
kNumBehaviorLevels <<- 3  # number of behavior levels.
source("../r-toolkit/checks.R")
CHECK_simplex <- function(x) {
  stopifnot(abs(sum(x)-1) < 1e-8 && all(x >= 0))
}

CHECK_params <- function(params) {
  CHECK_SETEQ(names(params), c("lambda", "psi", "phi", "Payoff_RowPlayer", "Payoff_ColPlayer"))
  CHECK_EQ(length(params$lambda), kNumBehaviorLevels, msg = "3 params for behavioral model")
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
param_beta0 <- function(params) {
  params$phi  
}
param_temporal <- function(params) {
  params$psi
}

## Helper functions.
logitSimplex <- function(x) {
  if(is.na(sum(x))) {
    print(x)
  }
  CHECK_NEAR(sum(x), 1, tol = 1e-4)
  x = x / sum(x)
  x_max = max(x)
  n <- length(x)
  if(x[n] < 1e-8) {
    x[n] <- 1e-5
    x <- x / sum(x)
  }
  x_ = head(x, n-1)
  y = log(x_ / x[n])
  return(y)
}

invLogitSimplex <- function(y) {
  x_n <- 1 / (1+sum(exp(y)))
  x <- exp(y) * x_n
  x <- c(x, x_n)
  return(x / sum(x))
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
  
  return(matrix(c(b0, b1, b2), ncol=kNumBehaviorLevels))
}

Rapoport_Data = data.frame(Game = c(1, 1, 1, 1, 2, 2, 2, 2),
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



