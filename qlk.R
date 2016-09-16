## QLk model. 
# qlk.params = list(a0, a1, a2, lambda1, lambda12, lambda2)
# B = LIST(1=(a0, a1, a2), 2=(b0, b1, b2))
# lambda1 = precision of level-1 agents
# lambda12 = perceived precision of level-1 by level-2
# lambda2 = precision of level-2 agents.

check.qlkParams <- function(params) {
  stopifnot(abs(1 - sum(params$B[[1]])) < 1e-3)
  stopifnot(abs(1 - sum(params$B[[2]])) < 1e-3)
  stopifnot(with(params, lambda1 >= 0, lambda12 >= 0, lambda2 >= 0))
}

logistic <- function(v) {
  v[v > 100] <- 100
  v[v < -100] <- -100
  exp(v) / sum(exp(v))
}

correct.counts <- function(Na, N) {
  # Correct for slight inconsistencies.
  excess = N - sum(Na)
  if(excess < 0 & excess > -3) {
    warning("Correcting counts")
    excess = -excess
    for(j in 1:excess) {
      i = which(Na > 0)
      k = sample(1:length(i), size=1)
      Na[i[k]] <- Na[i[k]] - 1
    }
    return(Na)
  }
  if(excess < -3 | excess > length(Na)) {
    stop("Error in observed frequencies. sum(Na) > total trials.")
  } else if(excess != 0) {
    # distribute the excess counts randomly
    warning("Correcting counts. +")
    e = sample(c(rep(1, excess), rep(0, length(Na)-excess)))
    Na <- Na + e
  }
  return(Na)
}

predict.mixedStrategy <- function(game, qlk.params) {
  # Given a game and the parameters of the QLk model
  # it returns the predicted mixed strategy
  #
  # Args:
  #   game = a game object (see games.R)
  #   qlk.params = LIST of parameters for QLk model (see check.* function)
  #
  # Returns:
  #   (piA, piB, ...)   as vector
  # where piA = distribution over action space of A
  #
  check.qlkParams(qlk.params)
  nActions = nrow(game[[1]])
  num.agents = length(game)
  new.mixedStrategy = c()
  
  other.player = function(a) {
    if(a==1) return(2)
    if(a==2) return(1)
    stop("not correct agent id")
  }
  
  for(a in 1:num.agents) {
    # level-0 : random actions.
    pi0 = rep(1 / nActions, nActions)
    # level-1
    u1.bar = game[[a]] %*% matrix(pi0, ncol=1)
    pi1 = as.vector(logistic(qlk.params$lambda1 * u1.bar))
    
    # level-2
    u12.bar = game[[other.player(a)]] %*% matrix(pi0, ncol=1)
    pi12 = logistic(qlk.params$lambda12 * u12.bar)
    u2.bar = game[[a]] %*% matrix(pi12, ncol=1)
    pi2 = as.vector(logistic(qlk.params$lambda2 * u2.bar))
    
    behavior.mix = qlk.params$B[[a]]
    # print(behavior.mix)
    # print(pi0)
    # print(pi1)
    # print(pi2)
    new.mixedStrategy = c(new.mixedStrategy,
                          behavior.mix[1] * pi0 + behavior.mix[2] * pi1 + behavior.mix[3] * pi2)
    stopifnot(all(!is.na(new.mixedStrategy)))
  }
  
  return(new.mixedStrategy)
}

qlk.logLikelihood <- function(qlk.params, game, obs.mixedStrategy, N) {
  # Computes the log-likelihood for the QLk model
  #
  # Args:
  #   game = GAME object (see games.R)
  #   obs.actionFrequency = vector (piA, piB, ...) where piA = distribution for actions of A
  #   N = total #observations that the freq corresponds to.
  #   qlk.params = parameters for the QLk model. 
  # 
  num.agents = length(game)
  pi.star = predict.mixedStrategy(game, qlk.params)
  ll = 0
  for(a  in 1:num.agents) {
    # 1. Get observed mixed strategy for player a.
    pi.agent.obs = player.mixedStrategy(obs.mixedStrategy, game, agentId = a)
    # 2. Get observed counts for each action.
    Na = round(N * pi.agent.obs)
    Na = correct.counts(Na, N)
    # 3. Get equilibrium mixed strategy for agent a
    pi.agent.star = player.mixedStrategy(pi.star, game, a)
    # 4. Compute log-likelihood 
    ll <- ll + dmultinom(Na, size=N, prob=pi.agent.star, log=T)
  }
  return(ll)
}

random.qlkParams <- function(max.lambda=10) {
  a = exp(runif(3))
  a = a / sum(a)
  lam = runif(3, min=0, max=max.lambda)
  return(list(a0=a[1], a1=a[2], a2=a[3],
              lambda1=lam[1], lambda12=lam[2], lambda2=lam[3]))
}

qlk.mle = function(game, obs.mixedStrategy, N) {
  print("Fixing errors..")
  for(a in length(game)) {
    pi.agent.obs = player.mixedStrategy(obs.mixedStrategy, game, agentId = a)
    # 2. Get observed counts for each action.
    Na = round(N * pi.agent.obs)
    excess = N - sum(Na)
    if(excess < 0 | excess > length(Na)) {
      stop("Error in observed frequencies. sum(Na) > total trials.")
    } else if(excess != 0) {
      # distribute the excess counts randomly
      e = sample(c(rep(1, excess), rep(0, length(Na)-excess)))
      Na <- Na + e
      warning(sprintf("Adding %d obs", excess))
      pi.agent.obs = Na / N
      ind = get.player.index(obs.mixedStrategy, game, a)
      obs.mixedStrategy[ind$start:ind$end] <- pi.agent.obs
    }
  }
  par.to.qlk <- function(par) {
    S = sum(exp(par[1:3]))
    a0 = exp(par[1]) / S
    a1 = exp(par[2]) / S
    a2 = exp(par[3]) / S
    list(a0 = a0, a1=a1, a2=a2,
         lambda1=par[4], lambda12=par[5], lambda2=par[6])
  }
  
  f <- function(param) {
    # waiting for 6 params.
    qlk.par = par.to.qlk(param)
    qlk.logLikelihood(qlk.par, game, obs.mixedStrategy, N)
  }
  
  fit = optim(par=rep(0, 6), fn = f, control=list(fnscale=-1),
              method="L-BFGS-B", lower=c(-3, -3, -3, 0, 0, 0),
              upper=rep(2, 6))
  
  best.par = par.to.qlk(fit$par)
  print(sprintf("Best log-likelihood for QLk = %.3f", fit$value))
  print("Predicted mixed strategy under QLk")
  print(round(predict.mixedStrategy(game, best.par), 3))
  print("Actual observed")
  print(round(obs.mixedStrategy, 3))
  return(best.par)
}
