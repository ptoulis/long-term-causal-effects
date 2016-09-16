##
logistic <- function(v) {
  v[v > 100] <- 100
  v[v < -100] <- -100
  exp(v) / sum(exp(v))
}

new.mixedStrategy <- function(game, old.mixedStrategy, lambda, verbose=F) {
  # Expected payoff for player 1.
  num.agents = length(game)
  new.mixedStrategy = c()
  for(agent in 1:num.agents) {
    p.other = other.mixedStrategy(old.mixedStrategy, game, agentId=agent)
    u.bar = game[[agent]] %*% p.other
    new.mixedStrategy = c(new.mixedStrategy, logistic(lambda * u.bar))
    stopifnot(all(!is.na(new.mixedStrategy)))
    # verbose
    if(verbose) {
      print(sprintf("Expected utilities for Player %d", agent))
      print(u.bar)
    }
  }

  return(new.mixedStrategy)
}

random.mixedStrategy <- function(game) {
  x = c()
  nplayers = length(game)
  nactions = nrow(game[[1]])
  for(j in 1:nplayers) {
    u = runif(nactions)
    u = exp(u) / sum(exp(u))
    x = c(x, u)
  }
  return(x)
}

solve.QRE = function(game, lambda, verbose=F) {
  num.players = length(game)
  ms.old = random.mixedStrategy(game)
  if(verbose) {
    print("Initial mixed strategy")
    print(ms.old)
  }
  
  diff = c()
  for(j in 1:2000) {
    ms.new = new.mixedStrategy(game, ms.old, lambda)
    # difference in two vectors.
    diff[j] = mean(abs(ms.new - ms.old))
    if(tail(diff[j], 1) < 1e-8) {
      if(verbose) {
        plot(diff, type="l")
        print("Info about convergence")
        print(tail(diff, 10))
      }
      return(ms.old)
    } else {
      ms.old = ms.new
    }
  }
  
  warning("Exhausted iters.")
  return(ms.old)
}

qre.logLikelihood <- function(lambda, game, obs.mixedStrategy, N) {
  # Computes the log-likelihood for the QRE model
  #
  # Args:
  #   game = GAME object (see games.R)
  #   obs.actionFrequency = vector (piA, piB, ...) where piA = distribution for actions of A
  #   N = total #observations that the freq corresponds to.
  #   lambda = parameter lambda of the QRE model.
  # 
  num.agents = length(game)
  num.actions = nrow(game[[1]])
  pi.star = solve.QRE(game, lambda)
  
  ll = 0
  for(a  in 1:num.agents) {
    # 1. Get observed mixed strategy for player a.
    pi.agent = player.mixedStrategy(obs.mixedStrategy, game, agentId = a)
    # 2. Get observed counts for each action.
    Na = round(N * pi.agent)
    # Correct for slight inconsistencies.
    excess = N - sum(Na)
    if(excess < 0 | excess > length(Na)) {
      stop("Error in observed frequencies. sum(Na) > total trials.")
    } else if(excess != 0) {
      # distribute the excess counts randomly
      e = sample(c(rep(1, excess), rep(0, length(Na)-excess)))
      Na <- Na + e
    }
    # 3. Get equilibrium mixed strategy for agent a
    pa.star = player.mixedStrategy(pi.star, game, a)
    # 4. Compute log-likelihood 
    ll <- ll + dmultinom(Na, size=N, prob=pa.star, log=T)
  }
  return(ll)
}

qre.mle <- function(game, obs.mixedStrategy, N, nLambdas=100) {
  # Computes the MLE for the QRE model.
  lambdas = seq(0, 2, length.out=nLambdas)
  y = sapply(lambdas, function(lam) qre.logLikelihood(lambda=lam, game, obs.mixedStrategy, N))
  plot(lambdas, y, type="l")
  
  best.lambda = lambdas[which.max(y)]
  print(sprintf("Best log-likelihood QRE = %.3f", 
                qre.logLikelihood(best.lambda, game, obs.mixedStrategy, N)))
  print("Predicted mixed strategy under QRE")
  print(round(solve.QRE(game, best.lambda), 3))
  print("Actual observed mixed strategy")
  print(round(obs.mixedStrategy, 3))
  
  
  return(best.lambda)
}
##  Specific examples.
##
reproduce.Figure1 <- function() {
  ## Reproduce Figure 1 from QRE paper.
  game = game.mckelvey()
  lambdas = c(0, exp(seq(-6, log(80), length.out=100)))
  DR.line = c()
  UL.line = c()
  M.line = c()
  for(lam in lambdas) {
    p = solve.QRE(game, lam)
    stopifnot(abs(p[1]-p[4]) < 1e-3, abs(p[3] - p[6]) < 1e-3)
    DR.line = c(DR.line, p[3])
    UL.line = c(UL.line, p[1])
    M.line = c(M.line, p[2])
  }
  plot(DR.line, xaxt="n", type="l", main="QRE", ylim=c(0, 1))
  lam.index = as.integer(seq(1, length(lambdas), length.out=10))
  axis(1, at=lam.index, 
       tick=T, labels=round(lambdas[lam.index], 2))
  lines( UL.line, lty=2)
  lines( M.line, lty=4)
}
