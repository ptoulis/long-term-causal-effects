## Main script
rm(list=ls())
source("games.R")
source("qre.R")
source("qlk.R")
source("simplex.R")

# stop("fix dataset for N=120")

lieberman <- function() {
  ## Reproduces Table III in the QRE paper.
  game = game.lieberman()
  obs.ms =  c(0.167, 0.027, 0.806, 0.013, 0.227, 0.76)  ## at time 11-20
  obs.ms =  c(0.053, 0.047, 0.9, 0.027, 0.053, 0.92) # at time 151-160
  lambda.mle = qre.mle(game, obs.ms, N=300 * 10, nLambdas = 500)
  print(sprintf("MLE for lambda (QRE) = %.3f", lambda.mle))
  print(paste(rep("=", 50), collapse=""))
  # QLk
  best.qlk = qlk.mle(game, obs.ms, N = 300 * 10)
  print(unlist(best.qlk))
}

random.proportion <- function(p) {
  # random proportion of lengthp 
  fix.proportion(rdirichlet(1, alpha=rep(1/p, p)))
}


fix.proportion <- function(prop, tol=1e-2) {
  i = which(prop < tol)
  if(length(i) > 0) {
    prop[i] <- tol
    s = 1 - length(i) * tol
    prop[-i] = prop[-i] * s
  }
  return(as.vector(prop / sum(prop)))
}

# rapoport.mcmc <- function() {
load("rapoport.rda")
game.rap = game.rapoport(W=10, L=-6)
num.actions = ncol(game.rap[[1]])  ## no of actions.
# dataset = rapoport[1:3, ]
dataset = rapoport[5:7, ]

num.periods = nrow(dataset)
long.term = rapoport[4,]
A0t = as.matrix(dataset[, 5:14])  ## action frequencies
N = 120
kMix = 0.9
kProposal.sd = 0.02
# Diffusion for player A
diffusion.A = A0t[, 1:5]
diffusion.B= A0t[, 6:10]

print(sprintf("Rapoport(1992) No of actions %d", num.actions))

## x = state of the chain.
log.lik <- function(x) {
  # Log-likelihood using 
  ## qlk.parameters and diffusion parameters.
  #   qlk = list(At, lambda1, lambda12, lambda2)
  #   psi = list(beta, sigma, mu)
  #   data = A0t matrix action frequencies for both players.
  qlk = x$qlk
  stopifnot(ncol(data)==ncol(A0t))
  ll.qlk = sum(sapply(1:nrow(A0t), function(i) {
    obs.ms = A0t[i,] # observed mixed strategy at round i
    par = list(a0=qlk$At[i, 1], a1=qlk$At[i, 2], a2=qlk$At[i, 3],
               lambda1=qlk$lambda1, lambda12=qlk$lambda12, lambda2=qlk$lambda2)
    qlk.logLikelihood(par, game = game.rap, obs.ms, N)
  }))
  ## QlK with priors.
  
  ll.qlk.prior = dexp(qlk$lambda1, mean=2, sd=1.5, log=T) +
    dnorm(qlk$lambda12, mean=2, sd=1.5,  log=T) + 
    dnorm(qlk$lambda2, mean=2, sd=1.5,  log=T)
  ll.qlk <- ll.qlk + ll.qlk.prior
  psi = x$psi
  # Diffusion with prors.
  # psi$mu = colMeans(qlk$At)
  ll.diffusion = diffusion.logLikelihood(qlk$At, psi) + diffusion.logLikelihood(qlk$At, psi) 
  ll.diffusion.prior = -log(psi$sigma^2)
  
  ll.diffusion <- ll.diffusion + ll.diffusion.prior
  
  ll.qlk + ll.diffusion 
}

## Proposal functions
## 
propose.perturb.one <- function(z) {
  ## used for slight perturbations in the params
  stopifnot(z >= 0)
  w = log(z)
  w.new = rnorm(1, mean=w, sd = kProposal.sd)
  return(exp(w.new))
}

## Computes log.density(x.new -> x.old) - log.density(x.old -> x.new)
log.density.perturb <- function(z.old, z.new) {
  # dnorm(z.new, mean=z.old, sd=kSD, log=T)
  w.old = log(z.old)
  w.new = log(z.new)
  dnorm(w.new, mean=w.old, sd=kProposal.sd, log=T)
}

sample.x0 <- function() {
  # samples starting point for MCMC.
  x = list(qlk=list(At=matrix(1/3, nrow=num.periods,ncol=3), 
                    lambda1=.5, lambda12=.55, lambda2=.5),
           psi=list(beta=runif(1), sigma=0.5, mu=random.proportion(3)))
  return(x)
}

propose.new.x <- function(x.old) {
  # Proposes new state.
  #
  propose.perturb.matrix <- function(A) {
    j = sample(1:nrow(A), size=1)
    u = random.proportion(ncol(A))
    A.new = A
    A.new[j, ] = kMix * A.new[j, ] + (1-kMix) * u  ## new composition
    stopifnot(all(abs(1-rowSums(A.new)) < 1e-5))
    return(A.new)
  }
  
  At.old = x.old$qlk$At
  At.new = propose.perturb.matrix(At.old)
  
  # propose lambdas
  lambda1 = propose.perturb.one(x.old$qlk$lambda1)
  lambda12 = propose.perturb.one(x.old$qlk$lambda12)
  lambda2 = propose.perturb.one(x.old$qlk$lambda2)
  
  # Propose new diffusion param
  sigma = propose.perturb.one(x.old$psi$sigma)
  beta = kMix * x.old$psi$beta + (1-kMix) * rbeta(1, shape1=1, shape=2)
  
  mu = 0.95 * x.old$psi$mu + (0.05) * colMeans(At.new) # there are three behaviors.
 
  
  x.new = x.old
  x.new$qlk$At = At.new
  x.new$qlk$lambda1 = lambda1
  x.new$qlk$lambda12 = lambda12
  x.new$qlk$lambda2 = lambda2
  x.new$psi$beta = beta
  x.new$psi$sigma = sigma
  # x.new$psi$mu = mu
  return(x.new)
}

log.density.proposal.diff <- function(x.new, x.old) {
  # Only lambda1, lambda12, lambda2, sigma, mu are changed:
  #
  # Compute ll(x' -> x) - ll(x -> x')
  ll = 0
  get.diff <- function(parname, el) {
    ll.new.old = log.density.perturb(x.new[[parname]][[el]], x.old[[parname]][[el]]) 
    ll.old.new = log.density.perturb(x.old[[parname]][[el]], x.new[[parname]][[el]]) 
    ll.new.old - ll.old.new
  }
  for(coeff in c("lambda1", "lambda12", "lambda2")) {
    ll <- ll + get.diff("qlk", coeff)
  }
  # sigma.
  ll <- ll + get.diff("psi", "sigma")

  return(ll)
}

mh.ratio.prob <- function(x.old, x.new) {
  dl = log.lik(x.new) - log.lik(x.old) + log.density.proposal.diff(x.new, x.old)
  stopifnot(!is.na(dl))
  return(min(1, exp(dl)))
}

# }

explore.likelihood <- function() { 
  x0 = x
  print(log.lik(x))
  par(mfrow =c(3, 2))
  l = seq(1e-2, 10, length.out=100)
  for(j in c("lambda1", "lambda12", "lambda2")) {
    y = sapply(l, function(i) { x$qlk[[j]] <- i; log.lik(x) })
    plot(l, y, type="l", main=j)
    abline(v=x0$qlk[[j]], col="red")
  }
  for(j in c("sigma", "beta")) {
    if(j=="beta") l = seq(0, 1-1e-3, length.out=100)
    y = sapply(l, function(i) { x$psi[[j]] <- i; log.lik(x) })
    plot(l, y, type="l", main=j)
  }
  
}

check.proposal <- function(niters=100) {
  mcmc_x.old = sample.x0()
  
  diff.ll  = c()
  mh.ratio = c()
  acc = c()
  nacc = 0
  for(i in 1:niters) {
    # print(unlist(mcmc_x.old))
    mcmc_x.new = propose.new.x(mcmc_x.old)
    diff.ll[i] = log.lik(mcmc_x.new) - log.lik(mcmc_x.old)
    mh.ratio[i] = mh.ratio.prob(mcmc_x.old, mcmc_x.new)
    acc[i] = min(1, exp(log.density.proposal.diff(mcmc_x.new, mcmc_x.old)))
    if(runif(1) < mh.ratio[i]) {
      mcmc_x.old = mcmc_x.new
      nacc <- nacc + 1
    }
  }
  par(mfrow=c(3, 1))
  hist(diff.ll)
  hist(mh.ratio)
  hist(acc)
  print(sprintf("Total accepts %.1f%%", 100 * nacc / niters))
}

run.mcmc <- function(niters.mcmc) {
  mcmc_x.old = sample.x0()
  nacc <- 0
  mcmc.out = list()
  mcmc.out[[1]] = mcmc_x.old
  checkpoints = seq(100, niters.mcmc, by=1500)
  
  for(j in 2:niters.mcmc) {
    mcmc_x.new = propose.new.x(mcmc_x.old)
    acc = mh.ratio.prob(mcmc_x.old, mcmc_x.new)
    if(j %in% checkpoints) {
      print(sprintf("(%d/%d) Acceptance %.1f%%", j, niters.mcmc, 100 * nacc / j))
      plot.mcmc(mcmc.out)
    }
    if(runif(1) < acc) {
      # accept move.
      mcmc_x.old <- mcmc_x.new
      nacc = nacc + 1
      # print(sprintf("%d/%d", nacc, niters.mcmc))
    }
    
    mcmc.out[[j]] <- mcmc_x.old
  }
  # plot results.
  plot.mcmc(mcmc.out, add.pop=T)
  return(mcmc.out)
}

plot.mcmc <- function(mcmc.obj, burnin.f=0.6, add.pop=F) {
  niters.mcmc = length(mcmc.obj)
  burnin = niters.mcmc * burnin.f
  ## plotting!
  par(mfrow=c(2,2))
 # xlab = c("iteration")
  for(coeff in c("lambda1", "lambda12")) {#}, "lambda2")) {
    y = unlist(lapply(mcmc.obj, function(x) x$qlk[[coeff]]))[burnin:niters.mcmc]
    plot(y, xlab="iteration", ylab=coeff, type="l")
  }
  
  for(coeff in c("beta")) {#}, "sigma")) {
    y = unlist(lapply(mcmc.obj, function(x) x$psi[[coeff]]))[burnin:niters.mcmc]
    plot(y, ylab="psi0", xlab="iteration", type="l")
  }
  if(add.pop) {
    i = sample(seq(1, niters.mcmc-burnin), size=1500, replace=F)
    
    for(t in 1:3) {
      tcol = rgb(1, 0, 0, alpha = .1)
      if(t==2) tcol = rgb(0, 1, 0, alpha = .1)
      if(t==3) tcol = rgb(0, 0, 1, alpha = .1)
      lims = c(0, .8)
  
      a0 = unlist(lapply(mcmc.obj, function(x) x$qlk$At[t,1]))[burnin:niters.mcmc]
      a1 = unlist(lapply(mcmc.obj, function(x) x$qlk$At[t,2]))[burnin:niters.mcmc]
      a0 = a0[i]
      a1 = a1[i]
    
      if(t==1) plot(a0 - 0.05 * t, a1 + 0.04 * t, xlab="b0", ylab="b1", pch=t+2, xlim=lims, ylim=lims, col=tcol)
      else  points(a0 - 0.05 * t, a1 + 0.04 * t, pch=t+2, xlab="b0", ylab="b1",  xlim=lims, ylim=lims, col=tcol)
    }
  }
  
}

posterior.table <- function(burnin.f=0.5) {
  load(file="mcmc.out1e5.rda")
  niters.mcmc = length(out)
  burnin = niters.mcmc * burnin.f
  group = c("qlk", "qlk", "qlk", "psi", "psi")
  par = c("lambda1", "lambda12", "lambda2", "beta", "sigma")
  
  drop = function(x) print(paste(x, collapse=" & "))
  row = c("parameter", "mean (sd)", "Q1-Q3")
  drop(row)
  
  for(i in 1:5) {
    x = sapply(burnin:niters.mcmc, function(j) as.numeric(out[[j]][[group[i]]][[par[i]]]))
    y = c(par[i], sprintf("%.1f (%.1f)", mean(x), sd(x)), 
          sprintf("[%.1f, %.1f]", round(quantile(x, 0.25),1), round(quantile(x, 0.75),1)))
    drop(y)
  }
}

posterior.inference <- function(burnin.f=0.5) {
  load(file="mcmc.out1e5.rda")
  load(file="mcmc.out1e5-game2.rda")
  ## out has the samples.
  niters.mcmc = length(out)
  burnin = niters.mcmc * burnin.f

  get.game.estimand <- function(game.id, j, verbose=F) {
    x = out[[j]]
    if(game.id==2) {
      x = out2[[j]]
    }
    bt = x$qlk$At[3,]
    mu = colMeans(x$qlk$At)
    Xt = tail(logit(bt), 2)
    mu.X = tail(logit(mu), 2)
    beta = x$psi$beta
    sigma = x$psi$sigma
  #  Xt.next = beta * Xt + mu.X + sigma * rmvnorm(1, mean=, sd=1)
    
    bt.next = logistic(c(0, Xt.next))
    #     print(sprintf("mu=%.2f, Xt=%s, mu.X=%s beta=%.3f sigma=%.2f bt.next=%s",
    #                   mu, paste(round(Xt, 3), collapse=", "), paste(round(mu.X, 3), collapse=","),
    #                   beta, sigma, paste(round(bt.next, 3), collapse=", ")))
    if(verbose) {
      print(sprintf("mu=%s, ....  Xt=%s,  ... mu.X=%s .... beta=%.3f ...... sigma=%.2f ... bt.next=%s",
                    paste(round(mu, 3), collapse=", "), paste(round(Xt, 3), collapse=", "), paste(round(mu.X, 3), collapse=","),
                    beta, sigma, paste(round(bt.next, 3), collapse=", ")))
    }
    aj = predict.mixedStrategy(game.rap, qlk.params = list(a0=bt.next[1], a1=bt.next[2], a2=bt.next[3],
                                                           lambda1=x$qlk$lambda1, 
                                                           lambda12=x$qlk$lambda12, 
                                                           lambda2=x$qlk$lambda2))
    return(aj)
  }
  
  revenue = c()
  for(j in burnin:niters.mcmc) {
    revenue.factor = 1 * c(0, 1, 0, 0, 1, 0, 0, 0, 1, 1)
    a0 = get.game.estimand(1, j)
    a1 = get.game.estimand(2, j)
    revenue = c(revenue, sum(revenue.factor * (a1 - a0)))
    # print(a)
  }
  par(mfrow=c(1, 1))
  hist(revenue)
  return(revenue)
}