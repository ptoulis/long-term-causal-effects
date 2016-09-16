## 
rm(list=ls())
rapoport = data.frame(Game = c(1, 1, 1, 1, 2, 2, 2, 2),
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

strategy.list = setdiff(names(rapoport), c("Game", "Session", "Period", "Obs"))
revenue.factor = 1 * c(0, 1, 0, 0, 1, 0, 0, 0, 1, 1)

get.revenue <- function(arg_game, arg_index, rev_factor) {
  arg_session = ifelse(arg_index > 2, 2, 1)
  arg_period = ifelse(arg_index %% 2==0, "late", "early")
  print(sprintf("Revenue for Game= %d, session= %d period= %s", arg_game, arg_session, arg_period))
  x = subset(rapoport, Game==arg_game & Session==arg_session & Period == arg_period,
              select=strategy.list)
  sum(x * rev_factor)
}

diff.revenue <- function(arg_index) {
  get.revenue(2, arg_index, revenue.factor) - get.revenue(1, arg_index, revenue.factor)
}

sweet.spot <- function() {
  # Finds a vector that exemplifies the long-term effect on the revenue.
  f <- function(x, verbose=F) {
    diff.early  = get.revenue(1, 2, "early", x) - get.revenue(2, 2, "early", x) # >0
    diff.late  =  get.revenue(1, 2, "late", x) - get.revenue(2, 2, "late", x)  # should be < 0
    if(verbose) {
      print(sprintf("Early revenue %.3f Late = %.3f", diff.early, diff.late))
    }
    exp(-diff.early * diff.late) + diff.early - 1.5 * diff.late
  }
  k = length(strategy.list)
  fit = optim(par=rep(0, k), fn=function(x) f(x, F), control=list(fnscale=-1),
          method="L-BFGS-B", lower=rep(0, k), upper=rep(1, k))
  
  f(fit$par, verbose=T)
  print(paste(fit$par, collapse=", "))
  return(fit$par)
}

estimand <- function() diff.revenue(4)

neyman.estimate <- function() {
  diff.revenue(3)
}

did.estimate <- function() {
  # Approach (2.) = consider difference-in-differences estimator.
  R1 = get.revenue(1, 3, revenue.factor) - get.revenue(1, 1, revenue.factor)
  R2 = get.revenue(2, 3, revenue.factor) - get.revenue(2, 1, revenue.factor)
  return(R2 - R1)
}

dag.estimate <- function() {
  # Approach (3.) = consider a DAG model for this data?
  
}

po.estimate <- function() {
  # Approach (4.) = consider potential outcomes, and game-theoretic model.
  #
  
}