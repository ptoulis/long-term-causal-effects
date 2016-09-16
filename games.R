# GAME = list(payoff matrix for agent A, payoff for agent B, ...)
# mixedStrategy = (piA, piB, ...) where piX = action distribution for agent X.
#
game.mckelvey <- function(A=5, B=5) {
  Payoff.A <- matrix(c(1, 0, 1,
                       0, 0, 0, 
                       1, A, 1), nrow=3, byrow=T)
  
  Payoff.B <- matrix(c(1, 0, 1,
                       0, 0, 0, 
                       1, B, 1), nrow=3, byrow=T)
  
  return(list(Payoff.A, Payoff.B))
}
game.lieberman <- function() {
  Payoff.A <- 3.373 * matrix(c(15, 0, -2,
                               0, -15, -1, 
                               1, 2, 0), nrow=3, byrow=T)
  
  Payoff.B <-  3.373 * matrix(c(-15, 0, -1,
                                0, 15, -2, 
                                2, 1, 0), nrow=3, byrow=T)
  return(list(A=Payoff.A, B=Payoff.B))
}

game.rapoport <- function(W, L) {
  Payoff.A <- 1 * matrix(c(W, L, L, L, L,
                               L, L, W, W, W,
                               L, W, L, L, W, 
                               L, W, L, W, L, 
                               L, W, W, L, L), nrow=5, byrow=T)
  Payoff.B <- -t(Payoff.A)
  return(list(A=Payoff.A, B=Payoff.B))
}

get.player.index <- function(mixed.strategy, game, agentId) {
  num.actions = nrow(game[[1]])
  start = 1 + (agentId-1) * num.actions
  end = start + num.actions - 1
  stopifnot(end <= length(mixed.strategy), start > 0)
  return(list(start=start, end=end))
}

player.mixedStrategy <- function(mixed.strategy, game, agentId) {
  # Gets the player strategy from the entire vector of mixed strategy profile
  #
  # Returns:
  #   distribution over actions of player agentId
  ind = get.player.index(mixed.strategy, game, agentId)
  st = matrix(mixed.strategy[ind$start:ind$end], ncol=1)
  stopifnot(abs(1-sum(st)) < 1e-3)
  return(st)
}

other.mixedStrategy <- function(mixed.strategy, game, agentId) {
  # 
  # Returns the mixed strategies of all other agents flattened out 
  # in a column vector.
  ind =  get.player.index(mixed.strategy, game, agentId)
  return(matrix(mixed.strategy[-seq(ind$start, ind$end)], ncol=1))
}