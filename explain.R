### picture for CA.

T = 200
d = 10
sigma=0.8

N = seq(1, T)
f = 2 * log(10) / T
y = -d * (1-exp(-N * f)) + rnorm(N, sd=sigma)
y1 = d * (1 - exp(-N * f)) + rnorm(N, sd=sigma)

t0 = as.integer(0.1 * T)

y <- c(rep(0, t0), y)
y1 <- c(rep(0, t0), 5 + y1)
par(mar=rep(2.5, 4))

to = as.integer(2 * T/3)

plot(head(N, to), head(y, to), ylim=c(-1.8*d, 1.8*d), type="l",
     xlim=c(0, T))
lines(tail(N, T-to), tail(y, T-to), lty=2, col="red")

lines(tail(head(N, to), to-t0), tail(head(y1, to), to-t0))
lines(tail(N, T-to), tail(y1, T-to), lty=3, col="red")

abline(v=0.1 * T)
# text(40, 15, "assigment Z")
# text(40, -10, "assigment Z=j1")
# text(180, -15, "β_jT(j1)")
# text(28, -15, "t=0")
# text(140, -15, "t=to")
abline(v=to, lty=3)
text(0, -2, "β^(0)")


