library(pbdDMAT, quiet=TRUE)
init.grid()

x <- matrix(rnorm(9),ncol=1)
y <- matrix(rnorm(9),ncol=1)

dx <- as.ddmatrix(x)
dy <- as.ddmatrix(y)

fit <- lm(dy~dx)
fit

summary(fit)

finalize()
