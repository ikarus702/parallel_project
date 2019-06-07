library(pbdMPI, quietly=TRUE)

.comm.rank <- comm.rank()

N <- 100
x <- split((1:N)+N*.comm.rank,rep(1:10,each=10))
y <- pbdLapply(x, sum, pbd.mode=\"spmd\")
comm.print(unlist(y))

finalize()
