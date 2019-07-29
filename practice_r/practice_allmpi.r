library(pbdDMAT, quiet=TRUE)
library(matrixcalc)

args = commandArgs(trailingOnly = TRUE)

Nnode = as.numeric(args[1])

init.grid()

.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 5

x <- (1:N) +N *.comm.rank
y <- allgather(matrix(x,nrow=1))
comm.print(y)


y <- allreduce(matrix(x, nrow=1))
comm.print(y)

finalize()
