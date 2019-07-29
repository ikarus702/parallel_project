library(pbdDMAT, quiet=TRUE)
init.grid()

dx <- ddmatrix(1:30, nrow=10)

x <- as.block(dx)
x
comm.print(submatrix(x))

x <- as.rowblock(dx)
x
comm.print(submatrix(x))


x <- as.colblock(dx)
x
comm.print(submatrix(x))


x <- as.rowcyclic(dx)
x 
comm.print(submatrix(x))


x <- as.colcyclic(dx)
x
comm.print(submatrix(x))


x <- as.blockcyclic(dx)
x
comm.print(submatrix(x))


finalize()
