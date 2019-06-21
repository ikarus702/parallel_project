suppressPackageStartupMessages( library( pbdDMAT ) )
library(psych)
init.grid()

## Construct on 0 and convert to distributed (communicates)
if ( comm.rank() == 0 ){
  x <- matrix(1:4,2,2)
  y <- matrix(3:6,2,2)
} else {
  x <- NULL
  y <- NULL
}
x.dmat <- as.ddmatrix(x, bldim=c(2,2))
y.dmat <- as.ddmatrix(y, bldim=c(2,2))

a <- 0.5*sum(diag(x.dmat %*% y.dmat))


print(a)
class(a)

z.dmat <- ddmatrix( 1:2,nrow=1,bldim=c(2,2))

z.mat <- as.matrix(z.dmat, proc.dest="all")

#x.dmat
#comm.print( submatrix( x.dmat ), all.rank = TRUE )

## More conveniently (no communication)
#y.dmat <- ddmatrix( 1:100, nrow = 10, ncol = 10, bldim = c( 2, 2 ) )
#comm.print( submatrix( y.dmat ), all.rank = TRUE )

cat("z.dmat is \n")

#z.dmat <- ddmatrix(c(x.dmat,y.dmat), nrow=20,ncol=10, bldim=c(2,2))
comm.print((z.mat), all.rank=TRUE)




## Some ddmatrix defaults

finalize()
