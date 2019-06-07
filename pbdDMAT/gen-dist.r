suppressPackageStartupMessages( library( pbdDMAT ) )
init.grid()

## Construct on 0 and convert to distributed (communicates)
if ( comm.rank() == 0 ){
  x <- matrix( 1:100, nrow = 10, ncol = 10 )
} else {
  x <- NULL
}
x.dmat <- as.ddmatrix( x, bldim = c( 2, 2 ) )

x.dmat
comm.print( submatrix( x.dmat ), all.rank = TRUE )

## More conveniently (no communication)
y.dmat <- ddmatrix( 1:100, nrow = 10, ncol = 10, bldim = c( 2, 2 ) )
comm.print( submatrix( y.dmat ), all.rank = TRUE )

all.equal( x.dmat, y.dmat )

## Some ddmatrix defaults
comm.cat(".pbd_env$BLDIM =", .pbd_env$BLDIM, "\n")
comm.cat(".pbd_env$ICTXT =", .pbd_env$ICTXT, "\n")

finalize()
