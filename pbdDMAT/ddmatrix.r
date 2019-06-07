suppressPackageStartupMessages( library( pbdDMAT ) )
init.grid()

## construct local piece of a row-block distributed matrix
x <- matrix( 100*comm.rank() + 1:(3*5), ncol = 5 )
comm.print( x, all.rank = TRUE )

## get correct row-block dimensions and glue into ddmatrix
gdim <- c( allreduce( nrow( x ) ), ncol( x ) )
bd <- c( allreduce( nrow( x ), op = "max" ), ncol( x ) )
xnew <- new( "ddmatrix", Data = x, dim = gdim, bldim = bd, ldim = dim( x ), ICTXT = 2 )
print( xnew )
comm.print( submatrix( xnew ), all.rank = TRUE )

xt.cb <- as.colblock( xnew )
print( xt.cb )
comm.print( submatrix( xt.cb ), all.rank = TRUE )

xt.bc <- as.blockcyclic( xnew, bldim = c( 2, 2 ) )
print( xt.bc )
comm.print( submatrix( xt.bc ), all.rank = TRUE )

xt <- as.rowblock( t( xnew ) )
print( xt )
comm.print( submatrix( xt ), all.rank = TRUE )
finalize()
