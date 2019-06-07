suppressPackageStartupMessages( library( pbdDMAT ) )
init.grid()

x.dmat <- ddmatrix( "rnorm", nrow = 10, ncol = 3 )
comm.print( x.dmat )

colsd <- apply( X = x.dmat, FUN = sd, MARGIN = 2 )
colsd.loc <- as.matrix( colsd )

comm.print( colsd )
comm.print( colsd.loc )

finalize ()
