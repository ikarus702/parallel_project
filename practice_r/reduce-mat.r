library( pbdMPI, quiet = TRUE )

x <- matrix( 10*comm.rank() + (1:6), nrow = 2 )

comm.print( x, all.rank = TRUE )

z <- gather( x ) # knows it's a matrix

comm.print( z, all.rank = TRUE )



my.rank <- comm.rank()
.comm.size <- comm.size()



finalize()

cat("The result of X is as follows.\n") 
print(x)

cat("The result of z is as follows.\n")
print(z)

print(unlist(z))







