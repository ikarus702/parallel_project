load("Time_Check_Kernel_Prod_Block.RData")

ls()

names(result.summary)

prod.summary <- result.summary

load("Time_Check_Kernel_Chol_Block.RData")
chol.summary <- result.summary

plot(prod.summary$mpi.time[,3]~prod.summary$block.size, type="b", lwd=2,xlim=c(0,33), ylim=c(0, ceiling(max(result.summary$single.time[3]))),
     main="Comp.Time for Kernel Matrix", xlab="block size", ylab="Time (sec.)")
lines(chol.summary$mpi.time[,3]~chol.summary$block.size,type="b",lwd=2,lty=2,col="blue")
points(0,result.summary$single.time[3],type="b",pch=2,lwd=2,col="red",lty=2)

legend("topright",c("MPI;Crossprod","MPI;Chol","Single"),lty=c(1,2,NA),pch=c(1,1,2),col=c("black","blue","red"),lwd=c(2,2,2))




plot(prod.summary$mpi.time.mat[301:400,3]~c(1:100),type="l", lwd=2,xlim=c(0,100), ylim=c(0, ceiling(max(prod.summary$mpi.time.mat[301:400,3]))),
     main="Comp.Time for Kernel Matrix", xlab="Iteration", ylab="Time (sec.)")
lines(chol.summary$mpi.time.mat[5:104,3]~c(1:100),type="l",lwd=2,lty=2,col="blue")
points(1,result.summary$single.time[3],type="b",pch=2,lwd=2,col="red",lty=2)

legend("topright",c("MPI;Crossprod","MPI;Chol","Single"),lty=c(1,2,NA),pch=c(NA,NA,2),col=c("black","blue","red"),lwd=c(2,2,2))



rm(list=ls())

load("Time_Check_Kernel_Prod_Row.RData")

ls()

names(result.summary)

prod.summary <- result.summary

load("Time_Check_Kernel_Chol_Row.RData")
chol.summary <- result.summary

plot(prod.summary$mpi.time[,3]~prod.summary$row.size, type="b", lwd=2, ylim=c(0, ceiling(max(prod.summary$single.time[,3]))),
     main="Comp.Time for Kernel Matrix", xlab="row size (n)", ylab="Time (sec.)")
lines(chol.summary$mpi.time[,3]~chol.summary$nrow,type="b",lwd=2,lty=2,col="blue")
lines(prod.summary$single.time[,3]~prod.summary$row.size,type="b",pch=2,lwd=2,col="red",lty=3)
legend("topleft",c("MPI;Crossprod","MPI;Chol","Single"),lty=c(1,2,3),pch=c(1,1,2),col=c("black","blue","red"),lwd=c(2,2,2))



rm(list=ls())

load("Time_Check_Kernel_Prod_Col.RData")

ls()

names(result.summary)

prod.summary <- result.summary

load("Time_Check_Kernel_Chol_Col.RData")
chol.summary <- result.summary

plot(prod.summary$mpi.time[,3]~chol.summary$ncol, type="b", lwd=2, ylim=c(0, ceiling(max(prod.summary$single.time[,3]))),
     main="Comp.Time for Kernel Matrix", xlab="column size (p)", ylab="Time (sec.)")
lines(chol.summary$mpi.time[,3]~chol.summary$ncol,type="b",lwd=2,lty=2,col="blue")
lines(prod.summary$single.time[,3]~chol.summary$ncol,type="b",pch=2,lwd=2,col="red",lty=3)
legend("topleft",c("MPI;Crossprod","MPI;Chol","Single"),lty=c(1,2,3),pch=c(1,1,2),col=c("black","blue","red"),lwd=c(2,2,2))
