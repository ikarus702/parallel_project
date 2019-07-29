load("Time_Check_Matmult_Block.RData")

ls()

names(result.summary)

plot(result.summary$mpi.time[,3]~result.summary$block.size, type="b", lwd=2, ylim=c(0, c(max(result.summary$mpi.time[,3])+0.1)),
     main="Comp.Time for OLS", xlab="block size", ylab="Time (sec.)")
abline(h=result.summary$single.time[3],lty=2,lwd=2,col="red")
legend("topright",c("MPI","Single"),lty=c(1,2),pch=c(1,NA),col=c("black","red"),lwd=c(2,2))

rm(list=ls())

load("Time_Check_Matmult_Matsize.RData")

plot(result.summary$mpi.time[,3]~result.summary$mat.size, type="b", lwd=2, ylim=c(0, ceiling(max(result.summary$single.time[,3]))),
     main="Comp.Time for OLS", xlab="row (column) size (n)", ylab="Time (sec.)")
lines(result.summary$single.time[,3]~result.summary$mat.size,type="b",pch=2,lwd=2,col="red",lty=2)
legend("topleft",c("MPI","Single"),lty=c(1,2),pch=c(1,2),col=c("black","red"),lwd=c(2,2))

rm(list=ls())

load("Time_Check_Matmult_Node.RData")

plot(result.summary$mpi.time[1:5,3]~result.summary$node, type="b", lwd=2, ylim=c(0, ceiling(max(result.summary$mpi.time[,3])+5)),
     main="Comp.Time for OLS", xlab="No. of Nodes", ylab="Time (sec.)")
lines(result.summary$mpi.time[6:10,3]~result.summary$node, type="b", lwd=2, lty=2)

abline(h=result.summary$single.time[1,3],lty=1,col="red")
abline(h=result.summary$single.time[2,3],lty=2,col="red",lwd=2)


legend("topright",c("MPI; n=1000","MPI; n=5000", "Single; n=1000", "Single; n=5000"),lty=c(1,2,1,2),pch=c(1,1,NA,NA),col=c(rep("black",2),rep("red",2)),lwd=c(rep(2,4)))
