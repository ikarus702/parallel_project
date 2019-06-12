load("Time_Check_Kernel_Col.RData")

ls()

names(result.summary)

plot(result.summary$mpi.time[,3]~result.summary$col.size, type="b", lwd=2, ylim=c(0, ceiling(max(result.summary$mpi.time[,3]))),
     main="Comp.Time for Kernel Matrix", xlab="column size (p)", ylab="Time (sec.)")
lines(result.summary$single.time[,3]~result.summary$col.size,type="b",pch=2,lwd=2,col="red",lty=2)
legend("topleft",c("MPI","Single"),lty=c(1,2),pch=c(1,2),col=c("black","red"),lwd=c(2,2))

rm(list=ls())

load("Time_Check_Kernel_Row.RData")

plot(result.summary$mpi.time[,3]~result.summary$nrows, type="b", lwd=2, ylim=c(0, ceiling(max(result.summary$mpi.time[,3]))),
     main="Comp.Time for Kernel Matrix", xlab="row size (n)", ylab="Time (sec.)")
lines(result.summary$single.time[,3]~result.summary$nrows,type="b",pch=2,lwd=2,col="red",lty=2)
legend("topleft",c("MPI","Single"),lty=c(1,2),pch=c(1,2),col=c("black","red"),lwd=c(2,2))

rm(list=ls())

load("Time_Check_Kernel_Prod_Col.RData")

plot(result.summary$mpi.time[,3]~result.summary$ncol.size, type="b", lwd=2, ylim=c(0, ceiling(max(result.summary$mpi.time[,3]))),
     main="Comp.Time for Kernel Matrix", xlab="column size (p)", ylab="Time (sec.)")
lines(result.summary$single.time[,3]~result.summary$ncol.size,type="b",pch=2,lwd=2,col="red",lty=2)
legend("topleft",c("MPI","Single"),lty=c(1,2),pch=c(1,2),col=c("black","red"),lwd=c(2,2))

rm(list=ls())

load("Time_Check_Kernel_Prod_Row.RData")

plot(result.summary$mpi.time[,3]~result.summary$nrow.size, type="b", lwd=2, ylim=c(0, ceiling(max(result.summary$mpi.time[,3]))),
     main="Comp.Time for Kernel Matrix", xlab="row size (n)", ylab="Time (sec.)")
lines(result.summary$single.time[,3]~result.summary$nrow.size,type="b",pch=2,lwd=2,col="red",lty=2)
legend("topleft",c("MPI","Single"),lty=c(1,2),pch=c(1,2),col=c("black","red"),lwd=c(2,2))
