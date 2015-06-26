name <- "report_09.dat"
name.clean <- paste(name,".clean",sep="")

r<-read.table(name,header=TRUE)

# Find duplicates
dpl<-duplicated(r[,1])
# Find index of duplicates (remove the first  keep the second)
removeindex<-which(dpl)-1 
rclean<-r[-removeindex,]

# check
a<-read.table("AvgLoss09.dat")
na <- nrow(a)
nr <- nrow(rclean)

cat(na,nr,"\n")

if (na==nr) {
 write.table(file=name.clean,rclean,quote=FALSE,row.names=FALSE)
 cat("Written to ",name.clean,"\n")
} else {
 cat("something went wrong... different row numbers\n")
}

