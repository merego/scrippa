pdb<-read.pdb("slc3a2.pdb")
tor<-torsion.pdb(pdb)

bb<-atom.select(pdb,elety=c("N","CA","C"))

first.res<-2
last.res<-21

nres <- last.res-first.res+1

phi.tmp <- bb$atom[3:62]
psi.tmp <- bb$atom[4:63]

first.index<-3 * (first.res-1) - 1
last.index<-3 * (last.res-1) -1 

resid<-first.res-1
Table <- matrix(0,nrow=((last.res-1)*2),ncol=10)
ti<-1
for (i in seq(first.index,last.index,3)) {
 resid <- resid + 1
 blkphi <- c(1:4)+i
 Table[ti,c(1:4)] <- bb$atom[blkphi]
 Table[ti,6] <- tor$phi[resid]
 blkpsi <- c(1:4)+i+1
 Table[ti+1,c(1:4)] <- bb$atom[blkpsi]
 Table[ti+1,6] <- tor$psi[resid]
 cat("phi",tor$phi[resid],"\n")
 ti <-ti + 2 
}
Table[,5]<-1
Table[,7]<-0
Table[,8]<-100

write.table(Table,"restraints.dat")




