dir <- "/home/HAO/JLM_Rpackage"
setwd(dir)

## load genotype data
geno <- read.table("gen10pop_in.csv",sep=",",header=T)
snp <- geno[,2:ncol(geno)]
parent <- 20
bin_num <- nrow(snp)
line_num <- ncol(snp)
kin = matrix(0,nrow=line_num,ncol=line_num)

##calculate kinship matrix
for(n in 1:bin_num) {
	k<-matrix(0,line_num,parent)
	a<-snp[n,]

	for(i in 1:parent){
		k[,i] <- 2*(a==i)
	}
	m = k %*% t(k)
	kin = kin + m	
}
diagnoal <- sum(diag(kin))
options(digits=10)
norm_factor <- diagnoal/line_num
kin_norm <- kin/norm_factor

ind <- c(1:line_num)
line <- paste("col",ind,sep="")
result1 <- rbind(line,kin_norm)
parm <- c("parm",rep(1,line_num))
row <- c("row",c(1:line_num))
result2 <- cbind(parm,row,result1)
write.table(result2,"kinship_out.csv",sep=",",quote=FALSE,row.names=FALSE,col.name=FALSE)

## The first column is eigen value, and the remaining are eigen vectors.
eig <- eigen(kin_norm)
save(eig,file="eigen_out.RData")


