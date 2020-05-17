dir <- "/home/HAO/JLM_Rpackage"
setwd(dir)

data <- read.table("gen2popC10_in.csv",sep=",",header=T)
pos <- data[,2]
snp <- data[,3:ncol(data)]
breakp <- c(rep(0,nrow(data)))
binnum <- c(rep(0,nrow(data)))
breaknum <- 1
num <- 1
for (i in 1:ncol(snp)){
	print(i)
	gc()
	for (j in 1:(nrow(snp)-1)){
		if (snp[j,i] != snp[j+1,i]){
			breakp[j+1] <- pos[j+1]
		}	
	}
}
for (m in 1:length(breakp)){
	if (breakp[m] != 0){	
		repnum = m - breaknum
		binnum[breaknum:(m-1)] = c(rep(num,repnum))
		breaknum = m
		num = num + 1
	}
	
}
breakp[1] <- pos[1]
lastrepnum <- length(breakp)-breaknum+1
binnum[breaknum:length(breakp)] = c(rep(num,lastrepnum))
result = cbind(binnum,breakp,pos)
write.table(result,"gen2popC10_out.csv",sep=",",quote=FALSE,row.names=FALSE,col.name=FALSE)
