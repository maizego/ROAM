dir <- "/home/HAO/JLM_Rpackage"
setwd(dir)

data <- read.table("genP1C1_in.csv",sep=",",header=T)
pos <- data[,3]
snp <- data[,4:ncol(data)]
###### project data to 300 markers, similar to the program of imputation.
## The first snp is missing.
for (i in 1:ncol(snp)){
	if (snp[1,i] == 0){
		for (j in 2:nrow(snp)){
			if (snp[j,i] != 0){			
				snp[1:j-1,i] = rep(snp[j,i],j-1)
				break
			}
		}	
	}
}
## The last snp is missing.
for (i in 1:ncol(snp)){
	if (snp[nrow(snp),i] == 0){
		for (j in 2:nrow(snp)){
			if (snp[nrow(snp)-j+1,i] != 0){			
				snp[(nrow(snp)-j+2):nrow(snp),i] = rep(snp[nrow(snp)-j+1,i],j-1)
				break
			}
		}	
	}
}
## missing snp is between two markers, imputat accoding to position.
num <- 1
for (i in 1:ncol(snp)){
	for (j in 2:nrow(snp)-1){
		if (snp[j,i] == 0 && snp[j+1,i] == 0){
			num <- num + 1
		}
		if (snp[j,i] == 0 && snp[j+1,i] != 0){
			snp_row <- j-num+1
			for (k in snp_row:j){
				if (pos[k]-pos[snp_row-1] >= pos[j+1]-pos[k]){
					snp[k,i] <- snp[j+1,i]	
				}
				else {
					snp[k,i] <- snp[snp_row-1,i]
				}
			}			
			num <- 1
		}
		
	}
}
snpr = cbind(data[,1:3],snp)
write.table(snpr,"genP1C1_out.csv",sep=",",quote=FALSE,row.names=FALSE,col.name=TRUE)

