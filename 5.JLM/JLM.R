dir <- "/home/HAO/JLM_Rpackage/"
setwd(dir)

library(regress)
library(stats4)
library(MASS)

######################################################################
####get the lambda 
k <- read.table(file=paste(dir,"kinship_in.csv",sep=""),sep=",",header=T)
line_num <- nrow(k)
k <- k[,3:(line_num+2)]
colnames(k) <- c(1:line_num)
kk <- as.matrix(k)
e <- diag(1,line_num)

phe <- read.table(file=paste(dir,"phenorm_in.csv",sep=""),sep=",",header=T)
y <- phe[,2]

result <- regress(y~1,~kk+e,identity=FALSE)
lin1 <- result$sigma[1]
lin2 <- result$sigma[2]
lin1n <- as.numeric(lin1)
lin2n <- as.numeric(lin2)
lambda <- lin1n/lin2n
lambda_result <- c(lin1n,lin2n,lambda)
write.table(lambda_result,file=paste(dir,"lambda_out.csv",sep=""),sep=",",row.names=F)

######################################################################
####mixed linear model 
gen<-read.table(file=paste(dir,"gen10pop_in.csv",sep=""),sep=",",header=TRUE)[,-1]
qq<-get(load(file=paste(dir,"eigen_in.RData",sep="")))
delta<-qq[[1]]
uu<-qq[[2]]
h<-1/(delta*lambda+1)

n<-ncol(gen)
m<-nrow(gen)
s<-1
r<-20 ###the number of parents

x<-matrix(1,n,1)
xu<-t(uu)%*%x
yu<-t(uu)%*%y
yy<-sum(yu*h*yu)
yx<-matrix(0,s,1)
xx<-matrix(0,s,s)

for(i in 1:s){
    yx[i]<-sum(yu*h*xu[,i])
    for(j in 1:s){
		xx[i,j]<-sum(xu[,i]*h*xu[,j])
    }
}

loglike<-function(theta){
   xi<-exp(theta)
   tmp0<-xi*zz+diag(r)
   tmp<-solve(tmp0)
   yHy<-yy-xi*t(zy)%*%tmp%*%zy
   yHx<-yx-xi*zx%*%tmp%*%zy
   xHx<-xx-xi*zx%*%tmp%*%t(zx)
   logdt2<-log(det(tmp0))
   loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
   return(-loglike)
}

fixed<-function(xi){
   tmp<-solve(xi*zz+diag(r))
   yHy<-yy-xi*t(zy)%*%tmp%*%zy
   yHx<-yx-xi*zx%*%tmp%*%zy
   xHx<-xx-xi*zx%*%tmp%*%t(zx)
   zHy<-zy-xi*zz%*%tmp%*%zy
   zHx<-zx-xi*zx%*%tmp%*%zz
   zHz<-zz-xi*zz%*%tmp%*%zz
   beta<-solve(xHx,yHx)
   tmp<-solve(xHx)
   sigma2<-(yHy-t(yHx)%*%tmp%*%yHx)/(n-s)
   gamma<-xi*zHy-xi*t(zHx)%*%tmp%*%yHx
   var<-abs((xi*diag(r)-xi^2*zHz)*as.numeric(sigma2))
   stderr<-sqrt(diag(var))
   return(c(gamma,stderr,beta,sigma2))
}
start<-date()
parr<-numeric()
blupp<-numeric()

for(k in 1:m) {
   z<-matrix(0,n,r)
   a<-gen[k,]

for(i in 1:r){
z[,i]<-2*(a==i)
}
zu<-t(uu)%*%z

   zx<-matrix(0,s,r)
   zy<-matrix(0,r,1)
   zz<-matrix(0,r,r)
   for(i in 1:s){
      for(j in 1:r){
         zx[i,j]<-sum(xu[,i]*h*zu[,j])
      }
   } 
   for(i in 1:r){
      zy[i]<-sum(yu*h*zu[,i])
      for(j in 1:r){
         zz[i,j]<-sum(zu[,i]*h*zu[,j])
      }
   }


theta<--0
parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
xi<-exp(parm$par)
conv<-parm$convergence
fn1<-parm$value
fn0<-loglike(c(-Inf))
lrt<-2*(fn0-fn1)
hess<-parm$hessian
parmfix<-fixed(xi)
gamma<-parmfix[1:r]
stderr<-parmfix[(r+1):(2*r)]
beta<-parmfix[2*r+1]
sigma2<-parmfix[2*r+2]
tau_k<-xi*sigma2
initial<-exp(theta)

par<-data.frame(k,initial,conv,fn1,fn0,lrt,beta,sigma2,xi,tau_k)
blup<-c(k,gamma,stderr)
parr<-rbind(parr,par)
blupp<-rbind(blupp,blup)
}

gamma<-NULL
stderr<-NULL
for(i in 1:r){
   gamma<-c(gamma,paste("gamma",i,sep=""))
   stderr<-c(stderr,paste("stderr",i,sep=""))
}
varnames<-c("bin",gamma,stderr)
colnames(blupp)<-varnames
colnames(parr)[1]<-"bin"

## parms file comprises the LRT and blup file comprises gamma and SE for each parent. 
write.table(parr,file=paste(dir,"parms_out.csv",sep=""),sep=",",row.names=F)
write.table(blupp,file=paste(dir,"blup_out.csv",sep=""),sep=",",row.names=F)





