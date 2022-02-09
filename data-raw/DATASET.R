## code to prepare `DATASET` dataset goes here
library(dplyr)
set.seed(123)
b0=b1=b2=1
#b3=2
narea=30
nsub=3
sigmav=8
sigmau=8
sigmae=8
n=narea*nsub
X1=runif(n,0,1)
X2=rnorm(n,10,1)
wr=runif(n,10,20)
codearea=c()
for(i in 1:narea){
  codearea=c(codearea,rep(i,nsub))
}
vardir=rep(sigmae,n)
wdat=data.frame(w=wr,code=codearea)
wmod=wdat%>%group_by(code)%>%summarise(weight=(w/sum(w)))
w=wmod$weight
u=c()
for(i in 1:narea){u=c(u,rnorm(nsub,0,sd=sqrt(sigmau)))}
v=rnorm(narea,0,sd=sqrt(sigmav))
vim=c()
for(i in 1:n){vim[i]=v[codearea[i]]}
e=rnorm(n,0,sqrt(vardir))
#simulated dataset of Y (SUBAREA)
Teta_sub=b0+b1*X1+b2*X2+vim+u
y_sub=Teta_sub+e
dataTwofold=data.frame(cbind(y_sub,X1,X2,codearea,w,vardir))
names(dataTwofold)=c("y","x1","x2","codearea","w","vardir")
#simulated dataset of Y (SUBAREA NONSAMPLED)
del_index=sample(1:nrow(dataTwofold), 10, replace=F)
dataTwofoldNS=dataTwofold
dataTwofoldNS[c(del_index),c(1,6)]=NA
dataTwofoldNS

usethis::use_data(dataTwofold,overwrite = TRUE)
usethis::use_data(dataTwofoldNS,overwrite = TRUE)
