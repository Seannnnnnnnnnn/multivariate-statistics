X = as.matrix(read.table("data-sets/Xtrainphoneme.txt",sep = ","))[, -257]
XXnew = as.matrix(read.table("data-sets/Xtestphoneme.txt",sep = ","))[, -257]
Z = as.numeric(read.table("data-sets/Ztrainphoneme.txt",sep = ","))
Znew = as.numeric(read.table("data-sets/Ztestphoneme.txt",sep = ","))
n=300
nnew=1417
Xcent=scale(X,scale=F)
Zcent=Z-mean(Z)
#1 Prediction by LS
betahat= solve(t(Xcent)%*%Xcent)%*%t(Xcent)%*%Zcent
repbarX=matrix(rep(colMeans(X),nnew),nrow=nnew,byrow=T)
ZpredLS=mean(Z)+ (XXnew-repbarX)%*%betahat
#2 Prediction by PC
compo=prcomp(X,retx=T)
phi=(compo$rotation)
Y=compo$x
Zpred=c()
for(q in 2:10)
{
  YDATA=Y[,1:q]
  betahatPC= solve(t(YDATA)%*%YDATA)%*%t(YDATA)%*%Zcent
  YNew=(XXnew-repbarX)%*%phi
  Zpred=cbind(Zpred,mean(Z)+YNew[,1:q]%*%betahatPC)
}
#3 Plot prediction errors
log_sq_error_list = list()
error_list = list()
for(i in 1:9){
  log_sq_error_list[[i]] = log((Zpred[,i]-Znew)^2)
  error_list[[i]] = (Zpred[,i]-Znew)
}
log_sq_error_list[[10]] = log((ZpredLS-Znew)^2)
error_list[[10]] = ZpredLS-Znew
par(mfrow=c(1,2))
boxplot(log_sq_error_list, ylim=c(-20,20),
        names = c(2:10,"LS"))
boxplot(error_list,ylim=c(-100,100),
        names = c(2:10,"LS"))
# Best PC predictor obtained when $q=10$:
plot(Znew,ZpredLS)
abline(0,1)
plot(Znew,Zpred[,9])
abline(0,1)
  
