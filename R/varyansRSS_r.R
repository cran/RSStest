varyansRSS_r= function(x) {
d=dim(x) 
r=d[1]
m=d[2]
xsm=colMeans(x)
MSE=sum((x-matrix(rep(xsm,r),nrow = r,byrow=TRUE))^2)/(m*(r-1))
MST=sum((x-mean(x))^2)/(m-1)-sum((x-matrix(rep(xsm,r),nrow = r,byrow=TRUE))^2)/(m-1)
return(((m-1)*MST+(m*r-m+1)*MSE)/(m*r))
}