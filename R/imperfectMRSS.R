#' Imperfect Median Ranked Set Sampling Data Generation from Finite Population
#'
#' This function chooses Median Ranked Set samples from specific finite population using 
#' auxiliary variable with cycle sizes r1 and r2 and set sizes m1 and m2. 
#'  
#' @export
#' @param df : dataframe of the finite population
#' @param cat  :  the indicator variable that shows the group of units
#' @param catname  :  the group names 
#' @param aux : auxilary variable
#' @param var : variable of interest
#' 
#' @param m1  : Set size of first group
#' @param m2  : Set size of second group
#'   
#' @param r1  : Cycle size of first group
#' @param r2 : Cycle size of second group
#'
#' @return two median ranked set sample matrix with sizes r1xm1 and r2xm2 from finite population. Each row indicates a cycle.
#' 
#' @examples 
#' data(otolith)
#' imperfectMRSS(otolith,"sex",c("F","M"),"fish.length","age",3,3,4,3)
#'  
#' @seealso \code{\link[RSStest]{datagen_RSS}},  \code{\link[RSStest]{teststat_RSS}}
#'  \code{\link[RSStest]{teststat_MRSS}},\code{\link[RSStest]{imperfectRSS}}
#' @references MacEachern, S. N., Öztürk, Ö., Wolfe, A. D. (2002). A new ranked set sample estimator of variance. Journal of the Royal Statistical Society: Series B., 64, Part 2 177–188.
#'
#' Özturk, Ö., Balakrishnan N (2009) Exact two-sample nonparametric test for quantile difference between two populations based on ranked set samples. Ann Inst Stat Math 61(1):235–249
#' 
#' Özdemir, Y. A., Ebegil, M., & Gökpinar, F. (2017). A test statistic based on ranked set sampling for two normal means. Communications in Statistics-Simulation and Computation, 46(10), 8077-8085.
#' 
#' Özdemir, Y. A., Ebegil, M., & Gökpinar, F. (2019). A test statistic for two normal means with median ranked set sampling. Iranian Journal of Science and Technology, Transactions A: Science, 43(3), 1109-1126.
#' 
imperfectMRSS=function (df,cat,catname,aux,var,r1,r2,m1,m2){
G1=df[df[,cat]==catname[1],]
G2=df[df[,cat]==catname[2],]

dx1=c();dx2=c()

for (i in 1:r1) {
  u1=sample(1:nrow(G1),size=m1^2,replace=FALSE)
  X1=matrix(G1[u1,aux],nrow=m1);
  Y1=matrix(G1[u1,var],nrow=m1);
  xs1=t(apply(X1,1,order))
  if (m1%%2==1) {
    YS1=diag(Y1[1:m1,xs1[,(m1+1)/2]])
    dx1=c( dx1, YS1)
  }
  else {
    YS1=diag(rbind(Y1[1:(m1/2),xs1[,m1/2]],Y1[((m1/2)+1):m1,xs1[,m1/2+1]]))
    dx1=c( dx1, YS1)
  }
}
  for (i in 1:r2) {
    u2=sample(1:nrow(G2),size=m2^2,replace=FALSE)
    X2=matrix(G2[u2,aux],nrow=m2);
    Y2=matrix(G2[u2,var],nrow=m2);
    xs2=t(apply(X2,1,order))
    if (m2%%2==1) {
      YS2=diag(Y2[1:m2,xs2[,(m2+1)/2]])
      dx2=c( dx2, YS2)
    }
    else {
      YS2=diag(rbind(Y2[1:(m2/2),xs2[,m2/2]],Y2[((m2/2)+1):m2,xs2[,m2/2+1]]))
      dx2=c( dx2, YS2)
    }
  }
y1=t(matrix(dx1,ncol=r1))
y2=t(matrix(dx2,ncol=r2))
out=list(y1,y2)
return(out)
}