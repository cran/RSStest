#' Ranked Set Sampling Test 
#'
#' This function tests for the difference of two population means using ranked set sampling given in Özdemir, Ebegil and Gökpınar (2017). 
#'  
#' @export
#' @import huxtable
#' @importFrom stats integrate median pchisq pf rchisq rnorm var
#' @importFrom graphics boxplot par
#' @importFrom utils data
#' @param x1 A (non-empty) numeric matrix (m1xr1) of ranked set sample for Group 1 with set size m1 and cycle size r1.
#' 
#' @param x2 A (non-empty) numeric matrix (m2xr2) of ranked set sample for Group 2 with set size m2 and cycle size r2.
#'  
#' @param alpha A scalar value of the significance level for hypothesis testing used in the table. Default is 0.05.
#'  
#' @param table A logical value that shows table gives the results of the hypothesis test are printed out. Default is TRUE.
#' 
#' @param alternative  A character string specifying the alternative hypothesis, must be one of "two-sided", "right" or "left". Can be abbreviated. Default is "two-sided". 
#'
#'@return If table is TRUE the hypothesis test results table includes sample sizes, test statistics, critical values and test results are printed out.
#' 
#' @examples 
#' x1=matrix(c(1,2.3, 3.4,4.5,5.6,4 ),nrow=3)
#' x2=matrix(c(2,3.2, 4.2,6.5,4.6,6 ),nrow=3)
#' teststat_RSS(x1,x2)
#'  
#' @references MacEachern, S. N., Öztürk, Ö., Wolfe, A. D. (2002). A new ranked set sample estimator of variance. Journal of the Royal Statistical Society: Series B., 64, Part 2 177–188.
#'
#' Özturk, Ö., Balakrishnan N (2009) Exact two-sample nonparametric test for quantile difference between two populations based on ranked set samples. Ann Inst Stat Math 61(1):235–249
#' 
#' Özdemir, Y. A., Ebegil, M., & Gökpinar, F. (2017). A test statistic based on ranked set sampling for two normal means. Communications in Statistics-Simulation and Computation, 46(10), 8077-8085.
#' 
#' Özdemir, Y. A., Ebegil, M., & Gökpinar, F. (2019). A test statistic for two normal means with median ranked set sampling. Iranian Journal of Science and Technology, Transactions A: Science, 43(3), 1109-1126.
#' 
#'  @seealso \code{\link[RSStest]{datagen_MRSS}}, \code{\link[RSStest]{datagen_RSS}}, \code{\link[RSStest]{teststat_MRSS}}



teststat_RSS=function (x1,x2,alpha=0.05,alternative="two-tailed",table=TRUE)
{

#T=data.matrix(read.csv("data/CVT.csv"))
#load(file="data/CVT.rda")
  #data(CVT)
  V=c()
a1=dim(x1);r11=a1[1];m11=a1[2]
a2=dim(x2);r22=a2[1];m22=a2[2]
if ((alternative=="right" | alternative=="left") & alpha==0.01){
  critical=CVT[,6][CVT[,1]==r11 & CVT[,3]==r22 & CVT[,2]==m11 & CVT[,4]==m22] 
} else if ((alternative=="right" | alternative=="left") & alpha==0.05){
  critical=CVT[,8][CVT[,1]==r11 & CVT[,3]==r22 & CVT[,2]==m11 & CVT[,4]==m22] 
} else if ((alternative=="right" | alternative=="left") & alpha==0.10){
  critical=CVT[,9][CVT[,1]==r11 & CVT[,3]==r22 & CVT[,2]==m11 & CVT[,4]==m22] 
}  else if (alternative=="two-tailed" & alpha==0.01){
  critical=CVT[,5][CVT[,1]==r11 & CVT[,3]==r22 & CVT[,2]==m11 & CVT[,4]==m22] 
}  else if (alternative=="two-tailed" & alpha==0.05){
  critical=CVT[,7][CVT[,1]==r11 & CVT[,3]==r22 & CVT[,2]==m11 & CVT[,4]==m22] 
} else if (alternative=="two-tailed" & alpha==0.10){
  critical=CVT[,8][CVT[,1]==r11 & CVT[,3]==r22 & CVT[,2]==m11 & CVT[,4]==m22] 
}
ssum=c(NA, 1.363380, 1.567605, 1.704340, 1.804940, 1.883436,1.947193,
       2.000506,2.046085,2.085728,2.120640,2.151820,2.179930,2.205440 ,
       2.228750,2.250180,2.270040,2.288420,2.305630,2.321740)
u1=varyansRSS_r(x1)
u2=varyansRSS_r(x2)
m1=ncol(x1);r1=nrow(x1)
m2=ncol(x2);r2=nrow(x2)
V=u1*ssum[m1]/(m1^2*r1)+u2*ssum[m2]/(m2^2*r2)
Th=(mean(x1)-mean(x2))/V^0.5
U=c(Th,mean(x1),mean(x2),u1, u2)
if (alternative=="left") {
  r=Th<critical
  decision=paste0("*Mean of Group 1  is ", ifelse(r==0," not",""), " smaller than Mean of Group 2 at significance level ",alpha,".")
}
else if  (alternative=="right") {
  r=Th>critical
  decision=paste0("*Mean of Group 1  is ", ifelse(r==0," not",""), " greater than Mean of Group 2 at significance level ",alpha,".")
}
else if (alternative=="two-tailed") {
  r=abs(Th)>critical
  decision=paste0("*Group Means are", ifelse(r==1," not",""), " equal at significance level ",alpha,".")
}
if (table==FALSE) {
  return(c(r,Th))
}
else if (table==TRUE)
{
 
  #C1=1:2
  #C2=c(ne)
  #C3=c(mut)
  #C4=c(s2tt)
  #C5=c(rep(NA,ifelse(k%%2==0,ceiling(k/2),floor(k/2))),T,rep(NA,floor(k/2-1)))
  #C6=c(rep(NA,ifelse(k%%2==0,ceiling(k/2),floor(k/2))),p,rep(NA,floor(k/2-1)))
  #HT<-hux(Group_No=C1,Sample_Size=C2,Sample_Mean=C3,Sample_Var=C4,Test_Stat=C5,p_value=C6)
  #align(HT)[,1:5]          <- 'center'
  #number_format(HT)[,3:6]      <- 2
  #number_format(HT)[,2]      <- 0
  uu1=c(critical,NA)
  uu2=c(-critical,NA)
  uu3=c(-critical,critical)
  a=c()
  a[1]=ifelse(alternative=="right",uu1[1],ifelse(alternative=="left",uu2[1],uu3[1]))
  a[2]=ifelse(alternative=="right",uu1[2],ifelse(alternative=="left",uu2[2],uu3[2]))
  print_md(hux(Group_No=1:2,Set_Size=c(m1,m2),C.Size=c(r1,r2),S.Size=c(m1,m2),S.Mean=c(mean(x1),mean(x2)),S.Var=c(u1,u2),T.Stat=c(Th,NA),C.Value=a))
  #return(HT)
  
  
   #
  #, C_Size=c(r1,r2), S_Size=c(r1*m1,r2*m2),S_Mean=c(mean(c(x1)),mean(c(x2))),Sample_Var=c(u1,u2)
  #        ,Test_Statistic=Th,c("Critical Value",ifelse(alternative=="right",critical,ifelse(alternative=="left",-critical,-critical)),ifelse(alternative=="right",NA,ifelse(alternative=="left",NA,crit))),add_colnames = FALSE)
  #bold(HT)[1,]<-TRUE
  #bold(HT)[,1]<-TRUE
  #header_cols(HT)=c(NA,NA,"Descriptives",NA)
  #bottom_border(HT)[1,]  <- 0.2
  #top_border(HT)[1,]  <- 0.2
  #right_border(HT)[,1] <- 0.2
  #left_border(HT)[,1] <- 0.2
  #right_border(HT)[,4] <- 0.2
  #right_border(HT)[,6] <- 0.2
  #right_border(HT)[,8] <- 0.2
  #align(HT)[,1:7]          <- 'center'
  #align(HT)[,8]          <- 'center'
  #number_format(HT)[,4:8]      <- 3
  #number_format(HT)[,2:4]      <- 0
  # HT <-add_footnote(HT,decision)
  # HT
  #return(HT)
}
}