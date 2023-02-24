#' Median Ranked Set Sampling Test 
#'
#' This function tests for the difference of two population means using ranked set sampling given in Özdemir, Ebegil and Gökpınar (2019). 
#'  
#' @export
#'
#' @param x1 A (non-empty) numeric matrix (m1xr1) of median ranked set sample for Group 1 with set size m1 and cycle size r1.
#' 
#' @param x2 A (non-empty) numeric matrix (m2xr2) of median ranked set sample for Group 2 with set size m2 and cycle size r2.
#'  
#' @param alpha A scalar value of the significance level for hypothesis testing used in the table. Default is 0.05.
#'  
#' @param tn A scalar value of the number of repetitions of Monte Carlo simulation. Default is 2000.
#' 
#' @param table A logical value that shows table gives the results of the hypothesis test are printed out. Default is TRUE.
#' 
#' @param alternative  A character string specifying the alternative hypothesis, must be one of "two-sided", "right" or "left". Can be abbreviated. Default is "two-sided". 
#'
#' @return If table is TRUE the hypothesis test results table includes sample sizes, test statistics, p values and test results are printed out.
#' 
#' @examples 
#' x1=matrix(c(1,2.3, 3.4,4.5,5.6,4 ),nrow=3)
#' x2=matrix(c(2,3.2, 4.2,6.5,4.6,6 ),nrow=3)
#' teststat_MRSS(x1,x2,tn=1000)


#' @references MacEachern, S. N., Öztürk, Ö., Wolfe, A. D. (2002). A new ranked set sample estimator of variance. Journal of the Royal Statistical Society: Series B., 64, Part 2 177–188.
#'
#' Özturk, Ö., Balakrishnan N (2009) Exact two-sample nonparametric test for quantile difference between two populations based on ranked set samples. Ann Inst Stat Math 61(1):235–249
#' 
#' Özdemir, Y. A., Ebegil, M., & Gökpinar, F. (2017). A test statistic based on ranked set sampling for two normal means. Communications in Statistics-Simulation and Computation, 46(10), 8077-8085.
#' 
#' Özdemir, Y. A., Ebegil, M., & Gökpinar, F. (2019). A test statistic for two normal means with median ranked set sampling. Iranian Journal of Science and Technology, Transactions A: Science, 43(3), 1109-1126.
#' 
#' @seealso \code{\link[RSStest]{datagen_MRSS}}, \code{\link[RSStest]{datagen_RSS}},
#'  \code{\link[RSStest]{teststat_RSS}}


teststat_MRSS=function(x1,x2,alpha=0.05,alternative="two-tailed",tn=2000,table=TRUE){
m11=ncol(x1);r11=nrow(x1)
m22=ncol(x2);r22=nrow(x2)
s1_2=ifelse(m11%%2==1,var(c(x1)),(var(c(x1[,1:m11/2]))+var(c(x1[,(m11/2+1):m11])))/2)
s2_2=ifelse(m22%%2==1,var(c(x2)),(var(c(x2[,1:m22/2]))+var(c(x2[,(m22/2+1):m22])))/2)
V=s1_2/(m11*r11)+s2_2/(m22*r22)
T=(mean(c(x1))-mean(c(x2)))/sqrt(V)
Tb=c()
for (i in 1:tn) {
x1b=datagen_MRSS(0,1,m11,r11)
x2b=datagen_MRSS(0,1,m22,r22)
s1_2b=ifelse(m11%%2==1,var(c(x1b)),(var(c(x1b[,1:m11/2]))+var(c(x1b[,(m11/2+1):m11])))/2)
s2_2b=ifelse(m22%%2==1,var(c(x2b)),(var(c(x2b[,1:m22/2]))+var(c(x2b[,(m22/2+1):m22])))/2)
Vb=s1_2b/(m11*r11)+s2_2b/(m22*r22)
Tb[i]=(mean(c(x1b))-mean(c(x2b)))/sqrt(Vb)
}
px=mean(Tb>T)
if (alternative=="left" | alternative=="L"| alternative=="l") {
  p=1-px
  r=p<alpha
  A=c("  p-value      (left-sided)",p,NA)
  decision=paste0("*Mean of Group 1  is ", ifelse(r==0," not",""), " smaller than Mean of Group 2 at significance level ",alpha,".")
  
}
  else if  (alternative=="right" | alternative=="R"| alternative=="r") {
    p=px
    r=p<alpha
    A=c("  p-value      (right-sided)",p,NA)
    decision=paste0("*Mean of Group 1  is ", ifelse(r==0," not",""), " greater than Mean of Group 2 at significance level ",alpha,".")
  }
  else if (alternative=="two-tailed"| alternative=="T" | alternative=="t") {
    p=2*min(px,1-px)
    r=p<alpha
    A=c("  p-value      (two-sided)",p,NA)
    decision=paste0("*Group Means are", ifelse(r==1," not",""), " equal at significance level ",alpha,".")
  }

if (table==FALSE) {
  return(c(r,p,T))
}
else if (table==TRUE) {
 

  print_md(hux(Group_No=1:2,Set_Size=c(m11,m22),C.Size=c(r11,r22),S.Size=c(m11,m22),S.Mean=c(mean(x1),mean(x2)),S.Var=c(s1_2,s2_2),T.Stat=c(T,NA),p_palue=c(p,NA)))
  
  
  # HT<-hux(c( "Group", paste("Group",1:2)),c("Set Size",c(m11,m22)),c("Cycle Size",c(r11,r22)),c("Sample Size",c(r11*m11,r22*m22)),c("Sample Mean",c(mean(c(x1)),mean(c(x2)))),c("Sample Variance",c(s1_2,s2_2)),c("Test Statistic",T,NA),A,add_colnames = FALSE)
   # bold(HT)[1,]<-TRUE
  #bold(HT)[,1]<-TRUE
  #header_cols(HT)=c(NA,NA,"Descriptives",NA)
  #bottom_border(HT)[1,]  <- 0.2
  #top_border(HT)[1,]  <- 0.2
  #right_border(HT)[,1] <- 0.2
  #left_border(HT)[,1] <- 0.2
  #right_border(HT)[,4] <- 0.2
  #right_border(HT)[,6] <- 0.2
  #right_border(HT)[,8] <- 0.2
  #align(HT)[,1:8]          <- 'center'
  #number_format(HT)[,4:7]      <- 3
  #number_format(HT)[,2:4]      <- 0
  #HT <-add_footnote(HT,decision)
  #HT
  #return(HT)
}
}