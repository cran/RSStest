#'  Ranked Set Sampling Data Generation
#'
#' This function generates random samples from normal population using ranked set sampling with
#'  mean \eqn{\mu}  and standard deviation \eqn{\sigma}  using cycle size r and set size m. 
#'  
#' @export
#'
#' @param mu : Normal population mean \eqn{\mu} 
#' 
#' @param s  :   Normal population standard deviation \eqn{\sigma}
#'  
#' @param m  : Set size 
#'  
#' @param r  : Cycle size
#'
#' @return A sample matrix with size rxm generated from normal distribution using ranked set sampling. Each row indicates a cycle.
#' 
#' @examples 
#' datagen_RSS(0,1,2,3)
#'  
#' @seealso \code{\link[RSStest]{datagen_MRSS}},  \code{\link[RSStest]{teststat_RSS}}
#'  \code{\link[RSStest]{teststat_MRSS}}
#' @references MacEachern, S. N., Öztürk, Ö., Wolfe, A. D. (2002). A new ranked set sample estimator of variance. Journal of the Royal Statistical Society: Series B., 64, Part 2 177–188.
#'
#' Özturk, Ö., Balakrishnan N (2009) Exact two-sample nonparametric test for quantile difference between two populations based on ranked set samples. Ann Inst Stat Math 61(1):235–249
#' 
#' Özdemir, Y. A., Ebegil, M., & Gökpinar, F. (2017). A test statistic based on ranked set sampling for two normal means. Communications in Statistics-Simulation and Computation, 46(10), 8077-8085.
#' 
#' Özdemir, Y. A., Ebegil, M., & Gökpinar, F. (2019). A test statistic for two normal means with median ranked set sampling. Iranian Journal of Science and Technology, Transactions A: Science, 43(3), 1109-1126.
#' 

datagen_RSS=function (mu,s,m,r){
dx=c()
for (i in 1:r) {
(xs=t(apply(matrix(rnorm(m^2,mu,s),m),1,sort)))
dx=c( dx, diag(xs))
}
x=t(matrix(dx,ncol=r))
rownames(x)=paste("Cycle",1:r)
return(x)
}