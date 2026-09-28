#rm(list=ls())
#library(foreach)
#library(parallel)
#library(doParallel)
#library(bigstatsr)


###################################################################
#### part II: component function for conditional probabilities ####
###         and adaptive procedures                            ####
###################################################################

####################################
##############################################
#' Group sequential designs / Adaptive sequential designs: Final analysis, With one sample size change
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param snew Vector of new information fractions after the sample size is changed
#' @param c_bry_new Critical boundaries after the sample size is changed
#' @param theta_hat_inter Treatment effect estimate observed at the interim analysis
#' @param theta_hat_inter_sd Standard error of 'theta_hat_inter'
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based)
#' @param theta_hat_last_ad Treatment effect estimate observed at the last analysis after the sample size change
#' @param theta_hat_last_ad_sd Standard error of 'theta_hat_last_ad'
#' @param last_vt_ad Index of the last analysis after the sample size change
#' @param alpha_2 Type I error for the two-sided test
#'
#' @return A named list containing the following elements: 'p-value', 'confidence interval'.
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_ci_est_back(s = c(0.3,0.7,1), c_bry = c(4.347011, 2.845787, 2.380956), snew = c(0.5,1),
#' c_bry_new = c(3.310743, 2.781890), theta_hat_inter = 1.2, theta_hat_inter_sd = 1, inter_vt = 1,
#' theta_hat_last_ad = 2.2, theta_hat_last_ad_sd = 1, last_vt_ad = 2, alpha_2 = 0.05)
#'
#' @export
asd_ci_est_back<-function(s,c_bry,snew,c_bry_new,theta_hat_inter,theta_hat_inter_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,alpha_2){
  #this function calculates point estimatre and confidence interval after sample size change
  #s is the infromation fractions for the original design
  #c_bry is the critical boundaries with the original design
  #snew is the vector of information fractions after sample size change #
  #c_bry_new is the criotical boundaries after sample size change##
  #inter_vt is the interim analysis at twhich the sample size was changed
  #theta_hat_inter is the theta_hat at inter_vt'th interim analysis
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #theta_hat_last_ad the theta_hat at the last interim analysis after sample size change
  #theta_hat_last_ad_sd is the standard error of theta_hat_last_ad
  #last_vt_ad is the last interim/final analysis after sample size change
  #alpha_2 is for calcularing the two-sided (1-alpha_2)*100% ci
  #typically, alpha_2=0.05
  T_inter=theta_hat_inter_sd^(-2)
  T=T_inter/s[inter_vt]
  t=T*s
  T_last_ad=theta_hat_last_ad_sd^(-2)
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  #T is the total information time for the original design
  #T_new is the new total information time with sampe size change

  ############################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(s,length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  ####################################################
  #######################

  cal<-function(theta){
    c_p=asd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad)
    bk_img=asd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p)
    last_vt_theta=length(bk_img)
    bkimg_wald_theta=last(bk_img)
    ### Calculating exit probability with lst_wald_theta,lst_vt_theta###
    last_vt=last_vt_theta

    ####################Boundary for vector of Brownian  motions#########################
    if (last_vt==1) {
      u_bry_theta=bkimg_wald_theta-theta*sqrt(T*s[last_vt])
    }else{
      u_c_bry=c(c_bry[1:(last_vt-1)],bkimg_wald_theta)
      u_bry_theta=u_c_bry*sqrt(s[1:last_vt])-theta*sqrt(T)*s[1:last_vt]
    }

    ###################################################
    sigma0=cov_mtxf=cov_mtx[1:last_vt,1:last_vt]
    if (last_vt==1){
      tmp = 1-stats::pnorm(u_bry_theta)
    }else {
      tmp = 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,last_vt),sigma=sigma0,algorithm=Miwa())[1]
    }
    return(tmp)
  }

  ##########################################################
  theta_ll=as.numeric(stats::uniroot(function(x) cal(x)-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal(x)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul=as.numeric(stats::uniroot(function(x) cal(x)-1+alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  p_value= cal(0)
  out=list(p_value,c(theta_ll,theta_pt,theta_ul))
  names(out)=c("p-value","confidence interval")
  return(out)
  #theta_pt is point estimate
}
#############################################
#asd_ci_est_back(s,c_bry,snew,c_bry_new,theta_hat_inter,theta_hat_inter_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,alpha_2)
#asd_ci_est_back(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=c(0.5,1),c_bry_new=c(3.310743, 2.781890),theta_hat_inter=1.2,theta_hat_inter_sd=1,inter_vt=1,theta_hat_last_ad=2.2,theta_hat_last_ad_sd=1,last_vt_ad=2,alpha_2=0.05)
#asd_ci_est_back(s=c(0.5,1),c_bry=c(5, 1.96),snew=1,c_bry_new=5,theta_hat_inter=1.83,theta_hat_inter_sd=1.3,inter_vt=1,theta_hat_last_ad=2.38,theta_hat_last_ad_sd=0.765,last_vt_ad=1,alpha_2=0.05)
#asd_ci_est_back(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=c(0.5,1),c_bry_new=c(3.310743, 2.781890),theta_hat_inter=1.2,theta_hat_inter_sd=1,inter_vt=1,theta_hat_last_ad=2.2,theta_hat_last_ad_sd=0.8,last_vt_ad=2,alpha_2=0.05)

##########################################################################
