#rm(list=ls())
#library(dplyr)
#library(foreach)
#library(parallel)
#library(doParallel)
#library(bigstatsr)

###################################################
##########   Part I                    ############
### Component functions, to be sourced ############
###################################################

##########Final analysi, without sample size change ###########

##############################################
###############################################################

#' Phase 2/3 seamless combination Group sequential design: Final analysis, Without sample size change, last_vt == 1
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param theta_hat_0 Vector of naive point estimates of the treatment effect for the comparisons at the
#'   end of the phase 2 stage (the first interim analysis); its length is the number of comparisons
#' @param theta_hat_0_sd Common standard error of 'theta_hat_0' at the end of the phase 2 stage, shared by
#'   all comparisons; must be > 0
#' @param alpha_2 Type I error for the two-sided test
#' @param dist Distribution indicator of the endpoint: 1 = the endpoint is normally distributed and the
#'   comparison is between multiple doses and a common control, otherwise the non-normal branch is used
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{p-value}{The p-value for the hypothesis test}
#'   \item{confidence interval}{Numeric vector of length 3 giving the lower confidence limit, the median
#'     unbiased point estimate and the upper confidence limit; the lower limit is consistent with the
#'     p-value}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_ci_0(s = c(1/3, 2/3, 1), theta_hat_0 = 0.2, theta_hat_0_sd = 0.5,
#' alpha_2 = 0.05, dist = 1)
#'
#' @export
two_stage_ci_0<-function(s,theta_hat_0,theta_hat_0_sd,alpha_2,dist){
  #This function is used if the null hypothesis is rejected (c_bry crossed )
  #at the end of phase 2 stage#
  #This function provides the point estimate , lower and upper confidence limits for the
  #slected comparsion (mu_selected,i.e., with the largest score value)
  #s is the vector of information fractions#
  #c_bry is the critical boundary for the hypothesis testing
  #alpha_2 is the significance level of the two-sided (1-alpha_2)*100% CI .
  #theta_hat_0 is the vector of point estimates at the end of the phase 2 stage (first interim analysis)##
  #theta_hat_0_sd is the common standard deviations for theta_hat_0  at the end of phase 2 stage###
  #### the next section constructs the covariance matrix###

  ####setting up covaraince matrices#######
  comp=length(theta_hat_0)
  t_0=(mean(theta_hat_0_sd))^(-2)
  s_0=s[1] #information fraction for the end of pahse 2 stage analysis#
  T=t_0/s_0
  ########## non-normal distribution- phase 2 stage matrix###
  cov_mtx_0_nonn=diag(comp)
  cov_mtx_0_nonn=s_0*cov_mtx_0_nonn
  ########## normal distribution- phase 2 stage matrix###
  cov_mtx_0_n=diag(comp)
  cov_mtx_0_n[lower.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n[upper.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n=s_0*cov_mtx_0_n
  ############################## calculating the Wald statistics and scores ###
  wald_0=theta_hat_0/theta_hat_0_sd
  wald_0_max=max(theta_hat_0/theta_hat_0_sd)

  #######calculating exit boundaries####

  ext_bry=wald_0_max*sqrt(s[1])
  # ext_bry is the exit boundary at the end of phase 2 stage

  ##############################

  tp00 <-function(theta,sigma0){

    rs= 1-pmvnorm(lower=-Inf, upper=rep(ext_bry-theta*sqrt(T)*s_0,comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]
  }


  theta_ulc=max(theta_hat_0)+stats::qnorm(1-alpha_2/2)/sqrt(T*s[1])
  if (dist==1){
    theta_llc=as.numeric(stats::uniroot(function(x) tp00(x, cov_mtx_0_n)-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
    theta_point=as.numeric(stats::uniroot(function(x) tp00(x, cov_mtx_0_n)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    p_value=tp00(0, cov_mtx_0_n)
  }else {
    theta_llc=as.numeric(stats::uniroot(function(x) tp00(x, cov_mtx_0_nonn)-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
    theta_point=as.numeric(stats::uniroot(function(x) tp00(x, cov_mtx_0_nonn)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    p_value=tp00(0, cov_mtx_0_nonn)
  }


  theta_ci=c(theta_llc,theta_point,theta_ulc)

  out=list(p_value, theta_ci)
  names(out)=c("p-value","confidence interval")
  return(out)
  #theta_pt is point estimate
  #(theta_ll,theta_ul) is the confidence interval
}
##########################################
#two_stage_ci_0(s,theta_hat_0,theta_hat_0_sd,alpha_2,dist)
#two_stage_ci_0(s=c(1/3,2/3,1),theta_hat_0=c(0.3,0.5),theta_hat_0_sd=0.6,alpha_2=0.05,dist=1)
#two_stage_ci_0(s=c(1/3,2/3,1),theta_hat_0=c(0.1,-0.1,0.3),c(0.23,0.2,0.24),alpha_2=0.05,dist=1)
#two_stage_ci_0(s=c(1/3,2/3,1),theta_hat_0=c(0.1,-0.1,0.3),c(0.23,0.2,0.24),alpha_2=0.05,dist=2)
###############################################################


###########################################
### The function Two_stage_est_1 provides point estimate and the confidence intervals, and p-value ##
### after a trial ended in the phase 3 stage (i.e. at s[2], s[3],....)###

#########################################
#' Phase 2/3 seamless combination Group sequential design: Final analysis, Without sample size change, last_vt > 1
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be
#'   1; here the length of 's' must be 2 or 3, i.e. the phase 2 stage plus one or two phase 3 analyses
#' @param c_bry_two_stg Vector of critical boundaries of the two-stage design for the overall test, one per
#'   analysis; must have the same length as 's'
#' @param c_bry_p3 Vector of critical boundaries of the prior (phase 3) study, one per phase 3 analysis;
#'   must have 'length(s) - 1' elements
#' @param comp Number of comparisons in the phase 2 stage, e.g. the number of doses or the number of
#'   endpoints; must be at least 2
#' @param theta_hat_p3 Treatment effect estimate observed at the analysis at which the trial stopped, using
#'   phase 3 data only
#' @param theta_hat_p3_sd Standard error of 'theta_hat_p3'; must be > 0
#' @param theta_hat_last Treatment effect estimate observed at the analysis at which the trial stopped,
#'   based on all data accumulated by then
#' @param theta_hat_last_sd Standard error of 'theta_hat_last'; must be > 0
#' @param last_vt Index of the analysis at which the trial stopped; must be > 1 because the trial stopped in
#'   the phase 3 stage
#' @param alpha_2 Type I error for the two-sided test
#' @param dist Distribution indicator of the endpoint: 1 = the endpoint is normally distributed and the
#'   comparison is between multiple doses and a common control, otherwise the non-normal branch is used
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{p-value}{The p-value for the hypothesis test}
#'   \item{confidence interval}{Numeric vector of length 3 giving the lower confidence limit, the median
#'     unbiased point estimate and the upper confidence limit; the lower limit is consistent with the
#'     p-value, while the point estimate and the upper limit use phase 3 data only}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_ci_1(s = c(1/3,2/3,1), c_bry_two_stg = c(3.938325, 2.784816,2.273793),
#' c_bry_p3 = c(3.954862,1.977431), comp = 3, theta_hat_p3 = 0.45, theta_hat_p3_sd = 0.6,
#' theta_hat_last = 2.784816*0.12, theta_hat_last_sd = 0.12, last_vt = 2, alpha_2 = 0.05, dist = 1)
#'
#' @export
two_stage_ci_1<-function(s,c_bry_two_stg,c_bry_p3,comp,theta_hat_p3,theta_hat_p3_sd,theta_hat_last,theta_hat_last_sd,last_vt,alpha_2,dist){
  #this function calculates confidence interval if the trial ended in the phase 3 stage#
  # the lower limit of the CI is consistent with the hypothesis test
  #the point estimate and the upper limit of CI only use phase 3 data.
  # the point estimate is median ubiased
  # the upper limit of CI is exact
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #alpha_2 is the significance level of the two-sided (1-alpha_2)*100% CI
  #comp is the number of comparisons, it can be for dose, or endpoints#
  #c_bry_two_stg is the critical boundary for the overall test
  #c_bry_p3 is the critical boundary for pahse 3 stage
  #this function calculates the p-value is the trial ended in the phase 3 stage#
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #c_bry_two_stg is the critical boundary for the overall test
  #theta_hat_last is the theta_hat at the last analysis
  #theta_hat_last_sd is the standard error of theta_hat_last
  # theta_hat_p3 is the theta_hat(only from phase 3 data) at the last analysis
  #theta_hat_p3_sd is the standard error of theta_hat_p3
  #comp is the number of comparisons, it can be for dose, or endpoints#
  #last_vt numbers the interim/final analysis at which the trial ended
  #dist is the distribution of the variable. dist==1 means that the variable is normally distributed and the comparison is for multiple dose comparison#

  t_last=theta_hat_last_sd^(-2)
  T=t_last/s[last_vt]
  s_0=s[1] #information fraction at the end of phase 2 stage
  r=(s[-1]-s[1])
  #information fraction conditonal on phase 2 sthage
  s_1=(s[-1]-s[1])/(last(s)-s[1])
  ## s_1 is the scaled information fraction for selection the ##
  ## critical boundary for phase 3 stage, this boundary will be only used for point estiante and upper limit of confiedence intgerval##
  ######################## exit boundaries##########
  e_bry_two_stg=c_bry_two_stg*sqrt(s)
  e_bry_p3=c_bry_p3*sqrt(s_1)
  ########## non-normal distribution- phase 2 stage matrix###
  cov_mtx_0_nonn=diag(comp)
  cov_mtx_0_nonn=s_0*cov_mtx_0_nonn ### cov_mtx_0_nonn=sigma0

  ##############non-normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_nonn=diag(comp-1)
  cov_mtx_0c_nonn=s_0*cov_mtx_01_nonn #cov_mtx_0c_nonn=sigma1

  ########## normal distribution- phase 2 stage matrix###
  cov_mtx_0_n=diag(comp)
  cov_mtx_0_n[lower.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n[upper.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n=s_0*cov_mtx_0_n  ### cov_mtx_0_n=sigma0

  ##############normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_n=diag(comp-1)
  cov_mtx_01_n[lower.tri(cov_mtx_01_n)]=0.5
  cov_mtx_01_n[upper.tri(cov_mtx_01_n)]=0.5

  a_mtx_01=(1/4)*matrix(rep(1,(comp-1)*(comp-1)),(comp-1),(comp-1))

  con_mtx_01_n=cov_mtx_01_n-a_mtx_01
  cov_mtx_0c_n=s_0*con_mtx_01_n #cov_mtx_0c_n=sigma1

  ###### Matirx for phase 3 stage #######
  cov_mtx_p3= matrix(rep(c(r),length(r)),nrow=length(r), ncol=length(r))
  cov_mtx_p3[lower.tri(cov_mtx_p3)]=t(cov_mtx_p3)[lower.tri(cov_mtx_p3)]

    #cov_mtx_p3=sigma2 to be used for integrations involving both phase 2 anf phas3 stage data
  #######################Covariance matrix for phase 3 ##########
  cov_mtx_p3_ind= matrix(rep(c(s_1),length(s_1)),nrow=length(s_1), ncol=length(s_1))
  cov_mtx_p3_ind[lower.tri( cov_mtx_p3_ind)]=t( cov_mtx_p3_ind)[lower.tri( cov_mtx_p3_ind)]
  ####  cov_mtx_p3_ind=sigma3 is for point estiate only using phase 3 data, and the phase 3 stage critical boundaries####

  ############################## calculating the Wald statistics and scores ###

  wald_last=theta_hat_last/theta_hat_last_sd
  score_last=wald_last/theta_hat_last_sd
  wald_last_p3=theta_hat_p3/theta_hat_p3_sd
  score_last_p3=wald_last_p3/theta_hat_p3_sd

  #######calculating exit boundaries -- for combined phase estimation####

  ext_bry_tail=c(e_bry_two_stg[1:(last_vt-1)],max(wald_last)*sqrt(s[last_vt]))

  ###### exit boundary for estimtes using pase 3 data only###

  if (last_vt==2){
    e_bry_p3_tail=wald_last_p3*sqrt(s_1[1])
  }
  if (last_vt>2){
    e_bry_p3_tail=c(e_bry_p3[1:(last_vt-2)],wald_last_p3*sqrt(s_1[last_vt-1]))
  }
  ##############################

  ### Define functions for estimation ################
  ############## tp1_norm calcualtes exit probability for normal distributions compared to a common control ######
  tp1_norm <-function(theta,ext_bry_tail,sigma0,sigma1,sigma2){
    ####sigma0 is the covariance matrix for all the comparisons in phase 2##
    ### sigma1 is the conditional covaraince matrix for phase 2
    #### sigma2 is the covariance matrix for phase 3 srtage only###

    intg <- function(x){
      up_bry_0_cond=rep(x,comp-1) ## up_bry_0_cond is for conditional probability integration for phase 2 stage
      up_bry_p3=(ext_bry_tail[-1]-x-theta*sqrt(T)*sqrt(s[2:last_vt])) #up_bry_p3 is the upper boundary for the phase 3 integration
      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0)) * pmvnorm(lower=-Inf,upper=up_bry_0_cond,mean=rep(1/2*x,comp-1),sigma=sigma1,algorithm=Miwa())[1] *
        (1-pmvnorm(lower=-Inf,upper= up_bry_p3,mean=rep(0,last_vt-1),sigma=sigma2,algorithm=Miwa())[1])

      return(rs)
    }
    up_bry_0_all=rep( e_bry_two_stg[1]-theta*sqrt(T)*s_0,comp) ##up_bry_0_all is the upper limit for all comparisons for the phase 2 stage integration
    up_bry_0_single= e_bry_two_stg[1]-theta*sqrt(T)*s_0  ## up_bry_0_single is the upper limit for the selected comparison in the phase 2 stage
    rs= 1 - pmvnorm(lower=-Inf, upper=up_bry_0_all,mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1] +
        comp* stats::integrate(Vectorize(intg),lower= -Inf,upper =up_bry_0_single,abs.tol=1e-5)$value

    return(rs)
  }

  ############## tp1_nnorm calcualtes exit probability for non-normal distribiutions, or normal distributions compared to independent controls ######
  tp1_nnorm <-function(theta,ext_bry_tail,sigma0,sigma1,sigma2){

    intg <- function(x){
      up_bry_0_cond=rep(x,comp-1) ## up_bry_0_cond is for conditional probability integration for phase 2 stage
      up_bry_p3=(ext_bry_tail[-1]-x-theta*sqrt(T)*sqrt(s[2:last_vt])) #up_bry_p3 is the upper boundary for the phase 3 integration
      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=up_bry_0_cond,mean=rep(0,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
        (1-pmvnorm(lower=-Inf,upper= up_bry_p3,mean=rep(0,last_vt-1),sigma=sigma2,algorithm=Miwa())[1]  )

      return(rs)
    }
    up_bry_0_all=rep( e_bry_two_stg[1]-theta*sqrt(T)*s_0,comp) ##up_bry_0_all is the upper limit for all comparisons for the phase 2 stage integration
    up_bry_0_single= e_bry_two_stg[1]-theta*sqrt(T)*s_0  ## up_bry_0_single is the upper limit for the selected comparison in the phase 2 stage
    rs= 1-pmvnorm(lower=-Inf, upper=up_bry_0_all,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]+
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper =up_bry_0_single,abs.tol=1e-5)$value

    return(rs)
  }

  ###############################################

  tp_p3 <-function(theta,sigma3){ #this function obtains unbised point estimate using phase stage data only
    up_bry_p3_est=e_bry_p3_tail-theta*sqrt(T*(1-s[1]))*s_1[1:(last_vt-1)] #up_bry_p3_est is for unbiased point estimate using only phase 3 data
    rs=1-pmvnorm(lower=-Inf,upper=up_bry_p3_est,mean=rep(0,last_vt-1),sigma=sigma3,algorithm=Miwa())[1]
    return(rs)
  }

  ############################## Estimates ##########
  if (dist==1) {
    theta_llc=as.numeric(stats::uniroot(function(x) tp1_norm (x,ext_bry_tail,cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_p3[1:(last_vt-1),1:(last_vt-1)])-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
    p_value=tp1_norm (0,ext_bry_tail,cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_p3[1:(last_vt-1),1:(last_vt-1)])
    sigma3=cov_mtx_p3_ind[1:(last_vt-1),1:(last_vt-1)]
    theta2_ulc=as.numeric(stats::uniroot(function(x) tp_p3(x,sigma3)-1+alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
    theta2_point=as.numeric(stats::uniroot(function(x) tp_p3(x,sigma3)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    #theta_ci=c(theta_llc, theta2_point,theta2_ulc)
    theta_ci=c( theta_llc,theta2_point,theta2_ulc)
  }else{
    theta_llc=as.numeric(stats::uniroot(function(x) tp1_nnorm(x,ext_bry_tail,cov_mtx_0_nonn, cov_mtx_0c_nonn,cov_mtx_p3[1:(last_vt-1),1:(last_vt-1)])-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
    p_value=tp1_nnorm(0,ext_bry_tail,cov_mtx_0_nonn, cov_mtx_0c_nonn,cov_mtx_p3[1:(last_vt-1),1:(last_vt-1)])
    sigma3=cov_mtx_p3_ind[1:(last_vt-1),1:(last_vt-1)]
    theta2_ulc=as.numeric(stats::uniroot(function(x) tp_p3(x,sigma3)-1+alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
    theta2_point=as.numeric(stats::uniroot(function(x) tp_p3(x,sigma3)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    #theta_ci=c(theta_llc, theta2_point,theta2_ulc)
    theta_ci=c(theta_llc, theta2_point,theta2_ulc)
  }

  out=list(p_value, theta_ci)
  names(out)=c("p-value","confidence interval")
  return(out)
  #theta_pt is point estimate
  #(theta_ll,theta_ul) is the confidence interval
}
#########################################
#two_stage_ci_1(s,c_bry_two_stg,c_bry_p3,comp,theta_hat_p3,theta_hat_p3_sd,theta_hat_last,theta_hat_last_sd,last_vt,alpha_2,dist)
#two_stage_ci_1(s=c(1/3,2/3,1),c_bry_two_stg=c(3.938325, 2.784816,2.273793),c_bry_p3=c(3.954862,1.977431),comp=3,theta_hat_p3=0.45,theta_hat_p3_sd=0.6,theta_hat_last=2.784816*0.12,theta_hat_last_sd=0.12,last_vt=2,alpha_2=0.05,dist=1)
#two_stage_ci_1(s=c(1/3,2/3,1),c_bry_two_stg=c(3.938325, 2.784816,2.273793),c_bry_p3=c(3.954862,1.977431),comp=3,theta_hat_p3=0.45,theta_hat_p3_sd=0.6,theta_hat_last=2.273793*0.12,theta_hat_last_sd=0.12,last_vt=3,alpha_2=0.05,dist=1)
#two_stage_ci_1(s=c(1/3,2/3,1),c_bry_two_stg=c( 4.057782,2.869285,2.342762),c_bry_p3=c(3.954862,1.977431),comp=3,theta_hat_p3=0.45,theta_hat_p3_sd=0.6,theta_hat_last=2.342762*0.12,theta_hat_last_sd=0.12,last_vt=2,alpha_2=0.05,dist=2)
#two_stage_ci_1(s=c(1/3,2/3,1),c_bry_two_stg=c( 4.057782,2.869285,2.342762),c_bry_p3=c(3.954862,1.977431),comp=3,theta_hat_p3=0.45,theta_hat_p3_sd=0.6,theta_hat_last=2.342762*0.12,theta_hat_last_sd=0.12,last_vt=3,alpha_2=0.05,dist=2)
#two_stage_ci_1(s=c(1/3,1),c_bry_two_stg=c( 3.382574,2.230020),c_bry_p3=c(3.954862,1.977431),comp=3,theta_hat_p3=0.45,theta_hat_p3_sd=0.6,theta_hat_last=2.342762*0.12,theta_hat_last_sd=0.12,last_vt=2,alpha_2=0.05,dist=2)

###################################################################
#### part II: component function for conditional probabilities ####
###         and adaptive procedures                            ####
###################################################################


##############################################
two_stage_llc_back<-function(s,c_bry_two_stg,comp,snew,c_bry_new,theta_hat_inter,theta_hat_inter_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,alpha_2,dist){
  #this function calculates the lower confidence limit after sample size change
  #s is the infromation fractions for the original design
  #c_bry_two_stg is the critical bpundaries with the original design
  # comp is the number of comparisons, such as the number of doses, or number of endpoint
  # snew is the new relative information fraction for the segment of trial after sample size change
  #inter_vt numbers the interim analysis at which the sample size was changed
  #theta_hat_inter is the theta_hat at inter_vt'th interim analysis
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #theta_hat_last_ad  is the theta_hat at the final analysia after sample size change
  # theta_hat_last_ad_sd is the standard error of theta_hat_last_ad
  # last_vt_ad numbers the last analysis where the trial ended
  # alpha_2 is the significance level of the two-sided (1-alpha_2)*100% CI
  #typically, alpha_2=0.05
  # dist=1 if the trial is for dose comparison, with normal distribution, and a common contorl. Otherwisw, dist=2
  ##############
  T_inter=theta_hat_inter_sd^(-2)
  T=T_inter/s[inter_vt]
  t=T*s
  T_last_ad=theta_hat_last_ad_sd^(-2)
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  #T is the total information time for the original design
  #T_new is the new total information time with sampe size change
  ###############################################
  c_bry=c_bry_two_stg
  s_0=s[1] #information fraction at the end of phase 2 stage
  r=(s[-1]-s[1])
  #information fraction conditonal on phase 2 sthage
  ######################## exit boundaries##########
  e_bry_two_stg=c_bry_two_stg*sqrt(s)
  ########## non-normal distribution- phase 2 stage matrix###
  cov_mtx_0_nonn=diag(comp)
  cov_mtx_0_nonn=s_0*cov_mtx_0_nonn ### cov_mtx_0_nonn=sigma0 if dist=2

  ##############non-normal distribution- phase 2 stage, conditional distribution matrix###
  cov_mtx_01_nonn=diag(comp-1)
  cov_mtx_0c_nonn=s_0*cov_mtx_01_nonn #cov_mtx_0c_nonn=sigma1 if dist=2
  ########## normal distribution- phase 2 stage matrix###
  cov_mtx_0_n=diag(comp)
  cov_mtx_0_n[lower.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n[upper.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n=s_0*cov_mtx_0_n  ### cov_mtx_0_n=sigma0 if dist=1

  ##############normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_n=diag(comp-1)
  cov_mtx_01_n[lower.tri(cov_mtx_01_n)]=0.5
  cov_mtx_01_n[upper.tri(cov_mtx_01_n)]=0.5

  a_mtx_01=(1/4)*matrix(rep(1,(comp-1)*(comp-1)),(comp-1),(comp-1))

  con_mtx_01_n=cov_mtx_01_n-a_mtx_01
  cov_mtx_0c_n=s_0*con_mtx_01_n #cov_mtx_0c_n=sigma1 if dist=1

  ###### Matirx for phase 3 stage #######
  cov_mtx_p3= matrix(rep(c(r),length(r)),nrow=length(r), ncol=length(r))
  cov_mtx_p3[lower.tri(cov_mtx_p3)]=t(cov_mtx_p3)[lower.tri(cov_mtx_p3)]
  #cov_mtx_p3=sigma2 to be used for integrations involving both phase 2 anf phas3 stage data
  #######################Covariance matrix for phase 3 ##########
  ####################################################
  ### Define functions for estimation ################
  ####################################################

  ############## cal_norm calcualtes exit probability for normal distributions compared to a common control ######
  cal_norm<-function(theta,sigma0,sigma1,sigma2){

    #asd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad)
    c_p=asd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad)
    #asd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p)
    bk_img=asd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p)
    last_vt_theta=length(bk_img)
    bkimg_wald_theta=last(bk_img)
    ### Calculating exit probability with lst_wald_theta,lst_vt_theta###
    last_vt=last_vt_theta

    wald_last=bkimg_wald_theta
    ####################################
    #######calculating exit boundaries -- for combined phase estimation####
    ext_bry_tail=c(e_bry_two_stg[1:(last_vt-1)],wald_last*sqrt(s[last_vt]))
    ##############################
    intg <- function(x){
      up_bry_0_cond=rep(x,comp-1) ## up_bry_0_cond is for conditional probability integration for phase 2 stage
      up_bry_p3=(ext_bry_tail[-1]-x-theta*sqrt(T)*sqrt(s[2:last_vt])) #up_bry_p3 is the upper boundary for the phase 3 integration
      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0)) * pmvnorm(lower=-Inf,upper=up_bry_0_cond,mean=rep(1/2*x,comp-1),sigma=sigma1,algorithm=Miwa())[1] *
        (1-pmvnorm(lower=-Inf,upper= up_bry_p3,mean=rep(0,last_vt-1),sigma=sigma2[1:(last_vt-1),1:(last_vt-1)],algorithm=Miwa())[1] )

      return(rs)
    }
    up_bry_0_all=rep( e_bry_two_stg[1]-theta*sqrt(T)*s_0,comp) ##up_bry_0_all is the upper limit for all comparisons for the phase 2 stage integration
    up_bry_0_single= e_bry_two_stg[1]-theta*sqrt(T)*s_0  ## up_bry_0_single is the upper limit for the selected comparison in the phase 2 stage
    rs= 1 - pmvnorm(lower=-Inf, upper=up_bry_0_all,mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1] +
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper =up_bry_0_single,abs.tol=1e-5)$value
   return(rs)
  }
  ###################################################
  ############## cal_non_norm calcualtes exit probability for if the trial does not compare normal distributions to a common control ######
  cal_non_norm<-function(theta,sigma0,sigma1,sigma2){
    #asd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad)
    c_p=asd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad)
    #asd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p)
    bk_img=asd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p)
    last_vt_theta=length(bk_img)
    bkimg_wald_theta=last(bk_img)
    ### Calculating exit probability with lst_wald_theta,lst_vt_theta###
    last_vt=last_vt_theta

    wald_last=bkimg_wald_theta
    ####################################
    #######calculating exit boundaries -- for combined phase estimation####
    ext_bry_tail=c(e_bry_two_stg[1:(last_vt-1)],wald_last*sqrt(s[last_vt]))
    ##############################
    intg <- function(x){
      up_bry_0_cond=rep(x,comp-1) ## up_bry_0_cond is for conditional probability integration for phase 2 stage
      up_bry_p3=(ext_bry_tail[-1]-0-theta*sqrt(T)*sqrt(s[2:last_vt])) #up_bry_p3 is the upper boundary for the phase 3 integration
      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0)) * pmvnorm(lower=-Inf,upper=up_bry_0_cond,mean=rep(0,comp-1),sigma=sigma1,algorithm=Miwa())[1] *
         (1-pmvnorm(lower=-Inf,upper= up_bry_p3,mean=rep(0,last_vt-1),sigma=sigma2[1:(last_vt-1),1:(last_vt-1)],algorithm=Miwa())[1])

      return(rs)
    }
    up_bry_0_all=rep( e_bry_two_stg[1]-theta*sqrt(T)*s_0,comp) ##up_bry_0_all is the upper limit for all comparisons for the phase 2 stage integration
    up_bry_0_single= e_bry_two_stg[1]-theta*sqrt(T)*s_0  ## up_bry_0_single is the upper limit for the selected comparison in the phase 2 stage
    rs= 1 - pmvnorm(lower=-Inf, upper=up_bry_0_all,mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1] +
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper =up_bry_0_single,abs.tol=1e-5)$value
    return(rs)
  }
  ###################################################

  ############## finished defining functions #########
  ##############################Obtain  Estimates ##########
  if (dist==1) {
    theta_llc=as.numeric(stats::uniroot(function(x) cal_norm(x,cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_p3)-alpha_2/2,c(-2^8,2^8),tol=1e-5)[1])
    p_value=cal_norm(0,cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_p3)
  }else{
   theta_llc=as.numeric(stats::uniroot(function(x) cal_non_norm(x,cov_mtx_0_nonn, cov_mtx_0c_nonn,cov_mtx_p3)-alpha_2/2,c(-2^8,2^8),tol=1e-5)[1])
   p_value=cal_non_norm(0,cov_mtx_0_nonn, cov_mtx_0c_nonn,cov_mtx_p3)
  }
  return(c(p_value,theta_llc))
  #lower limit of confidence interval
}
#############################################
#two_stage_llc_back(s,c_bry_two_stg,comp,snew,c_bry_new,theta_hat_inter,theta_hat_inter_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,alpha_2,dist)
#two_stage_llc_back(s=c(1/3,2/3,1),c_bry_two_stg=c(3.938325, 2.784816,2.273793),comp=3,snew=c(0.5,1),c_bry_new=c(3.310743, 2.781890),theta_hat_inter=0.45,theta_hat_inter_sd=0.6,inter_vt=2,theta_hat_last_ad=2.784816*0.12,theta_hat_last_ad_sd=0.12,last_vt_ad=2,alpha_2=0.05,dist=1)
#two_stage_llc_back(s=c(1/3,2/3,1),c_bry_two_stg=c(3.938325, 2.784816,2.273793),comp=3,snew=c(0.5,1),c_bry_new=c(3.310743, 2.781890),theta_hat_inter=0.45,theta_hat_inter_sd=0.6,inter_vt=2,theta_hat_last_ad=3.025*0.12,theta_hat_last_ad_sd=0.12,last_vt_ad=2,alpha_2=0.05,dist=1)
##############################################

#########################################################################
ph3_ci_est_back<-function(s_p3,c_bry_p3,snew,c_bry_new_p3,theta_hat_inter_p3,theta_hat_inter_p3_sd,inter_vt_p3,theta_hat_last_ad_p3,theta_hat_last_ad_p3_sd,last_vt_ad,alpha_2){
  s=s_p3
  c_bry=c_bry_p3
  snew=snew
  c_bry_new=c_bry_new_p3
  theta_hat_inter=theta_hat_inter_p3
  theta_hat_inter_sd=theta_hat_inter_p3_sd
  inter_vt=inter_vt_p3
  theta_hat_last_ad=theta_hat_last_ad_p3
  theta_hat_last_ad_sd=theta_hat_last_ad_p3_sd

  p3_out_1=asd_ci_est_back(s,c_bry,snew,c_bry_new,theta_hat_inter,theta_hat_inter_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,alpha_2)[2]
  p3_out=as.numeric(unlist(p3_out_1))
  return(p3_out)
}
#########################################
#ph3_ci_est_back(s_p3,c_bry_p3,snew,c_bry_new_p3,theta_hat_inter_p3,theta_hat_inter_p3_sd,inter_vt_p3,theta_hat_last_ad_p3,theta_hat_last_ad_p3_sd,last_vt_ad,alpha_2){
#ph3_ci_est_back(s_p3=c(0.5,1),c_bry_p3=c(3.954862,1.977431),snew=1,c_bry_new_p3=3.31,theta_hat_inter_p3=0.37,theta_hat_inter_p3_sd=0.17,inter_vt_p3=1,theta_hat_last_ad_p3=0.365,theta_hat_last_ad_p3_sd=0.15,last_vt_ad=1,alpha_2=0.05)

######################################################################

#' Phase 2/3 seamless combination Group sequential design: Final analysis, With one sample size change
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry_two_stg Vector of critical boundaries of the two-stage design for the overall test, one per
#'   analysis; must have the same length as 's'
#' @param c_bry_p3 Vector of critical boundaries of the prior (phase 3) study, one per phase 3 analysis;
#'   must have 'length(s) - 1' elements
#' @param comp Number of comparisons in the phase 2 stage, e.g. the number of doses or the number of
#'   endpoints; must be at least 2
#' @param snew Vector of new information fractions after the sample size is changed, one per analysis after
#'   the change; must have at least 'last_vt_ad' elements
#' @param c_bry_new Vector of critical boundaries after the sample size is changed, one per analysis after
#'   the change; must have the same length as 'snew'
#' @param c_bry_new_p3 Vector of critical boundaries after the sample size is changed on the prior
#'   (phase 3) study scale, one per analysis after the change; must have the same length as 'snew'
#' @param theta_hat_inter Treatment effect estimate observed at the interim analysis at which the sample
#'   size is changed
#' @param theta_hat_inter_sd Standard error of 'theta_hat_inter'; must be > 0
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based); must be
#'   less than 'length(s)', so that at least one further analysis remains
#' @param theta_hat_last_ad Treatment effect estimate observed at the last analysis after the sample size
#'   change
#' @param theta_hat_last_ad_sd Standard error of 'theta_hat_last_ad'; must be > 0
#' @param last_vt_ad Index of the last analysis after the sample size change
#' @param theta_hat_inter_p3 Treatment effect estimate observed at the interim analysis at which the sample
#'   size is changed, using phase 3 data only
#' @param theta_hat_inter_p3_sd Standard error of 'theta_hat_inter_p3'; must be > 0
#' @param theta_hat_last_ad_p3 Treatment effect estimate observed at the last analysis after the sample size
#'   change, using phase 3 data only
#' @param theta_hat_last_ad_p3_sd Standard error of 'theta_hat_last_ad_p3'; must be > 0
#' @param alpha_2 Type I error for the two-sided test
#' @param dist Distribution indicator of the endpoint: 1 = the endpoint is normally distributed and the
#'   comparison is between multiple doses and a common control, otherwise the non-normal branch is used
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{p-value}{The p-value for the hypothesis test}
#'   \item{confidence interval}{Numeric vector of length 3 giving the lower confidence limit, the median
#'     unbiased point estimate and the upper confidence limit; the lower limit is consistent with the
#'     p-value, while the point estimate and the upper limit use phase 3 data only}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_est_back(s=c(1/3,2/3,1), c_bry_two_stg=c(3.776606,2.670463,2.180424),
#' c_bry_p3=c(3.954862,1.977431), comp=2, snew=c(0.5,1), c_bry_new=c(3,2.2),
#' c_bry_new_p3=c(3.2,2.4), theta_hat_inter=0.4, theta_hat_inter_sd=0.35, inter_vt=1,
#' theta_hat_last_ad=0.5, theta_hat_last_ad_sd=0.23, last_vt_ad=1, theta_hat_inter_p3=0.2,
#' theta_hat_inter_p3_sd=0.4, theta_hat_last_ad_p3=0.32, theta_hat_last_ad_p3_sd=0.4,
#' alpha_2=0.05, dist=1)
#'
#' @export
two_stage_est_back<-function(s,c_bry_two_stg,c_bry_p3,comp,snew,c_bry_new,c_bry_new_p3,theta_hat_inter,theta_hat_inter_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,
                             theta_hat_inter_p3,theta_hat_inter_p3_sd,theta_hat_last_ad_p3,theta_hat_last_ad_p3_sd,alpha_2,dist){
  s_p3=(s[-1]-s[1])/(1-s[1])

  two_stage_out=two_stage_llc_back(s,c_bry_two_stg,comp,snew,c_bry_new,theta_hat_inter,theta_hat_inter_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,alpha_2,dist)
  p_value=two_stage_out[1]
  theta_ll=two_stage_out[2]
  inter_vt_p3=inter_vt-1
  last_vt_ad_p3=last_vt_ad
  #  T_inter=theta_hat_inter_sd^(-2)
  #  T=T_inter/s[inter_vt]
  #  t=T*s
  #  T_last_ad=theta_hat_last_ad_sd^(-2)
  #  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter

  if ( inter_vt_p3==0){
    alpha=alpha_2/2
    p3_out=as.numeric(unlist(gsd_ci_est(snew,c_bry=c_bry_new_p3,theta_hat_last=theta_hat_last_ad_p3,theta_hat_last_sd=theta_hat_last_ad_p3_sd,last_vt=last_vt_ad,alpha)[2]))
  }else{
    p3_out=ph3_ci_est_back(s_p3,c_bry_p3,snew,c_bry_new_p3,theta_hat_inter_p3,theta_hat_inter_p3_sd,inter_vt_p3,theta_hat_last_ad_p3,theta_hat_last_ad_p3_sd,last_vt_ad_p3,alpha_2)
  }
  theta_pt=p3_out[2]
  theta_ul=p3_out[3]
  theta_ci=c(theta_ll,theta_pt,theta_ul)
  out=list(p_value,theta_ci)
  names(out)=c("p-value","confidence interval")
  return(out)
}

######################################################################
#two_stage_est_back(s,c_bry_two_stg,c_bry_p3,comp,snew,c_bry_new,c_bry_new_p3,theta_hat_inter,theta_hat_inter_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,
#                             theta_hat_inter_p3,theta_hat_inter_p3_sd,theta_hat_last_ad_p3,theta_hat_last_ad_p3_sd,last_vt_ad_p3,alpha_2,dist)
#two_stage_est_back(s=c(1/3,2/3,1),c_bry_two_stg=c(3.776606,2.670463,2.180424),c_bry_p3=c(3.954862,1.977431),comp=2,snew=c(0.5,1),c_bry_new=c(3,2.2),c_bry_new_p3=c(3.2,2.4),theta_hat_inter=0.4,theta_hat_inter_sd=0.35,inter_vt=1,theta_hat_last_ad=0.5,theta_hat_last_ad_sd=0.23,last_vt_ad=1,theta_hat_inter_p3=0.2,theta_hat_inter_p3_sd=0.4,theta_hat_last_ad_p3=0.32,theta_hat_last_ad_p3_sd=0.4,alpha_2=0.05,dist=1)
#two_stage_est_back(s=c(1/3,2/3,1),c_bry_two_stg=c(3.776606,2.670463,2.180424),c_bry_p3=c(3.954862,1.977431),comp=2,snew=c(0.5,1),c_bry_new=c(3,2.2),c_bry_new_p3=c(3.2,2.4),theta_hat_inter=0.4,theta_hat_inter_sd=0.35,inter_vt=2,theta_hat_last_ad=0.5,theta_hat_last_ad_sd=0.23,last_vt_ad=1,theta_hat_inter_p3=0.2,theta_hat_inter_p3_sd=0.4,theta_hat_last_ad_p3=0.32,theta_hat_last_ad_p3_sd=0.4,alpha_2=0.05,dist=1)
#two_stage_est_back(s=c(1/3,2/3,1),c_bry_two_stg=c(4.057782,2.869285,2.342762),c_bry_p3=c(3.954862,1.977431),comp=2,snew=c(1),c_bry_new=c(2.404),c_bry_new_p3=c(3.31),theta_hat_inter=0.4,theta_hat_inter_sd=0.17,inter_vt=2,theta_hat_last_ad=0.39,theta_hat_last_ad_sd=0.127,last_vt_ad=1,theta_hat_inter_p3=0.37,theta_hat_inter_p3_sd=0.17,theta_hat_last_ad_p3=0.365,theta_hat_last_ad_p3_sd=0.15,alpha_2=0.05,dist=1)


