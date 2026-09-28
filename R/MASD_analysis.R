#rm(list=ls())
#library(foreach)
#library(parallel)
#library(doParallel)
#library(bigstatsr)

###################################################
##########   Part I                    ############
### Component functions, to be sourced ############
###################################################


#############################################################
#' Multiple comparisons Group sequential design: Final analysis, Without sample size change
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param theta_hat_last_m1 'theta_hat_last_m1' as used by this function; see Examples.
#' @param theta_hat_last_m1_sd 'theta_hat_last_m1_sd' as used by this function; see Examples.
#' @param last_vt Index of the analysis at which the trial stopped
#' @param alpha_2 Type I error for the two-sided test
#' @param dist Indicator of the endpoint / distribution type
#'
#' @return A named list containing the following elements: 'p-value', 'confidence interval'.
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' masd_est_ci(s = c(1/3, 2/3, 1), c_bry = c(3.984298, 3.961146, 2.349682),
#' theta_hat_last_m1 = 0.35, theta_hat_last_m1_sd = 0.6,
#' last_vt = 1, alpha_2 = 0.05, dist = 1)
#'
#' @export
masd_est_ci<-function(s,c_bry,theta_hat_last_m1,theta_hat_last_m1_sd,last_vt,alpha_2,dist){
  #this function provide the point estimate and confidence interval for the dose with the best efficacy#
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #c_bry is the critical bpundary
  #alpha_2  is the significanc level of the (1-alpha_2)*100% CI#
  #typically, alpha_2=0.05
  # theta_hat_last is the vector of the last observed theta_hat###
  #theta_hat_last_sd is  the  standard error of theta_hat_last. length(theta_hat_last_sd)=1, or length(theta_hat_last_sd)=length(theta_hat_last)
  #last_vt numbers the interim analysis at which the last analysis was performed
  #last_vt is the last interim/final analysis where the trial was terminated
  #dist=1 if the distirbution is normal, and the dose arms share a common control arm

  comp=length(theta_hat_last_m1)
  #####use last_obs_wald_max for estimation #########
  wald_last=theta_hat_last_m1/theta_hat_last_m1_sd
  wald_last_max=max(wald_last)

  T_last=(mean(theta_hat_last_m1_sd))^(-2)
  T=T_last/s[last_vt]

  t=T*s
  ############################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)

  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }
  cov_mtx_ws=matrix(nrow=length(s)*comp, ncol=length(s)*comp)
  for (i in 1:length(s)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s)*comp)]=s[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
  #######################
  cov_mtx_ws_mf=cov_mtx_ws[1:(last_vt*comp),1:(last_vt*comp)]
  ################ For a single Brownian  motion####
  cov_mtx_s= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx_s[lower.tri(cov_mtx_s)]=t(cov_mtx_s)[lower.tri(cov_mtx_s)]
  cov_mtx_sf=cov_mtx_s[1:last_vt,1:last_vt]
  ####################################################

  ## the function for the vector of Brownian motions ###
  cal1<-function(theta,sigma){
    ####################Boundary for vector of Brownian  motions#########################
    if (last_vt==1) {
      u_bry_theta=rep(wald_last_max*sqrt(s[last_vt])-theta*sqrt(T)*s[last_vt],comp)
    }else{
      u_c_bry=c(c_bry[1:(last_vt-1)],wald_last_max)
      u_bry_theta_mat=matrix(nrow=last_vt,ncol=comp)
      for (i in (1:last_vt)) {
        u_bry_theta_mat[i,]=rep(u_c_bry[i]*sqrt(s[i])-theta*sqrt(T)*s[i],comp)
      }
      u_bry_theta=c(t( u_bry_theta_mat))
    }
    ###################################################
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,last_vt*comp),sigma=sigma,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }

  ### the function for single Brownian motion###
  cal2<-function(theta,sigma){
    ##################Boundary for single Brownian motion ######################
    if (last_vt==1){
      up_bry_sf_theta=wald_last_max-theta*sqrt(T*s[last_vt])
    }else if (last_vt>1){
      u_c_bry=c(c_bry[1:(last_vt-1)],wald_last_max)
      up_bry_sf_theta=u_c_bry*sqrt(s[1:last_vt])-theta*sqrt(T)*s[1:last_vt]
    }
    if (last_vt==1){
      tmp = 1-stats::pnorm(up_bry_sf_theta)
    }else{
      tmp = 1-pmvnorm(lower=-Inf,upper=up_bry_sf_theta,mean=rep(0,last_vt),sigma=sigma,algorithm=mvtnorm::GenzBretz())[1]
    }
    return(tmp)
  }
  #####################

  u_bry_theta_sf=cal2(1,cov_mtx_sf)
  u_bry_theta_mf=cal1(1,cov_mtx_ws_mf)
  theta_ll=as.numeric(stats::uniroot(function(x) cal1(x,cov_mtx_ws_mf)-alpha_2/2,c(-2^32,2^32),tol=1e-6)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal1(x,cov_mtx_ws_mf)-0.5,c(-2^32,2^32),tol=1e-6)[1])
  theta_ul=as.numeric(stats::uniroot(function(x) cal2(x,cov_mtx_sf)-1+alpha_2/2,c(-2^32,2^32),tol=1e-6)[1])
  p_value=cal1(0,cov_mtx_ws_mf)

  out=list( p_value,c(theta_ll,theta_pt,theta_ul))
  names(out)=c("p-value","confidence interval")
  return(out)
}
########################
#masd_est_ci(s,c_bry,theta_hat_last_m1,theta_hat_last_m1_sd,last_vt,alpha,dist)
#masd_est_ci(s=c(1/3,2/3,1),c_bry=c(4.115415, 2.910038, 2.376036),theta_hat_last_m1=c(2.376,2,1),theta_hat_last_m1_sd=0.3,last_vt=1,alpha_2=0.05,dist=1)
#masd_est_ci(s=c(1/3,2/3,1),c_bry=c(4.115415, 2.910038, 2.376036),theta_hat_last_m1=c(2.376,2,1),theta_hat_last_m1_sd=0.3,last_vt=2,alpha_2=0.05,dist=1)
#masd_est_ci(s=c(1/3,2/3,1),c_bry=c(4.115415, 2.910038, 2.376036),theta_hat_last_m1=c(2.376*0.3,2*0.3,1*0.3),theta_hat_last_m1_sd=0.3,last_vt=3,alpha_2=0.05,dist=1)
#masd_est_ci(s=c(1/3,2/3,1),c_bry=c(4.115415, 2.910038, 2.376036),theta_hat_last_m1=c(1,2.376,2),theta_hat_last_m1_sd=0.3,last_vt=3,alpha_2=0.05,dist=2)
#masd_est_ci(s=c(1/3,2/3,1),c_bry=c(4.115415, 2.910038, 2.376036),theta_hat_last_m1=c(2.376*0.3,2*0.3,1*0.3),theta_hat_last_m1_sd=c(0.3,0.3,0.3),last_vt=3,alpha_2=0.05,dist=1)
#masd_est_ci(s=c(1/3,2/3,1),c_bry=c(4.115415, 2.910038, 2.376036),theta_hat_last_m1=c(2.376*0.3,2*0.3,1*0.3),theta_hat_last_m1_sd=c(0.3,0.3,0.3),last_vt=3,alpha_2=0.05,dist=1)
#masd_est_ci(s=c(0.5,1),c_bry=c(3.179276 ,2.248088 ),theta_hat_last_m1=c(0.26,0.31),theta_hat_last_m1_sd=0.13,last_vt=2,alpha_2=0.05,dist=1)
#masd_est_ci(s=c(0.5,1),c_bry=c(3.179276 ,2.248088 ),theta_hat_last_m1=c(0.31,0.31),theta_hat_last_m1_sd=0.13,last_vt=2,alpha_2=0.05,dist=1)
#masd_est_ci(s=c(0.5,1),c_bry=c(3.179276 ,2.248088 ),theta_hat_last_m1=c(0.31,0.31),theta_hat_last_m1_sd=0.13,last_vt=2,alpha_2=0.05,dist=2)
#############################################################

###################################################################
#### part II: component function for conditional probabilities ####
###         and adaptive procedures                            ####
###################################################################

masd_c_power_theta<-function(theta,s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist){
  #this function calculates conditional power, given that theta_hat_inter,theta_hat_inter_sd have been observed at inter_vt#
  # theta is the effect size (it is a vector. length(theta)=length(theta_hat_inter)
  #s is the vector of information fractions#
  #c_bry is the critical boundary
  # theta_hat_inter is the vector of theta_hat at the inter_vt interim analysis###
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the  analysis was conducted at the inter_vt-thinterim analysis
  #dist=1 if the the test arms share one common contorl and the endpoint is normally distributed
  #### the next section constructs the covariance matrix ###

  s1=(s[(inter_vt+1):length(s)]-s[inter_vt])/(1-s[inter_vt])

  comp=length(theta_hat_inter_m1)
  theta_hat=theta_hat_inter_m1

  T_inter=1/(mean(theta_hat_inter_m1_sd))^2
  score_inter=theta_hat_inter_m1*T_inter
  inter_wald=theta_hat_inter_m1*sqrt(T_inter)
  T=T_inter/s[inter_vt]
  T_inter_dif=T-T_inter
  t=T*s
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)

  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  }else {
    cov_mtx=diag(comp)
  }

  cov_mtx_ws=matrix(nrow=length(s1)*comp, ncol=length(s1)*comp)
  for (i in 1:length(s1)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s1)*comp)]=s1[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
  #####################################################
  ### boundary for the  vector of Brownian  motions####

  u_bry_theta_mat=matrix(ncol=comp,nrow=length(s1))

  theta1 = rep(theta,length.out=comp)
  priv_bry = c_bry*sqrt(t)
  priv_sqrt = sqrt(T-T_inter)
  for (i in (1:length(s1))) {
     for (j in 1:comp) {
       u_bry_theta_mat[i,j]=(priv_bry[i+inter_vt]-score_inter[j])/priv_sqrt - theta1[j]*priv_sqrt*s1[i]
     }
  }
  u_bry_theta=c(t(u_bry_theta_mat))

  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e0#

  c_power= 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,(length(s1)*comp)),sigma=cov_mtx_ws)[1]
  return(c_power)
}
########################################
#masd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)
#masd_c_power_theta(theta=c(0,0.2,0.1),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=1,theta_hat_inter_m1=c(0.3,0.4,1),theta_hat_inter_m1_sd=0.5,dist=1)
#masd_c_power_theta(c(0,0.2,0.1),c(1/3,2/3,1),c(3.984298, 3.961146, 2.349682),2,c(0.3,0.4,1),0.5,1)
#masd_c_power_theta(c(0,0,0),c(1/3,2/3,1),c(3.984298, 3.961146, 2.349682),1,c(0.3,0.4,1),0.5,1)
###########################################################

###########################################################
masd_c_power<-function(s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist){
  #this function calculates conditional power, given that theta_hat_inter,theta_hat_inter_sd have been observed at inter_vt#
  # the effect size theta is assumed to be equal to theta_hat_inter
  #s is the vector of information fractions#
  #c_bry is the critical boundary
  # theta_hat_inter is the vector of theta_hat at the inter_vt interim analysis###
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the  analysis was conducted at the inter_vt-thinterim analysis
  #dist=1 if the the test arms share one common contorl and the endpoint is normally distributed
  s1=(s[(inter_vt+1):length(s)]-s[inter_vt])/(1-s[inter_vt])

  comp=length(theta_hat_inter_m1)
  theta=theta_hat_inter_m1
  theta_hat=theta_hat_inter_m1

  T_inter=1/(mean(theta_hat_inter_m1_sd))^2
  score_inter=theta_hat_inter_m1*T_inter
  inter_wald=theta_hat_inter_m1*sqrt(T_inter)
  T=T_inter/s[inter_vt]
  T_inter_dif=T-T_inter
  t=T*s
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  }else {
    cov_mtx=diag(comp)
  }

  cov_mtx_ws = matrix(nrow=length(s1)*comp, ncol=length(s1)*comp)
  for (i in 1:length(s1)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s1)*comp)]=s1[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)] = t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
  #####################################################
  ### boundary for the  vector of Brownian  motions####

  u_bry_theta_mat=matrix(ncol=comp,nrow=length(s1))
  priv_bry = c_bry*sqrt(t)
  priv_sqrt = sqrt(T-T_inter)
  for (i in (1:length(s1))) {
    for (j in 1:comp){
      u_bry_theta_mat[i,j]=(priv_bry[i+inter_vt]-score_inter[j])/priv_sqrt - theta[j]*priv_sqrt*s1[i]
    }
  }
  u_bry_theta=c(t(u_bry_theta_mat))

  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e0#
  ## the function for the vector of Brownian motions ###
  cal1<-function(theta,sigma){
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,(length(s1)*comp)),sigma=sigma,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }
  #################################################
  c_power= cal1(theta,cov_mtx_ws)
  return(c_power)
}
########################################
#masd_c_power(s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)
#masd_c_power(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=1,theta_hat_inter_m1=c(0.3,0.4,1),theta_hat_inter_m1_sd=0.5,dist=1)
#masd_c_power(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=2,theta_hat_inter_m1=c(0.3,0.4,1),theta_hat_inter_m1_sd=0.5,dist=1)
###############################################################

###########################################################
masd_c_tyIer<-function(s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist){
  #this function calculates the conditional type I error, given that theta_hat_inter,theta_hat_inter_sd have been observed at inter_vt#
  # under the null hypothesis, the effect size theta is assumed to be zero
  #s is the vector of information fractions#
  #c_bry is the critical boundary
  # theta_hat_inter is the vector of thetha_hat at the inter_vt interim analysis###
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the  analysis was conducted at the inter_vt-thinterim analysis
  #dist=1 if the the test arms share one common contorl and the endpoint is normally distributed

  s1=(s[(inter_vt+1):length(s)]-s[inter_vt])/(1-s[inter_vt])

  comp=length(theta_hat_inter_m1)

  theta=theta_hat_inter_m1
  theta_hat=theta_hat_inter_m1

  T_inter=1/(mean(theta_hat_inter_m1_sd))^2
  score_inter=theta_hat_inter_m1*T_inter
  inter_wald=theta_hat_inter_m1*sqrt(T_inter)
  T=T_inter/s[inter_vt]
  T_inter_dif=T-T_inter
  t=T*s
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  }else {
    cov_mtx=diag(comp)
  }

  cov_mtx_ws = matrix(nrow=length(s1)*comp, ncol=length(s1)*comp)
  for (i in 1:length(s1)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s1)*comp)]=s1[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]

  #####################################################
  ### boundary for the  vector of Brownian  motions####

  u_bry_mat=matrix(ncol=comp,nrow=length(s1))
  priv_bry = c_bry*sqrt(t)
  priv_sqrt = sqrt(T-T_inter)
  for (i in (1:length(s1))) {
    for (j in 1:comp){
      u_bry_mat[i,j]=(priv_bry[i+inter_vt]-score_inter[j])/priv_sqrt
    }
  }

  u_bry=c(t(u_bry_mat))

  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e0#
  cal1<-function(theta,sigma){
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry,mean=rep(0,(length(s1)*comp)),sigma=sigma,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }
  ########################################################
  c_alpha= cal1(theta,cov_mtx_ws)
  return(c_alpha)
}
########################################
#masd_c_tyIer(s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)
#masd_c_tyIer(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=1,theta_hat_inter_m1=c(0.3,0.4,1),theta_hat_inter_m1_sd=0.5,dist=1)
#masd_c_tyIer(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=2,theta_hat_inter_m1=c(0.3,0.4,1),theta_hat_inter_m1_sd=0.5,dist=1)
##########################################################

new_OF_bundry_s<-function(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist){
  #this function calculates the new critical bundaries, given the conditional type I error of c_alpha#
  #the boundary is O'Brien-Fleming type#
  #snew is the vector of new information fractions after the sample size change.
  #c_alpha is the vector (with the same length as snew) of alpha spending##
  #s is the vector of information fractions#
  # theta_hat_inter is the vector of the theta_hat at the inter_vt interim analysis###
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the  analysis was conducted at the inter_vt-thinterim analysis
  ## sample_r is the ratio of the new sample size over the original sample size##
  #dist is the distribution of the variable. dist=1 means that the variable is normally distributed and the comparison is for multiple dose comparison#
  s_inter=s[inter_vt]

  T_inter=1/(max(theta_hat_inter_m2_sd))^2
  T_old=T_inter/s_inter
  if (sample_r<=s_inter){
    T_new=T_old*(s_inter*1.001)
  }else{
    T_new=T_old*sample_r
  }
  t_new=(T_new-T_inter)*snew+T_inter
  comp=length(theta_hat_inter_m2)
  score_inter=theta_hat_inter_m2*theta_hat_inter_m2_sd^(-2)
  ###################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)

  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }
  cov_mtx_ws = matrix(nrow=length(snew)*comp, ncol=length(snew)*comp)
  for (i in 1:length(snew)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(snew)*comp)]=snew[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
  #############################

  #### The functions c_cal calcuates the  exit probablities for the given exit bounadries e0_1
  c_cal<-function(e_1){
    u_bry_mat=matrix(ncol=comp,nrow=length(snew))
    for (i in (1:length(snew))) {
      for (j in 1:comp){
        u_bry_mat[i,j]=(e_1-score_inter[j])/sqrt(T_new-T_inter)
      }
    }
    u_bry=c(t(u_bry_mat))
    ########################################
    sigma2=cov_mtx_ws
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry,mean=rep(0,(length(snew)*comp)), sigma=sigma2,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }

  if (T_inter>=T_new){return(NA)}

  bry=as.numeric(stats::uniroot(function(x) c_cal(x)-c_alpha,c(-2^32,2^32),tol=1e-7)[1])
  #bry are exit boundaries. c_bry are the ctritical bundaries#
  c_bry=bry/sqrt(t_new)
  return(c_bry)
}
##############################
#new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
#new_OF_bundry_s(s=c(0.3,0.7,1),inter_vt=2,snew=1,c_alpha=0.025,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1.5,sample_r=1.5,dist=1)
#new_OF_bundry_s(s=c(0.3,0.7,1),inter_vt=2,snew=c(0.5,1),c_alpha=0.025,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1.5,sample_r=1.5,dist=1)
#new_OF_bundry_s(s=c(0.5,1),inter_vt=1,snew=1,c_alpha=0.1210272,theta_hat_inter_m2=0.25,theta_hat_inter_m2_sd=0.15,sample_r=2.25,dist=1)
##########################################################


######################################
masd_c_power_new<-function(s,inter_vt,snew,new_c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist){
  #this function calculates the conditional power with the new critical boundary#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #the effect size is taken as theta, which is a vector, length(theta)=length(theta_hat_inter)
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the analysis is conducted at the inter_vt-th  interim analysis
  #snew is the relative information fraction for the new schedule of interim analysis
  #new_c_bry is the new critical boundary after the sample size change
  # theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #sample_r=T_new/T=N_new/N, where T_new is the new information time, and T is the planned information time
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm
  s_inter=s[inter_vt]
  theta=theta_hat_inter_m2
  T_inter=1/(mean(theta_hat_inter_m2_sd))^2
  T_old=T_inter/s_inter
  T=T_old
  comp=length(theta_hat_inter_m2)
  T_new=T_old*sample_r

  t_new=(T_new-T_inter)*snew+T_inter
  score_inter=theta_hat_inter_m2*theta_hat_inter_m2_sd^(-2)

  new_e_bry=new_c_bry*sqrt(t_new)
  #################################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }

  cov_mtx_ws = matrix(nrow=length(snew)*comp, ncol=length(snew)*comp)
  for (i in 1:length(snew)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(snew)*comp)]=snew[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]

  #############################
  u_bry_mat=matrix(ncol=comp,nrow=length(snew))
  priv_sqrt = sqrt(T_new-T_inter)
  for (i in (1:length(snew))) {
    for (j in 1:comp){
      u_bry_mat[i,j]=(new_e_bry[i]-score_inter[j])/priv_sqrt-theta[j]*priv_sqrt*snew[i]
    }
  }

  u_bry=c(t(u_bry_mat))

  ########################################
  ## the function for the vector of Brownian motions ###
  cal1<-function(theta,sigma){
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry,mean=rep(0,(length(snew)*comp)),sigma=sigma,abseps=1e-8,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }

  ### the function for single Brownian motion###
  if (T_inter>=T_new){
    c_power_new=0
  } else {
    c_power_new= cal1(theta,cov_mtx_ws)
  }
  return(c_power_new)
}
##################################
#masd_c_power_new(s,inter_vt,snew,new_c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
#masd_c_power_new(s=c(0.3,0.7,1),inter_vt=2,snew=1, new_c_bry=2.781890,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.4,sample_r=2,dist=1)
#masd_c_power_new(s=c(0.3,0.7,1),inter_vt=2,snew=1,new_c_bry=c(3.310743, 2.781890),theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1,sample_r=2,dist=1)
#masd_c_power_new(s=c(0.3,0.7,1),inter_vt=2,snew=c(0.5,1),new_c_bry=c(3.310743, 2.781890),theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.4,sample_r=2,dist=1)
##################################

######################################
masd_c_power_new_theta<-function(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist){
  #this function calculates the conditional power with the new critical boundary#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #the effect size is taken as theta, which is a vector, length(theta)=length(theta_hat_inter)
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the analysis is conducted at the inter_vt-th  interim analysis
  #snew is the relative information fraction for the new schedule of interim analysis
  #new_c_bry is the new critical boundary after the sample size change
  # theta_hat_inter_m2 is the vector ob obsered theta_hat at the interim analysis ##
  #length(theta)=length(theta_hat_inter_m2)
  #theta_hat_inter_m2_sd is the standard error of theta_hat_inter
  #sample_r=T_new/T=N_new/N, where T_new is the new information time, and T is the planned information time
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm
  s_inter=s[inter_vt]

  T_inter=1/(mean(theta_hat_inter_m2_sd))^2
  T_old=T_inter/s_inter
  T=T_old
  comp=length(theta_hat_inter_m2)
  T_new=T_old*sample_r

  t_new=(T_new-T_inter)*snew+T_inter
  score_inter=theta_hat_inter_m2*theta_hat_inter_m2_sd^(-2)

  new_e_bry=new_c_bry*sqrt(t_new)
  #################################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }

  cov_mtx_ws <- matrix(nrow=length(snew)*comp, ncol=length(snew)*comp)
  for (i in 1:length(snew)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(snew)*comp)]=snew[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]

  #############################
  u_bry_mat=matrix(ncol=comp,nrow=length(snew))
  priv_sqrt = sqrt(T_new-T_inter)
  for (i in (1:length(snew))) {
    for (j in 1:comp){
      u_bry_mat[i,j]=(new_e_bry[i]-score_inter[j])/priv_sqrt - theta[j]*priv_sqrt*snew[i]
    }
  }
  u_bry=c(t(u_bry_mat))
  ########################################

  ## the function for the vector of Brownian motions ###
  cal1<-function(theta,sigma){
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry,mean=rep(0,(length(snew)*comp)),sigma=sigma,abseps=1e-8,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }
  ### the function for single Brownian motion###
  if (T_inter>=T_new){
    c_power_new=0
  } else {
    c_power_new= cal1(theta,cov_mtx_ws)
  }
  return(c_power_new)
}
##################################
#masd_c_power_new_theta(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist){
#masd_c_power_new_theta(theta=c(0,0.2),s=c(0.3,0.7,1),inter_vt=2,snew=1, new_c_bry=2.781890,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.4,sample_r=2,dist=1)
#masd_c_power_new_theta(theta=c(0.3,0.2),s=c(0.3,0.7,1),inter_vt=2,snew=1,new_c_bry=c(3.310743, 2.781890),theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1,sample_r=2,dist=1)
#masd_c_power_new_theta(theta=c(0.3,0.2),s=c(0.3,0.7,1),inter_vt=2,snew=c(0.5,1),new_c_bry=c(3.310743, 2.781890),theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.4,sample_r=2,dist=1)
##################################

masd_new_info_time<-function(theta,s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,dist){
  #this function calculates the new information time to achieve desired conditional power of 1-beta#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #the effect size is taken as theta, which is a vector, length(theta)=length(theta_hat_inter)
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the analysis is conducted at the inter_vt'th  interim analysis
  #c_alpha  is the target one-sided type I error#
  # theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  # for normal distribution, dist=1###
  # 1-beta is the desired power
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm
  s_inter=s[inter_vt]

  comp=length(theta)
  #sample_r=s_old_lst
  T_inter=(mean(theta_hat_inter_m2_sd))^(-2)
  T=T_inter/s[inter_vt]

  power<-function(sample_r){
    c_bry=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    tmp=masd_c_power_new_theta(theta,s,inter_vt,snew,c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    return(tmp)
  }

  new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.0001,2^32),tol=1e-5)[1])
  c_bry_f=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,new_ratio,dist)
  T_new=T*new_ratio

  #  out=list(T_new,c_bry_f)
  #  names(out)=c("New information fraction","new critical boundaries")
  # return(out)
  return(T_new)
}
#################################
#masd_new_info_time(theta,s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,dist){
#masd_new_info_time(theta=c(0,0.2),s=c(0.3,0.7,1),inter_vt=2,snew=1,c_alpha=0.4,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.2,beta=0.1,dist=1)
####################################
##################################
#masd_c_power_new(s,inter_vt,snew,new_c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)

#####################################
masd_new_samsz_theta_1<-function(theta,s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,n_control_inter,dist){
  #this function calculates the new sample size and new critical boundary to achieve desired conditional power of 1-beta#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #the effect size is taken as theta, which is a vector, length(theta)=length(theta_hat_inter_m1)
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the analysis is conducted at the inter_vt'th  interim analysis
  #c_alpha  is the target one-sided type I error#
  # theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  # for normal distribution, dist=1###
  # 1-beta is the desired power
  #n_control_inter is the sample size for the comtrol arm at the interim analysis
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm
  s_inter=s[inter_vt]

  comp=length(theta_hat_inter_m1)
  mu_keep=match(theta_hat_inter_m2,theta_hat_inter_m1)

  theta1= rep(theta,length.out=comp)
  theta_keep=theta1[mu_keep]
  T_inter=(mean(theta_hat_inter_m1_sd))^(-2)
  T=T_inter/s_inter
  wald_inter=theta_hat_inter_m1/theta_hat_inter_m1_sd
  wald_inter_max=max(wald_inter)
  test= wald_inter_max-c_bry[inter_vt]

  c_power_m1=masd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  c_alpha=masd_c_power_theta(rep(0,comp),s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  power<-function(sample_r){
    c_bry=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    tmp=masd_c_power_new_theta(theta_keep,s,inter_vt,snew,c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    return(tmp)
  }

  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    N_control_new="NA"
    c_bry_f="NA"
    c_power_m2="NA"
    new_ratio="NA"
    c_power_m2_beta="NA"
  }else {
    outcome="Null hypothesis not rejected at interim analysis"

    if (max(theta)<10e-4){
      N_control_new="NA"
      c_bry_f="NA"
      c_power_m2="NA"
    }else{
      c_power_m2=as.numeric(power(1))
      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.001,2^32),tol=1e-5)[1])
      c_bry_f=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,new_ratio,dist)
      T_new=T*new_ratio
      c_power_m2_beta=as.numeric(power(new_ratio))
      N_control_new=ceiling(n_control_inter*(T_new/T_inter))
    }
  }
  # out=list(outcome,wald_inter,c_bry[inter_vt],theta,theta_keep,c_alpha, c_power_m1,c_power_m2,N_control_new,snew,c_bry_f)
  # names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Assumed effect size: theta","remaining effect size","c_alpha","c_power_m1","c_power_m2","N_control_new","new information fraction","New critical boundary")
  # out=list(test,c_power_m2,c_power_m2_beta,new_ratio)
  # names(out)=c("test","c_power_m2","c_power_m2_beta","new_ratio")
  out=list(test,c_power_m2)
  names(out)=c("test","c_power_m2")
  return(out)
}
#################################
#masd_new_samsz_theta<-function(theta,s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,n_control_inter,dist)
#masd_new_samsz_theta_1(theta=c(0,0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(1.2,1.5,1.1),theta_hat_inter_m1_sd=0.2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz_theta_1(theta=c(0,0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(0.2,0.5,0.8),theta_hat_inter_m1_sd=0.5,theta_hat_inter_m2=c(0.5,0.8),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz_theta_1(theta=c(0,-0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(0.2,0.5,0.8),theta_hat_inter_m1_sd=0.5,theta_hat_inter_m2=c(0.5,0.8),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz_theta_1(theta=c(0.2,0.1,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(0.2,0.5,0.8),theta_hat_inter_m1_sd=0.5,theta_hat_inter_m2=c(0.5,0.8),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)

###################################
#' Multiple comparisons Group sequential design: Interim analysis, Non-survival analysis, Input theta
#'
#' @param theta Assumed effect size
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param snew Vector of new information fractions after the sample size is changed
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based)
#' @param theta_hat_inter_m1 'theta_hat_inter_m1' as used by this function; see Examples.
#' @param theta_hat_inter_m1_sd 'theta_hat_inter_m1_sd' as used by this function; see Examples.
#' @param theta_hat_inter_m2 'theta_hat_inter_m2' as used by this function; see Examples.
#' @param theta_hat_inter_m2_sd 'theta_hat_inter_m2_sd' as used by this function; see Examples.
#' @param beta Type II error rate, i.e. 1 - power
#' @param n_control_inter Sample size of the control arm at the interim analysis
#' @param dist Indicator of the endpoint / distribution type
#'
#' @return A named list containing the following elements: 'Interim analysis hypothesis test', 'Wald statistics', 'Critical boundary',
#' 'Assumed effect size: theta', 'remaining effect size', 'c_alpha', 'c_power_m1', 'c_power_m2', 'N_control_new', 'new information fraction',
#' 'New critical boundary'.
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' masd_new_samsz_theta(theta = c(0,0.2,0), s = c(1/3,2/3,1),
#' c_bry = c(3.984298, 3.961146, 2.349682), snew = 1, inter_vt = 2,
#' theta_hat_inter_m1 = c(1.2,1.5,1.1), theta_hat_inter_m1_sd = 0.2,
#' theta_hat_inter_m2 = c(1.2,1.5), theta_hat_inter_m2_sd = 0.2,
#' beta = 0.1, n_control_inter = 65, dist = 1)
#'
#' @export
masd_new_samsz_theta<-function(theta,s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,n_control_inter,dist){
  #this function calculates the new sample size and new critical boundary to achieve desired conditional power of 1-beta#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #the effect size is taken as theta, which is a vector, length(theta)=length(theta_hat_inter_m1)
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the analysis is conducted at the inter_vt'th  interim analysis
  #c_alpha  is the target one-sided type I error#
  # theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  # for normal distribution, dist=1###
  # 1-beta is the desired power
  #n_control_inter is the sample size for the comtrol arm at the interim analysis
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm
  s_inter=s[inter_vt]

  comp=length(theta_hat_inter_m1)
  mu_keep=match(theta_hat_inter_m2,theta_hat_inter_m1)

  theta1= rep(theta,length.out=comp)
  theta_keep=theta1[mu_keep]
  T_inter=(mean(theta_hat_inter_m1_sd))^(-2)
  T=T_inter/s_inter
  wald_inter=theta_hat_inter_m1/theta_hat_inter_m1_sd
  wald_inter_max=max(wald_inter)
  test= wald_inter_max-c_bry[inter_vt]

  c_power_m1=masd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  c_alpha=masd_c_power_theta(rep(0,comp),s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  power<-function(sample_r){
    c_bry=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    tmp=masd_c_power_new_theta(theta_keep,s,inter_vt,snew,c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    return(tmp)
  }

  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    N_control_new="NA"
    c_bry_f="NA"
    c_power_m2="NA"
  }else {
    outcome="Null hypothesis not rejected at interim analysis"
    if (max(theta)<10e-4){
      N_control_new="NA"
      c_bry_f="NA"
      c_power_m2="NA"
    }else{
      c_power_m2=as.numeric(power(1))
      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.001,2^32),tol=1e-5)[1])
      c_bry_f=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,new_ratio,dist)
      T_new=T*new_ratio
      N_control_new=ceiling(n_control_inter*(T_new/T_inter))
    }
  }
  out=list(outcome,wald_inter,c_bry[inter_vt],theta,theta_keep,c_alpha, c_power_m1,c_power_m2,N_control_new,snew,c_bry_f)
  names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Assumed effect size: theta","remaining effect size","c_alpha","c_power_m1","c_power_m2","N_control_new","new information fraction","New critical boundary")
  return(out)
}
#################################
#masd_new_samsz_theta<-function(theta,s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,n_control_inter,dist)
#masd_new_samsz_theta(theta=c(0,0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(1.2,1.5,1.1),theta_hat_inter_m1_sd=0.2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz_theta(theta=c(0,0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(0.2,0.5,0.8),theta_hat_inter_m1_sd=0.5,theta_hat_inter_m2=c(0.5,0.8),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz_theta(theta=c(0,-0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(0.2,0.5,0.8),theta_hat_inter_m1_sd=0.5,theta_hat_inter_m2=c(0.5,0.8),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz_theta(theta=c(0.2,0.1,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(0.2,0.5,0.8),theta_hat_inter_m1_sd=0.5,theta_hat_inter_m2=c(0.5,0.8),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)

###################################
##################################

#' Multiple comparisons Group sequential design: Interim analysis, Survival analysis, Input theta
#'
#' @param theta Assumed effect size
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param snew Vector of new information fractions after the sample size is changed
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based)
#' @param theta_hat_inter_m1 'theta_hat_inter_m1' as used by this function; see Examples.
#' @param theta_hat_inter_m1_sd 'theta_hat_inter_m1_sd' as used by this function; see Examples.
#' @param theta_hat_inter_m2 'theta_hat_inter_m2' as used by this function; see Examples.
#' @param theta_hat_inter_m2_sd 'theta_hat_inter_m2_sd' as used by this function; see Examples.
#' @param beta Type II error rate, i.e. 1 - power
#' @param tot_events_inter_m1 'tot_events_inter_m1' as used by this function; see Examples.
#' @param dist Indicator of the endpoint / distribution type
#'
#' @return A named list containing the following elements: 'Interim analysis hypothesis test', 'Wald statistics',
#' 'Critical boundary', 'Assumed effect size: theta', 'remaining effect size', 'c_alpha', 'c_power_m1', 'c_power_m2', 'new events between a test arm and control', 'new total events', 'new information fraction', 'New critical boundary'.
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' masd_new_samsz_theta_survival(theta=c(0,0.2,0), s=c(1/3,2/3,1),
#' c_bry=c(3.984298, 3.961146, 2.349682), snew=1, inter_vt=2,
#' theta_hat_inter_m1=c(1.2,1.5,1.1), theta_hat_inter_m1_sd=0.2, theta_hat_inter_m2=c(1.2,1.5),
#' theta_hat_inter_m2_sd=0.2, beta=0.1, tot_events_inter_m1=65, dist=1)
#'
#' @export
masd_new_samsz_theta_survival<-function(theta,s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,tot_events_inter_m1,dist){
  #this function calculates the new sample size and new critical boundary to achieve desired conditional power of 1-beta#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #the effect size is taken as theta, which is a vector, length(theta)=length(theta_hat_inter)
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the analysis is conducted at the inter_vt'th  interim analysis
  #c_alpha  is the target one-sided type I error#
  # theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  # for normal distribution, dist=1###
  # 1-beta is the desired power
  #n_control_inter is the sample size for the comtrol arm at the interim analysis
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm
  s_inter=s[inter_vt]

  comp=length(theta_hat_inter_m1)
  new_comp=length(theta_hat_inter_m2)
  mu_keep=match(theta_hat_inter_m2,theta_hat_inter_m1)
  theta1=rep(theta,length.out=comp)
  theta_keep=theta1[mu_keep]
  T_inter=(mean(theta_hat_inter_m2_sd))^(-2)
  T=T_inter/s_inter
  wald_inter=theta_hat_inter_m1/theta_hat_inter_m1_sd
  wald_inter_max=max(wald_inter)
  test= wald_inter_max-c_bry[inter_vt]

  c_power_m1=masd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  c_alpha=masd_c_power_theta(rep(0,comp),s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)


  power<-function(sample_r){
    c_bry=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    tmp=masd_c_power_new_theta(theta_keep,s,inter_vt,snew,c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    return(tmp)
  }

  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    two_arm_events_inter_m2="NA"
    tot_events_new_m2="NA"
    c_bry_f="NA"
    c_power_m2="NA"
    events_two_arm_new_m2="NA"
  }else {
    outcome="Null hypothesis not rejected at interim analysis"
    if (max(theta)<10e-4){
      two_arm_events_inter_m2="NA"
      tot_events_new_m2="NA"
      c_bry_f="NA"
      c_power_m2="NA"
      events_two_arm_new_m2="NA"
    }else{
      c_power_m2=as.numeric(power(1))
      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.001,2^32),tol=1e-5)[1])
      c_bry_f=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,new_ratio,dist)
      T_new=T*new_ratio
      two_arm_events_inter_m2=2*tot_events_inter_m1/(comp+1)
      two_arm_events_m2_planned=two_arm_events_inter_m2/s_inter

      events_two_arm_new_m2=ceiling(two_arm_events_m2_planned*(T_new/T))
      tot_events_new_m2=ceiling((new_comp+1)*events_two_arm_new_m2/2)
    }
  }

  out=list(outcome,wald_inter,c_bry[inter_vt],theta,theta_keep,c_alpha, c_power_m1,c_power_m2,two_arm_events_inter_m2,tot_events_new_m2,snew,c_bry_f)
  names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Assumed effect size: theta","remaining effect size","c_alpha","c_power_m1","c_power_m2","new events between a test arm and control","new total events","new information fraction","New critical boundary")
  # out=list(new_comp,new_ratio,outcome,wald_inter,c_bry[inter_vt],theta,theta_keep,c_alpha, c_power_m1,  c_power_m2,events_two_arm_new_m2,tot_events_new_m2,snew,c_bry_f)
  # names(out)=c("new_comp","new_ratio","Interim analysis hypothesis test","Wald statistics","Critical boundary","Assumed effect size: theta","remaining effect size","c_alpha","c_power_m1","c_power_m2","new events between a test arm and control","new total events","new information fraction","New critical boundary")
  # out=list(test,new_ratio,c_power_m2,c_power_m2_new_ratio)
  return(out)
}
#################################
#masd_new_samsz_theta_survival(theta,s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,tot_events_inter_m1,dist){
#masd_new_samsz_theta_survival(theta=c(0,0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(1.2,1.5,1.1),theta_hat_inter_m1_sd=0.2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.2,beta=0.1,tot_events_inter_m1=65,dist=1)
#masd_new_samsz_theta_survival(theta=c(0,0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(1.2,1.5,1.1),theta_hat_inter_m1_sd=1.2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1.2,beta=0.1,tot_events_inter_m1=65,dist=1)
#masd_new_samsz_theta_survival(theta=c(0,-0.2,0),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(1.2,1.5,1.1),theta_hat_inter_m1_sd=1.2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1.2,beta=0.1,tot_events_inter_m1=65,dist=1)

###################################
##################################
#####################################
#' Multiple comparisons Group sequential design: Interim analysis, Non-survival analysis, Use theta_interim
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param snew Vector of new information fractions after the sample size is changed
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based)
#' @param theta_hat_inter_m1 'theta_hat_inter_m1' as used by this function; see Examples.
#' @param theta_hat_inter_m1_sd 'theta_hat_inter_m1_sd' as used by this function; see Examples.
#' @param theta_hat_inter_m2 'theta_hat_inter_m2' as used by this function; see Examples.
#' @param theta_hat_inter_m2_sd 'theta_hat_inter_m2_sd' as used by this function; see Examples.
#' @param beta Type II error rate, i.e. 1 - power
#' @param n_control_inter Sample size of the control arm at the interim analysis
#' @param dist Indicator of the endpoint / distribution type
#'
#' @return A named list containing the following elements: 'Interim analysis hypothesis test', 'Wald statistics',
#' 'Critical boundary', 'Estimated effect size: theta', 'remaining effect size', 'c_alpha', 'c_power_m1', 'c_power_m2',
#' 'N_control_new', 'new information fraction', 'New critical boundary'.
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' masd_new_samsz(s=c(1/3,2/3,1), c_bry=c(3.984298, 3.961146, 2.349682), snew=1, inter_vt=2,
#' theta_hat_inter_m1=c(1.2,1.5,1.1), theta_hat_inter_m1_sd=0.2, theta_hat_inter_m2=c(1.2,1.5),
#' theta_hat_inter_m2_sd=0.2, beta=0.1, n_control_inter=65, dist=1)
#'
#' @export
masd_new_samsz<-function(s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,n_control_inter,dist){
  #this function calculates the new sample size and new critical boundary to achieve desired conditional power of 1-beta#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #the effect size is taken as theta, which is a vector, length(theta)=length(theta_hat_inter)
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the analysis is conducted at the inter_vt'th  interim analysis
  #c_alpha  is the target one-sided type I error#
  # theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  # for normal distribution, dist=1###
  # 1-beta is the desired power
  #n_control_inter is the sample size for the comtrol arm at the interim analysis
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm
  s_inter=s[inter_vt]
  theta=theta_hat_inter_m1
  comp=length(theta_hat_inter_m1)
  mu_keep=match(theta_hat_inter_m2,theta_hat_inter_m1)
  theta_keep=theta[mu_keep]
  T_inter=(mean(theta_hat_inter_m1_sd))^(-2)
  T=T_inter/s_inter
  wald_inter=theta_hat_inter_m1/theta_hat_inter_m1_sd
  wald_inter_max=max(wald_inter)
  test= wald_inter_max-c_bry[inter_vt]

  c_power_m1=masd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  c_alpha=masd_c_power_theta(rep(0,comp),s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  power<-function(sample_r){
    c_bry=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    tmp=masd_c_power_new_theta(theta_keep,s,inter_vt,snew,c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    return(tmp)
  }

  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    N_control_new="NA"
    c_bry_f="NA"
    c_power_m2="NA"
    new_ratio="NA"
    c_power_m2_new_ratio="NA"
  }else {
    outcome="Null hypothesis not rejected at interim analysis"
    if (max(theta)<10e-4){
      N_control_new="NA"
      c_bry_f="NA"
      c_power_m2="NA"
      new_ratio="NA"
      c_power_m2_new_ratio="NA"
    }else{
      c_power_m2=as.numeric(power(1))
      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.001,2^32),tol=1e-5)[1])
      c_bry_f=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,new_ratio,dist)
      T_new=T*new_ratio
      N_control_new=ceiling(n_control_inter*(T_new/T_inter))
      c_power_m2_new_ratio=as.numeric(power(new_ratio))
    }
  }
  out=list(outcome,wald_inter,c_bry[inter_vt],theta,theta_keep,c_alpha, c_power_m1,c_power_m2,N_control_new,snew,c_bry_f)
  names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Estimated effect size: theta","remaining effect size","c_alpha","c_power_m1","c_power_m2","N_control_new","new information fraction","New critical boundary")
  return(out)
}
#################################
#masd_new_samsz(theta,s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,n_control_inter,dist)
#masd_new_samsz(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(1.2,1.5,1.1),theta_hat_inter_m1_sd=0.2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(0.2,0.5,0.8),theta_hat_inter_m1_sd=0.5,theta_hat_inter_m2=c(0.5,0.8),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(-0.2,-0.5,-0.8),theta_hat_inter_m1_sd=0.5,theta_hat_inter_m2=c(0.5,0.8),theta_hat_inter_m2_sd=0.2,beta=0.1,n_control_inter=65,dist=1)
#masd_new_samsz(s=c(0.5,1),c_bry=c(3.179276, 2.248088 ),snew=1,inter_vt=1,theta_hat_inter_m1=c(0.09,0.25),theta_hat_inter_m1_sd=0.15,theta_hat_inter_m2=0.25,theta_hat_inter_m2_sd=0.15,beta=0.1,n_control_inter=100,dist=1)
#masd_new_samsz(s=c(0.5,1),c_bry=c(3.179276, 2.248088 ),snew=1,inter_vt=1,theta_hat_inter_m1=c(0.09,0.3),theta_hat_inter_m1_sd=0.15,theta_hat_inter_m2=0.3,theta_hat_inter_m2_sd=0.15,beta=0.1,n_control_inter=100,dist=1)

##################################
##################################

#' Multiple comparisons Group sequential design: Interim analysis, Survival analysis, Use theta_interim
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param snew Vector of new information fractions after the sample size is changed
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based)
#' @param theta_hat_inter_m1 'theta_hat_inter_m1' as used by this function; see Examples.
#' @param theta_hat_inter_m1_sd 'theta_hat_inter_m1_sd' as used by this function; see Examples.
#' @param theta_hat_inter_m2 'theta_hat_inter_m2' as used by this function; see Examples.
#' @param theta_hat_inter_m2_sd 'theta_hat_inter_m2_sd' as used by this function; see Examples.
#' @param beta Type II error rate, i.e. 1 - power
#' @param tot_events_inter_m1 'tot_events_inter_m1' as used by this function; see Examples.
#' @param dist Indicator of the endpoint / distribution type
#'
#' @return A named list containing the following elements: 'Interim analysis hypothesis test', 'Wald statistics', 'Critical boundary',
#' 'Assumed effect size: theta', 'remaining effect size', 'c_alpha', 'c_power_m1', 'c_power_m2', 'new events between a test arm and control',
#' 'new total events', 'new information fraction', 'New critical boundary'.
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' masd_new_samsz_survival(s = c(1/3,2/3,1), c_bry = c(3.984298, 3.961146, 2.349682), snew = 1,
#' inter_vt = 2, theta_hat_inter_m1 = c(1.2,1.5,1.1),
#' theta_hat_inter_m1_sd = 1.2, theta_hat_inter_m2 = c(1.2,1.5),
#' theta_hat_inter_m2_sd = 1.2, beta = 0.1, tot_events_inter_m1 = 65, dist = 1)
#'
#' @export
masd_new_samsz_survival<-function(s,c_bry,snew,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,tot_events_inter_m1,dist){
  #this function calculates the new sample size and new critical boundary to achieve desired conditional power of 1-beta#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #the effect size is taken as theta, which is a vector, length(theta)=length(theta_hat_inter)
  #inter_vt numbers the interim analysis at which the analysis was performed
  # the analysis is conducted at the inter_vt'th  interim analysis
  #c_alpha  is the target one-sided type I error#
  # theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  # for normal distribution, dist=1###
  # 1-beta is the desired power
  #n_control_inter is the sample size for the comtrol arm at the interim analysis
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm
  s_inter=s[inter_vt]
  theta=theta_hat_inter_m1

  comp=length(theta_hat_inter_m1)
  new_comp=length(theta_hat_inter_m2)
  mu_keep=match(theta_hat_inter_m2,theta_hat_inter_m1)
  theta_keep=theta[mu_keep]
  T_inter=(mean(theta_hat_inter_m2_sd))^(-2)
  T=T_inter/s_inter
  wald_inter=theta_hat_inter_m1/theta_hat_inter_m1_sd
  wald_inter_max=max(wald_inter)
  test= wald_inter_max-c_bry[inter_vt]

  c_power_m1=masd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  c_alpha=masd_c_power_theta(rep(0,comp),s,c_bry,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  power<-function(sample_r){
    c_bry=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    tmp=masd_c_power_new_theta(theta_keep,s,inter_vt,snew,c_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
    return(tmp)
  }

  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    two_arm_events_inter_m2="NA"
    tot_events_new_m2="NA"
    c_bry_f="NA"
    c_power_m2="NA"
    new_ratio="NA"
    c_power_m2_new_ratio="NA"
  }else {
    outcome="Null hypothesis not rejected at interim analysis"
    if (max(theta)<10e-4){
      two_arm_events_inter_m2="NA"
      tot_events_new_m2="NA"
      c_bry_f="NA"
      c_power_m2="NA"
      new_ratio="NA"
      c_power_m2_new_ratio="NA"
    }else{
      c_power_m2=as.numeric(power(1))
      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.001,2^32),tol=1e-5)[1])
      c_bry_f=new_OF_bundry_s(s,inter_vt,snew,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,new_ratio,dist)
      T_new=T*new_ratio
      two_arm_events_inter_m2=2*tot_events_inter_m1/(comp+1)
      two_arm_events_m2_planned=two_arm_events_inter_m2/s_inter
      events_two_arm_new_m2=ceiling(two_arm_events_m2_planned*(T_new/T))
      tot_events_new_m2=ceiling((new_comp+1)*events_two_arm_new_m2/2)
      c_power_m2_new_ratio=as.numeric(power(new_ratio))
    }
  }

  out=list(outcome,wald_inter,c_bry[inter_vt],theta,theta_keep,c_alpha, c_power_m1,  c_power_m2,events_two_arm_new_m2,tot_events_new_m2,snew,c_bry_f)
  names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Assumed effect size: theta","remaining effect size","c_alpha","c_power_m1","c_power_m2","new events between a test arm and control","new total events","new information fraction","New critical boundary")
  #  out=list(new_comp,new_ratio,outcome,wald_inter,c_bry[inter_vt],theta,theta_keep,c_alpha, c_power_m1,  c_power_m2,events_two_arm_new_m2,tot_events_new_m2,snew,c_bry_f)
  #  names(out)=c("new_comp","new_ratio","Interim analysis hypothesis test","Wald statistics","Critical boundary","Assumed effect size: theta","remaining effect size","c_alpha","c_power_m1","c_power_m2","new events between a test arm and control","new total events","new information fraction","New critical boundary")
  #  out=list(test,new_ratio,c_power_m2,c_power_m2_new_ratio)
  #  names(out)=c("test","new_ratio","c_power_m2","c_power_m2_new_ratio")
  return(out)
}
#################################
#masd_new_samsz_survival(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),snew=1,inter_vt=2,theta_hat_inter_m1=c(1.2,1.5,1.1),theta_hat_inter_m1_sd=1.2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1.2,beta=0.1,tot_events_inter_m1=65,dist=1)
#masd_new_samsz_survival(s=c(0.5,1),c_bry=c(3.179276 ,2.248088 ),snew=1,inter_vt=1,theta_hat_inter_m1=c(0.25,0.35),theta_hat_inter_m1_sd=0.15,theta_hat_inter_m2=c(0.25,0.35),theta_hat_inter_m2_sd=0.15,beta=0.1,tot_events_inter_m1=150,dist=2)
#masd_new_samsz_survival(s=c(0.5,1),c_bry=c(3.179276 ,2.248088 ),snew=1,inter_vt=1,theta_hat_inter_m1=c(0.09,0.25),theta_hat_inter_m1_sd=0.15,theta_hat_inter_m2=0.25,theta_hat_inter_m2_sd=0.15,beta=0.1,tot_events_inter_m1=150,dist=2)
###################################

####################################
masd_cprob_new<-function(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter_m2,theta_hat_inter_m2_sd,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,dist){
  #this function calculates conditional probability
  #associated with the theta_hat_last_ad at the last_vt_ad'th interim analysis,
  #snew is the vector of information fractions after sample size change #
  #c_bry_new is the criotical bpundaries after sample size change##
  #theta is a vector of effect size, same length as theta_hat_inter####
  #theta_hat_last_ad is the vector of the last observed theta_hat afer sample size change###
  #inter_vt numbers the interim analysis at which the analysis was performed
  #s is the vector of information fractions  #
  #theta_hat_inter is the vector of the  observed theta_hat at the inter_vt'th interim analysis###
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #inter_vt numbers the interim analysis at which the analysis was performed
  # last_vt_ad numbers the interim analysis after sample size change at which theta_hat_last_ad was observed
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm

  T_last_ad=(mean(theta_hat_last_ad_m2_sd))^(-2)
  T_inter=1/(mean(theta_hat_inter_m2_sd))^2
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  s_new=((T_new-T_inter)*snew+T_inter)/T_new
  t_new=(T_new-T_inter)*snew+T_inter

  comp=length(theta_hat_inter_m2)

  score_last_ad=theta_hat_last_ad_m2*theta_hat_last_ad_m2_sd^(-2)
  wald_last_ad=theta_hat_last_ad_m2/theta_hat_last_ad_m2_sd
  score_inter=theta_hat_inter_m2*theta_hat_inter_m2_sd^(-2)

  new_e_bry=c_bry_new*sqrt(t_new)
  if (last_vt_ad==1){
    new_e_bry_f=max(score_last_ad)
  }else{
    new_e_bry_f=c(new_e_bry[1:(last_vt_ad-1)],max(score_last_ad))
  }
  #################################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)

  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }
  cov_mtx_ws = matrix(nrow=length(snew)*comp, ncol=length(snew)*comp)
  for (i in 1:length(snew)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(snew)*comp)]=snew[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]

  #############################

  u_bry_mat=matrix(ncol=comp,nrow=last_vt_ad)
  priv_sqrt = sqrt(T_new-T_inter)
  for (i in (1:last_vt_ad)) {
    for (j in 1:comp){
      u_bry_mat[i,j]=(new_e_bry_f[i]-score_inter[j])/priv_sqrt - theta[j]*priv_sqrt*snew[i]
    }
  }
  u_bry=c(t(u_bry_mat))

  ########################################
  ## the function for the vector of Brownian motions ###
  cal1<-function(theta,sigma){
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry,mean=rep(0,(last_vt_ad*comp)),sigma=sigma,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }
  ##################################################
  if (T_new<=T_inter){
    c_power_new=0
  } else {
    c_power_new= cal1(theta,cov_mtx_ws[1:(last_vt_ad*comp),1:(last_vt_ad*comp)])
  }
  #  out=list(T_inter,T_new,c_power_new)
  #  names(out)=c("T_inter","T_new","c_power_new")
  return(c_power_new)
  # return(out)
}
##################################
#masd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter_m2,theta_hat_inter_m2_sd,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,dist)
#masd_cprob_new(theta=c(0,0.3),s=c(1/3,2/3,1),snew=1,c_bry_new=2.545644,inter_vt=2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1,theta_hat_last_ad_m2=c(3,5),theta_hat_last_ad_m2_sd=0.8,last_vt_ad=1,dist=1)
#masd_cprob_new(theta=c(0.3,0.9),s=c(1/3,2/3,1),snew=c(0.5,1),c_bry_new=c(2.545644,2),inter_vt=2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1,theta_hat_last_ad_m2=c(3,5),theta_hat_last_ad_m2_sd=0.8,last_vt_ad=2,dist=1)
#masd_cprob_new(theta=c(0,0.3),s=c(1/3,2/3,1),snew=c(0.5,1),c_bry_new=c(2.545644,2),inter_vt=2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1,theta_hat_last_ad_m2=c(3,5),theta_hat_last_ad_m2_sd=0.8,last_vt_ad=1,dist=1)
#masd_cprob_new(theta=c(0,0.3),s=c(1/3,2/3,1),snew=1,c_bry_new=2.545644,inter_vt=2,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1,theta_hat_last_ad_m2=c(3,5),theta_hat_last_ad_m2_sd=0.8,last_vt_ad=1,dist=1)
#masd_cprob_new(theta=c(0.1,15),s=c(0.3,0.7,1),snew=c(0.5,1),c_bry_new=c(7.512163,7.428292),inter_vt=1,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.8,theta_hat_last_ad_m2=c(7.428292*0.7,7.428292*0.7),theta_hat_last_ad_m2_sd=0.7,last_vt_ad=2,dist=1)
####################################
###########################Calculating backward image#######

masd_bk_img_gen<-function(theta,s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,c_p,dist){
  #this function calculates the backward image given the conditional probability c_p
  #s is the vector of information fractions  #
  #c_bry is the critical bpundaries ##
  #theta is a vector of effect size, same length as theta_hat_inter####
  #theta_hat_inter is the vector of the  observed theta_hat at the inter_vt'th interim analysis###
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #inter_vt numbers the interim analysis at which the analysis was performed
  #c_p is a conditional probability
  # dist is the distribution. dist=1 means normal distribution with a commmon control arm

  score_inter=theta_hat_inter_m1*theta_hat_inter_m1_sd^(-2)
  prob=NULL
  for (q in (inter_vt+1):length(s)){
    prob_q=masd_c_power_theta(theta,s[1:q],c_bry[1:q],inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)
    prob=c(prob,prob_q)
  }
  if (min(prob-c_p)<0){
    p_diff=prob-c_p
    neg=p_diff[p_diff<0]
    cut=match(last(neg),p_diff)
  }else{
    cut=0
  }

  bk_prob_1<-function(c_lst){
    c_bry_new=c(c_bry[1:(cut+inter_vt)],c_lst)
    tmp=masd_c_power_theta(theta,s[1:(cut+inter_vt+1)],c_bry_new,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)

  }

  bk_prob_2<-function(c_lst){
    c_bry_new=c(c_bry[1:(length(s)-1)],c_lst)
    tmp=prob_q=masd_c_power_theta(theta,s,c_bry_new,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)
  }

  if ((cut+inter_vt)<length(s)){
    bk_img=as.numeric(stats::uniroot(function(x) bk_prob_1(x)-c_p,c(-2^32,2^32),tol=1e-5)[1])
    bry_f=c(c_bry[1:(cut+inter_vt)],bk_img)
  } else {
    bk_img=as.numeric(stats::uniroot(function(x) bk_prob_2(x)-c_p,c(-2^32,2^32),tol=1e-5)[1])
    bry_f=c(c_bry[1:(length(s)-1)],bk_img)
  }
  #  out=list(cut,c_bry,bry_f,c(prob,c_p))
  #  names(out)=c("cut","c_bry","bry_f","prob")
  return(bry_f)
}
########################################
#masd_bk_img_gen(theta,s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,c_p,dist)
#masd_bk_img_gen(theta=c(0,0.2),s=c(1/3,2/3,1),c_bry=c(3.888374, 3.864257, 2.212819),theta_hat_inter_m1=c(1.2,1.5),theta_hat_inter_m1_sd=1,inter_vt=1,c_p=0.4,dist=1)
#masd_bk_img_gen(theta=c(0,0.2),s=c(1/3,2/3,1),c_bry=c(3.888374, 3.864257, 2.212819),theta_hat_inter_m1=c(1.2,1.5),theta_hat_inter_m1_sd=1,inter_vt=2,c_p=0.4,dist=1)
#masd_bk_img_gen(theta=c(0,0.2),s=c(1/3,2/3,1),c_bry=c(3.888374, 3.864257, 2.212819),theta_hat_inter_m1=c(1.2,1.5),theta_hat_inter_m1_sd=1,inter_vt=2,c_p=0,dist=1)

##############################################


#' Multiple comparisons Group sequential design: Final analysis, With one sample size change, orgin_comp==comp2
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param snew Vector of new information fractions after the sample size is changed
#' @param c_bry_new Critical boundaries after the sample size is changed
#' @param theta_hat_inter_m1 'theta_hat_inter_m1' as used by this function; see Examples.
#' @param theta_hat_inter_m1_sd 'theta_hat_inter_m1_sd' as used by this function; see Examples.
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based)
#' @param theta_hat_last_ad_m1 'theta_hat_last_ad_m1' as used by this function; see Examples.
#' @param theta_hat_last_ad_m1_sd 'theta_hat_last_ad_m1_sd' as used by this function; see Examples.
#' @param last_vt_ad Index of the last analysis after the sample size change
#' @param alpha_2 Type I error for the two-sided test
#' @param dist Indicator of the endpoint / distribution type
#'
#' @return A named list containing the following elements: 'p-value', 'confidence interval'.
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' masd_est_back_no_drop(s = c(0.3,0.7,1), c_bry = c(4.347011, 2.845787, 2.380956), snew = c(0.5,1),
#' c_bry_new = c(3.310743,2.781890), theta_hat_inter_m1 = c(1.2,1.3,1.5), theta_hat_inter_m1_sd = 1,
#' inter_vt = 1, theta_hat_last_ad_m1 = c(2.2,1.4,1.8), theta_hat_last_ad_m1_sd = 1,
#' last_vt_ad = 2, alpha_2 = 0.05, dist = 1)
#'
#' @export
masd_est_back_no_drop<-function(s,c_bry,snew,c_bry_new,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,theta_hat_last_ad_m1,theta_hat_last_ad_m1_sd,last_vt_ad,alpha_2,dist){
  #this function calculates point estimate and lower confidence limit after sample size change
  # when no doses were droped at the interim analysis
  #s is the infromation fractions for the original design
  #c_bry is the critical bpundaries with the original design
  #snew is the vector of information fractions after sample size change #
  #c_bry_new is the critical boundaries after sample size change##
  #theta_hat_inter is the theta_hat at the interim analysis at which the sample size was changed##
  #theta_hat_inter_sd is the standard error
  #inter_vt numbers the interim analysis at which the sample size was changed
  #theta_hat_last_ad is the vector of the last observed theta_hat afer sample size change###
  #last_vt_ad is the  interim analysis after sample size change at which the trial was terminated
  #alpha_2 is the significance level for the (1-alpha_2)*100% confidence interval
  # typically, alpha_2=0.05
  # dist is the distribution. dist=1 means normal distribution
  #####################################

  comp=length(theta_hat_inter_m1)

  t=T*s
  T_inter=1/(mean(theta_hat_inter_m1_sd))^2
  T=T_inter/s[inter_vt]
  T_last_ad=(mean(theta_hat_last_ad_m1_sd))^(-2)
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  s_new=((T_new-T_inter)*snew+T_inter)/T_new
  t_new=(T_new-T_inter)*snew+T_inter

  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)

  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }

  cov_mtx_ws <- matrix(nrow=length(s)*comp, ncol=length(s)*comp)
  for (i in 1:length(s)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s)*comp)]=s[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
  #################################################

  ####################################################
  cal1<-function(theta){
    theta_vec=rep(theta,comp)
    #masd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,dist)
    c_p=masd_cprob_new(theta_vec,s,snew,c_bry_new,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_last_ad_m1,theta_hat_last_ad_m1_sd,last_vt_ad,dist)
    #masd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p,dist)
    bk_img=masd_bk_img_gen(theta_vec,s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,c_p,dist)
    last_vt_theta=length(bk_img)
    bkimg_wald_theta=last(bk_img)
    ### Calculating exit probability with lst_wald_theta,lst_vt_theta###
    last_vt=last_vt_theta
    wald_last_max=bkimg_wald_theta
    ####################Boundary for vector of Brownian  motions#########################
    if (last_vt_theta==1) {
      u_bry_theta=rep(wald_last_max*sqrt(s[last_vt])-theta*sqrt(T)*s[last_vt],comp)
    }else{
      u_c_bry=c(c_bry[1:(last_vt-1)],wald_last_max)

      u_bry_theta_mat=matrix(ncol=comp,nrow=last_vt)
      for (i in (1:last_vt)) {
        u_bry_theta_mat[i,]=rep(u_c_bry[i]*sqrt(s[i])-theta*sqrt(T)*s[i],comp)
      }
      u_bry_theta=c(t(u_bry_theta_mat))
    }
    ###################################################
    sigma0=cov_mtx_ws_mf=cov_mtx_ws[1:(last_vt*comp),1:(last_vt*comp)]
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,last_vt*comp),sigma=sigma0,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }
  theta_hat_last_ad_m2=theta_hat_last_ad_m1
  theta_hat_last_ad_m2_sd=theta_hat_last_ad_m1_sd
  m_select=match(max(theta_hat_last_ad_m1),theta_hat_last_ad_m1)
  theta_ll=as.numeric(stats::uniroot(function(x) cal1(x)-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal1(x)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul=masd_ucl_est_back(s,c_bry,snew,c_bry_new,m_select,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,alpha_2)
  p_value=cal1(0)
  ##########################################
  theta_ci=c(theta_ll,theta_pt,theta_ul)
  out=list(p_value,theta_ci)
  names(out)=c("p-value","confidence interval")
  return(out)
}
#############################################
#masd_est_back_no_drop(s,c_bry,snew,c_bry_new,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,theta_hat_last_ad_m1,theta_hat_last_ad_m1_sd,last_vt_ad,alpha,dist)
#masd_est_back_no_drop(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=c(0.5,1),c_bry_new=c(3.310743,2.781890),theta_hat_inter_m1=c(1.2,1.3,1.5),theta_hat_inter_m1_sd=1,inter_vt=1,theta_hat_last_ad_m1=c(2.2,1.4,1.8),theta_hat_last_ad_m1_sd=1,last_vt_ad=2,alpha_2=0.05,dist=1)
#masd_est_back_no_drop(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=c(0.5,1),c_bry_new=c(7.512163,7.428292),theta_hat_inter_m1=c(1.2,1.5),theta_hat_inter_m1_sd=0.8,inter_vt=1,theta_hat_last_ad_m1=c(7.428292*0.7,7.428292*0.7),theta_hat_last_ad_m1_sd=0.7,last_vt_ad=2,alpha_2=0.05,dist=1)
#masd_est_back_no_drop(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=c(0.5,1),c_bry_new=c(3.310743,2.781890),theta_hat_inter_m1=c(1.2,1.5),theta_hat_inter_m1_sd=1,inter_vt=1,theta_hat_last_ad_m1=c(3.310743,2.781890),theta_hat_last_ad_m1_sd=0.8,last_vt_ad=2,alpha_2=0.05,dist=1)
#输入masd_est_back_no_drop(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=c(0.5,1),c_bry_new=c(3.310743,2.781890),theta_hat_inter_m1=c(1.2,1.3,1.5),theta_hat_inter_m1_sd=1,inter_vt=1,theta_hat_last_ad_m1=c(2.2,1.4,1.8),theta_hat_last_ad_m1_sd=1,last_vt_ad=2,alpha_2=0.05,dist=1)输出1.9982  ???3.7553,6.307???
##########################################################################
##############################################


masd_est_back_no_drop_com<-function(s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,theta_hat_last_ad_m1,theta_hat_last_ad_m1_sd,alpha_2,dist){
  #this function calculates point estimate and lower confidence limit after sample size change
  # when no doses were droped at the interim analysis
  # this function applies to a common case, when there is only final analysis after sample size change
  # no more interim analysis
  #s is the infromation fractions for the original design
  #c_bry is the critical bpundaries with the original design
  #snew is the vector of information fractions after sample size change #
  # for the common case, snew=1
  #c_bry_new is the critical boundaries after sample size change##
  #theta_hat_inter is the theta_hat at the interim analysis at which the sample size was changed##
  #theta_hat_inter_sd is the standard error
  #for the common case,inter_vt=length(s)-1 numbers the interim analysis at which the sample size was changed
  #theta_hat_last_ad is the vector of the last observed theta_hat afer sample size change###
  #last_vt_ad is the  interim analysis after sample size change at which the trial was terminated
  #alpha_2 is the significance level for the (1-alpha_2)*100% confidence interval
  # typically, alpha_2=0.05
  # dist is the distribution. dist=1 means normal distribution
  #####################################

  comp=length(theta_hat_inter_m1)
  inter_vt=length(s)-1
  last_vt_ad=1
  snew=1
  c_bry_new=3
  t=T*s
  T_inter=1/(mean(theta_hat_inter_m1_sd))^2
  T=T_inter/s[inter_vt]
  T_last_ad=(mean(theta_hat_last_ad_m1_sd))^(-2)
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  s_new=((T_new-T_inter)*snew+T_inter)/T_new
  t_new=(T_new-T_inter)*snew+T_inter
  #############################
  score_last_ad=theta_hat_last_ad_m1*theta_hat_last_ad_m1_sd^(-2)
  wald_last_ad=theta_hat_last_ad_m1/theta_hat_last_ad_m1_sd
  mu_select=match(max(wald_last_ad),wald_last_ad)
  score_last_ad_slct=score_last_ad[mu_select]
  score_inter=theta_hat_inter_m1*theta_hat_inter_m1_sd^(-2)
  score_inter_slct=score_inter[mu_select]
  ##########################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }
  cov_mtx_ws <- matrix(nrow=length(s)*comp, ncol=length(s)*comp)
  for (i in 1:length(s)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s)*comp)]=s[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
  #################################################
  ################ For a single Brownian  motion####
  cov_mtx_s= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx_s[lower.tri(cov_mtx_s)]=t(cov_mtx_s)[lower.tri(cov_mtx_s)]
  #cov_mtx_sf=cov_mtx_s[1:last_vt,1:last_vt]
  #####################################################
  ####################################################
  ## the function cal_m applies to the vector of brownian motions

  cal_m<-function(theta){
    theta_vec=rep(theta,comp)
    #masd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,dist)
    c_p=masd_cprob_new(theta_vec,s,snew,c_bry_new,inter_vt,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_last_ad_m1,theta_hat_last_ad_m1_sd,last_vt_ad,dist)
    #masd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p,dist)
    bk_img=masd_bk_img_gen(theta_vec,s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,c_p,dist)
    last_vt_theta=length(bk_img)
    bkimg_wald_theta=last(bk_img)
    ### Calculating exit probability with lst_wald_theta,lst_vt_theta###
    last_vt=last_vt_theta
    wald_last_max=bkimg_wald_theta
    ####################Boundary for vector of Brownian  motions#########################
    if (last_vt_theta==1) {
      u_bry_theta=rep(wald_last_max*sqrt(s[last_vt])-theta*sqrt(T)*s[last_vt],comp)
    }else{
      u_c_bry=c(c_bry[1:(last_vt-1)],wald_last_max)

      u_bry_theta_mat=matrix(ncol=comp,nrow=last_vt)
      for (i in (1:last_vt)) {
        u_bry_theta_mat[i,]=rep(u_c_bry[i]*sqrt(s[i])-theta*sqrt(T)*s[i],comp)
      }
      u_bry_theta=c(t(u_bry_theta_mat))
    }
    ###################################################
    sigma0=cov_mtx_ws_mf=cov_mtx_ws[1:(last_vt*comp),1:(last_vt*comp)]
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,last_vt*comp),sigma=sigma0,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }

  ########################################################
  ## the function cal_m applies to a single  brownian motions
  cal_s<-function(theta){
    bkimg_score_theta=(sqrt(T-T_inter)/sqrt(T_new-T_inter))*(score_last_ad_slct-score_inter_slct-theta*(T_new-T_inter))+score_inter_slct+theta*(T-T_inter)
    #masd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,dist)
    bkimg_wald_theta=bkimg_score_theta/sqrt(T)
    last_vt=length(s)
    ### ###
    u_c_bry=c(c_bry[1:(last_vt-1)],bkimg_wald_theta)
    ####################Boundary for vector of Brownian  motions#########################
    u_e_bry_theta=u_c_bry*sqrt(s)-theta*sqrt(T)*s
    ###################################################
    tmp = 1-pmvnorm(lower=-Inf,upper=u_e_bry_theta,mean=rep(0,last_vt),sigma=cov_mtx_s,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }

  #########################################################
  theta_ll=as.numeric(stats::uniroot(function(x) cal_m(x)-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal_m(x)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul=as.numeric(stats::uniroot(function(x) cal_s(x)-1+alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  ##########################################
  return(c(theta_ll,theta_pt,theta_ul))
}
#############################################
#masd_est_back_no_drop_com(s,c_bry,c_bry_new,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,theta_hat_last_ad_m1,theta_hat_last_ad_m1_sd,alpha_2,dist){
#masd_est_back_no_drop_com(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),theta_hat_inter_m1=c(1.2,1.3,1.5),theta_hat_inter_m1_sd=1,theta_hat_last_ad_m1=c(2.2,1.4,1.8),theta_hat_last_ad_m1_sd=1,alpha_2=0.05,dist=1)
##########################################################################

##########################################################################
#' Multiple comparisons Group sequential design: Final analysis, With one sample size change, orgin_comp>comp2
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param snew Vector of new information fractions after the sample size is changed
#' @param c_bry_new Critical boundaries after the sample size is changed
#' @param theta_hat_inter_m1 'theta_hat_inter_m1' as used by this function; see Examples.
#' @param theta_hat_inter_m1_sd 'theta_hat_inter_m1_sd' as used by this function; see Examples.
#' @param theta_hat_inter_m2 'theta_hat_inter_m2' as used by this function; see Examples.
#' @param theta_hat_inter_m2_sd 'theta_hat_inter_m2_sd' as used by this function; see Examples.
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based)
#' @param theta_hat_last_ad_m2 'theta_hat_last_ad_m2' as used by this function; see Examples.
#' @param theta_hat_last_ad_m2_sd 'theta_hat_last_ad_m2_sd' as used by this function; see Examples.
#' @param last_vt_ad Index of the last analysis after the sample size change
#' @param alpha_2 Type I error for the two-sided test
#' @param dist Indicator of the endpoint / distribution type
#'
#' @return A named list containing the following elements: 'p-value', 'confidence interval'.
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' masd_est_back_with_drop(s = c(0.3,0.7,1), c_bry = c(4.347011,2.845787,2.380956), snew = c(0.5,1),
#' c_bry_new = c(3.310743,2.781890), theta_hat_inter_m1 = c(1.2,1.3,1.5), theta_hat_inter_m1_sd = 1,
#' theta_hat_inter_m2 = c(1.3,1.5), theta_hat_inter_m2_sd = 1, inter_vt = 1,
#' theta_hat_last_ad_m2 = c(2.2,1.4), theta_hat_last_ad_m2_sd = 1, last_vt_ad = 2,
#' alpha_2 = 0.05, dist = 1)
#'
#' @export
masd_est_back_with_drop<-function(s,c_bry,snew,c_bry_new,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,inter_vt,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,alpha_2,dist){
  #this function calculates point estimate and the lower limit of confidence interval after sample size change , with dose dropping
  #s is the infromation fractions for the original design
  #c_bry is the critical bpundaries with the original design
  #snew is the vector of information fractions after sample size change #
  #c_bry_new is the criotical boundaries after sample size change##
  #theta_hat_inter is the theta_hat at the interim analysis at which the sample size was changed##
  #reduced_theta_hat_inter is the theta_hat for the remaining doses after others were dripped at inter_vt for lack of efficacy
  #theta_hat_inter_sd is the standard error
  #inter_vt is the interim analysis at twhich the sample size was changed
  #theta_hat_last_ad is the vector of the last observed theta_hat afer sample size change###
  #last_vt_ad is the  interim analysis after sample size change at which the trial was terminated
  #alpha_2 is the significance level for the (1-alpha_2)*100% confidence interval
  # typically, alpha_2=0.05
  # dist is the distribution. dist=1 means normal distribution
  ##################################################
  #### high, orginal dimensional statistics #######
  t=T*s
  T_inter=1/(mean(theta_hat_inter_m1_sd))^2
  T=T_inter/s[inter_vt]
  T_last_ad=(mean(theta_hat_last_ad_m2_sd))^(-2)
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  s_new=((T_new-T_inter)*snew+T_inter)/T_new
  t_new=(T_new-T_inter)*snew+T_inter
  ##################################
  comp_h_dim=length(theta_hat_inter_m1)
  comp_l_dim=length(theta_hat_inter_m2)

  #########################################
  #### the next section constructs the covariance matrix###
  comp=comp_h_dim

  cov_mtx=diag(comp)
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }

  cov_mtx_ws <- matrix(nrow=length(s)*comp, ncol=length(s)*comp)
  for (i in 1:length(s)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s)*comp)]=s[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
  #######################
  ################ For a single Brownian  motion####
  ################ For a single Brownian  motion####
  cov_mtx_s= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx_s[lower.tri(cov_mtx_s)]=t(cov_mtx_s)[lower.tri(cov_mtx_s)]
  #######################

  cal1<-function(theta){
    theta_vec_l_dim=rep(theta,comp_l_dim)
    theta_vec_h_dim=rep(theta,comp_h_dim)
    s_inter=s[inter_vt]
    s_0_new=s_inter*(T/T_new)
    ####################################################
    #masd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,dist)

    c_p=masd_cprob_new(theta_vec_l_dim,s,snew,c_bry_new,inter_vt,theta_hat_inter_m2,theta_hat_inter_m2_sd,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,dist)
    #masd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p,dist)

    bk_img=masd_bk_img_gen(theta_vec_h_dim,s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,c_p,dist)
    last_vt_theta=length(bk_img)
    bkimg_wald_theta=last(bk_img)
    #############################
    last_vt=last_vt_theta
    wald_last_max=bkimg_wald_theta
    ####################Boundary for vector of Brownian  motions#########################
    if (last_vt_theta==1) {
      u_bry_theta=rep(wald_last_max*sqrt(s[last_vt])-theta*sqrt(T)*s[last_vt],comp)
    }else{
      u_c_bry=c(c_bry[1:(last_vt-1)],wald_last_max)

      u_bry_theta_mat=matrix(ncol=comp,nrow=last_vt)
      for (i in (1:last_vt)) {
        u_bry_theta_mat[i,]=rep(u_c_bry[i]*sqrt(s[i])-theta*sqrt(T)*s[i],comp)
      }
      u_bry_theta=c(t(u_bry_theta_mat))
    }
    ###################################################
    sigma0=cov_mtx_ws_mf=cov_mtx_ws[1:(last_vt*comp),1:(last_vt*comp)]
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,last_vt*comp),sigma=sigma0,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }

  #####################################
  cal2<-function(theta){
    theta_vec_l_dim=rep(theta,comp_l_dim)
    theta_vec_h_dim=rep(theta,comp_h_dim)
    s_inter=s[inter_vt]
    s_0_new=s_inter*(T/T_new)
    ###################################################
    c_p=masd_cprob_new(theta_vec_l_dim,s,snew,c_bry_new,inter_vt,theta_hat_inter_m2,theta_hat_inter_m2_sd,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,dist)
    #masd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,dist)

    bk_img=masd_bk_img_gen(theta_vec_h_dim,s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,c_p,dist)
    last_vt_theta=length(bk_img)
    bkimg_wald_theta=last(bk_img)
    #############################
    last_vt=last_vt_theta
    wald_last_max=bkimg_wald_theta
    ##################Boundary for single Brownian motion ######################
    if (last_vt==1){
      up_bry_sf_theta=wald_last_max-theta*sqrt(T)*s[last_vt]
    }else{
      u_c_bry=c(c_bry[1:(last_vt-1)],wald_last_max)
      up_bry_sf_theta=u_c_bry*sqrt(s[1:last_vt])-theta*sqrt(T)*s[1:last_vt]
    }

    sigma0=cov_mtx_sf=cov_mtx_s[1:last_vt,1:last_vt]
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry_sf_theta,mean=rep(0,last_vt),sigma=sigma0,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }

  ######################################
  ########################################
  mu_keep=match(theta_hat_inter_m2,theta_hat_inter_m1)
  last_select=match(max(theta_hat_last_ad_m2),theta_hat_last_ad_m2)
  m_select=mu_keep[last_select]

  theta_ll=as.numeric(stats::uniroot(function(x) cal1(x)-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal1(x)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul=masd_ucl_est_back(s,c_bry,snew,c_bry_new,m_select,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,alpha_2)
  p_value=cal1(0)

  theta_ci=c(theta_ll,theta_pt,theta_ul)
  out=list(p_value,theta_ci)
  names(out)=c("p-value","confidence interval")
  return(out)
}
#############################################
#masd_est_back_with_drop(s,c_bry,snew,c_bry_new,theta_hat_inter_m1,theta_hat_inter_m1_sd,theta_hat_inter_m2,theta_hat_inter_m2_sd,inter_vt,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,alpha_2,dist)
#masd_est_back_with_drop(s=c(0.3,0.7,1),c_bry=c(4.347011,2.845787,2.380956),snew=c(0.5,1),c_bry_new=c(3.310743,2.781890),theta_hat_inter_m1=c(1.2,1.3,1.5),theta_hat_inter_m1_sd=1,theta_hat_inter_m2=c(1.3,1.5),theta_hat_inter_m2_sd=1,inter_vt=1,theta_hat_last_ad_m2=c(2.2,1.4),theta_hat_last_ad_m2_sd=1,last_vt_ad=2,alpha_2=0.05,dist=1)
#masd_est_back_with_drop(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=1,c_bry_new=2.781890,theta_hat_inter_m1=c(1.2,1.3,1.5),theta_hat_inter_m1_sd=1,theta_hat_inter_m2=c(1.3,1.5),theta_hat_inter_m2_sd=1,inter_vt=1,theta_hat_last_ad_m2=c(2.2,1.4),theta_hat_last_ad_m2_sd=1,last_vt_ad=1,alpha_2=0.05,dist=1)
#masd_est_back_with_drop(s=c(0.5,1),c_bry=c(3.179276 ,2.248088),snew=1,c_bry_new= 2.170774,theta_hat_inter_m1=c(0.09,0.25),theta_hat_inter_m1_sd=0.15,theta_hat_inter_m2=0.25,theta_hat_inter_m2_sd=0.15,inter_vt=1,theta_hat_last_ad_m2=0.261,theta_hat_last_ad_m2_sd=0.12,last_vt_ad=1,alpha_2=0.05,dist=1)
#masd_est_back_with_drop(s=c(0.5,1),c_bry=c(3.179276 ,2.248088),snew=1,c_bry_new= 2.170774,theta_hat_inter_m1=c(0.09,0.25),theta_hat_inter_m1_sd=0.15,theta_hat_inter_m2=0.25,theta_hat_inter_m2_sd=0.15,inter_vt=1,theta_hat_last_ad_m2=0.205527,theta_hat_last_ad_m2_sd=0.09467914,last_vt_ad=1,alpha_2=0.05,dist=1)
#masd_est_back_with_drop(s=c(0.5,1),c_bry=c(3.179276 ,2.248088),snew=1,c_bry_new= 1.8174,theta_hat_inter_m1=c(0.09,0.25),theta_hat_inter_m1_sd=0.15,theta_hat_inter_m2=0.25,theta_hat_inter_m2_sd=0.15,inter_vt=1,theta_hat_last_ad_m2=0.18174,theta_hat_last_ad_m2_sd=0.1,last_vt_ad=1,alpha_2=0.05,dist=1)

##########################################


masd_est_back_with_drop_com<-function(s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,mu_m2,theta_hat_inter_m2,theta_hat_inter_m2_sd,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,alpha_2,dist){
  #this function calculates point estimate and the lower limit of confidence interval after sample size change , with dose dropping
  #s is the infromation fractions for the original design
  #c_bry is the critical bpundaries with the original design
  #snew is the vector of information fractions after sample size change #
  #c_bry_new is the criotical boundaries after sample size change##
  #theta_hat_inter is the theta_hat at the interim analysis at which the sample size was changed##
  #reduced_theta_hat_inter is the theta_hat for the remaining doses after others were dripped at inter_vt for lack of efficacy
  #theta_hat_inter_sd is the standard error
  #inter_vt is the interim analysis at twhich the sample size was changed
  #theta_hat_last_ad is the vector of the last observed theta_hat afer sample size change###
  #last_vt_ad is the  interim analysis after sample size change at which the trial was terminated
  #alpha_2 is the significance level for the (1-alpha_2)*100% confidence interval
  # typically, alpha_2=0.05
  # dist is the distribution. dist=1 means normal distribution
  #### high, orginal dimensional statistics #######

  inter_vt=length(s)-1
  snew=1
  last_vt_ad=1
  c_bry_new=3
  t=T*s
  T_inter=1/(mean(theta_hat_inter_m1_sd))^2
  T=T_inter/s[inter_vt]
  T_last_ad=(mean(theta_hat_last_ad_m2_sd))^(-2)
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  s_new=((T_new-T_inter)*snew+T_inter)/T_new
  t_new=(T_new-T_inter)*snew+T_inter
  #############################
  score_last_ad=theta_hat_last_ad_m2*theta_hat_last_ad_m2_sd^(-2)
  wald_last_ad=theta_hat_last_ad_m2/theta_hat_last_ad_m2_sd
  mu_select_ad=match(max(wald_last_ad),wald_last_ad)
  mu_select=mu_m2[mu_select_ad]
  score_last_ad_slct=score_last_ad[mu_select]
  score_inter=theta_hat_inter_m1*theta_hat_inter_m1_sd^(-2)
  score_inter_slct=score_inter[mu_select]
  ##########################################
  ##################################
  comp_h_dim=length(theta_hat_inter_m1)
  comp_l_dim=length(theta_hat_inter_m2)
  ########################################
  #### the next section constructs the covariance matrix###
  comp=comp_h_dim

  cov_mtx=diag(comp)
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }

  cov_mtx_ws <- matrix(nrow=length(s)*comp, ncol=length(s)*comp)
  for (i in 1:length(s)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s)*comp)]=s[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
  #######################
  ################ For a single Brownian  motion####
  cov_mtx_s= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx_s[lower.tri(cov_mtx_s)]=t(cov_mtx_s)[lower.tri(cov_mtx_s)]
  #######################

  cal_m<-function(theta){
    theta_vec_l_dim=rep(theta,comp_l_dim)
    theta_vec_h_dim=rep(theta,comp_h_dim)
    s_inter=s[inter_vt]
    s_0_new=s_inter*(T/T_new)
    ####################################################
    #masd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,dist)

    c_p=masd_cprob_new(theta_vec_l_dim,s,snew,c_bry_new,inter_vt,theta_hat_inter_m2,theta_hat_inter_m2_sd,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,dist)
    #masd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p,dist)

    bk_img=masd_bk_img_gen(theta_vec_h_dim,s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,c_p,dist)
    last_vt_theta=length(bk_img)
    bkimg_wald_theta=last(bk_img)
    #############################
    last_vt=last_vt_theta
    wald_last_max=bkimg_wald_theta
    ####################Boundary for vector of Brownian  motions#########################
    if (last_vt_theta==1) {
      u_bry_theta=rep(wald_last_max*sqrt(s[last_vt])-theta*sqrt(T)*s[last_vt],comp)
    }else{
      u_c_bry=c(c_bry[1:(last_vt-1)],wald_last_max)

      u_bry_theta_mat=matrix(ncol=comp,nrow=last_vt)
      for (i in (1:last_vt)) {
        u_bry_theta_mat[i,]=rep(u_c_bry[i]*sqrt(s[i])-theta*sqrt(T)*s[i],comp)
      }
      u_bry_theta=c(t(u_bry_theta_mat))
    }
    ###################################################
    sigma0=cov_mtx_ws_mf=cov_mtx_ws[1:(last_vt*comp),1:(last_vt*comp)]
    tmp = 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,last_vt*comp),sigma=sigma0,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }
  #####################################

  cal_s<-function(theta){

    bkimg_score_theta=(sqrt(T-T_inter)/sqrt(T_new-T_inter))*(score_last_ad_slct-score_inter_slct-theta*(T_new-T_inter))+score_inter_slct+theta*(T-T_inter)
    #masd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad,dist)
    bkimg_wald_theta=bkimg_score_theta/sqrt(T)
    last_vt=length(s)
    ### ###
    u_c_bry=c(c_bry[1:(last_vt-1)],bkimg_wald_theta)
    ####################Boundary for vector of Brownian  motions#########################
    u_e_bry_theta=u_c_bry*sqrt(s)-theta*sqrt(T)*s

    ###################################################
    tmp = 1-pmvnorm(lower=-Inf,upper=u_e_bry_theta,mean=rep(0,last_vt),sigma=cov_mtx_s,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)

  }
  ######################################
  theta_ll=as.numeric(stats::uniroot(function(x) cal_m(x)-alpha_2/2,c(-2^32,2^32),tol=1e-6)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal_m(x)-0.5,c(-2^32,2^32),tol=1e-6)[1])
  theta_ul=as.numeric(stats::uniroot(function(x) cal_s(x)-1+alpha_2/2,c(-2^32,2^32),tol=1e-6)[1])
  theta_ci=c( theta_ll,theta_pt,theta_ul)
  return(theta_ci)
}
#############################################
#masd_est_back_with_drop_com(s,c_bry,theta_hat_inter_m1,theta_hat_inter_m1_sd,mu_m2,theta_hat_inter_m2,theta_hat_inter_m2_sd,theta_hat_last_ad_m2,theta_hat_last_ad_sd_m2,alpha_2,dist){
#masd_est_back_with_drop_com(s=c(0.3,0.7,1),c_bry=c(4.347011,2.845787,2.380956),theta_hat_inter=c(1.2,1.6,1.5),theta_hat_inter_m1_sd=1,mu_m2=c(1,3),theta_hat_inter_m2=c(1.3,1.5),theta_hat_inter_m2_sd=1,theta_hat_last_ad_m2=c(2.2,1.4),theta_hat_last_ad_m2_sd=0.6,alpha_2=0.05,dist=1)
#masd_est_back_with_drop_com(s=c(0.3,0.7,1),c_bry=c(4.347011,2.845787,2.380956),theta_hat_inter=c(1.2,1.6,1.5),theta_hat_inter_m1_sd=1,mu_m2=c(2,3),theta_hat_inter_m2=c(1.6,1.5),theta_hat_inter_m2_sd=1,theta_hat_last_ad=c(2.2,1.4),theta_hat_last_ad_sd=0.6,alpha_2=0.05,dist=1)

##########################################


####################################################################

##############################################
masd_ucl_est_back<-function(s,c_bry,snew,c_bry_new,m_select,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,alpha_2){
  #this function calculates point estimatre and confidence interval after sample size change
  #s is the infromation fractions for the original design
  #c_bry is the critical bpundaries with the original design
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
  theta_hat_inter=theta_hat_inter_m1[m_select]
  if (length(theta_hat_inter_m1_sd)==1) {
    theta_hat_inter_sd=theta_hat_inter_m1_sd
  }else {
    theta_hat_inter_sd=theta_hat_inter_m1_sd[m_select]
  }
  score_last_ad_m2=theta_hat_last_ad_m2*theta_hat_last_ad_m2_sd^(-2)
  m2_slct=match(max(score_last_ad_m2),score_last_ad_m2)
  theta_hat_last_ad=theta_hat_last_ad_m2[m2_slct]
  if (length(theta_hat_last_ad_m2_sd)==1) {
    theta_hat_last_ad_sd=theta_hat_last_ad_m2_sd
  }else {
    theta_hat_last_ad_sd=mean(theta_hat_last_ad_m2_sd)
  }
  T_inter=theta_hat_inter_sd^(-2)
  T=T_inter/s[inter_vt]
  t=T*s
  T_last_ad=theta_hat_last_ad_sd^(-2)
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  #T is the total information time for the original design
  #T_new is the new total information time with sampe size change

  ############################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  ####################################################
  #######################
  #(theta,snew,c_bry_new,s_old_lst,inter_wald,lst_wald_new,lst_vt_new,comp,T_old,T_new,dist){

  cal<-function(theta){
    #asd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad)
    c_p=asd_cprob_new(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad)
    #asd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p)
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

      u_bry_theta=(u_c_bry*sqrt(s[1:last_vt]))-theta*sqrt(T)*s[1:last_vt]
    }

    ###################################################
    sigma0=cov_mtxf=cov_mtx[1:last_vt,1:last_vt]
    if (last_vt==1){
      tmp = 1-stats::pnorm(u_bry_theta)
    }else {
      tmp = 1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,last_vt),sigma=sigma0,abseps=1e-8,algorithm=mvtnorm::GenzBretz())[1]
    }
    return(tmp)
  }

  ##########################################################
  theta_ul=as.numeric(stats::uniroot(function(x) cal(x)-1+alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  ##########################################
  return(theta_ul)
}
#############################################
#masd_ucl_est_back(s,c_bry,snew,c_bry_new,m_select,theta_hat_inter_m1,theta_hat_inter_m1_sd,inter_vt,theta_hat_last_ad_m2,theta_hat_last_ad_m2_sd,last_vt_ad,alpha_2)
#masd_ucl_est_back(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=c(0.5,1),c_bry_new=c(3.310743, 2.781890),m_select=2,theta_hat_inter_m1=c(1.2,1.3,1.5),theta_hat_inter_m1_sd=1,inter_vt=1,theta_hat_last_ad_m2=c(2.2,1.4),theta_hat_last_ad_m2_sd=0.5,last_vt_ad=2,alpha_2=0.05)
#masd_ucl_est_back(s=c(0.3,0.7,1),c_bry=c(4.347011, 2.845787, 2.380956),snew=c(0.5,1),c_bry_new=c(3.310743,2.781890),theta_hat_inter_m1=c(1.2,1.3,1.5),theta_hat_inter_m1_sd=1,inter_vt=1,theta_hat_last_ad_m1=c(2.2,1.4,1.8),theta_hat_last_ad_m1_sd=1,last_vt_ad=2,alpha_2=0.05,dist=1)输出1.9982  ???3.7553,6.307???

##########################################################################
