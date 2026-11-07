asd_c_power_theta<-function(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd){
  #this function calculates conditional power, given that inter_vt,theta_hat_inter,theta_hat_inter_sd have been observed at interim analysis inter_vt#
  # if the effect size is given by theta
  #s is the vector of information fractions.
  #theta_hat_inter is the theta_hat at the interim analysis inter_vt
  # theta_hat_inter_sd is the standard error of theta_hat_inter
  s1=(s[(inter_vt+1):length(s)]-s[inter_vt])/(1-s[inter_vt])
  theta_hat=theta_hat_inter
  T_inter=1/theta_hat_inter_sd^2
  score_inter=theta_hat_inter*T_inter
  T=T_inter/s[inter_vt]
  t=T*s
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(s1,length(s1)),nrow=length(s1), ncol=length(s1))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################
  u_bry_theta=NULL

  for (i in (1:length(s1))) {

    u_bry_theta[i]=((c_bry*sqrt(t))[i+inter_vt]-score_inter)/sqrt(T-T_inter)-theta*sqrt(T-T_inter)*s1[i]
  }

  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries#

  cal1<-function(theta,sigma){

    tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,length(s1)),sigma=sigma,algorithm=mvtnorm::Miwa())[1]
    return(tmp)
  }
  ### the function for single Brownian motion###

  c_power= cal1(theta,cov_mtx)

  return(c_power)
  #  return(out)
}
########################################
#asd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)
#asd_c_power_theta(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=1,theta_hat_inter=0.3,theta_hat_inter_sd=0.5)
###########################################################


asd_c_power<-function(s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd){
  #this function calculates conditional power, given that inter_vt,theta_hat_inter,theta_hat_inter_sd have been observed at interim analysis inter_vt#
  # if the effect size is assumed to be equal to theta_hat_inter
  #s is the vector of information fractions.
  #theta_hat_inter is the theta_hat at the interim analysis inter_vt
  # theta_hat_inter_sd is the standard error of theta_hat_inter
  s1=(s[(inter_vt+1):length(s)]-s[inter_vt])/(1-s[inter_vt])
  theta=theta_hat_inter
  theta_hat=theta_hat_inter
  T_inter=1/theta_hat_inter_sd^2
  score_inter=theta_hat_inter*T_inter
  T=T_inter/s[inter_vt]
  t=T*s
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(s1,length(s1)),nrow=length(s1), ncol=length(s1))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################
  u_bry_theta=NULL

  for (i in (1:length(s1))) {
    u_bry_theta[i]=((c_bry*sqrt(t))[i+inter_vt]-score_inter)/sqrt(T-T_inter)-theta*sqrt(T-T_inter)*s1[i]
  }
  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries#

  cal1<-function(theta,sigma){

    tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,length(s1)),sigma=sigma,algorithm=mvtnorm::Miwa())[1]
    return(tmp)
  }
  ### the function for single Brownian motion###

  c_power= cal1(theta,cov_mtx)

  return(c_power)
  #  return(out)
}
########################################
#asd_c_power(s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)
#asd_c_power(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=1,theta_hat_inter=0.3,theta_hat_inter_sd=0.5)
###########################################################

###########################################################
new_OF_bundry<-function(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r){
  #this function calculates the new critical bundaries, given the conditional type I error of c_alpha#
  #the boundary is O'Brien-Fleming type#
  #s is the original information fraction
  #c_bry is the critical boundary
  #snew is the vector of new information fractions. The length of s can be either 2 or 3#
  ## sample_r is the ratio of the new sample size over the original sample size##
  #inter_vt is the where  the interim analysis is conducted
  # theta_hat_inter is the  last observed theta_hat at the interim analysis inter_vt###
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  s_inter=s[inter_vt]
  T_inter=theta_hat_inter_sd^(-2)
  T_old=T_inter/s_inter
  T_new=T_old*sample_r
  t_new=(T_new-T_inter)*snew+T_inter
  score_inter=theta_hat_inter*theta_hat_inter_sd^(-2)
  ###################################
  cov_mtx= matrix(rep(snew,length(snew)),nrow=length(snew), ncol=length(snew))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e0#
  c_cal<-function(e_0,sigma){
    #### set up boundary for given e_0
    if (length(snew)==1){
      u_bry=(e_0-score_inter)/sqrt(T_new-T_inter)
      tmp = 1-stats::pnorm(u_bry)
    }else{
      u_bry=c(rep((e_0-score_inter)/sqrt(T_new-T_inter),length(snew)))
      tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=u_bry,mean=rep(0,length(snew)),sigma=sigma,abseps=1e-8,algorithm=mvtnorm::Miwa())[1]
    }
    return(tmp)
  }
  e_bry=as.numeric(stats::uniroot(function(x) c_cal(x,cov_mtx)-c_alpha,c(-2^32,2^32),tol=1e-5)[1])
  #e_bry are exit boundaries
  c_bry=e_bry/sqrt(t_new)
  return(c_bry)
}
##############################
#new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
#new_OF_bundry(s=c(1/3,2/3,1), inter_vt=1,snew=1,c_alpha=0.025,theta_hat_inter=1.2,theta_hat_inter_sd=1.5,sample_r=1.5)
#new_OF_bundry(s=c(1/3,2/3,1), inter_vt=1,snew=1,c_alpha=0.025,theta_hat_inter=1.2,theta_hat_inter_sd=0.2,sample_r=1.5)

#new_OF_bundry(s=c(1/3,2/3,1), inter_vt=2,snew=c(0.5,1),c_alpha=0.025,theta_hat_inter=1.2,theta_hat_inter_sd=1.5,sample_r=1.5)
#new_OF_bundry(s=c(1/3,2/3,1), inter_vt=1,snew=c(0.5,1),c_alpha=0.025,theta_hat_inter=1.2,theta_hat_inter_sd=0.2,sample_r=1.5)

##########################################################

####################################################
asd_c_power_new<-function(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r){
  #this function calculates the conditional power with the new critical boundary#
  #theta =effet size
  #snew= new relative infromation fraction###
  # new_c_bry is the new critical bounddary
  #sample_r=T_new/T_old
  ## T=T_old
  #snew=c(0.5,1)
  #s_inter=0.4
  #c_bry=c(3.310743, 2.781890)
  s_inter=s[inter_vt]
  T_inter=theta_hat_inter_sd^(-2)
  T_old=T_inter/s_inter
  T_new=T_old*sample_r

  t_new=(T_new-T_inter)*snew+T_inter
  score_inter=theta_hat_inter*theta_hat_inter_sd^(-2)

  new_e_bry=new_c_bry*sqrt(t_new)
  #################################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(snew,length(snew)),nrow=length(snew), ncol=length(snew))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #############################

  u_bry=NULL
  if (length(snew)==1) {
    u_bry=(new_e_bry-score_inter)/sqrt(T_new-T_inter)-theta*sqrt(T_new-T_inter)

  }else{
    for (i in (1:length(snew))) {
      u_bry[i]=(new_e_bry[i]-score_inter)/sqrt(T_new-T_inter)-theta*sqrt(T_new-T_inter)*snew[i]
    }
  }

  ########################################

  cal1<-function(theta,sigma){

    tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=u_bry,mean=rep(0,(length(snew))),sigma=sigma,abseps=1e-8,algorithm=mvtnorm::Miwa())[1]
    return(tmp)
  }
  ###########################################
  if (T_inter>=T_new){c_power_new=0} else {

    if (length(snew)==1){
      c_power_new=1-stats::pnorm(u_bry)
    }else{
      c_power_new=cal1(theta,cov_mtx)
    }
  }

  return(c_power_new)

}
##################################
#asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
#asd_c_power_new(s=c(1/3,2/3,1), inter_vt=1,theta=0.2,snew=1, new_c_bry=2.781890,theta_hat_inter=1.2,theta_hat_inter_sd=0.4,sample_r=2)
#asd_c_power_new(s=c(1/3,2/3,1), inter_vt=1,theta=0.3,snew=c(0.5,1),new_c_bry=c(3.310743, 2.781890),theta_hat_inter=1.2,theta_hat_inter_sd=0.4,sample_r=2)
##################################

asd_new_info_time<-function(theta,s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,beta){
  #this function calculates the new information time to achieve desired conditional power of 1-beta#
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #c_alpha  is the target one-sided type I error#
  # theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  ##where sample size is being re-estimated statistics###
  s_inter=s[inter_vt]
  T_inter=1/theta_hat_inter_sd^2
  T=T_inter/s[inter_vt]

  power<-function(sample_r){
    #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    new_c_bry=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    #asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    tmp=asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    return(tmp)
  }

  new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.0001,2^32),tol=1e-5)[1])
  #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
  c_bry_f=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,new_ratio)
  T_new=T*new_ratio
  #  N_new=ceiling(T_new*2*sigma^2)
  out=list( T_new,c_bry_f)
  names(out)=c("T_new","New critical boundary")
  #names(out)=c("New information fraction","new critical boundaries")
  return(out)
}
#################################
#asd_new_info_time(theta,s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,beta)
#asd_new_info_time(theta=0.2,s=c(1/3,2/3,1), inter_vt=1,snew=1,c_alpha=0.4,theta_hat_inter=1.2,theta_hat_inter_sd=1.5,beta=0.1)
#asd_new_info_time(theta=0.2,s=c(1/3,2/3,1), inter_vt=1,snew=c(0.5,1),c_alpha=0.4,theta_hat_inter=1.2,theta_hat_inter_sd=0.2,beta=0.1)
####################################


##################################
#' Group sequential designs / Adaptive sequential designs: Interim analysis, Non-survival analysis, Input theta
#'
#' @param theta Assumed effect size used for the conditional power calculation; different values of
#'   theta can be entered to check the sensitivity of the conditional power to the effect size; if
#'   theta < 0.001 the new sample size and the new critical boundary are returned as 'NA'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based);
#'   must be less than length(s)
#' @param snew Vector of new information fractions after the sample size is changed, e.g. c(0.5, 1);
#'   the fractions are relative to the new total information and the last element must be 1
#' @param theta_hat_inter Treatment effect estimate observed at the interim analysis, e.g. the
#'   difference of the observed means or event rates between the test arm and the control arm
#' @param theta_hat_inter_sd Standard error of 'theta_hat_inter'
#' @param n_control_inter Sample size of the control arm at the interim analysis
#' @param beta Type II error rate, i.e. 1 - power
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Interim analysis hypothesis test}{Character string giving the outcome of the interim
#'     hypothesis test: whether the null hypothesis was rejected at the interim analysis}
#'   \item{Wald statistics}{Wald statistic of the treatment effect observed at the interim analysis}
#'   \item{Critical boundary}{Planned critical boundary at the interim analysis at which the sample
#'     size is re-estimated}
#'   \item{Assumed effect size: theta}{The assumed effect size theta used for the conditional power
#'     calculation}
#'   \item{c_alpha}{Conditional type I error at the interim analysis, calculated by setting theta = 0}
#'   \item{c_power_theta}{Conditional power at the interim analysis using the assumed effect size theta}
#'   \item{N_control_new}{New sample size for the control arm after the sample size re-estimation,
#'     or 'NA' if the null hypothesis was rejected at the interim analysis or if theta < 0.001}
#'   \item{new information fraction}{The new information fractions after the sample size change}
#'   \item{New critical boundary}{New critical boundaries for the new interim analyses
#'     (O'Brien-Fleming type), or 'NA' if they were not calculated}
#' }
#'
#' @examples
#' asd_new_design_theta(theta = 0.2, s = c(1/3, 2/3, 1), c_bry = c(3.984298, 3.961146, 2.349682),
#' inter_vt = 1, snew = c(0.5, 1), theta_hat_inter = 0.3, theta_hat_inter_sd = 0.5,
#' n_control_inter = 100, beta = 0.1)
#'
#' @export
asd_new_design_theta<-function(theta,s,c_bry,inter_vt,snew,theta_hat_inter,theta_hat_inter_sd,n_control_inter,beta){
  #this function calculates the new sample size to achieve desired conditional power of 1-beta#
  #theta is the effect size
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #c_alpha  is the target one-sided type I error#
  #theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #n_control_inter is the sample size of the control arm at the interim analysis
  ##where sample size is being re-estimated statistics###
  s_inter=s[inter_vt]
  T_inter=1/theta_hat_inter_sd^2
  T=T_inter/s[inter_vt]
  wald_inter=theta_hat_inter/theta_hat_inter_sd
  test= wald_inter-c_bry[inter_vt]

  c_power_theta=asd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)

  c_alpha=asd_c_power_theta(0,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)
  power<-function(sample_r){
    #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    new_c_bry=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    #asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    tmp=asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    return(tmp)
  }

  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    N_control_new="NA"
    c_bry_f="NA"

  }else {
    outcome="Null hypothesis not rejected at interim analysis"

    if (theta<10e-4){
      N_control_new="NA"
      c_bry_f="NA"
    }else{

      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.0001,2^32),tol=1e-5)[1])
      #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
      c_bry_f=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,new_ratio)
      T_new=T*new_ratio
      N_control_new=ceiling(n_control_inter*(T_new/T_inter))
      #  N_new=ceiling(T_new*2*sigma^2)
    }
  }
  out=list(outcome,wald_inter,c_bry[inter_vt],theta,c_alpha, c_power_theta,N_control_new,snew,c_bry_f)
  names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Assumed effect size: theta","c_alpha","c_power_theta","N_control_new","new information fraction","New critical boundary")

  return(out)

}
#################################
#asd_new_design_theta(theta,s,inter_vt,snew,theta_hat_inter,theta_hat_inter_sd,n_control_inter,beta)
#asd_new_design_theta(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=1,snew=1,theta_hat_inter=1.2,theta_hat_inter_sd=1.5,n_control_inter=100,beta=0.1)
#asd_new_design_theta(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=1.2,theta_hat_inter_sd=0.2,n_control_inter=100,beta=0.1)
#asd_new_design_theta(theta=-0.2,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=1.2,theta_hat_inter_sd=0.2,n_control_inter=100,beta=0.1)
#asd_new_design_theta(theta=-0.2,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=-1.2,theta_hat_inter_sd=0.2,n_control_inter=100,beta=0.1)

####################################
####################################
#' Group sequential designs / Adaptive sequential designs: Interim analysis, Non-survival analysis, Use theta_interim
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based);
#'   must be less than length(s)
#' @param snew Vector of new information fractions after the sample size is changed, e.g. c(0.5, 1);
#'   the fractions are relative to the new total information and the last element must be 1
#' @param theta_hat_inter Treatment effect estimate observed at the interim analysis, e.g. the
#'   difference of the observed means or event rates between the test arm and the control arm
#' @param theta_hat_inter_sd Standard error of 'theta_hat_inter'
#' @param n_control_inter Sample size of the control arm at the interim analysis
#' @param beta Type II error rate, i.e. 1 - power
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Interim analysis hypothesis test}{Character string giving the outcome of the interim
#'     hypothesis test: whether the null hypothesis was rejected at the interim analysis}
#'   \item{Wald statistics}{Wald statistic of the treatment effect observed at the interim analysis}
#'   \item{Critical boundary}{Planned critical boundary at the interim analysis at which the sample
#'     size is re-estimated}
#'   \item{Estimated effect size: theta}{The effect size used for the conditional power calculation,
#'     taken to be the treatment effect estimate observed at the interim analysis}
#'   \item{c_alpha}{Conditional type I error at the interim analysis, calculated by setting theta = 0}
#'   \item{c_power_theta}{Conditional power at the interim analysis using the estimated effect size
#'     theta}
#'   \item{N_control_new}{New sample size for the control arm after the sample size re-estimation,
#'     or 'NA' if the null hypothesis was rejected at the interim analysis or if the estimated effect
#'     size is < 0.001}
#'   \item{new information fraction}{The new information fractions after the sample size change}
#'   \item{New critical boundary}{New critical boundaries for the new interim analyses
#'     (O'Brien-Fleming type), or 'NA' if they were not calculated}
#' }
#'
#' @examples
#' asd_new_design(s = c(1/3, 2/3, 1), c_bry = c(3.984298, 3.961146, 2.349682), inter_vt = 1,
#' snew = c(0.5, 1), theta_hat_inter = 0.3, theta_hat_inter_sd = 0.5,
#' n_control_inter = 100, beta = 0.1)
#'
#' @export
asd_new_design<-function(s,c_bry,inter_vt,snew,theta_hat_inter,theta_hat_inter_sd,n_control_inter,beta){
  #this function calculates the new sample size to achieve desired conditional power of 1-beta#
  #theta is the effect size
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #c_alpha  is the target one-sided type I error#
  #theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #n_control_inter is the sample size of the control arm at the interim analysis
  ##where sample size is being re-estimated statistics###
  theta=theta_hat_inter
  s_inter=s[inter_vt]
  T_inter=1/theta_hat_inter_sd^2
  T=T_inter/s[inter_vt]
  wald_inter=theta_hat_inter/theta_hat_inter_sd
  test= wald_inter-c_bry[inter_vt]

  c_power_theta=asd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)

  c_alpha=asd_c_power_theta(0,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)
  power<-function(sample_r){
    #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    new_c_bry=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    #asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    tmp=asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    return(tmp)
  }
  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    N_control_new="NA"
    c_bry_f="NA"

  }else {
    outcome="Null hypothesis not rejected at interim analysis"

    if (theta<10e-4){
      N_control_new="NA"
      c_bry_f="NA"
    }else{
      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.0001,2^32),tol=1e-5)[1])
      #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
      c_bry_f=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,new_ratio)
      T_new=T*new_ratio
      N_control_new=ceiling(n_control_inter*(T_new/T_inter))
      #  N_new=ceiling(T_new*2*sigma^2)
    }
  }
  out=list(outcome,wald_inter,c_bry[inter_vt],theta,c_alpha, c_power_theta,N_control_new,snew,c_bry_f)
  names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Estimated effect size: theta","c_alpha","c_power_theta","N_control_new","new information fraction","New critical boundary")

  return(out)

}
#################################
#asd_new_design(s,inter_vt,snew,theta_hat_inter,theta_hat_inter_sd,n_control_inter,beta)
#asd_new_design(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),inter_vt=1,snew=1,theta_hat_inter=1.2,theta_hat_inter_sd=1.5,n_control_inter=100,beta=0.1)
#asd_new_design(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=1.2,theta_hat_inter_sd=0.2,n_control_inter=100,beta=0.1)
#asd_new_design(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=-1.2,theta_hat_inter_sd=0.2,n_control_inter=100,beta=0.1)

####################################
####################################
#' Group sequential designs / Adaptive sequential designs: Interim analysis, Survival analysis, Input theta
#'
#' @param theta Assumed effect size used for the conditional power calculation, i.e. minus the log
#'   of the hazard ratio; different values of theta can be entered to check the sensitivity of the
#'   conditional power to the effect size; if theta < 0.001 the new number of events and the new
#'   critical boundary are returned as 'NA'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based);
#'   must be less than length(s)
#' @param snew Vector of new information fractions after the sample size is changed, e.g. c(0.5, 1);
#'   the fractions are relative to the new total information and the last element must be 1
#' @param theta_hat_inter Treatment effect estimate observed at the interim analysis, i.e. minus the
#'   log of the observed hazard ratio
#' @param theta_hat_inter_sd Standard error of 'theta_hat_inter'
#' @param events_inter Number of events observed at the interim analysis
#' @param beta Type II error rate, i.e. 1 - power
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Interim analysis hypothesis test}{Character string giving the outcome of the interim
#'     hypothesis test: whether the null hypothesis was rejected at the interim analysis}
#'   \item{Wald statistics}{Wald statistic of the treatment effect observed at the interim analysis}
#'   \item{Critical boundary}{Planned critical boundary at the interim analysis at which the sample
#'     size is re-estimated}
#'   \item{Assumed effect size: theta}{The assumed effect size theta used for the conditional power
#'     calculation}
#'   \item{c_alpha}{Conditional type I error at the interim analysis, calculated by setting theta = 0}
#'   \item{c_power_theta}{Conditional power at the interim analysis using the assumed effect size theta}
#'   \item{events_new}{New total number of events required after the sample size re-estimation,
#'     or 'NA' if the null hypothesis was rejected at the interim analysis or if theta < 0.001}
#'   \item{New information fraction}{The new information fractions after the sample size change}
#'   \item{New critical boundary}{New critical boundaries for the new interim analyses
#'     (O'Brien-Fleming type), or 'NA' if they were not calculated}
#' }
#'
#' @examples
#' asd_new_design_theta_survival(theta = 0.2, s = c(1/3, 2/3, 1),
#' c_bry = c(3.984298, 3.961146, 2.349682), inter_vt = 1, snew = c(0.5, 1), theta_hat_inter = 0.3,
#' theta_hat_inter_sd = 0.5, events_inter = 100, beta = 0.1)
#'
#' @export
asd_new_design_theta_survival<-function(theta,s,c_bry,inter_vt,snew,theta_hat_inter,theta_hat_inter_sd,events_inter,beta){
  #this function calculates the new sample size to achieve desired conditional power of 1-beta#
  #theta is the effect size
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #c_alpha  is the target one-sided type I error#
  #theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #n_control_inter is the sample size of the control arm at the interim analysis
  ##where sample size is being re-estimated statistics###
  s_inter=s[inter_vt]
  T_inter=1/theta_hat_inter_sd^2
  T=T_inter/s[inter_vt]
  wald_inter=theta_hat_inter/theta_hat_inter_sd
  test= wald_inter-c_bry[inter_vt]

  c_power_theta=asd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)

  c_alpha=asd_c_power_theta(0,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)
  power<-function(sample_r){
    #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    new_c_bry=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    #asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    tmp=asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    return(tmp)
  }
  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    events_new="NA"
    c_bry_f="NA"

  }else {
    outcome="Null hypothesis not rejected at interim analysis"

    if (theta<10e-4){
      events_new="NA"
      c_bry_f="NA"
    }else{
      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.0001,2^32),tol=1e-5)[1])
      #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
      c_bry_f=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,new_ratio)
      T_new=T*new_ratio
      events_new=ceiling(events_inter*(T_new/T_inter))
    }
  }
  out=list(outcome,wald_inter,c_bry[inter_vt],theta,c_alpha, c_power_theta, events_new,snew,c_bry_f)
  names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Assumed effect size: theta","c_alpha","c_power_theta","events_new","New information fraction","New critical boundary")

  return(out)

}
#################################
#asd_new_design_theta_survival(theta,s,inter_vt,snew,theta_hat_inter,theta_hat_inter_sd,events_inter,beta)
#asd_new_design_theta_survival(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=1,theta_hat_inter=1.2,theta_hat_inter_sd=1.5,events_inter=100,beta=0.1)
#asd_new_design_theta_survival(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=1.2,theta_hat_inter_sd=0.2,events_inter=100,beta=0.1)
#asd_new_design_theta_survival(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=-1.2,theta_hat_inter_sd=0.2,events_inter=100,beta=0.1)

#####################################
####################################
#' Group sequential designs / Adaptive sequential designs: Interim analysis, Survival analysis, Use theta_interim
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param inter_vt Index of the interim analysis at which the sample size is re-estimated (1-based);
#'   must be less than length(s)
#' @param snew Vector of new information fractions after the sample size is changed, e.g. c(0.5, 1);
#'   the fractions are relative to the new total information and the last element must be 1
#' @param theta_hat_inter Treatment effect estimate observed at the interim analysis, i.e. minus the
#'   log of the observed hazard ratio
#' @param theta_hat_inter_sd Standard error of 'theta_hat_inter'
#' @param events_inter Number of events observed at the interim analysis
#' @param beta Type II error rate, i.e. 1 - power
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Interim analysis hypothesis test}{Character string giving the outcome of the interim
#'     hypothesis test: whether the null hypothesis was rejected at the interim analysis}
#'   \item{Wald statistics}{Wald statistic of the treatment effect observed at the interim analysis}
#'   \item{Critical boundary}{Planned critical boundary at the interim analysis at which the sample
#'     size is re-estimated}
#'   \item{Estimated effect size: theta}{The effect size used for the conditional power calculation,
#'     taken to be the treatment effect estimate observed at the interim analysis}
#'   \item{c_alpha}{Conditional type I error at the interim analysis, calculated by setting theta = 0}
#'   \item{c_power_theta}{Conditional power at the interim analysis using the estimated effect size
#'     theta}
#'   \item{events_new}{New total number of events required after the sample size re-estimation,
#'     or 'NA' if the null hypothesis was rejected at the interim analysis or if the estimated effect
#'     size is < 0.001}
#'   \item{New information fraction}{The new information fractions after the sample size change}
#'   \item{New critical boundary}{New critical boundaries for the new interim analyses
#'     (O'Brien-Fleming type), or 'NA' if they were not calculated}
#' }
#'
#' @examples
#' asd_new_design_survival(s = c(1/3, 2/3, 1), c_bry = c(3.984298, 3.961146, 2.349682),
#' inter_vt = 1, snew = c(0.5, 1), theta_hat_inter = 0.3, theta_hat_inter_sd = 0.5,
#' events_inter = 100, beta = 0.1)
#'
#' @export
asd_new_design_survival<-function(s,c_bry,inter_vt,snew,theta_hat_inter,theta_hat_inter_sd,events_inter,beta){
  #this function calculates the new sample size to achieve desired conditional power of 1-beta#
  #theta is the effect size
  #s is the vector of new information fractions. The length of s can be either 2 or 3#
  #c_alpha  is the target one-sided type I error#
  #theta_hat_inter is the vector ob obsered theta_hat at the interim analysis ##
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #n_control_inter is the sample size of the control arm at the interim analysis
  ##where sample size is being re-estimated statistics###
  theta=theta_hat_inter
  s_inter=s[inter_vt]
  T_inter=1/theta_hat_inter_sd^2
  T=T_inter/s[inter_vt]
  wald_inter=theta_hat_inter/theta_hat_inter_sd
  test= wald_inter-c_bry[inter_vt]

  c_power_theta=asd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)

  c_alpha=asd_c_power_theta(0,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)
  power<-function(sample_r){
    #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    new_c_bry=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
    #asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    tmp=asd_c_power_new(theta,s,inter_vt,snew,new_c_bry,theta_hat_inter,theta_hat_inter_sd,sample_r)
    return(tmp)
  }
  if (test>=0){
    outcome="Null hypothesis rejected at interim analysis"
    events_new="NA"
    c_bry_f="NA"

  }else {
    outcome="Null hypothesis not rejected at interim analysis"

    if (theta<10e-4){
      events_new="NA"
      c_bry_f="NA"
    }else{
      new_ratio=as.numeric(stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.0001,2^32),tol=1e-5)[1])
      #new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,sample_r)
      c_bry_f=new_OF_bundry(s,inter_vt,snew,c_alpha,theta_hat_inter,theta_hat_inter_sd,new_ratio)
      T_new=T*new_ratio
      events_new=ceiling(events_inter*(T_new/T_inter))
    }
  }
  out=list(outcome,wald_inter,c_bry[inter_vt],theta,c_alpha, c_power_theta, events_new,snew,c_bry_f)
  names(out)=c("Interim analysis hypothesis test","Wald statistics","Critical boundary","Estimated effect size: theta","c_alpha","c_power_theta","events_new","New information fraction","New critical boundary")

  return(out)

}
#################################
#asd_new_design_survival(s,inter_vt,snew,theta_hat_inter,theta_hat_inter_sd,events_inter,beta)
#asd_new_design_survival(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=1,theta_hat_inter=1.2,theta_hat_inter_sd=1.5,events_inter=100,beta=0.1)
#asd_new_design_survival(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=1.2,theta_hat_inter_sd=0.2,events_inter=100,beta=0.1)
#asd_new_design_survival(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682), inter_vt=1,snew=c(0.5,1),theta_hat_inter=-1.2,theta_hat_inter_sd=0.2,events_inter=100,beta=0.1)
#############################################################################
asd_cprob_new<-function(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad){
  #this function calculates conditional probability from the data after sample szie change
  #snew is the vector of information fractions after sample size change #
  #c_bry_new is the criotical boundaries after sample size change##
  #inter_vt is the interim analysis at which the sample size was changed
  #theta_hat_inter is the theta_hat at inter_vt'th interim analysis
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #theta is the  effect size####
  # theta_hat_last_ad is the last observed thet_hat afer sample size change###
  # theta_hat_last_ad_sd is the standard error of theta_hat_last_ad ###
  t=T*s
  T_inter=theta_hat_inter_sd^(-2)
  T=T_inter/s[inter_vt]
  T_last_ad=theta_hat_last_ad_sd^(-2)
  T_new=(T_last_ad-T_inter)/snew[last_vt_ad]+T_inter
  s_new=((T_new-T_inter)*snew+T_inter)/T_new
  t_new=(T_new-T_inter)*snew+T_inter
  score_last_ad=theta_hat_last_ad*theta_hat_last_ad_sd^(-2)
  score_inter=theta_hat_inter*theta_hat_inter_sd^(-2)
  T_inter=theta_hat_inter_sd^(-2)
  t_new=(T_new-T_inter)*snew+T_inter

  new_e_bry=c_bry_new*sqrt(snew)
  if (last_vt_ad==1){
    new_e_bry_f=(score_last_ad-score_inter)/sqrt(T_new-T_inter)-theta*sqrt(T_new-T_inter)*snew[last_vt_ad]
  }else{
    new_e_bry_tail=(score_last_ad-score_inter)/sqrt(T_new-T_inter)-theta*sqrt(T_new-T_inter)*snew[last_vt_ad]

    new_e_bry_base=new_e_bry[1:(last_vt_ad-1)]-theta*sqrt(T_new-T_inter)*snew[1:(last_vt_ad-1)]
    new_e_bry_f=c(new_e_bry_base,new_e_bry_tail)
  }
  #################################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(snew,length(snew)),nrow=length(snew), ncol=length(snew))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  ########################################
  ## the function for the vector of Brownian motions ###
  cal1<-function(theta,sigma){

    tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=new_e_bry_f,mean=rep(0,last_vt_ad),sigma=sigma,algorithm=mvtnorm::Miwa())[1]
    return(tmp)
  }
  ### the function for single Brownian motion###
  if (T_new<=T_inter){c_power_new=0} else {
    if (last_vt_ad==1){
      c_power_new=1-stats::pnorm(new_e_bry_f)
    }else{
      c_power_new= cal1(theta,cov_mtx[1:last_vt_ad,1:last_vt_ad])
    }
  }
  return(c_power_new)
}
##################################
#asd_cprob_new<-function(theta,s,snew,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd,theta_hat_last_ad,theta_hat_last_ad_sd,last_vt_ad)

#asd_cprob_new(theta=0.2,s=c(1/3,2/3,1),snew=1,c_bry_new=2.545644,inter_vt=2,theta_hat_inter=1.2,theta_hat_inter_sd=1,theta_hat_last_ad=3,theta_hat_last_ad_sd=0.8,last_vt_ad=1)
#asd_cprob_new(theta=0.1,s=c(1/3,2/3,1),snew=c(0.5,1),c_bry_new=c(2.545644,2),inter_vt=2,theta_hat_inter=1.2,theta_hat_inter_sd=1,theta_hat_last_ad=3,theta_hat_last_ad_sd=0.8,last_vt_ad=2)
#asd_cprob_new(theta=0.3,s=c(1/3,2/3,1),snew=1,c_bry_new=2.545644,inter_vt=2,theta_hat_inter=1.5,theta_hat_inter_sd=1,theta_hat_last_ad=1.4,theta_hat_last_ad_sd=0.8,last_vt_ad=1)


####################################
###########################Calculating backward image#######
asd_bk_img_gen<-function(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p){
  #this function calculates the backward image after sample size change , given aonditional probability c_p
  #theta is any effect size
  #s is the infromation fractions for the original design
  #c_bry is the critical bpundaries with the original design
  #inter_vt is the interim analysis at which the sample size was changed
  #theta_hat_inter is the theta_hat at inter_vt'th interim analysis
  #theta_hat_inter_sd is the standard error of theta_hat_inter
  #c_p is the conditional probability from the data after sample size change
  score_inter=theta_hat_inter*theta_hat_inter_sd^(-2)
  prob=NULL
  for (q in (inter_vt+1):length(s)){
    #  prob_q=asd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd){

    prob_q=asd_c_power_theta(theta,s[1:q],c_bry[1:q],inter_vt,theta_hat_inter,theta_hat_inter_sd)
    prob=c(prob,prob_q)
  }
  if (min(prob-c_p)<0){
    p_diff=prob-c_p
    neg=p_diff[p_diff<0]
    cut=match(dplyr::last(neg),p_diff)
  }else{ cut=0}
  bk_prob_1<-function(c_lst){
    c_bry_new=c(c_bry[1:(cut+inter_vt)],c_lst)
    tmp=asd_c_power_theta(theta,s[1:(cut+inter_vt+1)],c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd)

  }
  bk_prob_2<-function(c_lst){
    c_bry_new=c(c_bry[1:(length(s)-1)],c_lst)
    tmp=prob_q=asd_c_power_theta(theta,s,c_bry_new,inter_vt,theta_hat_inter,theta_hat_inter_sd)
  }
  if ((cut+inter_vt)<length(s)){
    bk_img=as.numeric(stats::uniroot(function(x) bk_prob_1(x)-c_p,c(-2^32,2^32),tol=1e-5)[1])
    bry_f=c(c_bry[1:(cut+inter_vt)],bk_img)
  } else {
    bk_img=as.numeric(stats::uniroot(function(x) bk_prob_2(x)-c_p,c(-2^32,2^32),tol=1e-5)[1])
    bry_f=c(c_bry[1:(length(s)-1)],bk_img)
  }

  out=list(cut,c_bry,bry_f,c(prob,c_p))
  names(out)=c("cut","c_bry","bry_f","prob")
  return(bry_f)
  #  return(out)
}

########################################
#asd_bk_img_gen(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,inter_vt,c_p)

#asd_bk_img_gen(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.888374, 3.864257, 2.212819),theta_hat_inter=1.2,theta_hat_inter_sd=1,inter_vt=1,c_p=0.4)
#asd_bk_img_gen(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.888374, 3.864257, 2.212819),theta_hat_inter=1.2,theta_hat_inter_sd=1,inter_vt=2,c_p=0.4)
###########################Calculating backward image#######

###################################################
##########   Part I                    ############
### Component functions, to be sourced ############
###################################################

##########Final analysi, without sample size change ###########

#' Group sequential designs / Adaptive sequential designs: Final analysis, Without sample size change
#'
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param c_bry Vector of critical boundaries for each analysis
#' @param theta_hat_last Treatment effect estimate observed at the final analysis, i.e. the naive
#'   point estimate of the parameter of interest at the analysis at which the trial terminated
#' @param theta_hat_last_sd Standard error of 'theta_hat_last'
#' @param last_vt Index of the analysis at which the trial stopped (1-based); this can be the final
#'   analysis, at which s = 1, or an interim analysis, and it must not exceed length(s)
#' @param alpha_2 Type I error for the two-sided test
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{p-value}{The p-value of the test of the null hypothesis that the treatment effect is 0}
#'   \item{confidence interval}{Numeric vector of length 3 giving the lower confidence limit, the
#'     median-unbiased point estimate and the upper confidence limit of the two-sided
#'     (1 - alpha_2) confidence interval}
#' }
#'
#' @examples
#' gsd_ci_est(s = c(1/3,2/3,1), c_bry = c(3.471095, 2.454434,2.004037), theta_hat_last = 2.004,
#' theta_hat_last_sd = 1, last_vt = 3, alpha_2 = 0.05)
#'
#' @export
gsd_ci_est<-function(s,c_bry,theta_hat_last,theta_hat_last_sd,last_vt,alpha_2){
  #This function calculates the confidence interval and the point estimate
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #alpha_2  is for the (1-alpha_2)*100% confidence interval #
  #typically, alpha_2=0.05
  #c_bry is the critical nboundary
  # theta_hat_last is  the last observed theta_hat###
  #theta_hat_last_sd is hte standard error of theta_hat_last
  #last_vt numbers the interim analysis where the trial terminated (it can be final analysis)
  T_last=theta_hat_last_sd^(-2)
  T=T_last/s[last_vt]

  wald_last=theta_hat_last/theta_hat_last_sd

  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(s,length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################
  e_bry=sqrt(s)*c_bry

  #############################
  ### exit boundary ####

  if (last_vt>1){
    e_bry_1=e_bry[1:(last_vt-1)]
    e_bry_f=c(e_bry_1,wald_last*sqrt(s[last_vt]))
  }else if (last_vt==1){
    e_bry_f=wald_last
  }

  ### the function calculates the exit probability###
  cal<-function(theta,e_bry_f,sigma,last_vt){
    if(last_vt==1){
      up_bry=e_bry_f-theta*sqrt(T*s[last_vt])
      tmp=1-stats::pnorm(up_bry)
    }else{
      sigma1=sigma[1:last_vt,1:last_vt]
      up_bry=e_bry_f-theta*sqrt(T)*s[1:last_vt]
      sigma1=sigma[1:last_vt,1:last_vt]
      tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,last_vt),sigma=sigma1,algorithm=mvtnorm::Miwa())[1]
    }
    return(tmp)
  }
  theta_ll=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-1+alpha_2/2,c(-2^32,2^32),tol=1e-5)[1])
  p_value= cal(0,e_bry_f,cov_mtx,last_vt)
  out=list(p_value,c(theta_ll,theta_pt,theta_ul))
  names(out)=c("p-value","confidence interval")
  return(out)
  #theta_pt is point estimate
  #(theta_ll,theta_ul) is the confidence interval
}
########################
#gsd_ci_est<-function(s,c_bry,theta_hat_last,theta_hat_last_sd,last_vt,alpha)
#gsd_ci_est(s=c(1/3,2/3,1),c_bry=c(3.471095, 2.454434,2.004037),theta_hat_last=2.004,theta_hat_last_sd=1,last_vt=3,alpha_2=0.05)
#gsd_ci_est(s=c(0.3,0.7,1),c_bry=c(2.293075, 2.293075, 2.293075),theta_hat_last=2.294,theta_hat_last_sd=1,last_vt=1,alpha_2=0.05)
##########################################

###########################Calculating backward image#######
bk_img_com_theta<-function(theta,s,theta_hat_inter,theta_hat_inter_sd,theta_hat_last,theta_hat_last_sd){
  #This function calculates the backward image for given theta
  #theta=0.2
  #s=c(1/3,2/3,1)
  #c_bry=c(3.888374, 3.864257, 2.212819)

  #s is the oiginal infromation fraction
  #c_bry is the original critical boundary
  #inter_vt is the order of the interim analysis
  #n_inter_c is the sample size of the control arm at the interim analysis
  # inter_theta_hat is the theta_hat observed at the interim analysis ###
  #inter_theta_hat_sd is the estimated variance at the interim analysis
  #rand_ratio is randomization ratio, sample size ratio of n_t/n_c
  inter_vt=length(s)-1

  T_inter=1/theta_hat_inter_sd^2
  score_inter=theta_hat_inter*T_inter
  wald_inter=theta_hat_inter*sqrt(T_inter)
  T_old=T_inter/s[inter_vt]


  T_last=1/theta_hat_last_sd^2
  score_last=theta_hat_last*T_last
  wald_last=theta_hat_last*sqrt(T_last)

  ##########backward image, given theta############

  bk_img_com_theta=(sqrt(T_old-T_inter)/sqrt(T_last-T_inter))*((score_last-score_inter)-theta*(T_last-T_inter))+score_inter+theta*(T_old-T_inter)

  return(bk_img_com_theta)
}
########################################
#bk_img_com_theta(theta,s,theta_hat_inter,theta_hat_inter_sd,theta_hat_last,theta_hat_last_sd)
#bk_img_com_theta(theta=0.2,s=c(0.5,1),theta_hat_inter=0.1,theta_hat_inter_sd=0.1,theta_hat_last=0.11,theta_hat_last_sd=0.05)

#########################################################################
c_power_1<-function(s,c_bry,inter_vt,inter_theta_hat,inter_theta_hat_sd){
  #this function calculates the conditional power after observing inter_theta_hat,inter_theta_hat_sd at inter_vt#
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  # inter_theta_hat is the theta_hat observed at the interim analysis ###
  #inter_theta_hat_sd is the estimated variance at the interim analysis
  #### the next section constructs the covariance matrix ###
  e_bry=c_bry*sqrt(s)
  #inter_vt=length(s)-1
  theta=inter_theta_hat


  T_inter=1/inter_theta_hat_sd^2
  inter_score=inter_theta_hat*T_inter
  T_total=T_inter/s[inter_vt]
  ###########################

  up_bry=(e_bry[length(s)]-inter_score)/sqrt(T_total-T_inter)-theta*sqrt(T_total-T_inter)
  c_power= 1-stats::pnorm(up_bry)

  return(c_power)
}
###################################################
#c_power_1(s,c_bry,inter_vt,inter_theta_hat,inter_theta_hat_sd)
#c_power_1(c(1/3,2/3,1),c(3.471095,2.454434,2.004037),2,0.05,0.1)
#c_power_1(c(1/3,2/3,1),c(3.471095,2.454434,2.004037),2,0.03,0.1)


########################################################################
cond_power_com<-function(s,c_bry_two_stg,theta_hat_inter,theta_hat_inter_sd){
  #this function calculates the conditional power after observing inter_theta_hat,inter_theta_hat_sd at inter_vt#
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  # inter_theta_hat is the theta_hat observed at the interim analysis ###
  #inter_theta_hat_sd is the estimated variance at the interim analysis
  #### the next section constructs the covariance matrix ###

  #inter_vt=length(s)-1
  theta=theta_hat_inter

  inter_vt=length(s)-1
  T_inter=1/theta_hat_inter_sd^2
  score_inter=theta_hat_inter*T_inter
  T_total=T_inter/s[inter_vt]
  t=T_total*s
  e_bry_two_stg=c_bry_two_stg*sqrt(t)
  ###########################

  up_bry=(e_bry_two_stg[length(s)]-score_inter)/sqrt(T_total-T_inter)-theta*sqrt(T_total-T_inter)
  c_power= 1-stats::pnorm(up_bry)


  return(c_power)
}
########################################
#comd_power_com(s,c_bry_two_stg,inter_theta_hat_inter,theta_hat_inter_sd)

#cond_power_com(s=c(1/3,2/3,1),c_bry_two_stg=c(3.471095,2.454434,2.004037),theta_hat_inter=0.3,theta_hat_inter_sd=0.03)


#################################################################
cond_power_com_theta<-function(theta,s,c_bry_two_stg,theta_hat_inter,theta_hat_inter_sd){
  #this function calculates the conditional power after observing inter_theta_hat,inter_theta_hat_sd at inter_vt#
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #theta_hat_inter is the theta_hat observed at the interim analysis ###
  #theta_hat_inter_sd is the estimated variance at the interim analysis
  #### the next section constructs the covariance matrix ###
  #Theta is the given effect size
  #inter_vt=length(s)-1

  inter_vt=length(s)-1
  T_inter=1/theta_hat_inter_sd^2
  score_inter=theta_hat_inter*T_inter
  T_total=T_inter/s[inter_vt]
  t=T_total*s
  e_bry_two_stg=c_bry_two_stg*sqrt(t)
  ###########################

  up_bry=(e_bry_two_stg[length(s)]-score_inter)/sqrt(T_total-T_inter)-theta*sqrt(T_total-T_inter)
  c_power= 1-stats::pnorm(up_bry)

  return(c_power)
}
########################################
#cond_power_com_theta(theta,s,c_bry_two_stg,inter_vt,theta_hat_inter,theta_hat_inter_sd){

#cond_power_com_theta(theta=0.2,s=c(1/3,2/3,1),c_bry_two_stg=c(3.471095,2.454434,2.004037),theta_hat_inter=0.3,theta_hat_inter_sd=0.03)

######################################
new_c_bry_com<-function(s,c_bry,theta_hat_inter,theta_hat_inter_sd,T_new){
  #this function selects new  critical boundaries for sample size change#
  #s is the oiginal infromation fraction
  #c_bry is the original critical boundary
  #inter_vt is the order of the interim analysis
  #theta_hat_inter is the theta_hat observed at the interim analysis ###
  #theta_hat_inbter_sd is the estimated variance at the interim analysis
  #T_new is the new information time after sampke size change

  inter_vt=length(s)-1

  theta_hat=theta_hat_inter
  T_inter=1/theta_hat_inter_sd^2
  score_inter=theta_hat_inter*T_inter
  wald_inter=theta_hat_inter*sqrt(T_inter)
  T_old=T_inter/s[inter_vt]
  t=T_old*s
  e_bry=c_bry*sqrt(t)

  ###################################
  new_e_bry=(sqrt(T_new-T_inter)/sqrt(T_old-T_inter))*(dplyr::last(e_bry)-score_inter)+score_inter
  new_c_bry=new_e_bry/sqrt(T_new)

  ############################################################

  return(new_c_bry)
}
##############################
#new_c_bry_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,T_new)
#new_c_bry_com(s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),theta_hat_inter=0.05,theta_hat_inter_sd=0.03,T_new=3408)



######################################
new_info_time_com<-function(s,c_bry,theta_hat_inter,theta_hat_inter_sd,beta){
  #this function selects new information time  for sample size change#
  #s is the oiginal infromation fraction
  #c_bry is the original critical boundary
  #inter_vt is the order of the interim analysis
  #theta_hat_inter is the theta_hat observed at the interim analysis ###
  #theta_hat_inbter_sd is the estimated variance at the interim analysis

  inter_vt=length(s)-1

  theta_hat=theta_hat_inter
  T_inter=1/theta_hat_inter_sd^2
  score_inter=theta_hat_inter*T_inter
  wald_inter=theta_hat_inter*sqrt(T_inter)
  T_old=T_inter/s[inter_vt]
  t=T_old*s
  e_bry=c_bry*sqrt(t)

  ###################################
  T_new=theta_hat^(-2)*((dplyr::last(e_bry)-score_inter)/(sqrt(T_old-T_inter))+stats::qnorm(1-beta))^2+T_inter

  ############################################################
  return(T_new)
}
##############################
#new_info_time_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,beta)
#new_info_time_com(s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),theta_hat_inter=0.05,theta_hat_inter_sd=0.03,beta=0.1)
##################################################################


##################################################################
new_info_time_com_theta<-function(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,beta){
  #this function selects new information time,for sample size change#
  #s is the oiginal infromation fraction
  #c_bry is the original critical boundary
  #inter_vt is the order of the interim analysis
  #theta_hat_inter is the theta_hat observed at the interim analysis ###
  #theta_hat_inbter_sd is the estimated variance at the interim analysis
  inter_vt=length(s)-1

  T_inter=1/theta_hat_inter_sd^2
  score_inter=theta_hat_inter*T_inter
  wald_inter=theta_hat_inter*sqrt(T_inter)
  T_old=T_inter/s[inter_vt]
  t=T_old*s
  e_bry=c_bry*sqrt(t)

  ###################################
  T_new=theta^(-2)*((dplyr::last(e_bry)-score_inter)/(sqrt(T_old-T_inter))+stats::qnorm(1-beta))^2+T_inter

  ############################################################
  return(T_new)
}
##############################
#new_info_time_com_theta(theta,s,c_bry,theta_hat_inter,theta_hat_inter_sd,beta){
#new_info_time_com_theta(theta=0.2,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),theta_hat_inter=0.05,theta_hat_inter_sd=0.03,beta=0.1)
######################################

