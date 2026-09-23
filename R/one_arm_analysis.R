#rm(list = ls())

#library(mvtnorm)

#########Interim analysis , sample size re-estimation ######

#' The interim analysis can be conducted for both two-stage and three-stage designs
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param r_final For two-stage designs, r_final=r2. For the three-stage design, r_final=r3.
#' @param n_inter Number of patients at the interim look. For two-stage designs, n_inter=n1. For the three-stage design, n_inter=n2.
#' @param n_final For two-stage designs, n_inter=n2. For the three-stage design, n_inter=n3.
#' @param x_inter Number of responses observed at the interim look.
#' @param beta    The conditional type II error. 1-beta is the desired conditional power.
#' @param N_max   Pre-selected maximum sample size.
#'
#' @return n_new: is the new sample size.
#' @return r_new (without continuity correction) and r_new (with continuity correction) are thresholds such that the null hypothesis will be
#' rejected at n_new if at least r_new responses are observed. Simulations on type I error will help to determine if continuity correction is needed.
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references P. Gao, J.H. Ware, C. Mehta, (2008), Sample size re-estimation for adaptive sequential designs.
#' Journal of Biopharmaceutical Statistics, 18: 1184–1196, 2008.
#'
#' @examples
#' interim_analysis(p_0=0.2,r_final=17,n_inter=22,n_final=56,x_inter=7,beta=0.1,N_max=120)
#' @export
interim_analysis<-function(p_0,r_final,n_inter,n_final,x_inter,beta,N_max){
  b=stats::qnorm(1-beta)
  ### convert r_2 to corresponding c_2 wihtout continuity correction ####
  c_f=(r_final/n_final-p_0)*sqrt(n_final^3/(r_final*(n_final-r_final)))

  ### convert r_2 to corresponding c_2 wiht continity correction ####
  r_f_0=r_final+0.5   #continuity correction
  c_f_0=(r_f_0/n_final-p_0)*sqrt(n_final^3/(r_f_0*(n_final-r_f_0)))  #

  ###### Interim analysis #####

  p_hat_inter=x_inter/n_inter
  var_hat_inter=p_hat_inter*(1-p_hat_inter)
  t=c(n_inter,n_final)/var_hat_inter
  t_inter=t[1]
  t_f=t[2]
  score_inter=p_hat_inter*t_inter
  ## new information time with continity correction ####
  T_new_0=(p_hat_inter-p_0)^(-2)*(b+(c_f_0*sqrt(t_f)-score_inter)/sqrt(t_f-t_inter))+t_inter
  ### new information time with continuity correction ###
  n_new_0=ceiling(T_new_0*p_hat_inter*(1-p_hat_inter))
  # T_new=n_new/2/sigma^2
  n_new=min(N_max,n_new_0)
  T_new=n_new/(p_hat_inter*(1-p_hat_inter))
  #####  new critical boundary without continuoty correction ####
  e_bry_new=(sqrt(T_new-t_inter)/sqrt(t_f-t_inter))*(c_f*sqrt(t_f)-score_inter)+score_inter

  c_bry_new=e_bry_new/sqrt(T_new)

  ############## convert c_bry_new to r_2_new ######################

  wald_f_new_sq=NULL  # Wald statistics
  for (i in 1:n_new){
    wald_f_new_sq[i]=(i/n_new-p_0)*sqrt(n_new^3/(i*(n_new-i)))
  }  # sq=sequence ###
  r_f_new_sq=cbind(1:n_new,wald_f_new_sq)

  ## new adjusted r_f_new without continuity correction  ####

  r_f_new_cut_sq=r_f_new_sq[r_f_new_sq[,2]>=c_bry_new,]
  r_f_new=r_f_new_cut_sq[1,1]  ### The number of response that corresponds to c_bry_new. This will be used as the rejection criterion ###

  ###################### new critical boundary with continuity correction ####
  e_bry_new_0=(sqrt(T_new-t[1])/sqrt(t[2]-t[1]))*(c_f_0*sqrt(t[2])-score_inter)+score_inter

  c_bry_new_0=e_bry_new_0/sqrt(T_new)

  ##############################################
  ## new adjusted r_f_new with continuity correction  ####

  r_f_new_cut_sq_0=r_f_new_sq[r_f_new_sq[,2]>=c_bry_new_0,]
  r_f_new_0=r_f_new_cut_sq_0[1,1]### The number of response that corresponds to c_bry_new_0. This will be used as the rejection criterion ###

  ###########################################################

  out=list(n_new,r_f_new,r_f_new_0)
  names(out)=c("new sample size","new discrete rejection boundary without continuity correction","new discrete rejection boundary with continuity correction")
  return(out)
}
###########################
#interim_analysis(p_0,r_final,n_inter,n_final,x_inter,beta,N_max)
#interim_analysis(p_0=0.2,r_final=17,n_inter=22,n_final=56,x_inter=7,beta=0.1,N_max=120)

#########Final analysis , with/without sample size change ######

#' Final analysis: two stage design, no sample size change
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param n_1   Number of patients at first look.
#' @param n_2   Number of patients at second look.
#' @param r_1_e   The trial would be stopped for superiority if at least r2e responses are observed
#' at the second look.
#' @param r_2    If no more than r_2 responses are observed at the second look, then the null hypothesis will not be rejected.
#' @param x_last Number of responses at the final visit.
#' @param last_vt 1.
#' @param alpha    One-sided type I error.
#'
#' @return UL: upper limit of confidence interval.
#' @return LL: lower limit of confidence interval.
#' @return The results with continuity correction are not necessarily different from those with continuity correction.
#' @return The need for continuity correction can be determined through simulations on type I error (see Gao, Zhang, 2004).
#'
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs.
#' Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' gsd_ci_est_2_stage(p_0=0.2,n_1=100,n_2=105,r_1_e=NA,r_2=29,x_last=31,last_vt=1,alpha=0.025)
#'
#' @export
gsd_ci_est_2_stage<-function(p_0,n_1,n_2,r_1_e,r_2,x_last,last_vt,alpha){
  #This function calculates the confidence interval and the point estimate

  #typically, alpha=0.025
  N=c(n_1,n_2)
  s=N/n_2

  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################

  # statistics for final analysis

  p_hat_last=x_last/N[last_vt]

  var_last=p_hat_last*(1-p_hat_last)
  t=N/var_last
  T=t[last_vt]

  wald_last=(p_hat_last-p_0)*sqrt(t[last_vt])

  ### convert r_2 to corresponding c_2 wihtout continuity correction ####
  # d_1=(r_1/n_1-p_0)*sqrt(n_1^3/(r_1*(n_1-r_1)))
  if (is.na(r_1_e)){r_1_e=n_1-0.1}
  r=c(r_1_e,r_2)
  c_bry=(r/N-p_0)*sqrt(N^3/(r*(N-r)))
  e_bry=sqrt(s)*c_bry
  ### convert r_2 to corresponding c_2 wiht continity correction ####

  r_1_0=min(r_1_e+0.5,n_1-0.1)
  r_2_0=r_2+0.5   #continuity correction
  r_0=c(r_1_0,r_2_0)
  c_bry_0=(r_0/N-p_0)*sqrt(N^3/(r_0*(N-r_0)))

  e_bry_0=sqrt(s)*c_bry_0

  #############################
  ### exit boundary ####

  e_bry_1=e_bry[1:(last_vt-1)]
  e_bry_1_0=e_bry_0[1:(last_vt-1)]

  if (last_vt==2){
    e_bry_f=c(e_bry[1],wald_last*sqrt(s[2]))
    e_bry_f_0=c(e_bry_0[1],wald_last*sqrt(s[2]))
  }else if (last_vt==1){
    e_bry_f=wald_last
    e_bry_f_0=wald_last
  }

  ### the function calculates the exit probability###

  cal<-function(theta,e_bry_f,sigma,last_vt){
    if(last_vt==1){
      up_bry=e_bry_f-theta*sqrt(T*s[last_vt])
      tmp=1-stats::pnorm(up_bry)
    }else{
      sigma1=sigma[1:last_vt,1:last_vt]
      up_bry=e_bry_f-theta*sqrt(T)*s[1:last_vt]
      # ll_bry=c(d_1*sqrt(s[1]),-Inf)-theta*sqrt(T)*s[1:2]
      tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,last_vt),sigma=sigma1,abseps=1e-8,algorithm=mvtnorm::Miwa())[1]
    }
    return(tmp)
  }

  ##### CI without continuity correction ###
  theta_ll=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-alpha,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
  ##########################
  p_no_correction=cal(0,e_bry_f,cov_mtx,last_vt)
  ##### CI with continuity correction ###
  theta_ll_0=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f_0,cov_mtx,last_vt)-alpha,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt_0=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f_0,cov_mtx,last_vt)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul_0=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f_0,cov_mtx,last_vt)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
  #####################################
  p_correction=cal(0,e_bry_f_0,cov_mtx,last_vt)

  CI_withot_correction=c(theta_ll+p_0,theta_pt+p_0,theta_ul+p_0)
  CI_with_correction=c(theta_ll_0+p_0,theta_pt_0+p_0,theta_ul_0+p_0)
  out=list(p_no_correction,p_correction,CI_withot_correction,CI_with_correction)
  names(out)=c(" p_no_correction","p_correction","CI_without_correction","CI_with_correction")
  # out_1=list(out,e_bry_f_0)
  return(out)
}
########################
#gsd_ci_est_2_stage(p_0,n_1,n_2,r_1_e,r_2,x_last,last_vt,alpha)
#gsd_ci_est_2_stage(p_0=0.2,n_1=100,n_2=105,r_1_e=33,r_2=29,x_last=28,last_vt=2,alpha=0.025)
#gsd_ci_est_2_stage(p_0=0.2,n_1=100,n_2=105,r_1_e=NA,r_2=29,x_last=28,last_vt=2,alpha=0.025)
#gsd_ci_est_2_stage(p_0=0.2,n_1=100,n_2=105,r_1_e=31,r_2=29,x_last=31,last_vt=1,alpha=0.025)
#gsd_ci_est_2_stage(p_0=0.2,n_1=100,n_2=105,r_1_e=NA,r_2=29,x_last=31,last_vt=1,alpha=0.025)
######################################################

#' Final analysis: three stage design, no sample size change
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param n_2   Number of patients at second look.
#' @param n_3   Number of patients at third look.
#' @param r_2_e   The trial would be stopped for superiority if at least r_2_e responses are observed
#' at the second look. If r_2_e was not used in the design, enter NA.
#' @param r_3    If no more than r_3 responses are observed at the third look, then the null hypothesis will not be rejected.
#' @param x_last  Number of responses at the final visit.
#' @param last_vt 1.
#' @param alpha    One-sided type I error.
#'
#' @return UL: upper limit of confidence interval.
#' @return LL: lower limit of confidence interval.
#' @return The results with continuity correction are not necessarily different from those with continuity correction.
#' @return The need for continuity correction can be determined through simulations on type I error (see Gao, Zhang, 2004).
#'
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs.
#' Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' gsd_ci_est_3_stage(p_0=0.2,n_2=100,n_3=105,r_2_e=31,r_3=29,x_last=31,last_vt=3,alpha=0.025)
#'
#' @export
gsd_ci_est_3_stage<-function(p_0,n_2,n_3,r_2_e,r_3,x_last,last_vt,alpha){
  #This function calculates the confidence interval and the point estimate
  #typically, alpha=0.025

  N=c(n_2,n_3)
  s=N/n_3

  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################

  # statistics for final analysis

  p_hat_last=x_last/N[last_vt-1]

  var_last=p_hat_last*(1-p_hat_last)
  t=N/var_last
  T=t[last_vt-1]

  wald_last=(p_hat_last-p_0)*sqrt(t[last_vt-1])

  ### convert r_2 to corresponding c_2 wihtout continuity correction ####
  if (is.na(r_2_e)){r_2_e=n_2-0.1}
  r=c(r_2_e,r_3)
  c_bry=(r/N-p_0)*sqrt(N^3/(r*(N-r)))

  e_bry=sqrt(s)*c_bry
  ### convert r_2 to corresponding c_2 wiht continity correction ####
  r_2_0=min(r_2_e+0.5,n_2-0.1)
  r_3_0=r_3+0.5   #continuity correction
  r_0=c(r_2_0,r_3_0)
  c_bry_0=(r_0/N-p_0)*sqrt(N^3/(r_0*(N-r_0)))

  e_bry_0=sqrt(s)*c_bry_0

  e_bry_0=sqrt(s)*c_bry_0

  #############################
  ### exit boundary ####
  if (last_vt==3){
    e_bry_f=c(e_bry[1],wald_last*sqrt(s[2]))
    e_bry_f_0=c(e_bry_0[1],wald_last*sqrt(s[2]))
  }else if (last_vt==2){
    e_bry_f=wald_last
    e_bry_f_0=wald_last
  }

  ### the function calculates the exit probability###

  cal<-function(theta,e_bry_f,sigma,last_vt){
    if(last_vt==2){
      up_bry=e_bry_f-theta*sqrt(T)*s[1]
      tmp=1-stats::pnorm(up_bry)
    }else if (last_vt==3){
      sigma1=sigma[1:2,1:2]
      up_bry=e_bry_f-theta*sqrt(T)*s[1:2]
      #    ll_bry=c(d_1*sqrt(s[1]),-Inf)-theta*sqrt(T)*s[1:2]
      tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,last_vt-1),sigma=sigma1,abseps=1e-8,algorithm=mvtnorm::Miwa())[1]
    }
    return(tmp)
  }

  ##### CI without continuity correction ###
  theta_ll=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-alpha,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f,cov_mtx,last_vt)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
  ##########################
  p_no_correction=cal(0,e_bry_f,cov_mtx,last_vt)
  ##### CI with continuity correction ###
  theta_ll_0=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f_0,cov_mtx,last_vt)-alpha,c(-2^32,2^32),tol=1e-5)[1])
  theta_pt_0=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f_0,cov_mtx,last_vt)-0.5,c(-2^32,2^32),tol=1e-5)[1])
  theta_ul_0=as.numeric(stats::uniroot(function(x) cal(x,e_bry_f_0,cov_mtx,last_vt)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
  #####################################
  p_correction=cal(0,e_bry_f_0,cov_mtx,last_vt)

  CI_withot_correction=c(theta_ll+p_0,theta_pt+p_0,theta_ul+p_0)
  CI_with_correction=c(theta_ll_0+p_0,theta_pt_0+p_0,theta_ul_0+p_0)
  out=list(p_no_correction,p_correction,CI_withot_correction,CI_with_correction)
  names(out)=c(" p_no_correction","p_correction","CI_withot_correction","CI_with_correction")
  return(out)
}
########################
##gsd_ci_est_3_stage(p_0,n_2,n_3,r_2_e,r_3,x_last,last_vt,alpha)
#gsd_ci_est_3_stage(p_0=0.2,n_2=100,n_3=105,r_2_e=31,r_3=29,x_last=31,last_vt=3,alpha=0.025)
#gsd_ci_est_3_stage(p_0=0.2,n_2=100,n_3=105,r_2_e=NA,r_3=29,x_last=31,last_vt=3,alpha=0.025)
#gsd_ci_est_3_stage(p_0=0.2,n_2=70,n_3=10,r_2_e=NA,r_3=29,x_last=31,last_vt=3,alpha=0.025)
######################################################
##############################################

#' Final analysis: two stage design, with sample size change
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param n_1   Number of patients at first look.
#' @param n_2   Number of patients at second look.
#' @param r_1_e   The trial would be stopped for superiority if at least r_1_e responses are observed
#' at the first look.
#' @param r_2    If no more than r_2 responses are observed at the second look, then the null hypothesis will not be rejected.
#' @param n_new  New sample size after the interim analysis.
#' @param x_1    Number of responses at the first look.
#' @param x_new  Number of responses at the final visit.
#' @param alpha    One-sided type I error.
#'
#' @return UL: upper limit of confidence interval.
#' @return LL: lower limit of confidence interval.
#' @return The results with continuity correction are not necessarily different from those with continuity correction.
#' @return The need for continuity correction can be determined through simulations on type I error (see Gao, Zhang, 2004).
#'
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs.
#' Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_ci_est_one_arm_2_stage( p_0=0.2,n_1=100,n_2=105,r_1_e=33,r_2=29,n_new=115,
#' x_1=28,x_new=34,alpha=0.025)
#'
#' @export
asd_ci_est_one_arm_2_stage<-function(p_0,n_1,n_2,r_1_e,r_2,n_new,x_1,x_new,alpha){
  #this function calculates point estimatre and confidence interval after sample size change

  #typically, alpha=0.05
  # x_1 is the number of responses at n_1
  #n_new is the new sample size
  # x_new is the number of respnses at n_new
  N=c(n_1,n_2)
  s=N/n_2
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################

  # statistics for final analysis

  p_hat_new=x_new/n_new
  p_hat_1=x_1/n_1

  var_new=p_hat_new*(1-p_hat_new)
  t_new=n_new/var_new
  t=N/var_new
  T=t[2]

  wald_new=(p_hat_new-p_0)*sqrt(t_new)
  score_new=(p_hat_new-p_0)*t_new
  score_1=(p_hat_1-p_0)*t[1]

  ### convert r_2 to corresponding c_2 wihtout continuity correction ####

  if( is.na(r_1_e)){

    c_2=(r_2/n_2-p_0)*sqrt(n_2^3/(r_2*(n_2-r_2)))
    ### convert r_2 to corresponding c_2 wiht continity correction ####
    r_2_0=r_2+0.5   #continuity correction
    c_2_0=(r_2_0/n_2-p_0)*sqrt(n_2^3/(r_2_0*(n_2-r_2_0)))
    ####################################################

    #######################
    cal1<-function(theta){
      score_back_theta=(sqrt(t[2]-t[1])/sqrt(t_new-t[1]))*(score_new-score_1-theta*(t_new-t[1]))+score_1+theta*(t[2]-t[1])
      wald_back_theta=score_back_theta/sqrt(t[2])
      u_bry_theta=wald_back_theta-theta*sqrt(T)
      tmp = 1-stats::pnorm(u_bry_theta)
      return(tmp)
    }

    ##########################################################
    ##### CI without continuity correction ###
    theta_ll=as.numeric(stats::uniroot(function(x) cal1(x)-alpha,c(-2^32,2^32),tol=1e-5)[1])
    theta_pt=as.numeric(stats::uniroot(function(x) cal1(x)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    theta_ul=as.numeric(stats::uniroot(function(x) cal1(x)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
    ##########################
    p_value=cal1(0)
    ##### CI with continuity correction ###

    ##########################################
    CI_1=c(theta_ll+p_0,theta_pt+p_0,theta_ul+p_0)

    out=list(p_value,CI_1)
    names(out)=c(" p_value","Confidence interval")
  }else{
    r=c(r_1_e,r_2)
    c_bry=(r/N-p_0)*sqrt(N^3/(r*(N-r)))
    e_bry=sqrt(s)*c_bry
    ### convert r_2 to corresponding c_2 wiht continity correction ####
    r_1_0=r_1_e+0.5
    r_2_0=r_2+0.5   #continuity correction
    r_0=c(r_1_0,r_2_0)
    c_bry_0=(r_0/N-p_0)*sqrt(N^3/(r_0*(N-r_0)))

    e_bry_0=sqrt(s)*c_bry_0

    ####################################################
    #######################

    cal2<-function(theta,c_bry){
      score_back_theta=(sqrt(t[2]-t[1])/sqrt(t_new-t[1]))*(score_new-score_1-theta*(t_new-t[1]))+score_1+theta*(t[2]-t[1])
      wald_back_theta=score_back_theta/sqrt(t[1])

      ####################Boundary for vector of Brownian  motions#########################
      u_c_bry=c(c_bry[1],wald_back_theta)
      u_bry_theta=u_c_bry*sqrt(s)-theta*sqrt(T)*s

      ###################################################

      tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,2),sigma=cov_mtx,abseps=1e-8,algorithm=mvtnorm::Miwa())[1]

      return(tmp)
    }

    ##########################################################
    ##### CI without continuity correction ###
    theta_ll=as.numeric(stats::uniroot(function(x) cal2(x,c_bry)-alpha,c(-2^32,2^32),tol=1e-5)[1])
    theta_pt=as.numeric(stats::uniroot(function(x) cal2(x,c_bry)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    theta_ul=as.numeric(stats::uniroot(function(x) cal2(x,c_bry)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
    ##########################
    p_no_correction=cal2(0,c_bry)
    ##### CI with continuity correction ###
    theta_ll_0=as.numeric(stats::uniroot(function(x) cal2(x,c_bry_0)-alpha,c(-2^32,2^32),tol=1e-5)[1])
    theta_pt_0=as.numeric(stats::uniroot(function(x) cal2(x,c_bry_0)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    theta_ul_0=as.numeric(stats::uniroot(function(x) cal2(x,c_bry_0)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
    ##########################
    p_correction=cal2(0,c_bry_0)
    ##########################################
    CI_1=c(theta_ll+p_0,theta_pt+p_0,theta_ul+p_0)
    CI_1_0=c(theta_ll_0+p_0,theta_pt_0+p_0,theta_ul_0+p_0)
    out=list(p_no_correction,p_correction,CI_1,CI_1_0)
    names(out)=c(" p_no_correction","p_correction","CI_withot_correction","CI_with_correction")
  }
  return(out)
}
#############################################
#asd_ci_est_one_arm_2_stage(p_0,n_1,n_2,r_1_e,r_2,n_new,x_1,x_new,alpha)
#asd_ci_est_one_arm_2_stage( p_0=0.2,n_1=100,n_2=105,r_1_e=33,r_2=29,n_new=115,x_1=28,x_new=34,alpha=0.025)
#asd_ci_est_one_arm_2_stage( p_0=0.2,n_1=100,n_2=105,r_1_e=NA,r_2=29,n_new=115,x_1=28,x_new=34,alpha=0.025)

##########################################################################
##############################################

#' Final analysis: three stage design, with sample size change
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param n_2   Number of patients at second look.
#' @param n_3   Number of patients at third look.
#' @param r_2_e   The trial would be stopped for superiority if at least r_2_e responses are observed
#' at the second look. If r_2_e was not used in the design, enter NA.
#' @param r_2    If no more than r_2 responses are observed at the second look, then the trial would be stopped for futility.
#' @param r_3    If no more than r_3 responses are observed at the third look, then the null hypothesis will not be rejected.
#' @param n_new  New sample size after the interim analysis (at the second look).
#' @param x_2    Number of responses at the second look.
#' @param x_new  Number of responses at the final visit.
#' @param alpha    One-sided type I error.
#'
#' @return UL: upper limit of confidence interval.
#' @return LL: lower limit of confidence interval.
#' @return The results with continuity correction are not necessarily different from those with continuity correction.
#' @return The need for continuity correction can be determined through simulations on type I error (see Gao, Zhang, 2004).
#'
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs.
#' Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_ci_est_one_arm_3_stage( p_0=0.2,n_2=105,n_3=110,r_2_e=33,r_3=29,n_new=115,
#' x_2=28,x_new=34,alpha=0.025)
#'
#' @export
asd_ci_est_one_arm_3_stage<-function(p_0,n_2,n_3,r_2_e,r_2,r_3,n_new,x_2,x_new,alpha){
  #this function calculates point estimatre and confidence interval after sample size change
  #typically, alpha=0.05
  # x_2 is the number of responses at n_2
  #n_new is the new sample size
  # x_new is the number of respnses at n_new
  N=c(n_2,n_3)
  s=N/n_3
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx=matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################

  # statistics for final analysis

  p_hat_new=x_new/n_new
  p_hat_2=x_2/n_2

  var_new=p_hat_new*(1-p_hat_new)
  t_new=n_new/var_new
  t=N/var_new
  T=t[2]

  wald_new=(p_hat_new-p_0)*sqrt(t_new)
  score_new=(p_hat_new-p_0)*t_new
  score_2=(p_hat_2-p_0)*t[1]

  #########################################

  if( is.na(r_2_e)){

    c_3=(r_3/n_3-p_0)*sqrt(n_3^3/(r_3*(n_3-r_3)))
    ### convert r_2 to corresponding c_2 wiht continity correction ####
    r_3_0=r_3+0.5   #continuity correction
    c_3_0=(r_3_0/n_3-p_0)*sqrt(n_3^3/(r_3_0*(n_3-r_3_0)))
    ####################################################
    #######################

    cal1<-function(theta){
      score_back_theta=(sqrt(t[2]-t[1])/sqrt(t_new-t[1]))*(score_new-score_2-theta*(t_new-t[1]))+score_2+theta*(t[2]-t[1])
      wald_back_theta=score_back_theta/sqrt(t[2])
      u_bry_theta=wald_back_theta-theta*sqrt(T)
      ###################################################
      tmp = 1-stats::pnorm(u_bry_theta)
      return(tmp)
    }

    ##########################################################
    ##### CI without continuity correction ###
    theta_ll=as.numeric(stats::uniroot(function(x) cal1(x)-alpha,c(-2^32,2^32),tol=1e-5)[1])
    theta_pt=as.numeric(stats::uniroot(function(x) cal1(x)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    theta_ul=as.numeric(stats::uniroot(function(x) cal1(x)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
    ##########################
    p_value=cal1(0)
    ##### CI with continuity correction ###

    ##########################################
    CI_1=c(theta_ll+p_0,theta_pt+p_0,theta_ul+p_0)

    out=list(p_value,CI_1)
    names(out)=c(" p_value","Confidence interval")
  }else{
    r=c(r_2_e,r_3)
    c_bry=(r/N-p_0)*sqrt(N^3/(r*(N-r)))
    e_bry=sqrt(s)*c_bry
    ### convert r_2 to corresponding c_2 wiht continity correction ####
    r_2_0=r_2_e+0.5
    r_3_0=r_3+0.5   #continuity correction
    r_0=c(r_2_0,r_3_0)
    c_bry_0=(r_0/N-p_0)*sqrt(N^3/(r_0*(N-r_0)))

    e_bry_0=sqrt(s)*c_bry_0

    ####################################################
    #######################

    cal2<-function(theta,c_bry){
      score_back_theta=(sqrt(t[2]-t[1])/sqrt(t_new-t[1]))*(score_new-score_2-theta*(t_new-t[1]))+score_2+theta*(t[2]-t[1])
      wald_back_theta=score_back_theta/sqrt(t[1])
      ####################Boundary for vector of Brownian  motions#########################
      u_c_bry=c(c_bry[1],wald_back_theta)
      u_bry_theta=u_c_bry*sqrt(s)-theta*sqrt(T)*s
      ###################################################
      tmp = 1-mvtnorm::pmvnorm(lower=-Inf,upper=u_bry_theta,mean=rep(0,2),sigma=cov_mtx,abseps=1e-8,algorithm=mvtnorm::Miwa())[1]
      return(tmp)
    }

    ##########################################################
    ##### CI without continuity correction ###
    theta_ll=as.numeric(stats::uniroot(function(x) cal2(x,c_bry)-alpha,c(-2^32,2^32),tol=1e-5)[1])
    theta_pt=as.numeric(stats::uniroot(function(x) cal2(x,c_bry)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    theta_ul=as.numeric(stats::uniroot(function(x) cal2(x,c_bry)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
    ##########################
    p_no_correction=cal2(0,c_bry)
    ##### CI with continuity correction ###
    theta_ll_0=as.numeric(stats::uniroot(function(x) cal2(x,c_bry_0)-alpha,c(-2^32,2^32),tol=1e-5)[1])
    theta_pt_0=as.numeric(stats::uniroot(function(x) cal2(x,c_bry_0)-0.5,c(-2^32,2^32),tol=1e-5)[1])
    theta_ul_0=as.numeric(stats::uniroot(function(x) cal2(x,c_bry_0)-1+alpha,c(-2^32,2^32),tol=1e-5)[1])
    ##########################
    p_correction=cal2(0,c_bry_0)
    ##########################################
    CI_1=c(theta_ll+p_0,theta_pt+p_0,theta_ul+p_0)
    CI_1_0=c(theta_ll_0+p_0,theta_pt_0+p_0,theta_ul_0+p_0)
    out=list(p_no_correction,p_correction,CI_1,CI_1_0)
    names(out)=c(" p_no_correction","p_correction","CI_withot_correction","CI_with_correction")
  }
  return(out)
}
#############################################
#asd_ci_est_one_arm_3_stage(p_0,n_2,n_3,r_2_e,r_2,r_3,n_new,x_2,x_new,alpha)
#asd_ci_est_one_arm_3_stage( p_0=0.2,n_2=105,n_3=110,r_2_e=33,r_3=29,n_new=115,x_2=28,x_new=34,alpha=0.025)
#asd_ci_est_one_arm_3_stage( p_0=0.2,n_2=105,n_3=110,r_2_e=NA,r_3=29,n_new=115,x_2=28,x_new=34,alpha=0.025)
#asd_ci_est_one_arm_3_stage( p_0=0.2,n_2=105,n_3=110,r_2_e=30,r_3=29,n_new=120,x_2=28,x_new=34,alpha=0.025)
#asd_ci_est_one_arm_3_stage( p_0=0.2,n_2=105,n_3=110,r_2_e=NA,r_3=29,n_new=115,x_2=28,x_new=34,alpha=0.025)
##########################################################################
