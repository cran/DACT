
## The code in this file contains the functions required for the hybrid sample size re-estimation procedure ###
#rm(list=ls())

#library(tidyverse)
#library("coin")
#library(sm)
##############################
#################
## Functions ####
#################
#############################


References1 = "Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: \
                     using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical \
                     Statistics,2024 Aug;34(5):737-752. doi: 10.1080/10543406.2024.2358796. Epub 2024 May 30. \
                     <a href='https://doi.org/10.1080/10543406.2024.2358796' target='_blank'>https://doi.org/10.1080/10543406.2024.2358796</a>"
References2 = "Mehta, C. R., Pocock, S. J. Adaptive increase in sample size when interim results are promising: A practical  \
                     guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284  \
                     <a href='https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.4102' target='_blank'> \
                      https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.4102</a>"
References3 = "Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian \
                    approach, Stat in Med, 2024"
References4 = "Zhu, H.,Lakkis, H.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387 February 2014"
References5 = "Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503."

##############################################
c_bry_adj_f=function(t_new,T_planned,c_bry,t_inter,score_inter){
  #  This function adjusts the critical boundary [c_bry] when there the final analysis is changed from T_planned to a new time t_new
  # t_new is the new information time
  # T_planned is the planned information time
  ## c_bry is the original critical boundary associated with T_planned
  # t_inter is the information time at the interim analysis
  # score_inter=t_inter*theta_hat_inter
  c_bry_new =(1/sqrt(t_new))*((sqrt(t_new-t_inter)/sqrt(T_planned-t_inter))*(c_bry*sqrt(T_planned)- score_inter)+score_inter)
  return(c_bry_new)
}
################################################

#c_bry_adj_f(t_new=150,T_planned=100,c_bry=1.96,t_inter=50,theta_hat_inter=0.3)

###################################################
### Conditional power and new sample size    ######
###################################################
##################################################

cond_power=function(T,c_bry,t_inter,theta_hat_inter,theta){
  # This function calculates the conditional power at the information time T_planned
  # T_planned is the information time at which the planned final analysis will be conducted
  # c_bry is the critical boundary associated with T_planned at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # theta ims the assumed effect size
  score_inter=theta_hat_inter*t_inter
  cond_pw=stats::pnorm(theta*sqrt(T-t_inter)-(c_bry*sqrt(T)-score_inter)/sqrt(T-t_inter))

  return(cond_pw)
}

#############################################################
#cond_power(T=100,c_bry=1.96,t_inter=50,theta_hat_inter=0.3,theta=0.3)
#cond_power(T=100,c_bry=1.96,t_inter=50,theta_hat_inter=0.2,theta=0.3)
#cond_power_1(T=100,c_bry=1.96,t_inter=50,theta_hat_inter=0.2,theta=0.3)
###########################
new_time_cp=function(T_planned,c_bry,t_inter,theta_hat_inter,beta){
  # This function calculates the required information time at which the mconditional power will be 1-beta
  # T_planned is the information time at which the planned final analysis will be conducted
  # c_bry is the critical boundary associated with T_planned at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # beta is the type II error level

  score_inter=theta_hat_inter*t_inter

  T_new_cp=theta_hat_inter^(-2)*((c_bry*sqrt(T_planned)-score_inter)/sqrt(T_planned-t_inter)+stats::qnorm(1-beta))^2+t_inter

  return(T_new_cp)
}

#new_time_cp(T_planned=100,c_bry=1.96,t_inter=50,theta_hat_inter=0.3,beta=0.1)
#new_time_cp(T_planned=100,c_bry=1.96,t_inter=20,theta_hat_inter=0.5,beta=0.1)

###################################################
####    predictive probability of success      ####
###################################################

#pp(T=52,c_bry=1.96,t_inter=20,theta_hat_inter=0.2)
#pp(T=1000000,c_bry=1.96,t_inter=20,theta_hat_inter=0.2)

pp=function(T,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter){

  # This function calculates the predictive power
  # T is the information time at which the final analysis will be conducted
  # c_bry is the critical boundary at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0

  score_inter=t_inter*theta_hat_inter

  Wald_inter=sqrt(t_inter)*theta_hat_inter
  t_1=t_prior+t_inter
#  tao_1=1/sqrt(t_1)
  theta_1=(t_prior*theta_prior+t_inter*theta_hat_inter)/t_1
  integ =function(x){
    temp_0=stats::pnorm((x/sqrt(t_1)+theta_1)*sqrt(T-t_inter)-(c_bry*sqrt(T)-score_inter)/sqrt(T-t_inter))*stats::dnorm(x)

  }
  #  a=1-stats::integrate( Vectorize(integ),lower= -Inf,upper = +Inf,abs.tol=1e-6)$value
  a=stats::integrate(integ,lower=-Inf,upper =Inf,abs.tol=1e-6)$value
  return(a)
}
#pp(T,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)
#########################################################
#pp(T=52,c_bry=1.96,t_prior=0,theta_prior=0,t_inter=20,theta_hat_inter=0.2)
#pp(T=52,c_bry=1.96,t_prior=10,theta_prior=0.2,t_inter=20,theta_hat_inter=0.2)

#pp(T=1000000,c_bry=1.96,t_inter=20,theta_hat_inter=0.2)

new_pp=function(t,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter){
  # This function calculates the predictive power at a new information time point t, instead of T_planned
  # t is the new information time at which the re-scheduled  final analysis will be conducted
  # T_planned is the information time at which the planned final analysis will be conducted
  # c_bry is the critical boundary associated with T_planned at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0

  c_bry_adj=c_bry_adj_f(t,T_planned,c_bry,t_inter,theta_hat_inter)
  # c_bry_adj is the adjusted crutical boundary after changing the final anlysis from T_planned to t
  new_pp=pp(t,c_bry=c_bry_adj,t_prior,theta_prior,t_inter,theta_hat_inter)
  #  out=cbind(new_pp)
  return(new_pp)
}
#new_pp(t,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)
##################################################

#new_pp(t=120,t_prior=20,theta_prior=0.3,T_planned=100,t_inter=50,c_bry=1.96,theta_hat_inter=0.3)
#new_pp(t=120,t_prior=0,theta_prior=0,T_planned=100,t_inter=50,c_bry=1.96,theta_hat_inter=0.3)
#new_pp(t=120,t_prior=0,theta_prior=0,T_planned=100,t_inter=50,c_bry=1.96,theta_hat_inter=0.3)
#########################################################################
###################

new_time_pp=function(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta){

  # This function calculates the required information time at which the modified predictive power will be 1-beta
  # T_planned is the information time at which the planned final analysis will be conducted
  # c_bry is the critical boundary associated with T_planned at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis

  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0
  ############################################################################
  pp_T_planned=new_pp(t=T_planned,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)
  pp_t_large=new_pp(t=2^32,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)
  if(pp_T_planned>=1-beta){
    T_new_pp=T_planned
  }else if (pp_t_large<1-beta){
    T_new_pp=9999

  }else {
    ##############  obtain new information time  #######################################
    T_new_pp=as.numeric(stats::uniroot(function(x) new_pp(x,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)-(1-beta),c(T_planned,2^32),tol=1e-5)[1])
  }
  ######################
#  new_pp(t=T_planned,t_prior,theta_prior,T_planned,t_inter,c_bry,theta_hat_inter)
#  new_pp(t=2^32,t_prior,theta_prior,T_planned,t_inter,c_bry,theta_hat_inter)

  #  new_pp(t=T_planned,t_prior=0,theta_prior=0,T_planned,t_inter,c_bry,theta_hat_inter)
  #  new_pp(t=2^32,t_prior=0,theta_prior=0,T_planned,t_inter,c_bry,theta_hat_inter)
  ###########
  return(T_new_pp)
}

#new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
#new_time_pp(t_prior=30,theta_prior=0.3,T_planned=100,c_bry=1.96,t_inter=50,theta_hat_inter=0.3,beta=0.1)
#new_time_pp(t_prior=0,theta_prior=0,T_planned=100,c_bry=1.96,t_inter=50,theta_hat_inter=0.3,beta=0.1)


##############################################################################

##################################################
#################################################################
#### modified predictive probability of success  ####
#################################################################

#########################################################
#mpp(T=52,c_bry=1.96,t_inter=20,theta_hat_inter=0.2)
#mpp(T=1000000,c_bry=1.96,t_inter=20,theta_hat_inter=0.2)
mpp=function(T,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter){

  # This function calculates the modified predictive power
  # T is the information time at which the final analysis will be conducted
  # c_bry is the critical boundary at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0

  score_inter=t_inter*theta_hat_inter

  Wald_inter=sqrt(t_inter)*theta_hat_inter
  t_1=t_prior+t_inter
  #  tao_1=1/sqrt(t_1)
  theta_1=(t_prior*theta_prior+t_inter*theta_hat_inter)/t_1
  integ =function(x){
    temp_0=stats::pnorm((x/sqrt(t_1)+theta_1)*sqrt(T-t_inter)-(c_bry*sqrt(T)-score_inter)/sqrt(T-t_inter))*stats::dnorm(x)

  }

  #  a=1-stats::integrate( Vectorize(integ),lower= -Inf,upper = +Inf,abs.tol=1e-6)$value
  a=stats::integrate(integ,lower=-Wald_inter,upper =Inf,abs.tol=1e-6)$value
  b=a/stats::pnorm(Wald_inter)
  return(b)
}
#########################################################
#mpp(T,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)
#mpp(T=52,c_bry=1.96,t_prior=0,theta_prior=0,t_inter=20,theta_hat_inter=0.2)
#mpp(T=52,c_bry=1.96,t_prior=10,theta_prior=0.2,t_inter=20,theta_hat_inter=0.2)

#mpp(T=1000000,c_bry=1.96,t_inter=20,theta_hat_inter=0.2)

########################
##################################################

new_mpp=function(t,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter){
  # This function calculates the modified predictive power at a new information time point t, instead of T_planned
  # t is the new information time at which the re-scheduled  final analysis will be conducted
  # T_planned is the information time at which the planned final analysis will be conducted
  # c_bry is the critical boundary associated with T_planned at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0

  c_bry_adj=c_bry_adj_f(t,T_planned,c_bry=c_bry,t_inter,theta_hat_inter)
  # c_bry_adj is the adjusted crutical boundary after changing the final anlysis from T_planned to t
  new_mpp=mpp(t,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)
  #  out=cbind(c_bry_adj,new_mpp)
  return(new_mpp)
}
#new_mpp(t,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)
##############
#new_mpp(t=120,t_prior=0,theta_prior=0,T_planned=100,t_inter=50,c_bry=1.96,theta_hat_inter=0.3)
#new_mpp(t=2^32,t_prior=0,theta_prior=0,T_planned=100,t_inter=50,c_bry=1.96,theta_hat_inter=0.3)

#new_mpp(t=120,t_prior=20,theta_prior=0.3,T_planned=100,t_inter=50,c_bry=1.96,theta_hat_inter=0.3)

#########################################################################

#new_time_mpp=function(T_planned,c_bry,t_inter,theta_hat_inter,beta){


###################
new_time_mpp=function(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta){

  # This function calculates the required information time at which the modified predictive power will be 1-beta
  # T_planned is the information time at which the planned final analysis will be conducted
  # c_bry is the critical boundary associated with T_planned at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # beta is type II error
  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0

  ##############################
  mpp_T_planned=new_mpp(t=T_planned,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)
  mpp_t_large=new_mpp(t=2^32,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)

  if(mpp_T_planned>=1-beta){
    T_new_mpp=T_planned
  }else if (mpp_t_large<1-beta){
    T_new_mpp=9999

  }else {
    ##############  obtain new information time using normalized pp #######################################
    T_new_mpp=as.numeric(stats::uniroot(function(x) new_mpp(x,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter)-(1-beta),c(T_planned,2^32),tol=1e-5)[1])
  }
  ######################
  ######################
  return(T_new_mpp)
}

###############################################

#new_time_mpp(T_planned=100,c_bry=1.96,t_prior=0,theta_prior=0,t_inter=50,theta_hat_inter=0.3,beta=0.1)
#new_time_mpp(T_planned=100,c_bry=1.96,t_prior=20,theta_prior=0.3,t_inter=50,theta_hat_inter=0.3,beta=0.1)
#new_time_mpp(T_planned=100,c_bry=1.96,t_prior=0,theta_prior=0,t_inter=50,theta_hat_inter=0.3,beta=0.1)

#new_time_mpp(T_planned=100,c_bry=1.96,t_inter=20,theta_hat_inter=0.5,beta=0.1)
##############################################################################

########################################################################################
#### modified predictive probability of success   --- using the (1-gamma) interval  ####
########################################################################################

ci_mpp=function(T,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma){

  # This function calculates the credible interval modified predictive power
  # T is the information time at which the final analysis will be conducted
  # c_bry is the critical boundary at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0

  score_inter=t_inter*theta_hat_inter

  Wald_inter=sqrt(t_inter)*theta_hat_inter

  a=stats::qnorm(1-gamma/2)
  a_upper=a
  a_lower=max(c(-Wald_inter,-a))
  a_dinom=stats::pnorm(a_upper)-stats::pnorm(a_lower)

  t_1=t_prior+t_inter
  theta_1=(t_prior*theta_prior+t_inter*theta_hat_inter)/t_1

  integ =function(x){
    temp_0=stats::pnorm((x/sqrt(t_1)+theta_1)*sqrt(T-t_inter)-(c_bry*sqrt(T)-score_inter)/sqrt(T-t_inter))*stats::dnorm(x)

  }


  #  a=1-stats::integrate( Vectorize(integ),lower= -Inf,upper = +Inf,abs.tol=1e-6)$value
  b_int=stats::integrate(integ,lower=a_lower,upper=a_upper,abs.tol=1e-6)$value
  b_out=b_int/a_dinom
  return(b_out)
}
#########################################################
#ci_mpp(T,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma)
#ci_mpp(T=52,c_bry=1.96,t_prior=0,theta_prior=0,t_inter=20,theta_hat_inter=0.2,gamma=0.5)
#ci_mpp(T=52,c_bry=1.96,t_prior=10,theta_prior=0.2,t_inter=20,theta_hat_inter=0.2,gamma=0.5)


#########################################################################
new_ci_mpp=function(t,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma){
  # This function calculates the credible interval modified predictive power at a new information time point t, instead of T_planned
  # t is the new information time at which the re-scheduled  final analysis will be conducted
  # T_planned is the information time at which the planned final analysis will be conducted
  # c_bry is the critical boundary associated with T_planned at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0

  c_bry_adj=c_bry_adj_f(t,T_planned,c_bry,t_inter,theta_hat_inter)
  # c_bry_adj is the adjusted crutical boundary after changing the final anlysis from T_planned to t
  new_mpp=ci_mpp(t,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma)
  #  out=cbind(c_bry_adj,new_mpp)
  return(new_mpp)
}
#new_ci_mpp(t,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma)
##################################################
#new_ci_mpp(t=120,T_planned=100,t_prior=0,theta_prior=0,t_inter=50,c_bry=1.96,theta_hat_inter=0.3,gamma=0.1)
#new_ci_mpp(t=120,T_planned=100,t_prior=20,theta_prior=0.3,t_inter=50,c_bry=1.96,theta_hat_inter=0.3,gamma=0.1)

#new_ci_mpp(t=120,T_planned=100,t_inter=50,c_bry=1.96,theta_hat_inter=0.3,gamma=0.5)
##############################################

#########################################################################

new_time_ci_mpp=function(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma, beta){
  # This function calculates the required information time at which the credible interval modified predictive power will be 1-beta
  # T_planned is the information time at which the planned final analysis will be conducted
  # c_bry is the critical boundary associated with T_planned at the final analysis
  # t_inter is the information time at the interim analysis
  # theta_hat_inter is the estimated theta at the interim analysis
  # t_prior is the information time from the Bayesian prior
  # theta_prior is the theta from the Bayesian prior
  # when non-informative prior is used, t_prior=0 and  theta_prior=0
  # gamma defines the credible interval
  # beta is type II error
  ##############################
  #ci_mpp_T_planned=ci_mpp(T=T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma)
  ci_mpp_T_planned=new_ci_mpp(t=T_planned,T_planned=T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma)

  ci_mpp_t_large=new_ci_mpp(t=2^32,T_planned,t_prior,theta_prior,t_inter,c_bry,theta_hat_inter,gamma)

  if(ci_mpp_T_planned>=1-beta){
    T_new_ci_mpp=T_planned
  }else if (ci_mpp_t_large<1-beta){
    T_new_ci_mpp=9999

  }else {
    ##############  obtain new information time using normalized pp #######################################
    T_new_ci_mpp=as.numeric(stats::uniroot(function(x) new_ci_mpp(x,T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma)-(1-beta),c(T_planned,2^32),tol=1e-5)[1])
  }
  ######################
  ######################
  return(T_new_ci_mpp)
}

#########################################################################
#new_time_ci_mpp(T_planned=100,c_bry=1.96,t_prior=0,theta_prior=0,t_inter=50,theta_hat_inter=0.3,gamma=0.5,beta=0.1)
#new_time_ci_mpp(T_planned=100,c_bry=1.96,t_prior=20,theta_prior=0.3,t_inter=50,theta_hat_inter=0.3,gamma=0.5,beta=0.1)
#new_time_ci_mpp( T_planned=22.96142, c_bry=1.959964,t_prior=10,  theta_prior=0.5108256,t_inter= 11.48071, theta_hat_inter=0.1202579, gamma= 0.5, beta= 0.1)

#new_time_ci_mpp(T_planned=100,c_bry=1.96,t_inter=20,theta_hat_inter=0.05,gamma=0.5,beta=0.1)
#######

#######################################################################
#################
#' Hybrid: New sample size estimation
#'
#' @param N_planned Planned sample size: the planned sample size at the final analysis
#' @param c_bry Critical boundary or threshold for rejecting the null hypothesis at the final analysis
#' @param t_prior Prior information time: the prior distribution of theta is assumed to have a
#'   normal distribution, and the Fisher's information time is the reciprocal of the variance of
#'   the prior distribution; set t_prior = 0 for a non-informative prior
#' @param theta_prior Mean effect size of prior: the mean of the normal prior distribution of
#'   theta; set theta_prior = 0 for a non-informative prior
#' @param N_inter Sample size at interim analysis: the number of subjects with data at the
#'   interim analysis
#' @param t_inter Information time at interim analysis: the Fisher information at the interim
#'   analysis, the reciprocal of the variance of the unbiased maximum likelihood estimate
#'   theta_hat_inter
#' @param theta_hat_inter Estimated effect size at interim analysis: the maximum likelihood
#'   estimate of theta at the interim analysis
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{New sample sizes}{A data frame with the columns "prior type", "SSR method" and "new sample size", giving the new sample size for every combination of prior type and sample size re-estimation (SSR) method. The prior types are the non-informative prior (mean of the prior zero and information time zero) and the informative prior entered by the user. Within each prior type the four SSR methods are the conditional power (cp), the predictive power (pp), the modified predictive power (mpp) and the credible interval modified predictive power (ci_mpp), each being the sample size that provides 1 - beta of that power; the sample size is NA when the target conditional power cannot be reached. For the informative prior the cp entry is based on the Bayesian posterior mean of theta}
#'   \item{Reference}{References for the methods}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' new_sample_size_hybrid(N_planned = 100, c_bry = 1.96, t_prior = 45, theta_prior = 0.5,
#' N_inter = 60, t_inter = 50, theta_hat_inter = 0.3, gamma = 0.5, beta = 0.1)
#'
#' @export
new_sample_size_hybrid<-function(N_planned,c_bry,t_prior,theta_prior,N_inter,t_inter,theta_hat_inter,gamma, beta){

  theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(theta_prior+theta_hat_inter)

  # wald_inter=theta_hat_inter*sqrt(t_inter)
  #  score_inter=theta_hat_inter*t_inter

  T_planned=t_inter*N_planned/N_inter
  ################# cp new sample size ####

  new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
  new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
  new_sampsz_cp=max(new_sampsz_0_cp,ceiling(N_inter*1.1))
  ################# cp_Bayesian new sample size ####

  new_time_cp_bayesian_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=theta_hat_bayesian,beta)
  new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_bayesian_0/t_inter)
  new_sampsz_cp_bayesian=max(new_sampsz_0_cp,ceiling(N_inter*1.1))
  ####### pp informative prior  new sample size ####
  new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)

  if (new_time_pp_inf_prior==9999){
    new_sampsz_pp_inf_prior=NA
  }else if (new_time_pp_inf_prior!=9999){

    new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
    new_sampsz_pp_inf_prior=max(new_sampsz_0_pp_inf_prior,ceiling(N_inter*1.1))
  }
  ####### pp non-informative prior  new sample size ####

  new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)

  if (new_time_pp_non_inf_prior==9999){

    new_sampsz_pp_non_inf_prior=NA

  }else if (new_time_pp_non_inf_prior!=9999){

    new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
    new_sampsz_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,ceiling(N_inter*1.1))
  }
  ####### mpp informative prior  new sample size ####
  new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
  if (new_time_mpp_inf_prior==9999){
    new_sampsz_mpp_inf_prior=NA
  }else if (new_time_mpp_inf_prior!=9999){

    new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
    new_sampsz_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,ceiling(N_inter*1.1))
  }
  ####### mpp non-informative prior  new sample size ####

  new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)

  if (new_time_mpp_non_inf_prior==9999){
    new_sampsz_mpp_non_inf_prior=NA
  }else if (new_time_mpp_non_inf_prior!=9999){

    new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
    new_sampsz_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,ceiling(N_inter*1.1))
  }
  ####### ci_mpp informative prior  new sample size ####
  new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)


  if (new_time_ci_mpp_inf_prior==9999){
    new_sampsz_ci_mpp_inf_prior=NA
  }else if (new_time_ci_mpp_inf_prior!=9999){

    new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
    new_sampsz_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,ceiling(N_inter*1.1))
  }
  ####### ci_mpp non-informative prior  new sample size ####
  new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)
  if (new_time_ci_mpp_non_inf_prior==9999){
    new_sampsz_ci_mpp_non_inf_prior=NA

  }else if (new_time_ci_mpp_non_inf_prior!=9999){
    new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
    new_sampsz_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,ceiling(N_inter*1.1))
  }
  #################

  names_2=rep(c("","cp","pp","mpp","ci_mpp"),2)
  names_1=c("non-infromative prior","","","","",
            "informatrive prior","","","","")
  new_samsz=c("",new_sampsz_cp,
              new_sampsz_pp_non_inf_prior,
              new_sampsz_mpp_non_inf_prior,
              new_sampsz_ci_mpp_non_inf_prior,
              "",
              new_sampsz_cp_bayesian,
              new_sampsz_pp_inf_prior,
              new_sampsz_mpp_inf_prior,
              new_sampsz_ci_mpp_inf_prior)

  out_0=cbind(names_1,names_2,new_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("prior type","SSR method","new sample size")
  out_2=References3

  out=list(out_1,out_2)
  names(out)=c("New sample sizes","Reference")
  return(out)
}
#############
#new_sample_size_hybrid(N_planned=100,c_bry=1.96,t_prior=45,theta_prior=0.5,N_inter=60,t_inter=50,theta_hat_inter=0.3,gamma=0.5, beta=0.1)

################
# function for new sample analysis with normal variable ####
new_samsz_test_normal=function(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz){
  ## This function is used in the simulation function for evaluating trhe operationg charateristics ##

  if (new_samsz==NF){
    test_new_last=Wald_last-c_bry
  }else {
  theta_hat_new_last= mean(x[1:new_samsz])-mean(x0[1:new_samsz])
  t_new_last=new_samsz/(stats::var(x[1:new_samsz])+stats::var(x0[1:new_samsz]))

  T_new=max(t_new_last,t_last)
  #t_inter_adjusted=T_new*event_inter/new_samsz
  # T_planned_adjusted=T_new*event_num/new_samsz
  #new_c_bry=c_bry_adj_f(t_new=T_new,T_planned= T_planned_adjusted,c_bry,t_inter=t_inter_adjusted,theta_hat_inter)
  #  new_c_bry=c_bry_adj_f(t_new=T_new,T_planned= T_planned,c_bry,t_inter=t_inter,theta_hat_inter)
  new_c_bry=c_bry_adj_f(t_new=T_new,T_planned=t_last,c_bry,t_inter,score_inter)

  wald_new_last=theta_hat_new_last*sqrt(t_new_last)
  score_new_last=theta_hat_new_last*t_new_last

  test_new_last=wald_new_last-max(new_c_bry,mean(c(new_c_bry,c_bry)))
  }
  return(test_new_last)
}
#new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz)
############################
################################
# function for new sample analysis with survival endpoint####
new_samsz_test_surv=function(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,event_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz){
  ## This function is used in the simulation function for evaluating trhe operationg charateristics ##
  ## see the function hybrid_pz_cp_rej_comp_survival_simu  ###

  evt_pos_new=evt_count[evt_count>new_samsz] #substring of  evt_count that are larger than new_samsz
  evt_cut_new=match(evt_pos_new[1],evt_count) #the first incidence that evt_count is larger than new_samsz

  #  vector of zeros, length of the vector is the same as the original smaple size
  time_new=c(time_x0[1:evt_cut_new],time_x[1:evt_cut_new]) #combined samples befroe total number of events exceed evt_num
  event_new=c(x0_status[1:evt_cut_new],x_status[1:evt_cut_new]) #censoring status before total number of events exceed evt_num
  grp_new=c(grp_x0[1:evt_cut_new],grp_x[1:evt_cut_new]) #vector of zero, with length=2*evt_cut


  test_new= list(time=time_new,status=event_new,grp_new)#survival data
  coxcf_new=summary(coxph(Surv(time,status)~grp_new,test_new)) #run the coxph function and get the estimates
  # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
  # test is the same as data=test
  theta_hat_new_last= coxcf_new$coefficients[1] #estimated log(HR)
  theta_hat_new_last_sd= coxcf_new$coefficients[3] #estimated standard error of log(HR)
  t_new_last=(theta_hat_new_last_sd)^(-2)
  T_new=max(t_new_last,T_planned)
  #t_inter_adjusted=T_new*event_inter/new_samsz
  # T_planned_adjusted=T_new*event_num/new_samsz
  #new_c_bry=c_bry_adj_f(t_new=T_new,T_planned= T_planned_adjusted,c_bry,t_inter=t_inter_adjusted,theta_hat_inter)
  #  new_c_bry=c_bry_adj_f(t_new=T_new,T_planned= T_planned,c_bry,t_inter=t_inter,theta_hat_inter)
  new_c_bry=c_bry_adj_f(t_new=T_new,T_planned,c_bry,t_inter,score_inter)

  wald_new_last=theta_hat_new_last/ theta_hat_new_last_sd
  score_new_last=theta_hat_new_last*t_new_last

  test_new_last=wald_new_last-max(new_c_bry,mean(new_c_bry,c_bry))
  return(test_new_last)
}
################
# function for new sample analysis with normal variable ####

new_samsz_test_binary=function(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz){
  ## This function is used in the simulation function for evaluating trhe operationg charateristics ##

  if (new_samsz==NF){
    test_new_last=Wald_last-c_bry
  }else {
    theta_hat_new_last= mean(x[1:new_samsz])-mean(x0[1:new_samsz])
    t_new_last=new_samsz/(stats::var(x[1:new_samsz])+stats::var(x0[1:new_samsz]))

    T_new=max(t_new_last,t_last)
    #t_inter_adjusted=T_new*event_inter/new_samsz
    # T_planned_adjusted=T_new*event_num/new_samsz
    #new_c_bry=c_bry_adj_f(t_new=T_new,T_planned= T_planned_adjusted,c_bry,t_inter=t_inter_adjusted,theta_hat_inter)
    #  new_c_bry=c_bry_adj_f(t_new=T_new,T_planned= T_planned,c_bry,t_inter=t_inter,theta_hat_inter)
    new_c_bry=c_bry_adj_f(t_new=T_new,T_planned=t_last,c_bry,t_inter,score_inter)

    wald_new_last=theta_hat_new_last*sqrt(t_new_last)
    score_new_last=theta_hat_new_last*t_new_last

    test_new_last=wald_new_last-max(new_c_bry,mean(c(new_c_bry,c_bry)))
  }
  return(test_new_last)
}

#new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz)

############################
new_samsz_test_poisson=function(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz){
  ## This function is used in the simulation function for evaluating trhe operationg charateristics ##

  if (new_samsz==NF){
    test_new_last=Wald_last-c_bry
  }else {
    theta_hat_new_last= mean(x[1:new_samsz])-mean(x0[1:new_samsz])
    t_new_last=new_samsz/( mean(x[1:new_samsz])+mean(x0[1:new_samsz]))

    T_new=max(t_new_last,t_last)
    #t_inter_adjusted=T_new*event_inter/new_samsz
    # T_planned_adjusted=T_new*event_num/new_samsz
    #new_c_bry=c_bry_adj_f(t_new=T_new,T_planned= T_planned_adjusted,c_bry,t_inter=t_inter_adjusted,theta_hat_inter)
    #  new_c_bry=c_bry_adj_f(t_new=T_new,T_planned= T_planned,c_bry,t_inter=t_inter,theta_hat_inter)
    new_c_bry=c_bry_adj_f(t_new=T_new,T_planned=t_last,c_bry,t_inter,score_inter)

    wald_new_last=theta_hat_new_last*sqrt(t_new_last)
    score_new_last=theta_hat_new_last*t_new_last

    test_new_last=wald_new_last-max(new_c_bry,mean(c(new_c_bry,c_bry)))
  }
  return(test_new_last)
}
################

############################
new_samsz_test_negbinom=function(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz){
  ## This function is used in the simulation function for evaluating trhe operationg charateristics ##

  v_fun=function(r_0,r_1,theta,nu_t,kappa){
    V_01=(1+theta)/(nu_t*theta*r_0)+(1+theta)*kappa/theta
    V_02=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
    V_03=(1+theta)**2/(nu_t*theta*(r_0+theta*r_1))+(1+theta)*kappa/theta
    V_max=max(V_01,V_02,V_02)
    return(V_max)
  }

  if (new_samsz==NF){
    test_new_last=Wald_last-c_bry
  }else {

    #########
    x_c_new_last=x_c[1:new_samsz,]
    x_t_new_last=x_t[1:new_samsz,]

    x_new_last=rbind(x_c_new_last,x_t_new_last)
####################
    negglm_new_last=glm.nb(as.numeric(x_new_last[,2])~as.factor(x_new_last[,1]),link="log")
    sum_new_last=summary(negglm_new_last)

    kappa_hat_new_last=sum_new_last$theta
    r_c_new_last=mean(x_c_new_last[,2])
    r_t_new_last=mean(x_t_new_last[,2])

    v_hat_new_last=v_fun(r_0=r_c_new_last,r_1=r_t_new_last,theta=rand_ratio,nu_t,kappa=kappa_hat_new_last)
    theta_hat_new_last=-log(r_t_new_last/r_c_new_last)
    t_new_last=new_samsz/v_hat_new_last

    T_new=max(t_new_last,t_last)

    new_c_bry=c_bry_adj_f(t_new=T_new,T_planned=t_last,c_bry,t_inter,score_inter)

    wald_new_last=theta_hat_new_last*sqrt(t_new_last)
  #  score_new_last=theta_hat_new_last*t_new_last

   # test_new_last=wald_new_last-max(new_c_bry,mean(c(new_c_bry,c_bry)))
    test_new_last=wald_new_last-new_c_bry

  }
  return(test_new_last)
}
###########
#new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,r_t_pf,r_c_pf,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz)
############################
###      Simulations  for operating characteristics ########
############################

############### Simulation function for negative binomial variable ####
########################

################################
#### Negative binomial distribution #######
################################

####### Fixed  prior ######

#' Hybrid: Simulations, Negative binomial distribution, Using given estimates of informative prior
#'
#' @param t_prior Prior information time: the prior distribution of theta = -log(r_t/r_c) is
#'   assumed to have a normal distribution, and the Fisher's information time is the reciprocal
#'   of the variance of the prior distribution; set t_prior = 0 for a non-informative prior
#' @param theta_prior Mean effect size of prior: the mean of the normal prior distribution of
#'   theta = -log(r_t/r_c); set theta_prior = 0 for a non-informative prior
#' @param r_t Event rate in the test arm; the effect size is theta = -log(r_t/r_c) and a smaller
#'   event rate is assumed to be associated with better efficacy. Enter r_t = r_c to verify type
#'   I error control and r_t < r_c for the power simulation
#' @param r_c Event rate in the control arm; the effect size is theta = -log(r_t/r_c)
#' @param kappa Dispersion parameter of the negative binomial distribution; must be > 0
#' @param nu_t Exposure time for the negative binomial distribution; must be > 0
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz Planned sample size: the planned sample size at the final analysis
#' @param N_max Maximum sample size: a user chosen, maximum allowable sample size for sample size modification
#' @param theta_cut Futility threshold: a cut-off for theta_hat. In the simulations, if
#'   theta_hat at the interim analysis is at or below theta_cut the trial is considered futile
#'   and proceeds to the planned sample size without any change
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods}
#'   \item{futility rate}{The rate at which the observed theta_hat at the interim analysis is at or below theta_cut, so that the trial proceeds without any sample size change}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_fixed_prior_rej_negbinom_simu(t_prior = 5, theta_prior = 0.2, r_t = 1, r_c = 1,
#' kappa = 1, nu_t = 1, rand_ratio = 1, s = c(0.5,1), sampsz = 100, N_max = 150,
#' theta_cut = 0.001, cp_min = 0.27, gamma = 0.5, alpha_2 = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_fixed_prior_rej_negbinom_simu<-function(t_prior,theta_prior,r_t,r_c,kappa,nu_t,rand_ratio,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number
  #####################################################
  ## varaince function per Zhu, harris, 2012 ####
  v_fun=function(r_0,r_1,theta,nu_t,kappa){
    V_01=(1+theta)/(nu_t*theta*r_0)+(1+theta)*kappa/theta
    V_02=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
    V_03=(1+theta)**2/(nu_t*theta*(r_0+theta*r_1))+(1+theta)*kappa/theta
    V_max=max(V_01,V_02,V_02)
    return(V_max)
  }

  ######################
  N=ceiling(sampsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  NF=N[length(s)] # is the sample size at the final analysis
  N_inter=N[1]

  c_bry=stats::qnorm(1-alpha_2/2)
  mu=-log(r_t/r_c)
  ############################################
  ####  beginning of the simulation loops  ####
  ############################################
  # rm(x)
  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation


    ###########################################
    ### geting random sample ################
    x_c=cbind(rep(1,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_c))
    x_t=cbind(rep(2,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_t))
    x=rbind(x_c,x_t)
    #### interim analysis parameters ###
    x_c_inter=x_c[1:N[1],]
    x_t_inter=x_t[1:N[1],]

    x_inter=rbind(x_c_inter,x_t_inter)

    negglm_inter=glm.nb(as.numeric(x_inter[,2])~as.factor(x_inter[,1]),link="log")
    sum_inter=summary(negglm_inter)
    # s_e_inter=as.numeric(sum_inter$coefficients[2,2])
    # t_inter=s_e_inter^(-2)
    # theta_nbglm_inter=-as.numeric(sum_inter$coefficients[2,1])
    kappa_hat_inter=sum_inter$theta
    r_c_inter=mean(x_c_inter[,2])
    r_t_inter=mean(x_t_inter[,2])

    v_hat_inter=v_fun(r_0=r_c_inter,r_1=r_t_inter,theta=rand_ratio,nu_t,kappa=kappa_hat_inter)
    #  v_true=v_fun(r_0=r_c_f,r_1=r_t_f,theta=rand_ratio,nu_t,kappa=kappa)
    theta_hat_inter=-log(r_t_inter/r_c_inter)
    t_inter=N[1]/v_hat_inter
    ######## Planned final analyss ##########
    x_c_pf=x_c[1:NF,]
    x_t_pf=x_t[1:NF,]

    x_pf=rbind(x_c_inter,x_t_inter)

    negglm_pf=glm.nb(as.numeric(x_pf[,2])~as.factor(x_pf[,1]),link="log")
    sum_pf=summary(negglm_pf)
    # s_e_inter=as.numeric(sum_inter$coefficients[2,2])
    # t_inter=s_e_inter^(-2)
    # theta_nbglm_inter=-as.numeric(sum_inter$coefficients[2,1])
    kappa_hat_pf=sum_inter$theta
    r_c_pf=mean(x_c_pf[,2])
    r_t_pf=mean(x_t_pf[,2])
    #  rm(r_t_pf)
    v_hat_pf=v_fun(r_0=r_c_pf,r_1=r_t_pf,theta=rand_ratio,nu_t,kappa=kappa_hat_pf)
    #  v_true=v_fun(r_0=r_c_f,r_1=r_t_f,theta=rand_ratio,nu_t,kappa=kappa)
    theta_hat_pf=-log(r_t_pf/r_c_pf)
    t_pf=NF/v_hat_pf
    t_last= t_pf

    z_zhu_pf= theta_hat_pf*sqrt(t_pf)
    Wald_last=z_zhu_pf
    ##### negative binomial glm ####

    if( Wald_last>=c_bry){
      rej_Ztest=1
    }else{
      rej_Ztest=0
    }
    ###############################################
    ###############################################
    ############## generate the parameters ########
    ###############################################

    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]


    ###### Defining promising zone ########

    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)

    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      samsz_simu_pz=NF

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(N[1]*new_time_pz_0/t_inter)
      new_sampsz_1_pz=max(new_sampsz_0_pz, NF)
      ############
      samsz_simu_pz=min(c(N_max,new_sampsz_1_pz))
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################
      test_new_last_pz=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_pz)

      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }
    }
    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      samsz_simu_cp=NF
      samsz_simu_cp_bayesian=NF
      samsz_simu_true_cp=NF
      samsz_simu_pp_non_inf_prior=NF
      samsz_simu_pp_inf_prior=NF
      samsz_simu_mpp_inf_prior=NF
      samsz_simu_mpp_non_inf_prior=NF
      samsz_simu_ci_mpp_inf_prior=NF
      samsz_simu_ci_mpp_non_inf_prior=NF


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
      new_sampsz_1_cp=max(new_sampsz_0_cp,NF)
      #######
      samsz_simu_cp=min(c(N_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp)

      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_0_bayesian/t_inter)
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian,NF)
      ##############
      samsz_simu_cp_bayesian=min(c(N_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp_bayesian=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp_bayesian)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(N_inter*new_time_true_cp_0/t_inter)
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, NF)
      #################
      samsz_simu_true_cp=min(c(N_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      test_new_last_true_cp=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################
      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_pp_inf_prior==9999){

        samsz_simu_pp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, NF)
        ###############
        samsz_simu_pp_inf_prior=min(c(N_max,new_sampsz_1_pp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1
        }else{
          rej_pp_inf_prior=0
        }
      }
      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){

        samsz_simu_pp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,NF)
        #########################
        samsz_simu_pp_non_inf_prior=min(c(N_max,new_sampsz_1_pp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_non_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){

        samsz_simu_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,NF)
        ###################
        samsz_simu_mpp_inf_prior=min(c(N_max,new_sampsz_1_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1
        }else{
          rej_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){

        samsz_simu_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,NF)
        #################
        samsz_simu_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_non_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1
        }else{
          rej_mpp_non_inf_prior=0
        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)

      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){

        samsz_simu_ci_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #     new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1
        }else{
          rej_ci_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){

        samsz_simu_ci_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #      new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1
        }else{
          rej_ci_mpp_non_inf_prior=0
        }
      }
      #######################################
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest,t_last)
    #out

  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  m_t_last=mean(t_last)
  #T_last=NF/2/sigma^2
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(NF,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""
  out_2_4=data.frame(References4)
  colnames(out_2_4)=""
  out_2=rbind(out_2_1,out_2_2,out_2_3,out_2_4)


  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
}


#hybrid_fixed_prior_rej_negbinom_simu(t_prior,theta_prior,r_t,r_c,kappa,nu_t,rand_ratio,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num)
#hybrid_fixed_prior_rej_negbinom_simu(t_prior=5,theta_prior=0.2,r_t=1,r_c=1,kappa=1,nu_t=1,rand_ratio=1,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_fixed_prior_rej_negbinom_simu(t_prior=5,theta_prior=0.2,r_t=1.3,r_c=1,kappa=1,nu_t=1,rand_ratio=1,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)
#hybrid_fixed_prior_rej_negbinom_simu(t_prior=5,theta_prior=0.2,r_t=1,r_c=1.3,kappa=1,nu_t=1,rand_ratio=1,s=c(0.5,1),sampsz=150,N_max=550,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)


###############################################


##################### Random Prior ####

#################
#' Hybrid: Simulations, Negative binomial distribution, Using random informative prior from previous trial
#'
#' @param r_t_0 Event rate in the test arm of the previous trial, r_t0; the previous-trial
#'   effect size is theta_0 = -log(r_t_0/r_c_0)
#' @param r_c_0 Event rate in the control arm of the previous trial, r_c0
#' @param sampsz_0 Sample size per arm from the previous trial, N_0; random samples of this size are drawn from the previous trial to form the random prior
#' @param r_t Event rate in the test arm of the current trial; the effect size is
#'   theta = -log(r_t/r_c) and a smaller event rate is assumed to be associated with better
#'   efficacy. Enter r_t = r_c to verify type I error control and r_t < r_c for the power
#'   simulation
#' @param r_c Event rate in the control arm of the current trial; the effect size is theta = -log(r_t/r_c)
#' @param kappa Dispersion parameter of the negative binomial distribution; must be > 0
#' @param nu_t Exposure time for the negative binomial distribution; must be > 0
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz Planned sample size: the planned sample size at the final analysis
#' @param N_max Maximum sample size: a user chosen, maximum allowable sample size for sample size modification
#' @param theta_cut Futility threshold: a cut-off for theta_hat. In the simulations, if
#'   theta_hat at the interim analysis is at or below theta_cut the trial is considered futile
#'   and proceeds to the planned sample size without any change
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods}
#'   \item{futility rate}{The rate at which the observed theta_hat at the interim analysis is at or below theta_cut, so that the trial proceeds without any sample size change}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_random_prior_rej_negbinom_simu(r_t_0 = 1.2, r_c_0 = 1, sampsz_0 = 50, r_t = 1, r_c = 1,
#' kappa = 1, nu_t = 1, rand_ratio = 1, s = c(0.5,1), sampsz = 100, N_max = 150,
#' theta_cut = 0.001, cp_min = 0.27, gamma = 0.5, alpha_2 = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_random_prior_rej_negbinom_simu<-function(r_t_0,r_c_0,sampsz_0,r_t,r_c,kappa,nu_t,rand_ratio,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number
  #####################################################
  ## varaince function per Zhu, harris, 2012 ####
  v_fun=function(r_0,r_1,theta,nu_t,kappa){
    V_01=(1+theta)/(nu_t*theta*r_0)+(1+theta)*kappa/theta
    V_02=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
    V_03=(1+theta)**2/(nu_t*theta*(r_0+theta*r_1))+(1+theta)*kappa/theta
    V_max=max(V_01,V_02,V_02)
    return(V_max)
  }

  ######################
  N=ceiling(sampsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  NF=N[length(s)] # is the sample size at the final analysis
  N_inter=N[1]

  c_bry=stats::qnorm(1-alpha_2/2)
  mu=-log(r_t/r_c)
  ############################################
  ####  beginning of the simulation loops  ####
  ############################################
  # rm(x)
  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation


    ################# Random prior ####
    x_c_0=cbind(rep(1,sampsz_0),stats::rnbinom(n=sampsz_0, size=kappa, mu=r_c_0))
    x_t_0=cbind(rep(2,sampsz_0),stats::rnbinom(n=sampsz_0, size=kappa, mu=r_t_0))
    x_0=rbind(x_c_0,x_t_0)


    negglm_0=glm.nb(as.numeric(x_0[,2])~as.factor(x_0[,1]),link="log")
    sum_0=summary(negglm_0)
    # s_e_inter=as.numeric(sum_inter$coefficients[2,2])
    # t_inter=s_e_inter^(-2)
    # theta_nbglm_inter=-as.numeric(sum_inter$coefficients[2,1])
    kappa_hat_0=sum_0$theta
    r_c_0=mean(x_c_0[,2])
    r_t_0=mean(x_t_0[,2])

    v_hat_0=v_fun(r_0=r_c_0,r_1=r_t_0,theta=rand_ratio,nu_t,kappa=kappa_hat_0)
    #  v_true=v_fun(r_0=r_c_f,r_1=r_t_f,theta=rand_ratio,nu_t,kappa=kappa)
    theta_hat_0=-log(r_t_0/r_c_0)
    t_0=sampsz_0/v_hat_0

    ##########


    if ( theta_hat_0>0.05){
      t_prior=t_0
      theta_prior=theta_hat_0
    }else{
      t_prior=0
      theta_prior=0
    }

    ###########################################
    ### geting random sample ################
    x_c=cbind(rep(1,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_c))
    x_t=cbind(rep(2,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_t))
    x=rbind(x_c,x_t)
    #### interim analysis parameters ###
    x_c_inter=x_c[1:N[1],]
    x_t_inter=x_t[1:N[1],]

    x_inter=rbind(x_c_inter,x_t_inter)

    negglm_inter=glm.nb(as.numeric(x_inter[,2])~as.factor(x_inter[,1]),link="log")
    sum_inter=summary(negglm_inter)
    # s_e_inter=as.numeric(sum_inter$coefficients[2,2])
    # t_inter=s_e_inter^(-2)
    # theta_nbglm_inter=-as.numeric(sum_inter$coefficients[2,1])
    kappa_hat_inter=sum_inter$theta
    r_c_inter=mean(x_c_inter[,2])
    r_t_inter=mean(x_t_inter[,2])

    v_hat_inter=v_fun(r_0=r_c_inter,r_1=r_t_inter,theta=rand_ratio,nu_t,kappa=kappa_hat_inter)
    #  v_true=v_fun(r_0=r_c_f,r_1=r_t_f,theta=rand_ratio,nu_t,kappa=kappa)
    theta_hat_inter=-log(r_t_inter/r_c_inter)
    t_inter=N[1]/v_hat_inter
    ######## Planned final analyss ##########
    ######## Planned final analyss ##########
    x_c_pf=x_c[1:NF,]
    x_t_pf=x_t[1:NF,]

    x_pf=rbind(x_c_inter,x_t_inter)

    negglm_pf=glm.nb(as.numeric(x_pf[,2])~as.factor(x_pf[,1]),link="log")
    sum_pf=summary(negglm_pf)
    # s_e_inter=as.numeric(sum_inter$coefficients[2,2])
    # t_inter=s_e_inter^(-2)
    # theta_nbglm_inter=-as.numeric(sum_inter$coefficients[2,1])
    kappa_hat_pf=sum_inter$theta
    r_c_pf=mean(x_c_pf[,2])
    r_t_pf=mean(x_t_pf[,2])
    #  rm(r_t_pf)
    v_hat_pf=v_fun(r_0=r_c_pf,r_1=r_t_pf,theta=rand_ratio,nu_t,kappa=kappa_hat_pf)
    #  v_true=v_fun(r_0=r_c_f,r_1=r_t_f,theta=rand_ratio,nu_t,kappa=kappa)
    theta_hat_pf=-log(r_t_pf/r_c_pf)
    t_pf=NF/v_hat_pf
    t_last= t_pf

    z_zhu_pf= theta_hat_pf*sqrt(t_pf)
    Wald_last=z_zhu_pf
    ##### negative binomial glm ####

    if( Wald_last>=c_bry){
      rej_Ztest=1
    }else{
      rej_Ztest=0
    }
    ###############################################
    ###############################################
    ###############################################
    ############## generate the parameters ########
    ###############################################

    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]


    ###### Defining promising zone ########

    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)

    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      samsz_simu_pz=NF

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(N[1]*new_time_pz_0/t_inter)
      new_sampsz_1_pz=max(new_sampsz_0_pz, NF)
      ############
      samsz_simu_pz=min(c(N_max,new_sampsz_1_pz))
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################
      test_new_last_pz=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_pz)

      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }
    }
    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      samsz_simu_cp=NF
      samsz_simu_cp_bayesian=NF
      samsz_simu_true_cp=NF
      samsz_simu_pp_non_inf_prior=NF
      samsz_simu_pp_inf_prior=NF
      samsz_simu_mpp_inf_prior=NF
      samsz_simu_mpp_non_inf_prior=NF
      samsz_simu_ci_mpp_inf_prior=NF
      samsz_simu_ci_mpp_non_inf_prior=NF


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
      new_sampsz_1_cp=max(new_sampsz_0_cp,NF)
      #######
      samsz_simu_cp=min(c(N_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp)

      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_0_bayesian/t_inter)
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian,NF)
      ##############
      samsz_simu_cp_bayesian=min(c(N_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp_bayesian=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp_bayesian)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(N_inter*new_time_true_cp_0/t_inter)
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, NF)
      #################
      samsz_simu_true_cp=min(c(N_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      test_new_last_true_cp=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################
      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_pp_inf_prior==9999){

        samsz_simu_pp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, NF)
        ###############
        samsz_simu_pp_inf_prior=min(c(N_max,new_sampsz_1_pp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1
        }else{
          rej_pp_inf_prior=0
        }
      }
      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){

        samsz_simu_pp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,NF)
        #########################
        samsz_simu_pp_non_inf_prior=min(c(N_max,new_sampsz_1_pp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_non_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){

        samsz_simu_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,NF)
        ###################
        samsz_simu_mpp_inf_prior=min(c(N_max,new_sampsz_1_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1
        }else{
          rej_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){

        samsz_simu_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,NF)
        #################
        samsz_simu_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_non_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1
        }else{
          rej_mpp_non_inf_prior=0
        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)

      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){

        samsz_simu_ci_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #     new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1
        }else{
          rej_ci_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){

        samsz_simu_ci_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_negbinom(x_c,x_t,t_inter,score_inter,T_planned,t_last,nu_t,rand_ratio,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #      new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1
        }else{
          rej_ci_mpp_non_inf_prior=0
        }
      }
      #######################################
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest,t_last)
    #out

  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  m_t_last=mean(t_last)
  #T_last=NF/2/sigma^2
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(NF,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""
  out_2_4=data.frame(References4)
  colnames(out_2_4)=""
  out_2=rbind(out_2_1,out_2_2,out_2_3,out_2_4)


  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
}
#####################
#hybrid_random_prior_rej_negbinom_simu(r_t_0,r_c_0,sampsz_0,r_t,r_c,kappa,nu_t,rand_ratio,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num)
#hybrid_random_prior_rej_negbinom_simu(r_t_0=1.2,r_c_0=1,sampsz_0=50,r_t=1,r_c=1,kappa=1,nu_t=1,rand_ratio=1,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_random_prior_rej_negbinom_simu(r_t_0=1.2,r_c_0=1,sampsz_0=50,r_t=1.3,r_c=1,kappa=1,nu_t=1,rand_ratio=1,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)
#hybrid_random_prior_rej_negbinom_simu(r_t_0=1.2,r_c_0=1,sampsz_0=50,r_t=1,r_c=1.3,kappa=1,nu_t=1,rand_ratio=1,s=c(0.5,1),sampsz=150,N_max=550,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)


#############################
#### Binary distribution ####
#############################

##### Fixed Prior  #########

#' Hybrid: Simulations, Binary distribution, Using given estimates of informative prior
#'
#' @param t_prior Prior information time: the prior distribution of theta = p_t - p_c is assumed
#'   to have a normal distribution, and the Fisher's information time is the reciprocal of the
#'   variance of the prior distribution; set t_prior = 0 for a non-informative prior
#' @param theta_prior Mean effect size of prior: the mean of the normal prior distribution of
#'   theta = p_t - p_c; set theta_prior = 0 for a non-informative prior
#' @param p_t Response rate (probability of response) in the test arm; 0 < p_t < 1. The effect
#'   size is theta = p_t - p_c. Enter p_t = p_c to verify type I error control and p_t > p_c for
#'   the power simulation
#' @param p_c Response rate (probability of response) in the control arm; 0 < p_c < 1. The
#'   effect size is theta = p_t - p_c
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz Planned sample size: the planned sample size at the final analysis
#' @param N_max Maximum sample size: a user chosen, maximum allowable sample size for sample size modification
#' @param theta_cut Futility threshold: a cut-off for theta_hat. In the simulations, if
#'   theta_hat at the interim analysis is at or below theta_cut the trial is considered futile
#'   and proceeds to the planned sample size without any change
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods}
#'   \item{futility rate}{The rate at which the observed theta_hat at the interim analysis is at or below theta_cut, so that the trial proceeds without any sample size change}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_fixed_prior_rej_binary_simu(t_prior = 5, theta_prior = 0.2, p_t = 0.3, p_c = 0.3,
#' s = c(0.5,1), sampsz = 100, N_max = 150, theta_cut = 0.001, cp_min = 0.27, gamma = 0.5,
#' alpha_2 = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_fixed_prior_rej_binary_simu<-function(t_prior,theta_prior,p_t,p_c,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number
  #####################################################

  N=ceiling(sampsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  N_inter=N[1]
  T_true<- N/(p_t*(1-p_t)+p_c*(1-p_c))  # t is the vector of information times at each interim anbalysis
  NF=sampsz  ## NF is the sample size at the end of the trial
  #t is the vector of information times at the interim analysis
  mu=p_t-p_c

  c_bry=stats::qnorm(1-alpha_2/2)
  ############################################
  ####  beginning of the simulation loops  ####
  ############################################
  # rm(x)
  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation


    ###########################################
    ### geting random sample ################

    x0=stats::rbinom(N_max,1,p_c)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rbinom(N_max,1,p_t)

    ###############################################
    ############## generate the parameters ########
    ###############################################
    theta_hat_inter=mean(x[1:N[1]])-mean(x0[1:N[1]])
    theta_hat_last=mean(x[1:NF])-mean(x0[1:NF])
    t_inter=N[1]/(stats::var(x[1:N[1]])+stats::var(x0[1:N[1]]))
    t_last=NF/(stats::var(x[1:NF])+stats::var(x0[1:NF]))

    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]
    Wald_inter=theta_hat_inter*sqrt(t_inter)
    Wald_last= theta_hat_last*sqrt(t_last)

    if(Wald_last>=c_bry){
      rej_Ztest=1
    }else{
      rej_Ztest=0
    }
    ###### Defining promising zone ########

    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)

    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      samsz_simu_pz=NF

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(N[1]*new_time_pz_0/t_inter)
      new_sampsz_1_pz=max(new_sampsz_0_pz, NF)
      ############
      samsz_simu_pz=min(c(N_max,new_sampsz_1_pz))
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################
      test_new_last_pz=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pz)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }
    }
    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      samsz_simu_cp=NF
      samsz_simu_cp_bayesian=NF
      samsz_simu_true_cp=NF
      samsz_simu_pp_non_inf_prior=NF
      samsz_simu_pp_inf_prior=NF
      samsz_simu_mpp_inf_prior=NF
      samsz_simu_mpp_non_inf_prior=NF
      samsz_simu_ci_mpp_inf_prior=NF
      samsz_simu_ci_mpp_non_inf_prior=NF


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
      new_sampsz_1_cp=max(new_sampsz_0_cp,NF)
      #######
      samsz_simu_cp=min(c(N_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_0_bayesian/t_inter)
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian,NF)
      ##############
      samsz_simu_cp_bayesian=min(c(N_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp_bayesian=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp_bayesian)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(N_inter*new_time_true_cp_0/t_inter)
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, NF)
      #################
      samsz_simu_true_cp=min(c(N_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      test_new_last_true_cp=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################
      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_pp_inf_prior==9999){

        samsz_simu_pp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, NF)
        ###############
        samsz_simu_pp_inf_prior=min(c(N_max,new_sampsz_1_pp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1
        }else{
          rej_pp_inf_prior=0
        }
      }
      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){

        samsz_simu_pp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,NF)
        #########################
        samsz_simu_pp_non_inf_prior=min(c(N_max,new_sampsz_1_pp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_non_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){

        samsz_simu_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,NF)
        ###################
        samsz_simu_mpp_inf_prior=min(c(N_max,new_sampsz_1_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1
        }else{
          rej_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){

        samsz_simu_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,NF)
        #################
        samsz_simu_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_non_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1
        }else{
          rej_mpp_non_inf_prior=0
        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)

      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){

        samsz_simu_ci_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1
        }else{
          rej_ci_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){

        samsz_simu_ci_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1
        }else{
          rej_ci_mpp_non_inf_prior=0
        }
      }
      #######################################
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest,t_last)
    #out

  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  m_t_last=mean(t_last)
  T_last=NF/(p_t*(1-p_t)+p_c*(1-p_c))
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(NF,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""

  out_2=rbind(out_2_1,out_2_2,out_2_3)

  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
}
######################################################################
#hybrid_fixed_prior_rej_binary_simu(t_prior,theta_prior,p_t,p_c,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num)
#################
#hybrid_fixed_prior_rej_binary_simu(t_prior=5,theta_prior=0.2,p_t=0.3,p_c=0.3,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_fixed_prior_rej_binary_simu(t_prior=5,theta_prior=0.2,p_t=0.5,p_c=0.3,s=c(0.5,1),sampsz=80,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_fixed_prior_rej_binary_simu(t_prior=5,theta_prior=0.2,p_t=0.45,p_c=0.3,s=c(0.5,1),sampsz=80,N_max=250,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)

##################### Random Prior ####
#' Hybrid: Simulations, Binary distribution, Using random informative prior from previous trial
#'
#' @param p_t_0 Response rate (probability of response) in the test arm of the previous trial,
#'   p_t0; the previous-trial effect size is theta_0 = p_t_0 - p_c_0
#' @param p_c_0 Response rate (probability of response) in the control arm of the previous trial
#' @param sampsz_0 Sample size per arm from the previous trial, N_0; random samples of this size are drawn from the previous trial to form the random prior
#' @param p_t Response rate (probability of response) in the test arm of the current trial;
#'   0 < p_t < 1. The effect size is theta = p_t - p_c. Enter p_t = p_c to verify type I error
#'   control and p_t > p_c for the power simulation
#' @param p_c Response rate (probability of response) in the control arm of the current trial;
#'   0 < p_c < 1. The effect size is theta = p_t - p_c
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz Planned sample size: the planned sample size at the final analysis
#' @param N_max Maximum sample size: a user chosen, maximum allowable sample size for sample size modification
#' @param theta_cut Futility threshold: a cut-off for theta_hat. In the simulations, if
#'   theta_hat at the interim analysis is at or below theta_cut the trial is considered futile
#'   and proceeds to the planned sample size without any change
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods}
#'   \item{futility rate}{The rate at which the observed theta_hat at the interim analysis is at or below theta_cut, so that the trial proceeds without any sample size change}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_random_prior_rej_binary_simu(p_t_0 = 0.4, p_c_0 = 0.3, sampsz_0 = 60, p_t = 0.3,
#' p_c = 0.3, s = c(0.5,1), sampsz = 100, N_max = 150, theta_cut = 0.001, cp_min = 0.27,
#' gamma = 0.5, alpha_2 = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_random_prior_rej_binary_simu<-function(p_t_0,p_c_0,sampsz_0,p_t,p_c,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number
  #####################################################

  N=ceiling(sampsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  N_inter=N[1]
  T_true<- N/(p_t*(1-p_t)+p_c*(1-p_c))  # t is the vector of information times at each interim anbalysis
  NF=sampsz  ## NF is the sample size at the end of the trial
  #t is the vector of information times at the interim analysis
  mu=p_t-p_c

  c_bry=stats::qnorm(1-alpha_2/2)
  ############################################
  ####  beginning of the simulation loops  ####
  ############################################
  # rm(x)
  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation

    ###########################################
    ### geting random sample for old trial ################

    x0_0=stats::rbinom(sampsz_0,1,p_c_0)
    #x0 is the random sample at the first interim analysis for the control
    x_0=stats::rbinom(sampsz_0,1,p_t_0)
    ###############################################
    ############## generate the parameters ########
    ###############################################

    #  inter_vt=length(s)-1
    theta_hat_0=mean(x_0)-mean(x0_0)
    t_0=sampsz_0/(stats::var(x_0)+stats::var(x0_0))

    if ( theta_hat_0>0.01){
      t_prior=t_0
      theta_prior=theta_hat_0
    }else{
      t_prior=0
      theta_prior=0
    }

    ###########################################
    ###########################################
    ### geting random sample for new trial ################

    x0=stats::rbinom(N_max,1,p_c)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rbinom(N_max,1,p_t)

    ###############################################
    ############## generate the parameters ########
    ###############################################
    theta_hat_inter=mean(x[1:N[1]])-mean(x0[1:N[1]])
    theta_hat_last=mean(x[1:NF])-mean(x0[1:NF])
    t_inter=N[1]/(stats::var(x[1:N[1]])+stats::var(x0[1:N[1]]))
    t_last=NF/(stats::var(x[1:NF])+stats::var(x0[1:NF]))

    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]
    Wald_inter=theta_hat_inter*sqrt(t_inter)
    Wald_last= theta_hat_last*sqrt(t_last)

    if(Wald_last>=c_bry){
      rej_Ztest=1
    }else{
      rej_Ztest=0
    }
    ###### Defining promising zone ########

    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)

    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      samsz_simu_pz=NF

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(N[1]*new_time_pz_0/t_inter)
      new_sampsz_1_pz=max(new_sampsz_0_pz, NF)
      ############
      samsz_simu_pz=min(c(N_max,new_sampsz_1_pz))
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################
      test_new_last_pz=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pz)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }
    }
    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      samsz_simu_cp=NF
      samsz_simu_cp_bayesian=NF
      samsz_simu_true_cp=NF
      samsz_simu_pp_non_inf_prior=NF
      samsz_simu_pp_inf_prior=NF
      samsz_simu_mpp_inf_prior=NF
      samsz_simu_mpp_non_inf_prior=NF
      samsz_simu_ci_mpp_inf_prior=NF
      samsz_simu_ci_mpp_non_inf_prior=NF


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
      new_sampsz_1_cp=max(new_sampsz_0_cp,NF)
      #######
      samsz_simu_cp=min(c(N_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_0_bayesian/t_inter)
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian,NF)
      ##############
      samsz_simu_cp_bayesian=min(c(N_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp_bayesian=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp_bayesian)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(N_inter*new_time_true_cp_0/t_inter)
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, NF)
      #################
      samsz_simu_true_cp=min(c(N_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      test_new_last_true_cp=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################
      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_pp_inf_prior==9999){

        samsz_simu_pp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, NF)
        ###############
        samsz_simu_pp_inf_prior=min(c(N_max,new_sampsz_1_pp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1
        }else{
          rej_pp_inf_prior=0
        }
      }
      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){

        samsz_simu_pp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,NF)
        #########################
        samsz_simu_pp_non_inf_prior=min(c(N_max,new_sampsz_1_pp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_non_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){

        samsz_simu_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,NF)
        ###################
        samsz_simu_mpp_inf_prior=min(c(N_max,new_sampsz_1_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1
        }else{
          rej_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){

        samsz_simu_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,NF)
        #################
        samsz_simu_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_non_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1
        }else{
          rej_mpp_non_inf_prior=0
        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)

      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){

        samsz_simu_ci_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1
        }else{
          rej_ci_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){

        samsz_simu_ci_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_binary(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1
        }else{
          rej_ci_mpp_non_inf_prior=0
        }
      }
      #######################################
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest,t_last)
    #out

  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  m_t_last=mean(t_last)
  T_last=NF/(p_t*(1-p_t)+p_c*(1-p_c))
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########

  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(NF,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""

  out_2=rbind(out_2_1,out_2_2,out_2_3)


  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
}
######################################################################
#hybrid_random_prior_rej_binary_simu(p_t_0,p_c_0,sampsz_0,p_t,p_c,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num)
#################
#hybrid_random_prior_rej_binary_simu(p_t_0=0.4,p_c_0=0.3,sampsz_0=60,p_t=0.3,p_c=0.3,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_random_prior_random_prior_rej_binary_simu(p_t_0=0.5,p_c_0=0.3,sampsz_0=60,p_t=0.5,p_c=0.3,s=c(0.5,1),sampsz=80,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_random_prior_rej_binary_simu(p_t_0=0.45,p_c_0=0.3,sampsz_0=60,p_t=0.45,p_c=0.3,s=c(0.5,1),sampsz=80,N_max=250,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)


################################
#### Normal distribution #######
################################

####### Fixed  prior ######
########################
#' Hybrid: Simulations, Normal distribution, Using given estimates of informative prior
#'
#' @param t_prior Prior information time: the prior distribution of theta = mu_t - mu_c is
#'   assumed to have a normal distribution, and the Fisher's information time is the reciprocal
#'   of the variance of the prior distribution; set t_prior = 0 for a non-informative prior
#' @param theta_prior Mean effect size of prior: the mean of the normal prior distribution of
#'   theta = mu_t - mu_c; set theta_prior = 0 for a non-informative prior
#' @param mu Effect size (difference in means) of the current trial: theta = mu_t - mu_c, where
#'   mu_t is the mean of the treatment arm and mu_c is the mean of the control arm. In the
#'   simulations the test arm is generated from a normal distribution with mean mu and the
#'   control arm from a normal distribution with mean 0, both with standard deviation sigma.
#'   Enter mu = 0 to verify type I error control and mu > 0 for the power simulation
#' @param sigma Common standard deviation shared by the test and control arms; must be > 0
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz Planned sample size: the planned sample size at the final analysis
#' @param N_max Maximum sample size: a user chosen, maximum allowable sample size for sample size modification
#' @param theta_cut Futility threshold: a cut-off for theta_hat. In the simulations, if
#'   theta_hat at the interim analysis is at or below theta_cut the trial is considered futile
#'   and proceeds to the planned sample size without any change
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods}
#'   \item{futility rate}{The rate at which the observed theta_hat at the interim analysis is at or below theta_cut, so that the trial proceeds without any sample size change}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_fixed_prior_rej_normal_simu(t_prior = 5, theta_prior = 0.2, mu = 0, sigma = 1,
#' s = c(0.5,1), sampsz = 100, N_max = 150, theta_cut = 0.001, cp_min = 0.27, gamma = 0.5,
#' alpha_2 = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_fixed_prior_rej_normal_simu<-function(t_prior,theta_prior,mu,sigma,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number
  #####################################################

  N=ceiling(sampsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  NF=N[length(s)] # is the sample size at the final analysis
  N_inter=N[1]

  c_bry=stats::qnorm(1-alpha_2/2)
  ############################################
  ####  beginning of the simulation loops  ####
  ############################################
  # rm(x)
  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation


    ###########################################
    ### geting random sample ################

    x0=stats::rnorm(N_max,0,sigma)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rnorm(N_max,mu,sigma)

    ###############################################
    ############## generate the parameters ########
    ###############################################

    #  inter_vt=length(s)-1
    theta_hat_inter=mean(x[1:N[1]])-mean(x0[1:N[1]])
    theta_hat_last=mean(x[1:NF])-mean(x0[1:NF])
    t_inter=N[1]/(stats::var(x[1:N[1]])+stats::var(x0[1:N[1]]))
    t_last=NF/(stats::var(x[1:NF])+stats::var(x0[1:NF]))

    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]
    Wald_inter=theta_hat_inter*sqrt(t_inter)
    Wald_last= theta_hat_last*sqrt(t_last)

    if(Wald_last>=c_bry){
      rej_Ztest=1
    }else{
      rej_Ztest=0
    }
    ###### Defining promising zone ########

    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)

    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      samsz_simu_pz=NF

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(N[1]*new_time_pz_0/t_inter)
      new_sampsz_1_pz=max(new_sampsz_0_pz, NF)
      ############
      samsz_simu_pz=min(c(N_max,new_sampsz_1_pz))
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################
      test_new_last_pz=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pz)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }
    }
    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      samsz_simu_cp=NF
      samsz_simu_cp_bayesian=NF
      samsz_simu_true_cp=NF
      samsz_simu_pp_non_inf_prior=NF
      samsz_simu_pp_inf_prior=NF
      samsz_simu_mpp_inf_prior=NF
      samsz_simu_mpp_non_inf_prior=NF
      samsz_simu_ci_mpp_inf_prior=NF
      samsz_simu_ci_mpp_non_inf_prior=NF


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
      new_sampsz_1_cp=max(new_sampsz_0_cp,NF)
      #######
      samsz_simu_cp=min(c(N_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_0_bayesian/t_inter)
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian,NF)
      ##############
      samsz_simu_cp_bayesian=min(c(N_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp_bayesian=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp_bayesian)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(N_inter*new_time_true_cp_0/t_inter)
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, NF)
      #################
      samsz_simu_true_cp=min(c(N_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      test_new_last_true_cp=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################
      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_pp_inf_prior==9999){

        samsz_simu_pp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, NF)
        ###############
        samsz_simu_pp_inf_prior=min(c(N_max,new_sampsz_1_pp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1
        }else{
          rej_pp_inf_prior=0
        }
      }
      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){

        samsz_simu_pp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,NF)
        #########################
        samsz_simu_pp_non_inf_prior=min(c(N_max,new_sampsz_1_pp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_non_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){

        samsz_simu_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,NF)
        ###################
        samsz_simu_mpp_inf_prior=min(c(N_max,new_sampsz_1_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1
        }else{
          rej_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){

        samsz_simu_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,NF)
        #################
        samsz_simu_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_non_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1
        }else{
          rej_mpp_non_inf_prior=0
        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)

      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){

        samsz_simu_ci_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1
        }else{
          rej_ci_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){

        samsz_simu_ci_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1
        }else{
          rej_ci_mpp_non_inf_prior=0
        }
      }
      #######################################
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest,t_last)
    #out

  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  m_t_last=mean(t_last)
  T_last=NF/2/sigma^2
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(NF,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)

  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""

  out_2=rbind(out_2_1,out_2_2,out_2_3)


  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
}
#####################

#hybrid_fixed_prior_rej_normal_simu(t_prior,theta_prior,mu,sigma,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num)
#hybrid_fixed_prior_rej_normal_simu(t_prior=5,theta_prior=0.2,mu=0,sigma=1,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)
#hybrid_fixed_prior_rej_normal_simu(t_prior=5,theta_prior=0.2,mu=0.3,sigma=1,s=c(0.5,1),sampsz=100,N_max=250,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)
###############################################

####### Random prior ######

#' Hybrid: Simulations, Normal distribution, Using random informative prior from previous trial
#'
#' @param mu_0 Effect size (difference in means) of the previous trial: theta_0 = mu_t0 - mu_c0;
#'   in the simulations the previous-trial test arm is generated from a normal distribution with
#'   mean mu_0 and the previous-trial control arm from a normal distribution with mean 0, both
#'   with standard deviation sigma_0
#' @param sigma_0 Common standard deviation of the test and control arms in the previous trial;
#'   must be > 0
#' @param sampsz_0 Sample size per arm from the previous trial, N_0; random samples of this size are drawn from the previous trial to form the random prior
#' @param mu Effect size (difference in means) of the current trial: theta = mu_t - mu_c, where
#'   mu_t is the mean of the treatment arm and mu_c is the mean of the control arm. In the
#'   simulations the test arm is generated from a normal distribution with mean mu and the
#'   control arm from a normal distribution with mean 0, both with standard deviation sigma.
#'   Enter mu = 0 to verify type I error control and mu > 0 for the power simulation
#' @param sigma Common standard deviation shared by the test and control arms; must be > 0
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz Planned sample size: the planned sample size at the final analysis
#' @param N_max Maximum sample size: a user chosen, maximum allowable sample size for sample size modification
#' @param theta_cut Futility threshold: a cut-off for theta_hat. In the simulations, if
#'   theta_hat at the interim analysis is at or below theta_cut the trial is considered futile
#'   and proceeds to the planned sample size without any change
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods}
#'   \item{futility rate}{The rate at which the observed theta_hat at the interim analysis is at or below theta_cut, so that the trial proceeds without any sample size change}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_random_prior_rej_normal_simu(mu_0 = 0.3, sigma_0 = 1, sampsz_0 = 60, mu = 0, sigma = 1,
#' s = c(0.5,1), sampsz = 100, N_max = 250, theta_cut = 0.01, cp_min = 0.27, gamma = 0.5,
#' alpha_2 = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_random_prior_rej_normal_simu<-function(mu_0,sigma_0,sampsz_0, mu,sigma,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number
  #####################################################

  N=ceiling(sampsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  NF=N[length(s)] # is the sample size at the final analysis
  N_inter=N[1]

  c_bry=stats::qnorm(1-alpha_2/2)
  ############################################
  ####  beginning of the simulation loops  ####
  ############################################
  # rm(x)
  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation

    ###########################################
    ### geting random sample for old trial ################

    x0_0=stats::rnorm(sampsz_0,0,sigma_0)
    #x0 is the random sample at the first interim analysis for the control
    x_0=stats::rnorm(sampsz_0,mu_0,sigma_0)

    ###############################################
    ############## generate the parameters ########
    ###############################################

    #  inter_vt=length(s)-1
    theta_hat_0=mean(x_0)-mean(x0_0)
    t_0=sampsz_0/(stats::var(x_0)+stats::var(x0_0))

    if ( theta_hat_0>0.05){
      t_prior=t_0
      theta_prior=theta_hat_0
    }else{
      t_prior=0
      theta_prior=0
    }

    ###########################################
    ### geting random sample for new trial ################

    x0=stats::rnorm(N_max,0,sigma)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rnorm(N_max,mu,sigma)

    ###############################################
    ############## generate the parameters ########
    ###############################################

    #  inter_vt=length(s)-1
    theta_hat_inter=mean(x[1:N[1]])-mean(x0[1:N[1]])
    theta_hat_last=mean(x[1:NF])-mean(x0[1:NF])
    t_inter=N[1]/(stats::var(x[1:N[1]])+stats::var(x0[1:N[1]]))
    t_last=NF/(stats::var(x[1:NF])+stats::var(x0[1:NF]))

    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]
    Wald_inter=theta_hat_inter*sqrt(t_inter)
    Wald_last= theta_hat_last*sqrt(t_last)

    if(Wald_last>=c_bry){
      rej_Ztest=1
    }else{
      rej_Ztest=0
    }
    ###### Defining promising zone ########

    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)

    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      samsz_simu_pz=NF

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(N[1]*new_time_pz_0/t_inter)
      new_sampsz_1_pz=max(new_sampsz_0_pz, NF)
      ############
      samsz_simu_pz=min(c(N_max,new_sampsz_1_pz))
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################
      test_new_last_pz=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pz)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }
    }
    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      samsz_simu_cp=NF
      samsz_simu_cp_bayesian=NF
      samsz_simu_true_cp=NF
      samsz_simu_pp_non_inf_prior=NF
      samsz_simu_pp_inf_prior=NF
      samsz_simu_mpp_inf_prior=NF
      samsz_simu_mpp_non_inf_prior=NF
      samsz_simu_ci_mpp_inf_prior=NF
      samsz_simu_ci_mpp_non_inf_prior=NF


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
      new_sampsz_1_cp=max(new_sampsz_0_cp,NF)
      #######
      samsz_simu_cp=min(c(N_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_0_bayesian/t_inter)
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian,NF)
      ##############
      samsz_simu_cp_bayesian=min(c(N_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp_bayesian=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp_bayesian)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(N_inter*new_time_true_cp_0/t_inter)
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, NF)
      #################
      samsz_simu_true_cp=min(c(N_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      test_new_last_true_cp=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################
      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_pp_inf_prior==9999){

        samsz_simu_pp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, NF)
        ###############
        samsz_simu_pp_inf_prior=min(c(N_max,new_sampsz_1_pp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1
        }else{
          rej_pp_inf_prior=0
        }
      }
      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){

        samsz_simu_pp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,NF)
        #########################
        samsz_simu_pp_non_inf_prior=min(c(N_max,new_sampsz_1_pp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_non_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){

        samsz_simu_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,NF)
        ###################
        samsz_simu_mpp_inf_prior=min(c(N_max,new_sampsz_1_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1
        }else{
          rej_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){

        samsz_simu_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,NF)
        #################
        samsz_simu_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_non_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1
        }else{
          rej_mpp_non_inf_prior=0
        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)

      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){

        samsz_simu_ci_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1
        }else{
          rej_ci_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){

        samsz_simu_ci_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_normal(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1
        }else{
          rej_ci_mpp_non_inf_prior=0
        }
      }
      #######################################
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest,t_last)
    #out

  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  m_t_last=mean(t_last)
  T_last=NF/2/sigma^2
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(NF,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""

  out_2=rbind(out_2_1,out_2_2,out_2_3)


  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
}

####################################################################################
#hybrid_random_prior_rej_normal_simu(mu_0,sigma_0,sampsz_0, mu,sigma,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num)
#hybrid_random_prior_rej_normal_simu(mu_0=0.3,sigma_0=1,sampsz_0=60,mu=0,sigma=1,s=c(0.5,1),sampsz=100,N_max=250,theta_cut=0.01,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_random_prior_rej_normal_simu(mu_0=0.3,sigma_0=1,sampsz_0=60,mu=0.3,sigma=1,s=c(0.5,1),sampsz=100,N_max=250,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)
#########################
#### Survival analysis ##
#########################

####### Fixed  prior ######
###################
#' Hybrid: Simulations, Survival analysis, Using given estimates of informative prior
#'
#' @param t_prior Prior information time: the prior distribution of theta = -log(HR) is assumed
#'   to have a normal distribution, and the Fisher's information time is the reciprocal of the
#'   variance of the prior distribution; set t_prior = 0 for a non-informative prior
#' @param theta_prior Mean effect size of prior: the mean of the normal prior distribution of
#'   theta = -log(HR); set theta_prior = 0 for a non-informative prior
#' @param HR Hazards ratio of test arm vs. the control arm. The effect size is theta = -log(HR)
#'   and a smaller HR indicates better efficacy. Enter HR = 1 to verify type I error control and
#'   HR < 1 for the power simulation
#' @param lambda Baseline hazard rate of the control arm for the exponential survival times: the
#'   control arm is generated as rexp(lambda) and the test arm as rexp(HR * lambda), so lambda
#'   is a rate (not a log-hazard) and must be > 0. It also sets the censoring time
#'   T = -log(surv_rate)/lambda
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param surv_rate Survival rate used to translate the number of events into a sample size; it
#'   is the survival probability at the (administrative) censoring time T = -log(surv_rate)/lambda, 0 < surv_rate < 1
#' @param evt_num Planned number of events: the planned number of events at the final analysis
#' @param evt_max Maximum number of events: a user chosen, maximum allowable number of events for sample size modification
#' @param HR_cut Futility hazards ratio threshold: a cut-off for HR_hat. In the simulations, if
#'   the estimated hazards ratio at the interim analysis is at or above HR_cut the trial is
#'   considered futile and proceeds to the planned number of events without any change; it is
#'   converted internally to theta_cut = -log(HR_cut)
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param samp_multip Multiplication factor applied when translating the maximum number of
#'   events into the sample size used in the simulations,
#'   sampsz = ceiling(samp_multip * evt_max / surv_rate); it must be large enough for the
#'   required number of events to be observed
#' @param sim_num Number of simulation replicates
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods. For this survival design the mean sample size is the mean number of events}
#'   \item{futility rate}{The rate at which the estimated hazards ratio at the interim analysis is at or above HR_cut, so that the trial proceeds without any change to the number of events}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_fixed_prior_rej_survival_simu(t_prior = 10, theta_prior = 0.8, HR = 1, lambda = 1,
#' s = c(0.5,1), alpha_2 = 0.05, surv_rate = 0.5, evt_num = 88, evt_max = 331, HR_cut = 0.95,
#' cp_min = 0.27, gamma = 0.5, samp_multip = 5, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_fixed_prior_rej_survival_simu<-function(t_prior,theta_prior,HR,lambda,s,alpha_2,surv_rate,evt_num,evt_max,HR_cut,cp_min, gamma,samp_multip,sim_num,beta){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number


  #####################################################

  evt_num_s=c(ceiling(evt_num*s),evt_max)
  sampsz=ceiling(samp_multip*evt_max/surv_rate)
  theta_cut=-log(HR_cut)
  mu=-log(HR)

  c_bry=stats::qnorm(1-alpha_2/2)


  ############################################
  ####  beginning of the simulation loops  ####
  ############################################

  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation


    ###########################################
    ### geting random sample ################

    x0=stats::rexp(sampsz,lambda) #generate exponential data for the control arm
    x=stats::rexp(sampsz,HR*lambda)  #generate exponential data for the test arm
    T=-log(surv_rate)/lambda  #Censoring time. Events after T will be censored
    time_x0=x0
    time_x=x
    x0_status=c(rep(1,length(x0)))
    x_status=c(rep(1,length(x)))

    ######## creating two long collum of censored data,  ###########
    for (i in (1:sampsz)){
      if (x0[i]>T){
        time_x0[i]=T  #x0[i] is censored after T
        x0_status[i]=0   #censoring indicator
      }
    }
    grp_x0=c(rep(1,length(time_x0)))#vector of zeros, length of the vector is the same as the original smaple size
    grp_x=c(rep(0,length(time_x))) #
    for (i in (1:sampsz)){
      if (x[i]>T){
        time_x[i]=T #x[i] is censored after T
        x_status[i]=0  #censoring indicator
      }
    }

    #### count the numnber of events in the
    evt_count=c(rep(0,sampsz))
    for (i in 1:sampsz){
      evt_count[i]=sum(x0_status[1:i])+sum(x_status[1:i]) #number of events after 2*i observations
    }

    #################################################################################################
    theta_hat=rep(0,length(evt_num_s))
    theta_hat_sd=rep(0,length(evt_num_s))
    wald=rep(0,length(evt_num_s))
    score=rep(0,length(evt_num_s))
    t=rep(0,length(evt_num_s))
    for (i in (1:length(evt_num_s))) {
      evt_pos=evt_count[evt_count>evt_num_s[i]]#substring of  evt_count that are larger than evt_num
      evt_cut=match(evt_pos[1],evt_count) #the first incidence that evt_count is larger than evt_num

      #  vector of zeros, length of the vector is the same as the original smaple size
      time=c(time_x0[1:evt_cut],time_x[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
      event=c(x0_status[1:evt_cut],x_status[1:evt_cut]) #censoring status before total number of events exceed evt_num
      grp=c(grp_x0[1:evt_cut],grp_x[1:evt_cut]) #vector of zero, with length=2*evt_cut


      test= list(time=time,status=event,grp)#survival data
      coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
      # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
      # test is the same as data=test
      theta_hat[i]= coxcf$coefficients[1] #estimated log(HR)
      theta_hat_sd[i]= coxcf$coefficients[3] #estimated standard error of log(HR)
      t[i]=(theta_hat_sd[i])^(-2)
      wald[i]=theta_hat[i]/ theta_hat_sd[i]
      score[i]=theta_hat[i]*t[i]
    }

    ###############################################
    ############## generate the parameters ########
    ###############################################

    inter_vt=length(s)-1
    theta_hat_inter=theta_hat[inter_vt]
    theta_hat_inter_sd=theta_hat_sd[inter_vt]
    event_inter=evt_num_s[inter_vt]

    t_inter=t[1]
    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]
    Wald_inter=wald[1]
    ###### Planned sample size test ####
    last_vt=length(s)
    theta_hat_last=theta_hat[last_vt]
    theta_hat_last_sd=theta_hat_sd[last_vt]
    t_last= theta_hat_last_sd^(-2)

    Wald_last= theta_hat_last/theta_hat_last_sd
    if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
      rej_Ztest=1

    }else { # the NULL hypothesis is NOT rejected at the final analysis
      rej_Ztest=0

    }
    ###### Defining promising zone ########
    #      cp_inter=cond_power(T=T_planned,cry=c_bry,t_inter,s_inter,theta_hat_inter,Wald_inter)
    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)
    #    cond_power_inter=asd_c_power_theta(theta=theta_hat_inter,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)
    #    cp_min=0.27
    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      #stg identifies the interim analysis at which the null hypothesis is rejected
      last_vt=length(s)
      theta_hat_last=theta_hat[last_vt]
      theta_hat_last_sd=theta_hat_sd[last_vt]
      t_last= theta_hat_last_sd^(-2)

      Wald_last= theta_hat_last/theta_hat_last_sd


      samsz_simu_pz=evt_num_s[length(s)]

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(event_inter*new_time_pz_0/t[inter_vt])
      new_sampsz_1_pz=max(new_sampsz_0_pz, ceiling(event_inter*1.1))
      new_sampsz_pz=min(c(evt_max,new_sampsz_1_pz))

      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      samsz_simu_pz=max(new_sampsz_pz,evt_num)

      ##### Generate addtional random samples #######
      inter_vt=length(s)-1
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################

      #     test_new_last_pz=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,T_planned,event_inter,t_inter,new_samsz=samsz_simu_pz)
      #     test_new_last_pz=test_new_last_pz=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_pz)
      test_new_last_pz=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_pz)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }

    }


    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    # ## new information estimated using theta_hat_inter  ####



    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      #stg identifies the interim analysis at which the null hypothesis is rejected
      last_vt=length(s)
      theta_hat_last=theta_hat[last_vt]
      theta_hat_last_sd=theta_hat_sd[last_vt]
      t_last= theta_hat_last_sd^(-2)

      Wald_last= theta_hat_last/theta_hat_last_sd

      samsz_simu_cp=evt_num_s[length(s)]
      samsz_simu_cp_bayesian=evt_num_s[length(s)]
      samsz_simu_true_cp=evt_num_s[length(s)]
      samsz_simu_pp_non_inf_prior=evt_num_s[length(s)]
      samsz_simu_pp_inf_prior=evt_num_s[length(s)]
      samsz_simu_mpp_inf_prior=evt_num_s[length(s)]
      samsz_simu_mpp_non_inf_prior=evt_num_s[length(s)]
      samsz_simu_ci_mpp_inf_prior=evt_num_s[length(s)]
      samsz_simu_ci_mpp_non_inf_prior=evt_num_s[length(s)]


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(event_inter*new_time_cp_0/t[inter_vt])
      new_sampsz_1_cp=max(new_sampsz_0_cp, ceiling(event_inter*1.1))
      new_sampsz_cp=min(c(evt_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      samsz_simu_cp=max(new_sampsz_cp,evt_num)

      ##### Generate addtional random samples #######
      inter_vt=length(s)-1
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################

      #  test_new_last_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_cp)

      test_new_last_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_cp)

      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(event_inter*new_time_cp_0_bayesian/t[inter_vt])
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian, evt_num)
      samsz_simu_cp_bayesian=min(c(evt_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      #  samsz_simu_cp_bayesian=max(new_sampsz_cp_bayesian,evt_num)

      ##### Generate addtional random samples #######
      inter_vt=length(s)-1
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################

      #  test_new_last_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_cp)

      test_new_last_cp_bayesian=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_cp_bayesian)

      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(event_inter*new_time_true_cp_0/t[inter_vt])
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, ceiling(event_inter*1.1))
      new_sampsz_true_cp=min(c(evt_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      samsz_simu_true_cp=new_sampsz_true_cp

      ##### Generate addtional random samples #######
      inter_vt=length(s)-1
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################

      #    test_new_last_true_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_true_cp)
      test_new_last_true_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################

      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      #                        new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)

      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)

      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)

      ### informative prior #####
      if (new_time_pp_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_pp_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(event_inter*new_time_pp_inf_prior/t[inter_vt])

        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_pp_inf_prior=min(c(evt_max,new_sampsz_1_pp_inf_prior))

        samsz_simu_pp_inf_prior=max(new_sampsz_pp_inf_prior,evt_num)

        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################

        #      test_new_last_pp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_pp)
        test_new_last_pp_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_pp_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1

        }else{
          rej_pp_inf_prior=0

        }
      }

      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_pp_non_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(event_inter*new_time_pp_non_inf_prior/t[inter_vt])

        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_pp_non_inf_prior=min(c(evt_max,new_sampsz_1_pp_non_inf_prior))

        samsz_simu_pp_non_inf_prior=max(new_sampsz_pp_non_inf_prior,evt_num)

        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################

        #      test_new_last_pp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_pp)
        test_new_last_pp_non_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_pp_non_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)

      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_mpp_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(event_inter*new_time_mpp_inf_prior/t[inter_vt])
        #      T_new=new_sampsz/2/sigma^2
        #     new_c_bry=new_c_bry_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,T_new)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_mpp_inf_prior=min(c(evt_max,new_sampsz_1_mpp_inf_prior))

        samsz_simu_mpp_inf_prior=max(new_sampsz_mpp_inf_prior,evt_num)


        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################



        #         test_new_last_mpp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_mpp)
        test_new_last_mpp_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_mpp_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1

        }else{
          rej_mpp_inf_prior=0

        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_mpp_non_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(event_inter*new_time_mpp_non_inf_prior/t[inter_vt])
        #      T_new=new_sampsz/2/sigma^2
        #     new_c_bry=new_c_bry_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,T_new)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_mpp_non_inf_prior=min(c(evt_max,new_sampsz_1_mpp_non_inf_prior))

        samsz_simu_mpp_non_inf_prior=max(new_sampsz_mpp_non_inf_prior,evt_num)


        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################



        #         test_new_last_mpp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_mpp)
        test_new_last_mpp_non_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_mpp_non_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1

        }else{
          rej_mpp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)


      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_ci_mpp_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(event_inter*new_time_ci_mpp_inf_prior/t[inter_vt])

        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_ci_mpp_inf_prior=min(c(evt_max,new_sampsz_1_ci_mpp_inf_prior))

        samsz_simu_ci_mpp_inf_prior=max(new_sampsz_ci_mpp_inf_prior,evt_num)
        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################



        #     test_new_last_ci_mpp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_ci_mpp)

        test_new_last_ci_mpp_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_ci_mpp_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1

        }else{
          rej_ci_mpp_inf_prior=0

        }
      }
      #######################################


      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_ci_mpp_non_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(event_inter*new_time_ci_mpp_non_inf_prior/t[inter_vt])

        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_ci_mpp_non_inf_prior=min(c(evt_max,new_sampsz_1_ci_mpp_non_inf_prior))

        samsz_simu_ci_mpp_non_inf_prior=max(new_sampsz_ci_mpp_non_inf_prior,evt_num)
        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################



        #     test_new_last_ci_mpp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_ci_mpp)

        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_ci_mpp_non_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1

        }else{
          rej_ci_mpp_non_inf_prior=0

        }
      }
      #######################################


    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################

    #out
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest)


  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  #m_t_last=mean(t_last)
  #T_last=NF/2/sigma^2
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(evt_num,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)
  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""
  out_2_4=data.frame(References5)
  colnames(out_2_4)=""
  out_2=rbind(out_2_1,out_2_2,out_2_3,out_2_4)

  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
  #return(out_3)
}
#hybrid_fixed_prior_rej_survival_simu(t_prior=10,theta_prior=0.8,HR=1,lambda=1,s=c(0.5,1),alpha_2=0.05,surv_rate=0.5,evt_num=88,evt_max=331,HR_cut=0.95,cp_min=0.27,gamma=0.5,samp_multip=5,beta=0.1,sim_num=1000)
#hybrid_fixed_prior_rej_survival_simu(t_prior=20,theta_prior=0.5,HR=0.66,lambda=1,s=c(0.5,1),alpha_2=0.05,surv_rate=0.5,evt_num=88,evt_max=331,HR_cut=0.95,cp_min=0.27,gamma=0.5,samp_multip=5,beta=0.1,sim_num=10000)
#hybrid_fixed_prior_rej_survival_simu(t_prior=20,theta_prior=0.5,HR=0.63,lambda=1,s=c(0.5,1),alpha_2=0.05,surv_rate=0.5,evt_num=88,evt_max=331,HR_cut=0.95,cp_min=0.27,gamma=0.5,samp_multip=5,beta=0.1,sim_num=10000)
#hybrid_fixed_prior_rej_survival_simu(t_prior=20,theta_prior=0.66,HR=0.63,lambda=1,s=c(0.5,1),alpha_2=0.05,surv_rate=0.5,evt_num=88,evt_max=331,HR_cut=0.95,cp_min=0.27,gamma=0.5,samp_multip=5,beta=0.1,sim_num=1000)
#hybrid_fixed_prior_rej_survival_simu(t_prior=10,theta_prior=0.8,HR=1,lambda=1,s=c(0.5,1),alpha_2=0.05,surv_rate=0.5,evt_num=88,evt_max=331,HR_cut=0.95,cp_min=0.27,gamma=0.5,samp_multip=5,beta=0.1,sim_num=1000)


####### Random prior ######

#' Hybrid: Simulations, Survival analysis, Using random informative prior from previous trial
#'
#' @param HR_0 Hazards ratio of test vs. control in the previous trial, HR_0; the
#'   previous-trial effect size is theta_0 = -log(HR_0)
#' @param evt_num_0 Number of events from the previous trial, N_0; random samples are drawn from
#'   the previous trial up to this number of events to form the random prior
#' @param HR Hazards ratio of test arm vs. the control arm for the current trial. The effect size
#'   is theta = -log(HR) and a smaller HR indicates better efficacy. Enter HR = 1 to verify type
#'   I error control and HR < 1 for the power simulation
#' @param lambda Baseline hazard rate of the control arm for the exponential survival times: the
#'   control arm is generated as rexp(lambda) and the test arm as rexp(HR * lambda), so lambda
#'   is a rate (not a log-hazard) and must be > 0. It also sets the censoring time
#'   T = -log(surv_rate)/lambda for both the previous and the current trial
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param surv_rate Survival rate used to translate the number of events into a sample size; it
#'   is the survival probability at the (administrative) censoring time T = -log(surv_rate)/lambda, 0 < surv_rate < 1
#' @param evt_num Planned number of events: the planned number of events at the final analysis of the current trial
#' @param evt_max Maximum number of events: a user chosen, maximum allowable number of events for sample size modification
#' @param HR_cut Futility hazards ratio threshold: a cut-off for HR_hat. In the simulations, if
#'   the estimated hazards ratio at the interim analysis is at or above HR_cut the trial is
#'   considered futile and proceeds to the planned number of events without any change; it is
#'   converted internally to theta_cut = -log(HR_cut)
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param samp_multip Multiplication factor applied when translating the number of events of
#'   the previous and the current trial into sample sizes, for example
#'   sampsz = ceiling(samp_multip * evt_max / surv_rate); it must be large enough for the
#'   required number of events to be observed
#' @param sim_num Number of simulation replicates
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods. For this survival design the mean sample size is the mean number of events}
#'   \item{futility rate}{The rate at which the estimated hazards ratio at the interim analysis is at or above HR_cut, so that the trial proceeds without any change to the number of events}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' \donttest{
#' hybrid_random_prior_rej_survival_simu(HR_0 = 0.63, evt_num_0 = 40, HR = 1, lambda = 1,
#' s = c(0.5,1), alpha_2 = 0.05, surv_rate = 0.5, evt_num = 88, evt_max = 331, HR_cut = 0.95,
#' cp_min = 0.27, gamma = 0.5, samp_multip = 5, beta = 0.1, sim_num = 20)
#' }
#'
#' @export
hybrid_random_prior_rej_survival_simu<-function(HR_0,evt_num_0,HR,lambda,s,alpha_2,surv_rate,evt_num,evt_max,HR_cut,cp_min, gamma,samp_multip,sim_num,beta){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number


  #####################################################

  sampsz_0=ceiling(samp_multip*evt_num_0/surv_rate)

  evt_num_s=c(ceiling(evt_num*s),evt_max)
  sampsz=ceiling(samp_multip*evt_max/surv_rate)
  theta_cut=-log(HR_cut)
  mu=-log(HR)

  c_bry=stats::qnorm(1-alpha_2/2)


  ############################################
  ####  beginning of the simulation loops  ####
  ############################################

  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation
    ######################
    ### getting random sample for old study ###################

    x0_0=stats::rexp(sampsz_0,lambda) #generate exponential data for the control arm
    x_0=stats::rexp(sampsz_0,HR_0*lambda)  #generate exponential data for the test arm

    T=-log(surv_rate)/lambda
    #T_t=-log(surv_rate)/(HR*lambda)#Censoring time. Events after T will be censored
    time_x0_0=x0_0
    time_x_0=x_0
    x0_0_status=c(rep(1,sampsz_0))
    x_0_status=c(rep(1,sampsz_0))

    #######################


    ######## creating two long collum of censored data,  ###########
    for (i in (1:sampsz_0)){
      if (x0_0[i]>T){
        time_x0_0[i]=T  #x0[i] is censored after T
        x0_0_status[i]=0   #censoring indicator
      }
    }
    grp_x0_0=c(rep(1,length(time_x0_0)))#vector of zeros, length of the vector is the same as the original smaple size
    grp_x_0=c(rep(0,length(time_x_0))) #
    for (i in (1:sampsz_0)){
      if (x_0[i]>T){
        time_x_0[i]=T #x[i] is censored after T
        x_0_status[i]=0  #censoring indicator
      }
    }

    #### count the numnber of events in the
    evt_0_count=c(rep(0,sampsz_0))
    for (i in 1:sampsz_0){
      evt_0_count[i]=sum(x0_0_status[1:i])+sum(x_0_status[1:i]) #number of events after 2*i observations
    }

    #################################################################################################

    evt_0_pos=evt_0_count[evt_0_count>evt_num_0]#substring of  evt_count that are larger than evt_num
    evt_0_cut=match(evt_0_pos[1],evt_0_count) #the first incidence that evt_count is larger than evt_num

    #  vector of zeros, length of the vector is the same as the original smaple size
    time_0=c(time_x0_0[1:evt_0_cut],time_x_0[1:evt_0_cut]) #combined samples befroe total number of events exceed evt_num
    event_0=c(x0_0_status[1:evt_0_cut],x_0_status[1:evt_0_cut]) #censoring status before total number of events exceed evt_num
    grp_0=c(grp_x0_0[1:evt_0_cut],grp_x_0[1:evt_0_cut]) #vector of zero, with length=2*evt_cut


    test_0= list(time=time_0,status=event_0,grp_0)#survival data
    ########################################################################################################
    ####################### genrate score function , Wald testr stratistics from Cox partial likelihood regression model #####
    ###################################################################################################################
    coxcf_0=summary(coxph(Surv(time,status)~grp_0,test_0)) #run the coxph function and get the estimates
    # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
    # test is the same as data=test
    theta_hat_0= coxcf_0$coefficients[1] #estimated log(HR)
    theta_hat_sd_0= coxcf_0$coefficients[3] #estimated standard error of log(HR)
   # wald_cox=theta_hat/theta_hat_sd

    ### infromative prior parameters ###
    if(theta_hat_0>0.05){
    t_prior=theta_hat_sd_0^(-2)
    theta_prior=theta_hat_0
    }else{
      t_prior=0
      theta_prior=0
    }
    ###########################################
    ### getting random sample for new study ###################

    x0=stats::rexp(sampsz,lambda) #generate exponential data for the control arm
    x=stats::rexp(sampsz,HR*lambda)  #generate exponential data for the test arm
    T=-log(surv_rate)/lambda  #Censoring time. Events after T will be censored
    time_x0=x0
    time_x=x
    x0_status=c(rep(1,length(x0)))
    x_status=c(rep(1,length(x)))

    ######## creating two long collum of censored data,  ###########
    for (i in (1:sampsz)){
      if (x0[i]>T){
        time_x0[i]=T  #x0[i] is censored after T
        x0_status[i]=0   #censoring indicator
      }
    }
    grp_x0=c(rep(1,length(time_x0)))#vector of zeros, length of the vector is the same as the original smaple size
    grp_x=c(rep(0,length(time_x))) #
    for (i in (1:sampsz)){
      if (x[i]>T){
        time_x[i]=T #x[i] is censored after T
        x_status[i]=0  #censoring indicator
      }
    }

    #### count the numnber of events in the
    evt_count=c(rep(0,sampsz))
    for (i in 1:sampsz){
      evt_count[i]=sum(x0_status[1:i])+sum(x_status[1:i]) #number of events after 2*i observations
    }

    #################################################################################################
    theta_hat=rep(0,length(evt_num_s))
    theta_hat_sd=rep(0,length(evt_num_s))
    wald=rep(0,length(evt_num_s))
    score=rep(0,length(evt_num_s))
    t=rep(0,length(evt_num_s))
    for (i in (1:length(evt_num_s))) {
      evt_pos=evt_count[evt_count>evt_num_s[i]]#substring of  evt_count that are larger than evt_num
      evt_cut=match(evt_pos[1],evt_count) #the first incidence that evt_count is larger than evt_num

      #  vector of zeros, length of the vector is the same as the original smaple size
      time=c(time_x0[1:evt_cut],time_x[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
      event=c(x0_status[1:evt_cut],x_status[1:evt_cut]) #censoring status before total number of events exceed evt_num
      grp=c(grp_x0[1:evt_cut],grp_x[1:evt_cut]) #vector of zero, with length=2*evt_cut


      test= list(time=time,status=event,grp)#survival data
      coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
      # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
      # test is the same as data=test
      theta_hat[i]= coxcf$coefficients[1] #estimated log(HR)
      theta_hat_sd[i]= coxcf$coefficients[3] #estimated standard error of log(HR)
      t[i]=(theta_hat_sd[i])^(-2)
      wald[i]=theta_hat[i]/ theta_hat_sd[i]
      score[i]=theta_hat[i]*t[i]
    }

    ###############################################
    ############## generate the parameters ########
    ###############################################

    inter_vt=length(s)-1
    theta_hat_inter=theta_hat[inter_vt]
    theta_hat_inter_sd=theta_hat_sd[inter_vt]
    event_inter=evt_num_s[inter_vt]

    t_inter=t[1]
    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]
    Wald_inter=wald[1]
    ###### Planned sample size test ####
    last_vt=length(s)
    theta_hat_last=theta_hat[last_vt]
    theta_hat_last_sd=theta_hat_sd[last_vt]
    t_last= theta_hat_last_sd^(-2)

    Wald_last= theta_hat_last/theta_hat_last_sd
    if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
      rej_Ztest=1

    }else { # the NULL hypothesis is NOT rejected at the final analysis
      rej_Ztest=0

    }
    ###### Defining promising zone ########
    #      cp_inter=cond_power(T=T_planned,cry=c_bry,t_inter,s_inter,theta_hat_inter,Wald_inter)
    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)
    #    cond_power_inter=asd_c_power_theta(theta=theta_hat_inter,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)
    #    cp_min=0.27
    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      #stg identifies the interim analysis at which the null hypothesis is rejected
      last_vt=length(s)
      theta_hat_last=theta_hat[last_vt]
      theta_hat_last_sd=theta_hat_sd[last_vt]
      t_last= theta_hat_last_sd^(-2)

      Wald_last= theta_hat_last/theta_hat_last_sd


      samsz_simu_pz=evt_num_s[length(s)]

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(event_inter*new_time_pz_0/t[inter_vt])
      new_sampsz_1_pz=max(new_sampsz_0_pz, ceiling(event_inter*1.1))
      new_sampsz_pz=min(c(evt_max,new_sampsz_1_pz))

      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      samsz_simu_pz=max(new_sampsz_pz,evt_num)

      ##### Generate addtional random samples #######
      inter_vt=length(s)-1
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################

      #     test_new_last_pz=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,T_planned,event_inter,t_inter,new_samsz=samsz_simu_pz)
      #     test_new_last_pz=test_new_last_pz=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_pz)
      test_new_last_pz=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_pz)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }

    }


    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    # ## new information estimated using theta_hat_inter  ####



    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      #stg identifies the interim analysis at which the null hypothesis is rejected
      last_vt=length(s)
      theta_hat_last=theta_hat[last_vt]
      theta_hat_last_sd=theta_hat_sd[last_vt]
      t_last= theta_hat_last_sd^(-2)

      Wald_last= theta_hat_last/theta_hat_last_sd

      samsz_simu_cp=evt_num_s[length(s)]
      samsz_simu_cp_bayesian=evt_num_s[length(s)]
      samsz_simu_true_cp=evt_num_s[length(s)]
      samsz_simu_pp_non_inf_prior=evt_num_s[length(s)]
      samsz_simu_pp_inf_prior=evt_num_s[length(s)]
      samsz_simu_mpp_inf_prior=evt_num_s[length(s)]
      samsz_simu_mpp_non_inf_prior=evt_num_s[length(s)]
      samsz_simu_ci_mpp_inf_prior=evt_num_s[length(s)]
      samsz_simu_ci_mpp_non_inf_prior=evt_num_s[length(s)]


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(event_inter*new_time_cp_0/t[inter_vt])
      new_sampsz_1_cp=max(new_sampsz_0_cp, ceiling(event_inter*1.1))
      new_sampsz_cp=min(c(evt_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      samsz_simu_cp=max(new_sampsz_cp,evt_num)

      ##### Generate addtional random samples #######
      inter_vt=length(s)-1
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################

      #  test_new_last_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_cp)

      test_new_last_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_cp)

      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(event_inter*new_time_cp_0_bayesian/t[inter_vt])
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian, evt_num)
      samsz_simu_cp_bayesian=min(c(evt_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      #  samsz_simu_cp_bayesian=max(new_sampsz_cp_bayesian,evt_num)

      ##### Generate addtional random samples #######
      inter_vt=length(s)-1
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################

      #  test_new_last_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_cp)

      test_new_last_cp_bayesian=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_cp_bayesian)

      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(event_inter*new_time_true_cp_0/t[inter_vt])
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, ceiling(event_inter*1.1))
      new_sampsz_true_cp=min(c(evt_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      samsz_simu_true_cp=new_sampsz_true_cp

      ##### Generate addtional random samples #######
      inter_vt=length(s)-1
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################

      #    test_new_last_true_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_true_cp)
      test_new_last_true_cp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################

      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      #                        new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)

      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)

      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)

      ### informative prior #####
      if (new_time_pp_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_pp_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(event_inter*new_time_pp_inf_prior/t[inter_vt])

        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_pp_inf_prior=min(c(evt_max,new_sampsz_1_pp_inf_prior))

        samsz_simu_pp_inf_prior=max(new_sampsz_pp_inf_prior,evt_num)

        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################

        #      test_new_last_pp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_pp)
        test_new_last_pp_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_pp_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1

        }else{
          rej_pp_inf_prior=0

        }
      }

      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_pp_non_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(event_inter*new_time_pp_non_inf_prior/t[inter_vt])

        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_pp_non_inf_prior=min(c(evt_max,new_sampsz_1_pp_non_inf_prior))

        samsz_simu_pp_non_inf_prior=max(new_sampsz_pp_non_inf_prior,evt_num)

        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################

        #      test_new_last_pp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_pp)
        test_new_last_pp_non_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_pp_non_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)

      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_mpp_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(event_inter*new_time_mpp_inf_prior/t[inter_vt])
        #      T_new=new_sampsz/2/sigma^2
        #     new_c_bry=new_c_bry_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,T_new)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_mpp_inf_prior=min(c(evt_max,new_sampsz_1_mpp_inf_prior))

        samsz_simu_mpp_inf_prior=max(new_sampsz_mpp_inf_prior,evt_num)


        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################



        #         test_new_last_mpp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_mpp)
        test_new_last_mpp_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_mpp_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1

        }else{
          rej_mpp_inf_prior=0

        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_mpp_non_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(event_inter*new_time_mpp_non_inf_prior/t[inter_vt])
        #      T_new=new_sampsz/2/sigma^2
        #     new_c_bry=new_c_bry_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,T_new)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_mpp_non_inf_prior=min(c(evt_max,new_sampsz_1_mpp_non_inf_prior))

        samsz_simu_mpp_non_inf_prior=max(new_sampsz_mpp_non_inf_prior,evt_num)


        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################



        #         test_new_last_mpp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_mpp)
        test_new_last_mpp_non_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_mpp_non_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1

        }else{
          rej_mpp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)


      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_ci_mpp_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(event_inter*new_time_ci_mpp_inf_prior/t[inter_vt])

        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_ci_mpp_inf_prior=min(c(evt_max,new_sampsz_1_ci_mpp_inf_prior))

        samsz_simu_ci_mpp_inf_prior=max(new_sampsz_ci_mpp_inf_prior,evt_num)
        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################



        #     test_new_last_ci_mpp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_ci_mpp)

        test_new_last_ci_mpp_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_ci_mpp_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1

        }else{
          rej_ci_mpp_inf_prior=0

        }
      }
      #######################################


      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=length(s)
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        t_last= theta_hat_last_sd^(-2)

        Wald_last= theta_hat_last/theta_hat_last_sd

        samsz_simu_ci_mpp_non_inf_prior=evt_num_s[length(s)]
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(event_inter*new_time_ci_mpp_non_inf_prior/t[inter_vt])

        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior, ceiling(event_inter*1.1))

        new_sampsz_ci_mpp_non_inf_prior=min(c(evt_max,new_sampsz_1_ci_mpp_non_inf_prior))

        samsz_simu_ci_mpp_non_inf_prior=max(new_sampsz_ci_mpp_non_inf_prior,evt_num)
        ##### Generate addtional random samples #######
        inter_vt=length(s)-1
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################



        #     test_new_last_ci_mpp=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,new_samsz=samsz_simu_ci_mpp)

        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_surv(evt_count,grp_x0,grp_x,time_x0,time_x,x0_status,x_status,evt_num,event_inter,t_inter,score_inter,T_planned,c_bry,new_samsz=samsz_simu_ci_mpp_non_inf_prior)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1

        }else{
          rej_ci_mpp_non_inf_prior=0

        }
      }
      #######################################


    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################

    #out
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest)


  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  #m_t_last=mean(t_last)
  #T_last=NF/2/sigma^2
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)

  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(evt_num,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""
  out_2_4=data.frame(References5)
  colnames(out_2_4)=""
  out_2=rbind(out_2_1,out_2_2,out_2_3,out_2_4)


  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
  #return(out_3)
}
#################
#hybrid_random_prior_rej_survival_simu(HR_0,evt_num_0,HR,lambda,s,alpha_2,surv_rate,evt_num,evt_max,HR_cut,cp_min, gamma,samp_multip,sim_num,beta)
#hybrid_random_prior_rej_survival_simu(HR_0=0.63,evt_num_0=40,HR=1,lambda=1,s=c(0.5,1),alpha_2=0.05,surv_rate=0.5,evt_num=88,evt_max=331,HR_cut=0.95,cp_min=0.27,gamma=0.5,samp_multip=5,beta=0.1,sim_num=10000)
#hybrid_random_prior_rej_survival_simu(HR_0=0.63,evt_num_0=40,HR=0.63,lambda=1,s=c(0.5,1),alpha_2=0.05,surv_rate=0.5,evt_num=88,evt_max=331,HR_cut=0.95,cp_min=0.27,gamma=0.5,samp_multip=5,beta=0.1,sim_num=10000)
################################
#### Poisson distribution ######
################################

####### Fixed  prior ######


#####################
#' Hybrid: Simulations, Poisson distribution, Using given estimates of informative prior
#'
#' @param t_prior Prior information time: the prior distribution of
#'   theta = lambda_t - lambda_c is assumed to have a normal distribution, and the Fisher's
#'   information time is the reciprocal of the variance of the prior distribution; set
#'   t_prior = 0 for a non-informative prior
#' @param theta_prior Mean effect size of prior: the mean of the normal prior distribution of
#'   theta = lambda_t - lambda_c; set theta_prior = 0 for a non-informative prior
#' @param lambda_t Event rate in the test arm of the current trial: the number of events per
#'   subject, generated as rpois(lambda_t). The effect size is theta = lambda_t - lambda_c.
#'   Enter lambda_t = lambda_c to verify type I error control and lambda_t > lambda_c for the
#'   power simulation
#' @param lambda_c Event rate in the control arm of the current trial, generated as
#'   rpois(lambda_c); the effect size is theta = lambda_t - lambda_c
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz Planned sample size: the planned sample size at the final analysis
#' @param N_max Maximum sample size: a user chosen, maximum allowable sample size for sample size modification
#' @param theta_cut Futility threshold: a cut-off for theta_hat. In the simulations, if
#'   theta_hat at the interim analysis is at or below theta_cut the trial is considered futile
#'   and proceeds to the planned sample size without any change
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods}
#'   \item{futility rate}{The rate at which the observed theta_hat at the interim analysis is at or below theta_cut, so that the trial proceeds without any sample size change}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_fixed_prior_rej_poisson_simu(t_prior = 5, theta_prior = 0.2, lambda_t = 0.5,
#' lambda_c = 0.5, s = c(0.5,1), sampsz = 100, N_max = 150, theta_cut = 0.001, cp_min = 0.27,
#' gamma = 0.5, alpha_2 = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_fixed_prior_rej_poisson_simu<-function(t_prior,theta_prior,lambda_t,lambda_c,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number
  #####################################################

  N=ceiling(sampsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  N_inter=N[1]
  T_true<- N/(lambda_t+lambda_c)  # t is the vector of information times at each interim anbalysis
  NF=sampsz  ## NF is the sample size at the end of the trial
  mu=lambda_t-lambda_c

  c_bry=stats::qnorm(1-alpha_2/2)
  ############################################
  ####  beginning of the simulation loops  ####
  ############################################
  # rm(x)
  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation

    ########################################
    x0=stats::rpois(N_max,lambda_c)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rpois(N_max,lambda_t)

    #####################################################################
    ### Step 2: calculate the Wald statistics for hypothesis testing ####
    #####################################################################
    lambda_c_inter=mean(x0[1:N_inter])
    lambda_t_inter=mean(x[1:N_inter])
    t_inter= N_inter/(lambda_t_inter+lambda_c_inter)

    lambda_c_last=mean(x0[1:NF])
    lambda_t_last=mean(x[1:NF])
    t_last= NF/(lambda_t_last+lambda_c_last)

    ###############################################
    ############## generate the parameters ########
    ###############################################
    theta_hat_inter=mean(x[1:N[1]])-mean(x0[1:N[1]])
    theta_hat_last=mean(x[1:NF])-mean(x0[1:NF])


    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]
    Wald_inter=theta_hat_inter*sqrt(t_inter)
    Wald_last= theta_hat_last*sqrt(t_last)

    if(Wald_last>=c_bry){
      rej_Ztest=1
    }else{
      rej_Ztest=0
    }
    ###### Defining promising zone ########

    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)

    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      samsz_simu_pz=NF

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(N[1]*new_time_pz_0/t_inter)
      new_sampsz_1_pz=max(new_sampsz_0_pz, NF)
      ############
      samsz_simu_pz=min(c(N_max,new_sampsz_1_pz))
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################
      test_new_last_pz=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pz)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }
    }
    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      samsz_simu_cp=NF
      samsz_simu_cp_bayesian=NF
      samsz_simu_true_cp=NF
      samsz_simu_pp_non_inf_prior=NF
      samsz_simu_pp_inf_prior=NF
      samsz_simu_mpp_inf_prior=NF
      samsz_simu_mpp_non_inf_prior=NF
      samsz_simu_ci_mpp_inf_prior=NF
      samsz_simu_ci_mpp_non_inf_prior=NF


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
      new_sampsz_1_cp=max(new_sampsz_0_cp,NF)
      #######
      samsz_simu_cp=min(c(N_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_0_bayesian/t_inter)
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian,NF)
      ##############
      samsz_simu_cp_bayesian=min(c(N_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp_bayesian=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp_bayesian)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(N_inter*new_time_true_cp_0/t_inter)
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, NF)
      #################
      samsz_simu_true_cp=min(c(N_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      test_new_last_true_cp=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################
      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_pp_inf_prior==9999){

        samsz_simu_pp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, NF)
        ###############
        samsz_simu_pp_inf_prior=min(c(N_max,new_sampsz_1_pp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1
        }else{
          rej_pp_inf_prior=0
        }
      }
      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){

        samsz_simu_pp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,NF)
        #########################
        samsz_simu_pp_non_inf_prior=min(c(N_max,new_sampsz_1_pp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_non_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){

        samsz_simu_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,NF)
        ###################
        samsz_simu_mpp_inf_prior=min(c(N_max,new_sampsz_1_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1
        }else{
          rej_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){

        samsz_simu_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,NF)
        #################
        samsz_simu_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_non_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1
        }else{
          rej_mpp_non_inf_prior=0
        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)

      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){

        samsz_simu_ci_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1
        }else{
          rej_ci_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){

        samsz_simu_ci_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1
        }else{
          rej_ci_mpp_non_inf_prior=0
        }
      }
      #######################################
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest,t_last)
    #out

  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  m_t_last=mean(t_last)
  T_last=NF/(lambda_t+lambda_c)
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)


  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(NF,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""

  out_2=rbind(out_2_1,out_2_2,out_2_3)


  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
}
######################################################################
#hybrid_fixed_prior_rej_poisson_simu(t_prior,theta_prior,lambda_t,lambda_c,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num)
########################
#hybrid_fixed_prior_rej_poisson_simu(t_prior=5,theta_prior=0.2,lambda_t=0.5,lambda_c=0.5,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_fixed_prior_rej_poisson_simu(t_prior=5,theta_prior=0.2,lambda_t=0.7,lambda_c=0.5,s=c(0.5,1),sampsz=100,N_max=250,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)

#####################

####### Random prior ######

#' Hybrid: Simulations, Poisson distribution, Using random informative prior from previous trial
#'
#' @param lambda_t_0 Event rate in the test arm of the previous trial, generated as
#'   rpois(lambda_t_0); the previous-trial effect size is theta_0 = lambda_t_0 - lambda_c_0
#' @param lambda_c_0 Event rate in the control arm of the previous trial, generated as
#'   rpois(lambda_c_0)
#' @param sampsz_0 Sample size per arm from the previous trial, N_0; random samples of this size are drawn from the previous trial to form the random prior
#' @param lambda_t Event rate in the test arm of the current trial: the number of events per
#'   subject, generated as rpois(lambda_t). The effect size is theta = lambda_t - lambda_c.
#'   Enter lambda_t = lambda_c to verify type I error control and lambda_t > lambda_c for the
#'   power simulation
#' @param lambda_c Event rate in the control arm of the current trial, generated as
#'   rpois(lambda_c); the effect size is theta = lambda_t - lambda_c
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz Planned sample size: the planned sample size at the final analysis
#' @param N_max Maximum sample size: a user chosen, maximum allowable sample size for sample size modification
#' @param theta_cut Futility threshold: a cut-off for theta_hat. In the simulations, if
#'   theta_hat at the interim analysis is at or below theta_cut the trial is considered futile
#'   and proceeds to the planned sample size without any change
#' @param cp_min Promising zone cp_min: CP_min in the promising zone sample size re-estimation;
#'   a conditional power below cp_min is the unfavorable zone. The user should select this
#'   value from Mehta, Pocock (2011)
#' @param gamma Parameter chosen for the credible interval modified predictive power (CI_MPP);
#'   it defines the (1 - gamma) credible interval of the posterior for theta, e.g. gamma = 0.5
#'   gives a 50 percent interval
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Conditional type II error; 1 - beta is the target conditional power used for the
#'   sample size re-estimation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate and mean sample size}{A data frame with the columns "Prior type", "SSR method", "Rejection rate" and "Mean sample size", with one row per sample size re-estimation (SSR) method. The methods compared are the planned sample size with no SSR, the promising zone SSR, the SSR using the true effect size with conditional power, and, for both the non-informative prior and the informative prior, the cp SSR, pp SSR, mpp SSR and ci_mpp SSR methods}
#'   \item{futility rate}{The rate at which the observed theta_hat at the interim analysis is at or below theta_cut, so that the trial proceeds without any sample size change}
#'   \item{References}{A data frame of the references for the methods compared}
#' }
#'
#' @references Gao, P. Improving the effectiveness of sample size re-estimation: an operating characteristic focused, hybrid frequentist-Bayesian approach, Statistics in Medicine, 2025; 44:e10310. To link to this article: \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.10310}.
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references C. R. Mehta, S. J. Pocock. Adaptive increase in sample size when interim results are promising: A practical guide with examples. Statistics in Medicine. Volume30, Issue28, 10 December 2011, Pages 3267-3284
#'
#' @examples
#' hybrid_random_prior_rej_poisson_simu(lambda_t_0 = 0.5, lambda_c_0 = 0.3, sampsz_0 = 60,
#' lambda_t = 0.5, lambda_c = 0.5, s = c(0.5,1), sampsz = 100, N_max = 150, theta_cut = 0.001,
#' cp_min = 0.27, gamma = 0.5, alpha_2 = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
hybrid_random_prior_rej_poisson_simu<-function(lambda_t_0,lambda_c_0,sampsz_0,lambda_t,lambda_c,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num){

  ## This function conducts simulations to evaluate the operating characteristics for the promising zone and the hybrid procedure ##
  #s is the vector of information fractions
  #sim_num=simulation number
  #####################################################

  N=ceiling(sampsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  N_inter=N[1]
  T_true<- N/(lambda_t+lambda_c)  # t is the vector of information times at each interim anbalysis
  NF=sampsz  ## NF is the sample size at the end of the trial
  mu=lambda_t-lambda_c

  c_bry=stats::qnorm(1-alpha_2/2)
  ############################################
  ####  beginning of the simulation loops  ####
  ############################################
  # rm(x)
  simr <- simuloop(sim_num, "HYBRID", function() {

    #the loop specifies what is generated for each simulation
    ###########################################
    ### geting random sample for old trial ################

    x0_0=stats::rpois(sampsz_0,lambda_c_0)
    #x0 is the random sample at the first interim analysis for the control
    x_0=stats::rpois(sampsz_0,lambda_t_0)
    ###############################################
    ############## generate the parameters ########
    ###############################################

    #  inter_vt=length(s)-1
    theta_hat_0=mean(x_0)-mean(x0_0)
    t_0=sampsz_0/(mean(x_0)+mean(x0_0))

    if ( theta_hat_0>0.05){
      t_prior=t_0
      theta_prior=theta_hat_0
    }else{
      t_prior=0
      theta_prior=0
    }

    ###########################################
    ########################################
    x0=stats::rpois(N_max,lambda_c)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rpois(N_max,lambda_t)

    #####################################################################
    ### Step 2: calculate the Wald statistics for hypothesis testing ####
    #####################################################################
    lambda_c_inter=mean(x0[1:N_inter])
    lambda_t_inter=mean(x[1:N_inter])
    t_inter= N_inter/(lambda_t_inter+lambda_c_inter)

    lambda_c_last=mean(x0[1:NF])
    lambda_t_last=mean(x[1:NF])
    t_last= NF/(lambda_t_last+lambda_c_last)

    ###############################################
    ############## generate the parameters ########
    ###############################################
    theta_hat_inter=mean(x[1:N[1]])-mean(x0[1:N[1]])
    theta_hat_last=mean(x[1:NF])-mean(x0[1:NF])


    s_inter=s[1]
    score_inter=theta_hat_inter*t_inter
    T_planned=t_inter/s[1]
    Wald_inter=theta_hat_inter*sqrt(t_inter)
    Wald_last= theta_hat_last*sqrt(t_last)

    if(Wald_last>=c_bry){
      rej_Ztest=1
    }else{
      rej_Ztest=0
    }
    ###### Defining promising zone ########

    cp_inter=cond_power(T=T_planned,c_bry,t_inter,theta_hat_inter,theta=theta_hat_inter)

    if (cp_inter<cp_min){
      zone=1  #unfavorable zone
    }else if (cp_inter>=1-beta){
      zone=3  #favorable zone
    }else if ((cp_inter<1-beta)&&(cp_inter>=cp_min)){
      zone=2  #promising zone
    }

    #### The promising zone procedure  ####
    if ((zone==1)||(zone==3)){

      samsz_simu_pz=NF

      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_pz=1

      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_pz=0

      }

    }else if (zone==2){

      ############################# new infromation time using cp #####
      new_time_pz_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_pz=ceiling(N[1]*new_time_pz_0/t_inter)
      new_sampsz_1_pz=max(new_sampsz_0_pz, NF)
      ############
      samsz_simu_pz=min(c(N_max,new_sampsz_1_pz))
      #################################################################################
      ###                   Step 4.5: generate new random samples                  ####
      ###                         and calculate statistics                         ####
      #################################################################################
      test_new_last_pz=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pz)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_pz >0){
        rej_pz=1
      }else{
        rej_pz=0
      }
    }
    ##########################################################################
    ### Step 3:obtain the new sample sizes using four different methods  ####
    ##########################################################################

    if (theta_hat_inter<=theta_cut){futi=1}else {futi=0}

    if (futi==1){

      samsz_simu_cp=NF
      samsz_simu_cp_bayesian=NF
      samsz_simu_true_cp=NF
      samsz_simu_pp_non_inf_prior=NF
      samsz_simu_pp_inf_prior=NF
      samsz_simu_mpp_inf_prior=NF
      samsz_simu_mpp_non_inf_prior=NF
      samsz_simu_ci_mpp_inf_prior=NF
      samsz_simu_ci_mpp_non_inf_prior=NF


      if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
        rej_cp=1
        rej_cp_bayesian=1
        rej_true_cp=1
        rej_pp_non_inf_prior=1
        rej_pp_inf_prior=1
        rej_mpp_non_inf_prior=1
        rej_mpp_inf_prior=1
        rej_ci_mpp_non_inf_prior=1
        rej_ci_mpp_inf_prior=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej_cp=0
        rej_cp_bayesian=0
        rej_true_cp=0
        rej_pp_non_inf_prior=0
        rej_pp_inf_prior=0
        rej_mpp_non_inf_prior=0
        rej_mpp_inf_prior=0
        rej_ci_mpp_non_inf_prior=0
        rej_ci_mpp_inf_prior=0
      }

    }else if (futi==0){

      ################################################################
      ##### The conditional power  procedure estimated CP ############
      ################################################################
      theta_hat_bayesian=( t_prior*theta_prior+t_inter*theta_hat_inter)/(t_prior+t_inter)
      ############################# new infromation time using cp #####
      new_time_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter,beta)
      new_sampsz_0_cp=ceiling(N_inter*new_time_cp_0/t_inter)
      new_sampsz_1_cp=max(new_sampsz_0_cp,NF)
      #######
      samsz_simu_cp=min(c(N_max,new_sampsz_1_cp))
      #################################################################################
      ###                   Step 4.4: Compare theta_hat with theta_cut             ####
      ###                     if new_sampsz<=N_max, proceed to new sample size     ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp >0){
        rej_cp=1
      }else{
        rej_cp=0
      }
      #######################################
      ############################# new infromation time using Bayesian cp #####
      new_time_cp_0_bayesian=new_time_cp(T_planned,c_bry,t_inter,theta_hat_bayesian,beta)
      new_sampsz_0_cp_bayesian=ceiling(N_inter*new_time_cp_0_bayesian/t_inter)
      new_sampsz_1_cp_bayesian=max(new_sampsz_0_cp_bayesian,NF)
      ##############
      samsz_simu_cp_bayesian=min(c(N_max,new_sampsz_1_cp_bayesian))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      test_new_last_cp_bayesian=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_cp_bayesian)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_cp_bayesian >0){
        rej_cp_bayesian=1
      }else{
        rej_cp_bayesian=0
      }
      #######################################################################
      ##### The conditional power  procedure , using true theta(mu)  ############
      #######################################################################

      ############################# new infromation time using cp #####
      new_time_true_cp_0=new_time_cp(T_planned,c_bry,t_inter,theta_hat_inter=mu,beta)
      new_sampsz_0_true_cp=ceiling(N_inter*new_time_true_cp_0/t_inter)
      new_sampsz_1_true_cp=max(new_sampsz_0_true_cp, NF)
      #################
      samsz_simu_true_cp=min(c(N_max,new_sampsz_1_true_cp))
      #################################################################################
      ###                   Step 4.3: Compare theta_hat with theta_cut             ####
      ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################

      test_new_last_true_cp=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_true_cp)
      #################################################################################
      ###                   Step 4.6: Hypothesis testing                           ####
      ########                   and calculate confidenc einterval                 ####
      #################################################################################
      if (test_new_last_true_cp >0){
        rej_true_cp=1
      }else{
        rej_true_cp=0
      }
      #######################################
      ##### The pp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_pp_inf_prior=new_time_pp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_pp_non_inf_prior=new_time_pp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_pp_inf_prior==9999){

        samsz_simu_pp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_inf_prior=0
        }
      }else if (new_time_pp_inf_prior!=9999){

        new_sampsz_0_pp_inf_prior=ceiling(N_inter*new_time_pp_inf_prior/t_inter)
        new_sampsz_1_pp_inf_prior=max(new_sampsz_0_pp_inf_prior, NF)
        ###############
        samsz_simu_pp_inf_prior=min(c(N_max,new_sampsz_1_pp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_pp_inf_prior>0){
          rej_pp_inf_prior=1
        }else{
          rej_pp_inf_prior=0
        }
      }
      ### non-informative prior #####
      if (new_time_pp_non_inf_prior==9999){

        samsz_simu_pp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_pp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_pp_non_inf_prior=0
        }
      }else if (new_time_pp_non_inf_prior!=9999){

        new_sampsz_0_pp_non_inf_prior=ceiling(N_inter*new_time_pp_non_inf_prior/t_inter)
        new_sampsz_1_pp_non_inf_prior=max(new_sampsz_0_pp_non_inf_prior,NF)
        #########################
        samsz_simu_pp_non_inf_prior=min(c(N_max,new_sampsz_1_pp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_pp_non_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_pp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (test_new_last_pp_non_inf_prior>0){
          rej_pp_non_inf_prior=1

        }else{
          rej_pp_non_inf_prior=0

        }
      }
      #######################################
      #######################################
      ##### The mpp procedure ############
      #######################################
      ############################# new information time using mpp #####
      new_time_mpp_inf_prior=new_time_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,beta)
      new_time_mpp_non_inf_prior=new_time_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,beta)
      ### informative prior #####
      if (new_time_mpp_inf_prior==9999){

        samsz_simu_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_inf_prior=0
        }
      }else if (new_time_mpp_inf_prior!=9999){

        new_sampsz_0_mpp_inf_prior=ceiling(N_inter*new_time_mpp_inf_prior/t_inter)
        new_sampsz_1_mpp_inf_prior=max(new_sampsz_0_mpp_inf_prior,NF)
        ###################
        samsz_simu_mpp_inf_prior=min(c(N_max,new_sampsz_1_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_inf_prior>0){
          rej_mpp_inf_prior=1
        }else{
          rej_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_mpp_non_inf_prior==9999){

        samsz_simu_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_mpp_non_inf_prior=0
        }
      }else if (new_time_mpp_non_inf_prior!=9999){

        new_sampsz_0_mpp_non_inf_prior=ceiling(N_inter*new_time_mpp_non_inf_prior/t_inter)
        new_sampsz_1_mpp_non_inf_prior=max(new_sampsz_0_mpp_non_inf_prior,NF)
        #################
        samsz_simu_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_mpp_non_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_mpp_non_inf_prior>0){
          rej_mpp_non_inf_prior=1
        }else{
          rej_mpp_non_inf_prior=0
        }
      }
      #######################################
      #######################################
      ##### The ci_mpp procedure #######
      #######################################
      ############################# new information time using mpp #####
      new_time_ci_mpp_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior,theta_prior,t_inter,theta_hat_inter,gamma,beta)
      new_time_ci_mpp_non_inf_prior=new_time_ci_mpp(T_planned,c_bry,t_prior=0,theta_prior=0,t_inter,theta_hat_inter,gamma,beta)

      ### informative prior #####
      if (new_time_ci_mpp_inf_prior==9999){

        samsz_simu_ci_mpp_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_inf_prior=0
        }
      }else if (new_time_ci_mpp_inf_prior!=9999){

        new_sampsz_0_ci_mpp_inf_prior=ceiling(N_inter*new_time_ci_mpp_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_inf_prior=max(new_sampsz_0_ci_mpp_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_inf_prior>0){
          rej_ci_mpp_inf_prior=1
        }else{
          rej_ci_mpp_inf_prior=0
        }
      }
      #######################################
      ### non-informative prior #####
      if (new_time_ci_mpp_non_inf_prior==9999){

        samsz_simu_ci_mpp_non_inf_prior=NF
        if (Wald_last>c_bry){ # the NULL hypothesis is rejected at the final analysis
          rej_ci_mpp_non_inf_prior=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej_ci_mpp_non_inf_prior=0
        }
      }else if (new_time_ci_mpp_non_inf_prior!=9999){

        new_sampsz_0_ci_mpp_non_inf_prior=ceiling(N_inter*new_time_ci_mpp_non_inf_prior/t_inter)
        new_sampsz_1_ci_mpp_non_inf_prior=max(new_sampsz_0_ci_mpp_non_inf_prior,NF)
        ##################
        samsz_simu_ci_mpp_non_inf_prior=min(c(N_max,new_sampsz_1_ci_mpp_non_inf_prior))
        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################
        test_new_last_ci_mpp_non_inf_prior=new_samsz_test_poisson(x,x0,t_inter,score_inter,T_planned,t_last,Wald_last,c_bry,NF,new_samsz=samsz_simu_ci_mpp_non_inf_prior)
        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (test_new_last_ci_mpp_non_inf_prior>0){
          rej_ci_mpp_non_inf_prior=1
        }else{
          rej_ci_mpp_non_inf_prior=0
        }
      }
      #######################################
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out_simu=c(rej_pz,rej_cp,rej_cp_bayesian,rej_true_cp,rej_pp_non_inf_prior,rej_pp_inf_prior,rej_mpp_non_inf_prior,rej_mpp_inf_prior,
               rej_ci_mpp_non_inf_prior, rej_ci_mpp_inf_prior,
               samsz_simu_pz,samsz_simu_cp,samsz_simu_cp_bayesian,samsz_simu_true_cp,
               samsz_simu_pp_inf_prior,samsz_simu_pp_non_inf_prior,
               samsz_simu_mpp_inf_prior,samsz_simu_mpp_non_inf_prior,samsz_simu_ci_mpp_inf_prior,
               samsz_simu_ci_mpp_non_inf_prior,
               futi,theta_hat_inter,zone,cp_inter,
               rej_Ztest,t_last)
    #out

  })#(end of simulation loops)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################

  rej_pz =as.numeric(simr[,1])
  rej_cp =as.numeric(simr[,2])
  rej_cp_bayesian =as.numeric(simr[,3])
  rej_true_cp =as.numeric(simr[,4])
  rej_pp_non_inf_prior =as.numeric(simr[,5])
  rej_pp_inf_prior =as.numeric(simr[,6])
  rej_mpp_non_inf_prior =as.numeric(simr[,7])
  rej_mpp_inf_prior =as.numeric(simr[,8])
  rej_ci_mpp_non_inf_prior =as.numeric(simr[,9])
  rej_ci_mpp_inf_prior =as.numeric(simr[,10])
  samsz_simu_pz =as.numeric(simr[,11])
  samsz_simu_cp =as.numeric(simr[,12])
  samsz_simu_cp_bayesian =as.numeric(simr[,13])
  samsz_simu_true_cp =as.numeric(simr[,14])
  samsz_simu_pp_inf_prior =as.numeric(simr[,15])
  samsz_simu_pp_non_inf_prior =as.numeric(simr[,16])
  samsz_simu_mpp_inf_prior =as.numeric(simr[,17])
  samsz_simu_mpp_non_inf_prior =as.numeric(simr[,18])
  samsz_simu_ci_mpp_inf_prior =as.numeric(simr[,19])
  samsz_simu_ci_mpp_non_inf_prior =as.numeric(simr[,20])
  futi =as.numeric(simr[,21])
  theta_hat_inter =as.numeric(simr[,22])
  zone =as.numeric(simr[,23])
  cp_inter =as.numeric(simr[,24])
  rej_Ztest =as.numeric(simr[,25])
  t_last =as.numeric(simr[,26])
  ################################################
  #min(samsz_simu_cp)
  ### promising zone summary ###
  rej_pz_rt=mean(rej_pz)
  m_samsz_pz=mean(samsz_simu_pz)

  ### All other SSR summary ######

  ### rejection rate ###
  rej_cp_rt=mean(rej_cp)
  rej_cp_bayesian_rt=mean(rej_cp_bayesian)
  rej_true_cp_rt=mean(rej_true_cp)
  rej_pp_non_inf_prior_rt=mean(rej_pp_non_inf_prior)
  rej_pp_inf_prior_rt=mean(rej_pp_inf_prior)
  rej_mpp_non_inf_prior_rt=mean(rej_mpp_non_inf_prior)
  rej_mpp_inf_prior_rt=mean(rej_mpp_inf_prior)
  rej_ci_mpp_non_inf_prior_rt=mean(rej_ci_mpp_non_inf_prior)
  rej_ci_mpp_inf_prior_rt=mean(rej_ci_mpp_inf_prior)
  #########
  rej_Ztest_rt=mean(rej_Ztest)
  m_t_last=mean(t_last)
  T_last=NF/(lambda_t+lambda_c)
  ### Sample size #####

  m_samsz_cp=mean(samsz_simu_cp)
  m_samsz_true_cp=mean(samsz_simu_true_cp)
  m_samsz_cp_bayesian=mean(samsz_simu_cp_bayesian)
  m_samsz_pp_non_inf_prior=mean(samsz_simu_pp_non_inf_prior)
  m_samsz_pp_inf_prior=mean(samsz_simu_pp_inf_prior)
  m_samsz_mpp_non_inf_prior=mean(samsz_simu_mpp_non_inf_prior)
  m_samsz_mpp_inf_prior=mean(samsz_simu_mpp_inf_prior)
  m_samsz_ci_mpp_non_inf_prior=mean(samsz_simu_ci_mpp_non_inf_prior)
  m_samsz_ci_mpp_inf_prior=mean(samsz_simu_ci_mpp_inf_prior)


  futi_rt=mean(futi)


  #######################
  names_2=c("Planned sample size","Promising zone SSR","True effect size cp SSR","",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"),"",c("cp SSR","pp SSR","mpp SSR","ci_mpp SSR"))
  names_1=c("","","","non-informative prior","","","","",
            "informative prior","","","","")
  rej_rt=c(rej_Ztest_rt,
           rej_pz_rt,
           rej_true_cp_rt,
           "",
           rej_cp_rt,
           rej_pp_non_inf_prior_rt,
           rej_mpp_non_inf_prior_rt,
           rej_ci_mpp_non_inf_prior_rt,
           "",
           rej_cp_bayesian_rt,
           rej_pp_inf_prior_rt,
           rej_mpp_inf_prior_rt,
           rej_ci_mpp_inf_prior_rt
  )
  m_samsz=c(NF,
            m_samsz_pz,
            m_samsz_true_cp,
            "",
            m_samsz_cp,
            m_samsz_pp_non_inf_prior,
            m_samsz_mpp_non_inf_prior,
            m_samsz_ci_mpp_non_inf_prior,
            "",
            m_samsz_cp_bayesian,
            m_samsz_pp_inf_prior,
            m_samsz_mpp_inf_prior,
            m_samsz_ci_mpp_inf_prior
  )
  # out_0=cbind(names_1,names_2)
  out_0=cbind(names_1,names_2, rej_rt,m_samsz)
  out_1=data.frame(out_0)
  names(out_1)=c("Prior type","SSR method","Rejection rate","Mean sample size")
  ###############################
  # out_2=c(rej_Ztest_rt, m_t_last,T_last)
  #  out_2=rbind(References3,
  #              References1)

  out_2_1=data.frame(References3)
  colnames(out_2_1)=""
  out_2_2=data.frame(References1)
  colnames(out_2_2)=""
  out_2_3=data.frame(References2)
  colnames(out_2_3)=""

  out_2=rbind(out_2_1,out_2_2,out_2_3)


  #  out_2.df=data.frame(out_2)
  #  rownames(out_2.df)=c("","")
  #  colnames(out_2.df)=""
  out=list(out_1,futi_rt,out_2)
  names(out)=c("rejection rate and mean sample size","futility rate","References")
  return(out)
}
######################################################################
#hybrid_random_prior_rej_poisson_sim(lambda_t_0,lambda_c_0,sampsz_0,lambda_t,lambda_c,s,sampsz,N_max,theta_cut,cp_min, gamma,alpha_2,beta,sim_num)
########################
#hybrid_random_prior_rej_poisson_simu(lambda_t_0=0.5,lambda_c_0=0.3,sampsz_0=60,lambda_t=0.5,lambda_c=0.5,s=c(0.5,1),sampsz=100,N_max=150,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=100000)
#hybrid_random_prior_rej_poisson_simu(lambda_t_0=0.5,lambda_c_0=0.3,sampsz_0=60,lambda_t=0.7,lambda_c=0.5,s=c(0.5,1),sampsz=100,N_max=350,theta_cut=0.001,cp_min=0.27, gamma=0.5,alpha_2=0.05,beta=0.1,sim_num=10000)

######



