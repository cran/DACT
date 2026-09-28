#rm(list=ls())
#library(dplyr)
#library(foreach)
#library(parallel)
#library(doParallel)
#library(bigstatsr)
#library(tidyverse)
#library(ggplot2)


################################################
##   part I. Component functions            ####
##  for group sequential trial designs      ####
################################################

################################################################
###  Part I.1: functions for critical boundary calculations ####
################################################################

## critical bundaries according to pre-selected alpha -spending ###
alsp_bundry<-function(s,alpha){
  #this function selects critical boundaries #
  #length(s)<=3
  #the boundary selection is based on pre-determined alpha spending at each interm analysis#
  cov_mtx= matrix(rep(s,length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #### The functions c_cal1, c_cal2, c_cal3 calcuate the  exit probablities for the given exit bounadries e0_1, e0_2,e0_3, bry1, bry2#

  c_cal1<-function(e0_1){
    #e0_1 is the  exit boundary for the first interim analysis#
    up_bry=e0_1
    tmp = 1-stats::pnorm(up_bry)
    return(tmp)
  }

  ##########################################
  c_cal2<-function(bry,e0_2,sigma1){
    #bry1 is the  exit boundary for the first interim analysis#
    #e0_2 is the  exit boundary for the second interim analysis#
    #sigma2 is the covariance matrix at the second interim analysis#
    up_bry_0=c(bry,e0_2)
    rank=length(up_bry_0)
    up_bry=up_bry_0*sqrt(s[1:rank])
    sigma2=sigma1[1:rank,1:rank]
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rank),sigma=sigma2,algorithm=Miwa())[1]
    return(tmp)
  }

  ###########################################
  c_cal3<-function(up_bry_0,rank,sigma){
    #e0 is the  exit boundary #
    #sigma is the covariance matrix, it has the rank of =rep_time x reptime#
    sigma1=sigma[1:rank,1:rank]
    up_bry=up_bry_0*sqrt(s[1:rank])
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rank),sigma=sigma1,algorithm=Miwa())[1]
    return(tmp)
  }

  if(length(s)==1){
    bry=as.numeric(stats::uniroot(function(x) c_cal1(x)-alpha,c(-2^8,2^8))[1])
    #bry1,bry2, bry3 are exit boundaries. c_bry are the ctritical bundaries#
    c_bry=bry
    cum_alpha_0=1-stats::pnorm(bry)
  } else if(length(s)>1){
    bry=as.numeric(stats::uniroot(function(x) c_cal1(x)-alpha[1],c(-2^8,2^8))[1])
    for (i in (2:length(s))){
      bry_i=as.numeric(stats::uniroot(function(x) c_cal2(bry,x,cov_mtx)-sum(alpha[1:i]),c(-2^8,2^8),tol=1e-8)[1])
      bry=c(bry,bry_i)
    }
  }
  c_bry=bry
  ## the next loop canculate the exit bpundaries bry1, bry2,bry3 and critical boundaries c_bry, when there are a to of three interim analyses##
  #####################################

  cumu_alpha=c(rep(0,length(s)))
  for (i in 1:length(s)){
    cumu_alpha[i]=c_cal3(c_bry[1:i],i,cov_mtx)
  }
  alpha_s=c(rep(0,length(s)))
  if(length(s)==1){
    alpha_s=cumu_alpha
  } else {
    alpha_s[1]=cumu_alpha[1]
    for (i in 2:length(s)){
      alpha_s[i]=cumu_alpha[i]-cumu_alpha[i-1]
    }
  }

  nom_alpha=1-stats::pnorm(c_bry)
  sum_nom_alpha=sum(nom_alpha)
  out=list(c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
##############################
##############################
#alsp_bundry(1,0.025)
#alsp_bundry(c(0.25,1),c(0.0001,0.0249))
#alsp_bundry(c(0.5,1),c(0.0001,0.0249))
#alsp_bundry(c(0.75,1),c(0.0001,0.0249))
#alsp_bundry(c(1/3,2/3,1),c(0.0001,0.0001,0.0248))

###############################

### This function calculates the O'Brien-Fleming boundary #####
OF_bundry<-function(s,alpha){
  #this function selects critical boundaries #
  #the boundary selection is O'Brien-Fleming type#
  #This function also calculates the actual alpha spending at each interim analysis########
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #alpha  is the target one-sided type I error#
  cov_mtx= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]

  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e0#
  c_cal1<-function(e0,sigma){
    #e0 is the  exit boundary #
    #sigma is the covariance matrix, it has the rank of =rep_time x reptime#
    rank=length(s)
    up_bry=rep(e0,rank)
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rank),sigma=sigma,algorithm=Miwa())[1]
    return(tmp)
  }

  #########################
  c_cal2<-function(up_bry_0,rank,sigma){
    #e0 is the  exit boundary #
    #sigma is the covariance matrix, it has the rank of =rep_time x reptime#
    sigma1=sigma[1:rank,1:rank]
    up_bry=up_bry_0*sqrt(s[1:rank])
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rank),sigma=sigma1,algorithm=Miwa())[1]
    return(tmp)
  }

  ########################################
  bry=as.numeric(stats::uniroot(function(x) c_cal1(x,cov_mtx)-alpha,c(-2^8,2^8),tol=1e-8)[1])
  #bry are exit boundaries. c_bry are the ctritical bundaries#
  c_bry=bry/sqrt(s)

  ############### alpha spending######
  cumu_alpha=c(rep(0,length(s)))
  for (i in 1:length(s)){
    cumu_alpha[i]=c_cal2(c_bry[1:i],i,cov_mtx)
  }
  alpha_s=c(rep(0,length(s)))
  if(length(s)==1){
    alpha_s=cumu_alpha
  } else {
    alpha_s[1]=cumu_alpha[1]
    for (i in 2:length(s)){
      alpha_s[i]=cumu_alpha[i]-cumu_alpha[i-1]
    }
  }

  nom_alpha=1-stats::pnorm(c_bry)
  #################################################
  sum_nom_alpha=sum( nom_alpha)

  out=list(c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
##############################
#OF_bundry(1,0.025)
#OF_bundry(c(0.25,1),0.025)
#OF_bundry(c(0.5,1),0.025)
#OF_bundry(c(0.75,1),0.025)
#OF_bundry(c(1/3,2/3,1),0.025)
#OF_bundry(c(0.3,0.7,1),0.025)
#OF_bundry(c(1/4,1/3,2/3,1),0.025)
###############################

#### This function calculates the Pocock boundary ######
pocock_bundry<-function(s,alpha){
  #this function selects critical boundaries #
  #the boundary selection is O'Brien-Fleming type#
  #This function also calculates the actual alpha spending at each interim analysis########
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #alpha  is the target one-sided type I error#
  cov_mtx= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]

  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e0#
  c_cal<-function(e0,sigma){
    #e0 is the  exit boundary #
    rank=length(s)
    #sigma is the covariance matrix, it has the rank of =rep_time x reptime#
    up_bry=e0*sqrt(s[1:rank])
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rank),sigma=sigma,algorithm=Miwa())[1]
    return(tmp)
  }

  c_bry=as.numeric(stats::uniroot(function(x) c_cal(x,cov_mtx)-alpha,c(-2^8,2^8),tol=1e-8)[1])
  #c_bry are the ctritical bundaries#

  ############### alpha spending######
  c_cal3<-function(up_bry_0,rank,sigma){
    #e0 is the  exit boundary #
    #sigma is the covariance matrix, it has the rank of =rep_time x reptime#
    sigma1=sigma[1:rank,1:rank]
    up_bry=up_bry_0*sqrt(s[1:rank])
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rank),sigma=sigma1,algorithm=Miwa())[1]
    return(tmp)
  }

  ######################################
  cumu_alpha=c(rep(0,length(s)))
  for (i in 1:length(s)){
    cumu_alpha[i]=c_cal3(rep(c_bry,i),i,cov_mtx)
  }
  alpha_s=c(rep(0,length(s)))
  if(length(s)==1){
    alpha_s=cumu_alpha
  } else {
    alpha_s[1]=cumu_alpha[1]
    for (i in 2:length(s)){
      alpha_s[i]=cumu_alpha[i]-cumu_alpha[i-1]
    }
  }

  nom_alpha=1-stats::pnorm(c_bry)
  #################################################
  sum_nom_alpha=sum( nom_alpha)

  out=list(rep(c_bry,length(s)),nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
##############################
#pocock_bundry(1,0.025)
#pocock_bundry(c(0.5,1),0.025)
#pocock_bundry(c(1/3,2/3,1),0.025)
#pocock_bundry(c(0.3,0.7,1),0.025)
#pocock_bundry(c(1/4,1/3,2/3,1),0.025)
############################


#####################################################################
####             Part I.2: component fimctions                   ####
####          for Sample size and power calculation              ####
#####################################################################

###################################################
# This function calculates the required information time for the desired power
info_time<-function(mu,s,c_bry,beta,alpha){
  #This function calculates the infromation required for 1-beta power
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #c_bry is the critial booundary
  #1-beta is the desired power
  #mu is the effct size
  #alpha  is the target one-sided type I error#
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################

  ### exit boundary ####
  e_bry=sqrt(s)*c_bry

  ### the function calculates the exit probability###
  cal1<-function(time){
    up_bry=c_bry-sqrt(time)*mu
    tmp = 1-stats::pnorm(up_bry)
    return(tmp)
  }

  cal2<-function(time,sigma1){
    up_bry=sqrt(s)*c_bry-sqrt(time)*mu*s
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,length(s)),sigma=sigma1,abseps=1e-8,algorithm=Miwa())[1]
    return(tmp)
  }

  if (length(s)==1){
    if (mu==0) {total_time=NA}else{
      total_time=as.numeric(stats::uniroot(function(x) cal1(x)-1+beta,c(0,2^32),tol=1e-5)[1])
    }
  } else {
    if (mu==0) {total_time=NA}else{
      total_time=as.numeric(stats::uniroot(function(x) cal2(x,cov_mtx)-1+beta,c(0,2^32),tol=1e-5)[1])
    }
  }
  return(total_time)
}
########################
#info_time(mu,s,c_bry,beta,alpha)
#info_time(0.2,c(0.3,0.7,1),c(2.293075, 2.293075, 2.293075),0.1,0.025)
#info_time(0.2,1,1.96,0.1,0.025)

######################################


######################################
## import c_bry from OF_boundry, pococok_bunadry, or alsp_bundry##
##This function convets total reqired inofmation time to sample size for each arm##
## if the efficacy endpoint is a  binary variable##

####################################
time_power<-function(mu,s,c_bry,T,alpha){
  #This function calculates the  power for given information time T
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #c_bry is the critial booundary
  #1-beta is the desired power
  #mu is the effct size
  #alpha  is the target one-sided type I error#
  #### the next section constructs the covariance matrix ###
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx= matrix(rep(c(s),length(s)),nrow=length(s), ncol=length(s))
  cov_mtx[lower.tri(cov_mtx)]=t(cov_mtx)[lower.tri(cov_mtx)]
  #####################################################

  ### exit boundary ####
  e_bry=sqrt(s)*c_bry

  ### the function calculates the exit probability###
  cal1<-function(time,sigma1){
    up_bry=sqrt(s)*c_bry-sqrt(time)*mu*s
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,length(s)),sigma=sigma1,abseps=1e-8,algorithm=Miwa())[1]
    return(tmp)
  }

  cal2<-function(time) {
    up_bry=c_bry-sqrt(time)*mu*s
    tmp=1-stats::pnorm(up_bry)
    return(tmp)
  }

  if (length(s)==1) {
    power=cal2(T)
  } else {
    power=cal1(T,cov_mtx)
  }

  return(power)
}
########################
#power(mu,s,c_bry,T, alpha)
#time_power(0.2,c(0.3,0.7,1),c(2.293075, 2.293075, 2.293075),100,0.025)
#time_power(0.2,1,1.96,100,0.025)
##################################

#####################################################################
####            Part II sample size and power calculation      ####
#####################################################################

##############################################################################
#### Part II.1. Sample size and power clculation: Normal distributions   ####
##############################################################################


####################################
####     supeiority trials   #######
####################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################

#' Group sequential designs / Adaptive sequential designs: Normal distribution, Superiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_normal(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' s = c(1/3, 2/3, 1), alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_OF_boundary_normal<-function(mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is O'Brien-Fleming boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #mu=abs(mu_t-mu_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  mu=abs(mu_t-mu_c)
  if (mu==0){return(NA)}

  alpha=alpha_2/2
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  #######################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_normal(mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio){
#Sample_size_OF_boundary_normal(0.2,0,1,1,c(1/3,2/3,1),0.05,0.1,1)
#Sample_size_OF_boundary_normal(0.2,0,1,1,c(1/3,2/3,1),0.05,0.1,3)
#Sample_size_OF_boundary_normal(0.2,0.2,1,1,c(1/3,2/3,1),0.05,0.1,1)
##################################################################

#' Group sequential designs / Adaptive sequential designs: Normal distribution, Superiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param alpha_2 Type I error for the two-sided test
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_normal(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86, s = c(1/3, 2/3, 1),
#' sampsz_control = 100, rand_ratio = 1, alpha_2 = 0.05)
#'
#' @export
power_OF_boundary_normal<-function(mu_t,mu_c,std_t,std_c,s,sampsz_control,rand_ratio,alpha_2){
  ### This function calculates the critical boundary and power for given mu_t,mu_c,std_t,std_c,s,alpha_2,sampsz_control,rand_ratio
  # the critical boundary is O'Brien-Fleming boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #mu=abs(mu_t-mu_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta -s type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(mu_t-mu_c)

  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_OF_boundary_normal(mu_t,mu_c,std_t,std_c,s,sampsz_control,rand_ratio,alpha_2)
#power_OF_boundary_normal(0.2,0,1,1,c(1/3,2/3,1),534,1,0.05)
#power_OF_boundary_normal(0.2,0,1,1,c(1/3,2/3,1),356,3,0.05)
#################################################################

#######################
### Pocock boundary ###
#######################

###########################################################
#' Group sequential designs / Adaptive sequential designs: Normal distribution, Superiority trial, Sample size calculation, Pocock boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_normal(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' s = c(1/3, 2/3, 1), alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_pocock_boundary_normal<-function(mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is Pocock boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #mu=abs(mu_t-mu_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  mu=abs(mu_t-mu_c)
  if (mu==0){return(NA)}

  alpha=alpha_2/2
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ############################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_normal(mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio){
#Sample_size_pocock_boundary_normal(0.2,0,1,1,c(1/3,2/3,1),0.05,0.1,1)
#Sample_size_pocock_boundary_normal(0.2,0.2,1,1,c(1/3,2/3,1),0.05,0.1,1)
#####################################################################

#' Group sequential designs / Adaptive sequential designs: Normal distribution, Superiority trial, Power calculation, Pocock boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param alpha_2 Type I error for the two-sided test
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_normal(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' s = c(1/3, 2/3, 1), sampsz_control = 100, rand_ratio = 1, alpha_2 = 0.05)
#'
#' @export
power_pocock_boundary_normal<-function(mu_t,mu_c,std_t,std_c,s,sampsz_control,rand_ratio,alpha_2){
  ### This function calculates the critical boundary and power for given mu_t,mu_c,std_t,std_c,s,alpha_2,sampsz_control,rand_ratio
  # the critical boundary is Pocock boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #mu=abs(mu_t-mu_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta -s type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(mu_t-mu_c)

  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_pocock_boundary_normal(mu_t,mu_c,std_t,std_c,s,sampsz_test,rand_ratio,alpha_2)
#power_pocock_boundary_normal(0.2,0,1,1,c(1/3,2/3,1),605,1,0.05)
#################################################################

###############################
### Alpha spending boundary ###
###############################

####################################################
#' Group sequential designs / Adaptive sequential designs: Normal distribution, Superiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_normal(mu_t = 0.2, mu_c = 0, std_t = 1, std_c = 1, s = c(1/3,2/3,1),
#' alpha_2 = c(0.001,0.001,0.048), beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_alsp_boundary_normal<-function(mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is alpha-spending boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #mu=abs(mu_t-mu_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  mu=abs(mu_t-mu_c)
  if (mu==0){return(NA)}

  alpha=alpha_2/2
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ############################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}

####################################################################
#Sample_size_alsp_boundary_normal(mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio){
##Sample_size_alsp_boundary_normal(0.2,0,1,1,c(1/3,2/3,1),c(0.001,0.001,0.048),0.1,1)
##Sample_size_alsp_boundary_normal(0.2,0.2,1,1,c(1/3,2/3,1),c(0.001,0.001,0.048),0.1,1)
#####################################################################

#' Group sequential designs / Adaptive sequential designs: Normal distribution, Superiority trial, Power calculation, Alpha-spending boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param alpha_2 Type I error for the two-sided test
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_normal(mu_t = 0.2, mu_c = 0, std_t = 1, std_c = 1, s = c(1/3,2/3,1),
#' sampsz_control = 527, rand_ratio = 1, alpha_2 = c(0.001,0.001,0.048))
#'
#' @export
power_alsp_boundary_normal<-function(mu_t,mu_c,std_t,std_c,s,sampsz_control,rand_ratio,alpha_2){
  ### This function calculates the critical boundary and power for given mu_t,mu_c,std_t,std_c,s,alpha_2,sampsz_control,rand_ratio
  # the critical boundary is alpha-spending boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #mu=abs(mu_t-mu_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta -s type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(mu_t-mu_c)

  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_alsp_boundary_normal(mu_t,mu_c,std_t,std_c,s,sampsz_test,rand_ratio,alpha_2)
#power_alsp_boundary_normal(0.2,0,1,1,c(1/3,2/3,1),527,1,c(0.001,0.001,0.048))
#################################################################

#########################################
####     non-inferiority trials   #######
#########################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################
#' Group sequential designs / Adaptive sequential designs: Normal distribution, Non-inferiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_normal_NI(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' margin = 0.1, direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_OF_boundary_normal_NI<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio
  # the critical boundary is O'Brien-Fleming boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm  #mu_t is the mean of the tesrt arm
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ##############################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_normal_NI(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio){
#Sample_size_OF_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_OF_boundary_normal_NI(0.2,0,1,2,0.15,1,c(1/3,2/3,1),0.025,0.1,1)

#################################################################
#' Group sequential designs / Adaptive sequential designs: Normal distribution, Non-inferiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param alpha One-sided significance level (type I error rate)
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_normal_NI(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), sampsz_control = 100, rand_ratio = 1, alpha = 0.025)
#'
#' @export
power_OF_boundary_normal_NI<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,sampsz_control,rand_ratio,alpha){
  ### This function calculates the critical boundary and powerfor given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio for a non-inferiority trial
  # the critical boundary is O'Brien-Fleming boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm  #mu_t is the mean of the tesrt arm
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better

  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }

  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_OF_boundary_normal_NI(mu_t,mu_c,std_t,std_c,margin,direction,s,sampsz_control,rand_ratio,alpha)
#power_OF_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/3,2/3,1),175,1,0.025)
#power_OF_boundary_normal_NI(0.2,0,1,1,0.15,2,c(1/3,2/3,1),534,1,0.025)
################################################################


#######################
### Pocock boundary ###
#######################

###########################################################

#' Group sequential designs / Adaptive sequential designs: Normal distribution, Non-inferiority trial, Sample size calculation, Pocock boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_normal_NI(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' margin = 0.1, direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_pocock_boundary_normal_NI<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size  for given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio
  # the critical boundary is Pocock boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm  #mu_t is the mean of the tesrt arm
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  #########################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_normal_NI(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha_2,beta,rand_ratio)
#Sample_size_pocock_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_pocock_boundary_normal_NI(0.2,0,1,1,0.15,2,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_pocock_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/3,1),0.025,0.1,1)
#Sample_size_pocock_boundary_normal_NI(0.2,0,1,1,0.15,2,c(1/3,1),0.025,0.1,1)
#Sample_size_pocock_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/4,1/3,2/3,1),0.025,0.1,1)
##########################################################################

#' Group sequential designs / Adaptive sequential designs: Normal distribution, Non-inferiority trial, Power calculation, Pocock boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param alpha One-sided significance level (type I error rate)
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_normal_NI(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' margin = 0.1, direction = 0, s = c(1/3, 2/3, 1), sampsz_control = 100,
#' rand_ratio = 1, alpha = 0.025)
#'
#' @export
power_pocock_boundary_normal_NI<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,sampsz_control,rand_ratio,alpha){
  ### This function calculates the critical boundary and powerfor given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio for a non-inferiority trial
  # the critical boundary is Pocock boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm  #mu_t is the mean of the tesrt arm
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }

  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_pocock_boundary_normal_NI(mu_t,mu_c,std_t,std_c,margin,direction,s,sampsz_control,rand_ratio,alpha)
#power_pocock_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/3,2/3,1),198,1,0.025)
########################################################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################

#' Group sequential designs / Adaptive sequential designs: Normal distribution, Non-inferiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_normal_NI(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' margin = 0.1, direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_alsp_boundary_normal_NI<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given mu_t,mu_c,std_t,std_c,s,alpha_2,beta,rand_ratio for non-inferiority trial
  # the critical boundary is alpha-spending boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm  #mu_t is the mean of the tesrt arm
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  #############################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_normal_NI(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio)
#Sample_size_alsp_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1)
#Sample_size_alsp_boundary_normal_NI(0.2,0,1,1,0.15,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1)
#Sample_size_alsp_boundary_normal_NI(0,0.2,1,1,0.15,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1)
#Sample_size_alsp_boundary_normal_NI(0,0.2,1,1,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1)
################################################################

#' Group sequential designs / Adaptive sequential designs: Normal distribution, Non-inferiority trial, Power calculation, Alpha-spending boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param alpha One-sided significance level (type I error rate)
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_normal_NI(mu_t = 0.2, mu_c = 0, std_t = 1, std_c = 1, margin = 0.15,
#' direction = 1, s = c(1/3,2/3,1), sampsz_control = 534, rand_ratio = 1, alpha = c(0.001,0.001,0.023))
#'
#' @export
power_alsp_boundary_normal_NI<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,sampsz_control,rand_ratio,alpha){
  ### This function calculates the critical boundary and powerfor given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio for a non-inferiority trial
  # the critical boundary is alpha-spending boundary
  #mu_t is the mean of the test arm
  #mu_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #std_c>0 is the standard deviation for control arm
  #std_t>0 is the standard deviation for test arm  #mu_t is the mean of the tesrt arm
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_alsp_boundary_normal_NI(mu_t,mu_c,std_t,std_c,margin,direction,s,sampsz_control,rand_ratio,alpha_2)
#power_alsp_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/3,2/3,1),534,1,c(0.001,0.001,0.023))
#power_alsp_boundary_normal_NI(0.2,0,1,1,0.15,2,c(1/3,2/3,1),534,1,c(0.001,0.001,0.023))
#power_alsp_boundary_normal_NI(0.2,0,1,1,0.15,1,c(1/3,2/3,1),173,1,c(0.001,0.001,0.023))

########################################################################

###############################################################################
####     Part II.2 Power and Sample size clculation:Binary distributions    ####
###############################################################################
####################################
####     supeiority trials   #######
####################################

################################
### O'Brien_Fleming boundary ###
################################

######################################################

#' Group sequential designs / Adaptive sequential designs: Binary distribution, Superiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_binary_diff(p_t = 0.4, p_c = 0.3, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_OF_boundary_binary_diff<-function(p_t,p_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given p_t,p_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the rate of the test arm
  #p_c is the rate of the control arm
  #mu=abs(p_t-p_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(p_t-p_c)
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ##########################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_binary_diff(p_t,p_c,s,alpha_2,beta,rand_ratio)
#Sample_size_OF_boundary_binary_diff(0.2,0.6,c(1/3,2/3,1),0.05,0.1,1)
#Sample_size_OF_boundary_binary_diff(0.2,0.6,c(1/3,2/3,1),0.05,0.1,3)

#Sample_size_OF_boundary_binary_diff(0.2,0.2,c(1/3,2/3,1),0.05,0.1,1)
######################################################################
#' Group sequential designs / Adaptive sequential designs: Binary distribution, Superiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_binary_diff(p_t = 0.4, p_c = 0.3, s = c(1/3, 2/3, 1), alpha_2 = 0.05,
#' sampsz_control = 100, rand_ratio = 1)
#'
#' @export
power_OF_boundary_binary_diff<-function(p_t,p_c,s,alpha_2,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given p_t,p_c,s,alpha_2,sampsz_control,rand_ratio
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the rate of the test arm
  #p_c is the rate of the control arm
  #mu=abs(p_t-p_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #sampsz_control is the sample size for the control arm
  alpha=alpha_2/2
  mu=abs(p_t-p_c)
  r=rand_ratio
  T=sampsz_control/(p_t*(1-p_t)/r+p_c*(1-p_c))

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################

  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_OF_boundary_binary_diff(p_t,p_c,s,alpha_2,sampsz_control,rand_ratio)
#power_OF_boundary_binary_diff(0.2,0.6,c(1/3,2/3,1),0.05,27,1)
#power_OF_boundary_binary_diff(0.2,0.6,c(1/3,2/3,1),0.05,20,3)

################################################################

#######################
### Pocock boundary ###
#######################

########################################################
#' Group sequential designs / Adaptive sequential designs: Binary distribution, Superiority trial, Sample size calculation, Pocock boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_binary_diff(p_t = 0.4, p_c = 0.3, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_pocock_boundary_binary_diff<-function(p_t,p_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given p_t,p_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is Pocock boundary
  #p_t is the rate of the test arm
  #p_c is the rate of the control arm
  #mu=abs(p_t-p_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(p_t-p_c)
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ##########################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_binary_diff(p_t,p_c,s,alpha_2,beta,rand_ratio)
#Sample_size_pocock_boundary_binary_diff(0.2,0.6,c(1/3,2/3,1),0.05,0.1,1)
#Sample_size_pocock_boundary_binary_diff(0.2,0.2,c(1/3,2/3,1),0.05,0.1,1)
############################################################################

#' Group sequential designs / Adaptive sequential designs: Binary distribution, Superiority trial, Power calculation, Pocock boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_binary_diff(p_t = 0.4, p_c = 0.3, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, sampsz_control = 100, rand_ratio = 1)
#'
#' @export
power_pocock_boundary_binary_diff<-function(p_t,p_c,s,alpha_2,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given p_t,p_c,s,alpha_2,sampsz_control,rand_ratio
  # the critical boundary is Pocock boundary
  #p_t is the rate of the test arm
  #p_c is the rate of the control arm
  #mu=abs(p_t-p_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #sampsz_control is the sample size for the control arm
  alpha=alpha_2/2
  mu=abs(p_t-p_c)
  r=rand_ratio
  T=sampsz_control/(p_t*(1-p_t)/r+p_c*(1-p_c))

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################

  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_pocock_boundary_binary_diff(p_t,p_c,s,alpha_2,sampsz_control,rand_ratio)
#power_pocock_boundary_binary_diff(0.2,0.6,c(1/3,2/3,1),0.05,31,1)
########################################################################

###############################
### Alpha spending boundary ###
###############################

##########################################################
#' Group sequential designs / Adaptive sequential designs: Binary distribution, Superiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_binary_diff(p_t = 0.2, p_c = 0.6, s = c(1/3,2/3,1),
#' alpha_2 = c(0.001,0.001,0.048), beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_alsp_boundary_binary_diff<-function(p_t,p_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given p_t,p_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is alpha-spending boundary
  #p_t is the rate of the test arm
  #p_c is the rate of the control arm
  #mu=abs(p_t-p_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(p_t-p_c)
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ##########################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_binary_diff(p_t,p_c,s,alpha_2,beta,rand_ratio)
#Sample_size_alsp_boundary_binary_diff(0.2,0.6,c(1/3,2/3,1),c(0.001,0.001,0.048),0.1,1)
#Sample_size_alsp_boundary_binary_diff(0.2,0.2,c(1/3,2/3,1),c(0.001,0.001,0.048),0.1,1)
#######################################################################

#' Group sequential designs / Adaptive sequential designs: Binary distribution, Superiority trial, Power calculation, Alpha-spending boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_binary_diff(p_t = 0.2, p_c = 0.6, s = c(1/3,2/3,1),
#' alpha_2 = c(0.001,0.001,0.048), sampsz_control = 27, rand_ratio = 1)
#'
#' @export
power_alsp_boundary_binary_diff<-function(p_t,p_c,s,alpha_2,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given p_t,p_c,s,alpha_2,sampsz_control,rand_ratio
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the rate of the test arm
  #p_c is the rate of the control arm
  #mu=abs(p_t-p_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #sampsz_control is the sample size for the control arm
  alpha=alpha_2/2
  mu=abs(p_t-p_c)
  r=rand_ratio
  T=sampsz_control/(p_t*(1-p_t)/r+p_c*(1-p_c))

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_alsp_boundary_binary_diff(p_t,p_c,s,alpha_2,sampsz_control,rand_ratio)
#power_alsp_boundary_binary_diff(0.2,0.6,c(1/3,2/3,1),c(0.001,0.001,0.048),27,1)
########################################################################

#########################################
####     non-inferiority trials   #######
#########################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################
#' Group sequential designs / Adaptive sequential designs: Binary distribution, Non-inferiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger p_t is better, 0 = smaller p_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_binary_diff_NI(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_OF_boundary_binary_diff_NI<-function(p_t,p_c,margin, direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given p_t,p_c,margin, direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the mean of the test arm
  #p_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better

  r=rand_ratio
  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###########################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_binary_diff_NI(p_t,p_c,margin, direction,s,alpha_2,beta,rand_ratio)
#Sample_size_OF_boundary_binary_diff_NI(0.2,0.6,0.15,1,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_OF_boundary_binary_diff_NI(0.2,0.6,0.15,2,c(1/3,2/3,1),0.025,0.1,1)
#########################################################

#' Group sequential designs / Adaptive sequential designs: Binary distribution, Non-inferiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger p_t is better, 0 = smaller p_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_binary_diff_NI(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1)
#'
#' @export
power_OF_boundary_binary_diff_NI<-function(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given p_t,p_c,margin,direction,s,sampsz_control,rand_ratio,alpha for non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the mean of the test arm
  #p_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger mu is better
  #direction ne 1 if smaller p_t is better

  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  r=rand_ratio
  T=sampsz_control/(p_t*(1-p_t)/r+p_c*(1-p_c))

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################

  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_OF_boundary_binary_diff_NI(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio)
#power_OF_boundary_binary_diff_NI(0.2,0.6,0.15,1,c(1/3,2/3,1),0.025,27,1)
#power_OF_boundary_binary_diff_NI(0.2,0.6,0.15,2,c(1/3,2/3,1),0.025,15,1)


################################################

#######################
### Pocock boundary ###
#######################

########################################################

#' Group sequential designs / Adaptive sequential designs: Binary distribution, Non-inferiority trial, Sample size calculation, Pocock boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger p_t is better, 0 = smaller p_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_binary_diff_NI(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_pocock_boundary_binary_diff_NI<-function(p_t,p_c,margin, direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given p_t,p_c,margin, direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #p_t is the mean of the test arm
  #p_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better

  r=rand_ratio
  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###########################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_binary_diff_NI(p_t,p_c,margin, direction,s,alpha,beta,rand_ratio)
#Sample_size_pocock_boundary_binary_diff_NI(0.2,0.6,0.15,1,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_pocock_boundary_binary_diff_NI(0.2,0.6,0.15,2,c(1/3,2/3,1),0.025,0.1,1)
##################################################

#' Group sequential designs / Adaptive sequential designs: Binary distribution, Non-inferiority trial, Power calculation, Pocock boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger p_t is better, 0 = smaller p_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_binary_diff_NI(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1)
#'
#' @export
power_pocock_boundary_binary_diff_NI<-function(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given p_t,p_c,margin,direction,s,sampsz_control,rand_ratio,alpha for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #p_t is the mean of the test arm
  #p_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger mu is better
  #direction ne 1 if smaller p_t is better

  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }

  r=rand_ratio
  T=sampsz_control/(p_t*(1-p_t)/r+p_c*(1-p_c))

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################

  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_pocock_boundary_binary_diff_NI(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio)
#power_pocock_boundary_binary_diff_NI(0.2,0.6,0.15,1,c(1/3,2/3,1),0.025,16,1)
#power_pocock_boundary_binary_diff_NI(0.2,0.6,0.15,2,c(1/3,2/3,1),0.025,16,1)


##################################################################


#################################
### Alpha -spending boundary#####
#################################

###########################################################

#' Group sequential designs / Adaptive sequential designs: Binary distribution, Non-inferiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger p_t is better, 0 = smaller p_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_binary_diff_NI(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_alsp_boundary_binary_diff_NI<-function(p_t,p_c,margin, direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given p_t,p_c,margin, direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #p_t is the mean of the test arm
  #p_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  r=rand_ratio

  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###########################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_binary_diff_NI(p_t,p_c,margin, direction,s,alpha,beta,rand_ratio)
#Sample_size_alsp_boundary_binary_diff_NI(0.2,0.6,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1)
#Sample_size_alsp_boundary_binary_diff_NI(0.2,0.6,0.15,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1)
################################################################

#' Group sequential designs / Adaptive sequential designs: Binary distribution, Non-inferiority trial, Power calculation, Alpha-spending boundary
#'
#' @param p_t Mean of the test arm
#' @param p_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger p_t is better, 0 = smaller p_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_binary_diff_NI(p_t = 0.2, p_c = 0.6, margin = 0.15, direction = 1,
#' s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), sampsz_control = 14, rand_ratio = 1)
#'
#' @export
power_alsp_boundary_binary_diff_NI<-function(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given p_t,p_c,margin,direction,s,sampsz_control,rand_ratio,alpha for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #p_t is the mean of the test arm
  #p_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger mu is better
  #direction ne 1 if smaller p_t is better
  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }

  r=rand_ratio
  T=sampsz_control/(p_t*(1-p_t)/r+p_c*(1-p_c))

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_alsp_boundary_binary_diff_NI(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio){
#power_alsp_boundary_binary_diff_NI(0.2,0.6,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),14,1)
#power_alsp_boundary_binary_diff_NI(0.2,0.6,0.15,2,c(1/3,2/3,1),c(0.001,0.001,0.023),14,1)

##############################################################################
#### Part II.3. Sample size and power clculation: Survivial analysis     ####
##############################################################################

####################################
####     supeiority trials   #######
####################################

################################
### O'Brien_Fleming boundary ###
################################

####################################################################
#' Group sequential designs / Adaptive sequential designs: Survival analysis, Superiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_survival_Schoenfeld(HR = 0.7, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_OF_boundary_survival_Schoenfeld<-function(HR,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given HR,s,alpha_2,beta,rand_ratio
  # the critical boundary is O'Brien-Fleming boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(log(HR))
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ####################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### sample size for specified distribution ####
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_tot=ceiling(info_time/(p_c*p_t))

  ######################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(n_tot,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("total number of events","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_survival_Schoenfeld(HR,s,alpha_2,beta,rand_ratio)
#Sample_size_OF_boundary_survival_Schoenfeld(0.8,c(1/3,2/3,1),0.05,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld(0.6,c(0.5,1),0.05,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld(1,c(0.5,1),0.05,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld(0.8,1,0.05,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),0.05,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),0.05,0.1,3)

#Sample_size_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),(1-0.54),0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),0.05,0.1,3)
#Sample_size_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),(1-0.4897),0.1,3)
#Sample_size_OF_boundary_survival_Schoenfeld(HR=0.62,s=c(0.6,1),0.05,0.1,1)
#############################################################################

#' Group sequential designs / Adaptive sequential designs: Survival analysis, Superiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param events_tot Total number of events
#' @param alpha_2 Type I error for the two-sided test
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_survival_Schoenfeld(HR = 0.7, s = c(1/3, 2/3, 1), events_tot = 100,
#' alpha_2 = 0.05, rand_ratio = 1)
#'
#' @export
power_OF_boundary_survival_Schoenfeld<-function(HR,s,events_tot,alpha_2,rand_ratio){
  ### This function calculates the critical boundary and power for given HR,s,events_tot,alpha_2,rand_ratio
  # the critical boundary is O'btien-Fleming boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #events_tot is the combined number of events in both the test arm and the control arm
  alpha=alpha_2/2
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  events_control=events_tot/(1+rand_ratio)
  T=events_tot*(p_c*p_t)
  mu=abs(log(HR))

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
########################
#power_OF_boundary_survival_Schoenfeld(HR,s,events_tot,alpha_2,rand_ratio)
#power_OF_boundary_survival_Schoenfeld(0.8,c(1/3,2/3,1),858,0.05,1)
#power_OF_boundary_survival_Schoenfeld(0.8,1,845,0.05,1)
#power_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),138,0.05,1)
#power_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),276,0.05,3)
#power_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),335,0.05,1)
#power_OF_boundary_survival_Schoenfeld(HR=0.7,s=c(0.6,1),118,0.05,1)
#power_OF_boundary_survival_Schoenfeld(HR=0.62,s=c(0.6,1),376,0.05,3)

#######################
### Pocock boundary ###
#######################

###############################################
#' Group sequential designs / Adaptive sequential designs: Survival analysis, Superiority trial, Sample size calculation, Pocock boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_survival_Schoenfeld(HR = 0.7, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_pocock_boundary_survival_Schoenfeld<-function(HR,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given HR,s,alpha_2,beta,rand_ratio
  # the critical boundary is Pocock boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(log(HR))
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  #####################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### sample size for specified distribution ####
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_tot=ceiling(info_time/(p_c*p_t))

  ######################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(n_tot,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("total number of events","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_survival_Schoenfeld(HR,s,alpha_2,beta,rand_ratio)
#Sample_size_pocock_boundary_survival_Schoenfeld(0.8,1,0.05,0.1,1)
#Sample_size_pocock_boundary_survival_Schoenfeld(0.8,c(1/3,2/3,1),0.05,0.1,1)
#Sample_size_pocock_boundary_survival_Schoenfeld(1,c(1/3,2/3,1),0.05,0.1,1)
#############################################################################

#' Group sequential designs / Adaptive sequential designs: Survival analysis, Superiority trial, Power calculation, Pocock boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param events_tot Total number of events
#' @param alpha_2 Type I error for the two-sided test
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_survival_Schoenfeld(HR = 0.7, s = c(1/3, 2/3, 1), events_tot = 100,
#' alpha_2 = 0.05, rand_ratio = 1)
#'
#' @export
power_pocock_boundary_survival_Schoenfeld<-function(HR,s,events_tot,alpha_2,rand_ratio){
  ### This function calculates the critical boundary and power for given HR,s,sampsz_tot,alpha_2,rand_ratio
  # the critical boundary is Pocock boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #events_tot is the combined number of events in both the test arm and the control arm
  alpha=alpha_2/2
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  T=events_tot*(p_c*p_t)
  mu=abs(log(HR))
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
########################
#power_pocock_boundary_survival_Schoenfeld(HR,s,events_tot,alpha_2,rand_ratio)
#power_pocock_boundary_survival_Schoenfeld(0.8,c(1/3,2/3,1),972,0.05,1)
#power_pocock_boundary_survival_Schoenfeld(0.8,1,845,0.05,1)
################################################################

###############################
### Alpha spending boundary ###
###############################

#######################################################
#' Group sequential designs / Adaptive sequential designs: Survival analysis, Superiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_survival_Schoenfeld(HR = 0.59, s = c(0.596,1),
#' alpha_2 = c(0.001,0.0249), beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_alsp_boundary_survival_Schoenfeld<-function(HR,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given HR,s,alpha_2,beta,rand_ratio
  # the critical boundary is alpha-spending boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(log(HR))
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ######################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_tot=ceiling(info_time/(p_c*p_t))
  ######################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(n_tot,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("total number of events","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_survival_Schoenfeld(HR,s,alpha_2,beta,rand_ratio)
#Sample_size_alsp_boundary_survival_Schoenfeld(0.59,c(0.596,1),c(0.001,0.0249),0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld(0.8,1,0.05,0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld(0.8,c(0.5,1),c(0.001,0.049),0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld(1,c(0.5,1),c(0.001,0.049),0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld(HR=0.62,s=c(0.6,1),alpha_2=c(0.001,0.049),0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld(HR=0.62,s=c(0.6,1),alpha_2=c(0.008,0.042),0.1,1)

#####################################################
#' Group sequential designs / Adaptive sequential designs: Survival analysis, Superiority trial, Power calculation, Alpha-spending boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param events_tot Total number of events
#' @param alpha_2 Type I error for the two-sided test
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_survival_Schoenfeld(HR = 0.8, s = c(0.5,1), events_tot = 841,
#' alpha_2 = c(0.001,0.0499), rand_ratio = 1)
#'
#' @export
power_alsp_boundary_survival_Schoenfeld<-function(HR,s,events_tot,alpha_2,rand_ratio){
  ### This function calculates the critical boundary and power for given HR,s,events_tot,alpha_2,rand_ratio
  # the critical boundary is alpha-spending boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #events_tot is the combined number of events in both the test arm and the control arm
  alpha=alpha_2/2
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  T=events_tot*(p_c*p_t)
  mu=abs(log(HR))

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
########################
#power_alsp_boundary_survival_Schoenfeld(HR,s,events_tot,alpha_2,rand_ratio)
#power_alsp_boundary_survival_Schoenfeld(0.8,c(0.5,1),841,c(0.001,0.0499),1)
#power_alsp_boundary_survival_Schoenfeld(0.8,1,845,0.05,1)
################################################################

#########################################
####     non-inferiority trials   #######
#########################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################
#' Group sequential designs / Adaptive sequential designs: Survival analysis, Non-inferiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_survival_Schoenfeld_NI(HR = 0.7, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_OF_boundary_survival_Schoenfeld_NI<-function(HR,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given HR,margin,direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is O'brien-Fleming boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better

  mu_0=log(HR/margin)

  if ((direction==1)&&(mu_0<=0)){return(NA)}
  if ((direction !=1)&&(mu_0>=0)){return(NA)}

  mu=abs(mu_0)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ####################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### sample size for specified distribution ####
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_tot=ceiling(info_time/(p_c*p_t))

  ######################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(n_tot,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("total number of events","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_survival_Schoenfeld_NI(HR,margin,direction,s,alpha,beta,rand_ratio)
#Sample_size_OF_boundary_survival_Schoenfeld_NI(0.8,1.1,1,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld_NI(0.8,1.1,2,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld_NI(0.6,1.1,1,c(0.5,1),0.025,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld_NI(0.8,1.1,1,1,0.025,0.1,1)
#Sample_size_OF_boundary_survival_Schoenfeld_NI(0.8,1.1,2,1,0.025,0.1,1)

###################################################
#' Group sequential designs / Adaptive sequential designs: Survival analysis, Non-inferiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param events_tot Total number of events
#' @param alpha One-sided significance level (type I error rate)
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_survival_Schoenfeld_NI(HR = 0.7, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), events_tot = 100, alpha = 0.025, rand_ratio = 1)
#'
#' @export
power_OF_boundary_survival_Schoenfeld_NI<-function(HR,margin,direction,s,events_tot,alpha,rand_ratio){
  ### This function calculates the critical boundary and power for given HR,margin,direction,s,sampsz_tot,alpha,rand_ratio for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>0
  #events_tot is the combined number of events in both the test arm and the control arm
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  T=events_tot*(p_c*p_t)
  mu_0=log(HR/margin)
  if (direction==1){
    mu=mu_0
  }else{
    mu=-mu_0
  }

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)

}
########################
#power_OF_boundary_survival_Schoenfeld_NI(HR,margin,direction,s,events_tot,alpha,rand_ratio)
#power_OF_boundary_survival_Schoenfeld_NI(0.8,1.1,1,c(1/3,2/3,1),422,0.025,1)
#power_OF_boundary_survival_Schoenfeld_NI(0.8,1.1,2,c(1/3,2/3,1),422,0.025,1)
#power_OF_boundary_survival_Schoenfeld_NI(0.8,1.1,1,1,423,0.025,1)
#power_OF_boundary_survival_Schoenfeld_NI(0.8,1.1,2,1,423,0.025,1)
#power_OF_boundary_survival_Schoenfeld_NI(0.8,1,2,1,845,0.025,1)
################################################################

#######################
### Pocock boundary ###
#######################

###############################################
#' Group sequential designs / Adaptive sequential designs: Survival analysis, Non-inferiority trial, Sample size calculation, Pocock boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_survival_Schoenfeld_NI(HR = 0.7, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_pocock_boundary_survival_Schoenfeld_NI<-function(HR,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given HR,margin,direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  mu_0=log(HR/margin)

  if ((direction==1)&&(mu_0<=0)){return(NA)}
  if ((direction!=1)&&(mu_0>=0)){return(NA)}

  mu=abs(mu_0)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ##############################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_tot=ceiling(info_time/(p_c*p_t))

  ######################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(n_tot,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("total number of events","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_survival_Schoenfeld_NI(HR,margin,direction,s,alpha,beta,rand_ratio)
#Sample_size_pocock_boundary_survival_Schoenfeld_NI(0.8,1.1,1,1,0.025,0.1,1)
#Sample_size_pocock_boundary_survival_Schoenfeld_NI(0.8,1.1,2,1,0.025,0.1,1)
#Sample_size_pocock_boundary_survival_Schoenfeld_NI(0.8,1.1,1,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_pocock_boundary_survival_Schoenfeld_NI(0.8,1.1,2,c(1/3,2/3,1),0.025,0.1,1)
##################################################

#' Group sequential designs / Adaptive sequential designs: Survival analysis, Non-inferiority trial, Power calculation, Pocock boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param events_tot Total number of events
#' @param alpha One-sided significance level (type I error rate)
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_survival_Schoenfeld_NI(HR = 0.7, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), events_tot = 100, alpha = 0.025, rand_ratio = 1)
#'
#' @export
power_pocock_boundary_survival_Schoenfeld_NI<-function(HR,margin,direction,s,events_tot,alpha,rand_ratio){
  ### This function calculates the critical boundary and power for given HR,margin,direction,s,sampsz_tot,alpha,rand_ratio for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>0
  #events_tot is the combined number of events in both the test arm and the control arm
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  T=events_tot*(p_c*p_t)
  ####################
  mu_0=log(HR/margin)
  if (direction==1){
    mu=mu_0
  }else{
    mu=-mu_0
  }

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
########################
#power_pocock_boundary_survival_Schoenfeld_NI(HR,margin,diredtion,s,events_tot,alpha,rand_ratio)
#power_pocock_boundary_survival_Schoenfeld_NI(0.8,1.1,1,c(1/3,2/3,1),477,0.025,1)
#power_pocock_boundary_survival_Schoenfeld_NI(0.8,1.1,2,c(1/3,2/3,1),477,0.025,1)
#power_pocock_boundary_survival_Schoenfeld_NI(0.8,1.1,1,1,415,0.025,1)
#power_pocock_boundary_survival_Schoenfeld_NI(0.8,1.1,2,1,415,0.025,1)

###################################################


#################################
### Alpha -spending boundary#####
#################################

###########################################################

#' Group sequential designs / Adaptive sequential designs: Survival analysis, Non-inferiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_survival_Schoenfeld_NI(HR = 0.7, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_alsp_boundary_survival_Schoenfeld_NI<-function(HR,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given HR,margin,direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  mu_0=log(HR/margin)

  if ((direction==1)&&(mu_0<=0)){return(NA)}
  if ((direction!=1)&&(mu_0>=0)){return(NA)}

  mu=abs(mu_0)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ######################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_tot=ceiling(info_time/(p_c*p_t))
  ######################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(n_tot,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("total number of events","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_survival_Schoenfeld_NI(HR,margin,s,alpha,beta,rand_ratio)
#Sample_size_alsp_boundary_survival_Schoenfeld_NI(0.59,1,1,c(0.596,1),c(0.005,0.02),0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld_NI(0.59,1,2,c(0.596,1),c(0.005,0.02),0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld_NI(0.8,1.1,1,1,0.025,0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld_NI(0.8,1.1,2,1,0.025,0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld_NI(0.8,1.1,1,c(0.5,1),c(0.001,0.024),0.1,1)
#Sample_size_alsp_boundary_survival_Schoenfeld_NI(0.8,1.1,2,c(0.5,1),c(0.001,0.024),0.1,1)

################################################################
#' Group sequential designs / Adaptive sequential designs: Survival analysis, Non-inferiority trial, Power calculation, Alpha-spending boundary
#'
#' @param HR Hazards ratio of test arm vs. the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param events_tot Total number of events
#' @param alpha One-sided significance level (type I error rate)
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_survival_Schoenfeld_NI(HR = 0.8, margin = 1.1, direction = 1, s = c(0.5,1),
#' events_tot = 416, alpha = c(0.001,0.024), rand_ratio = 1)
#'
#' @export
power_alsp_boundary_survival_Schoenfeld_NI<-function(HR,margin,direction,s,events_tot,alpha,rand_ratio){
  ### This function calculates the critical boundary and power for given HR,margin,direction,s,sampsz_tot,alpha,rand_ratio for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #HR is the hazards ratio of test arm vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>0
  #events_tot is the combined number of events in both the test arm and the control arm
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  T=events_tot*(p_c*p_t)
  ####################
  mu_0=log(HR/margin)
  if (direction==1){
    mu=mu_0
  }else{
    mu=-mu_0
  }
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","reference")
  return(out)
}
########################
#power_alsp_boundary_survival_Schoenfeld_NI(HR,margin,direction,s,events_tot,alpha,rand_ratio)
#power_alsp_boundary_survival_Schoenfeld_NI(0.8,1.1,1,c(0.5,1),416,c(0.001,0.024),1)
#power_alsp_boundary_survival_Schoenfeld_NI(0.8,1.1,2,c(0.5,1),416,c(0.001,0.024),1)
#power_alsp_boundary_survival_Schoenfeld_NI(0.8,1,1,1,845,0.025,1)
#power_alsp_boundary_survival_Schoenfeld_NI(0.8,1,2,1,415,0.025,1)



################################################################

############################################################
####                        Part III.4.                  ###
##                     Sample size and power calculation ###
####                     Poisson distributions           ###
############################################################

###################################
####    supeiority trials   #######
###################################

################################
### O'Brien_Fleming boundary ###
################################

###################################################
#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Superiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_poisson(lambda_t = 1.1, lambda_c = 1.0, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_OF_boundary_poisson<-function(lambda_t,lambda_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given lambda_t,lambda_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is O'Brien-Fleming boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #mu=abs(lambda_t -lambda_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(lambda_t-lambda_c)
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_poisson(lambda_t,lambda_c,s,alpha_2,beta,rand_ratio)
#Sample_size_OF_boundary_poisson(0.2,0.4,c(1/3,2/3,1),0.05,0.1,1)
#Sample_size_OF_boundary_poisson(0.2,0.2,c(1/3,2/3,1),0.05,0.1,1)
###################################################################

#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Superiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_poisson(lambda_t = 1.1, lambda_c = 1.0, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, sampsz_control = 100, rand_ratio = 1)
#'
#' @export
power_OF_boundary_poisson<-function(lambda_t,lambda_c,s,alpha_2,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given lambda_t,lambda_c,s,sampsz_control,rand_ratio,alpha_2
  # the critical boundary is O'brien-Fleming boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #mu=abs(lambda_t -lambda_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #sampsz_control is the sample size of the control arm
  alpha=alpha_2/2
  mu=abs(lambda_t-lambda_c)
  r=rand_ratio
  T=sampsz_control/(lambda_t/r+lambda_c)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_OF_boundary_poisson(lambda_t,lambda_c,s,alpha_2,sampsz_control,rand_ratio){
#power_OF_boundary_poisson(0.2,0.4,c(1/3,2/3,1),0.05,161,1)
################################################################

#######################
### Pocock boundary ###
#######################

#######################################################
#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Superiority trial, Sample size calculation, Pocock boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_poisson(lambda_t = 0.2, lambda_c = 0.4, s = c(1/3,2/3,1),
#' alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_pocock_boundary_poisson<-function(lambda_t,lambda_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given lambda_t,lambda_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is Pocock boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #mu=abs(lambda_t -lambda_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(lambda_t-lambda_c)
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_pcocok_boundary_poisson(lambda_t,lambda_c,s,alpha_2,beta,rand_ratio)
#Sample_size_pocock_boundary_poisson(0.2,0.4,c(1/3,2/3,1),0.05,0.1,1)
#Sample_size_pocock_boundary_poisson(0.2,0.2,c(1/3,2/3,1),0.05,0.1,1)
###################################################################

#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Superiority trial, Power calculation, Pocock boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_poisson(lambda_t = 1.1, lambda_c = 1.0, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, sampsz_control = 100, rand_ratio = 1)
#'
#' @export
power_pocock_boundary_poisson<-function(lambda_t,lambda_c,s,alpha_2,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given lambda_t,lambda_c,s,sampsz_control,rand_ratio,alpha_2
  # the critical boundary is Pocock boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #mu=abs(lambda_t -lambda_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #sampsz_control is the sample size of the control arm
  alpha=alpha_2/2
  mu=abs(lambda_t-lambda_c)
  r=rand_ratio
  T=sampsz_control/(lambda_t/r+lambda_c)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_pocock_boundary_poisson(lambda_t,lambda_c,s,alpha_2,sampsz_control,rand_ratio)
#power_pocock_boundary_poisson(0.2,0.4,c(1/3,2/3,1),0.05,182,1)
########################################################################

###############################
### Alpha spending boundary ###
###############################

################################################################
#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Superiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_poisson(lambda_t = 0.2, lambda_c = 0.4, s = c(1/3,2/3,1),
#' alpha_2 = c(0.001,0.001,0.048), beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_alsp_boundary_poisson<-function(lambda_t,lambda_c,s,alpha_2,beta,rand_ratio){
  ### This function calculates the critical boundary and sample size for given lambda_t,lambda_c,s,alpha_2,beta,rand_ratio
  # the critical boundary is alpha-spending boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #mu=abs(lambda_t -lambda_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  alpha=alpha_2/2
  mu=abs(lambda_t-lambda_c)
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_poisson(lambda_t,lambda_c,s,alpha_2,beta,rand_ratio)
#Sample_size_alsp_boundary_poisson(0.2,0.4,c(1/3,2/3,1),c(0.001,0.001,0.048),0.1,1)
#Sample_size_alsp_boundary_poisson(0.2,0.2,c(1/3,2/3,1),c(0.001,0.001,0.048),0.1,1)
###################################################################

#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Superiority trial, Power calculation, Alpha-spending boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_poisson(lambda_t = 0.2, lambda_c = 0.4, s = c(1/3,2/3,1),
#' alpha_2 = c(0.001,0.001,0.048), sampsz_control = 158, rand_ratio = 1)
#'
#' @export
power_alsp_boundary_poisson<-function(lambda_t,lambda_c,s,alpha_2,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given lambda_t,lambda_c,s,sampsz_control,rand_ratio,alpha_2
  # the critical boundary is alpha-spending  boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #mu=abs(lambda_t -lambda_c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #sampsz_control is the sample size of the control arm
  alpha=alpha_2/2
  mu=abs(lambda_t-lambda_c)
  r=rand_ratio
  T=sampsz_control/(lambda_t/r+lambda_c)

  ##### Step 1.1: boundary calcualtion -alsp boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_alsp_boundary_poisson(lambda_t,lambda_c,s,alpha_2,sampsz_control,rand_ratio)
#power_alsp_boundary_poisson(0.2,0.4,c(1/3,2/3,1),c(0.001,0.001,0.048),158,1)
########################################################################

#########################################
####     non-inferiority trials   #######
#########################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################
#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Non-inferiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_poisson_NI(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_OF_boundary_poisson_NI<-function(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is O'brien-Fleming boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better

  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_poisson_NI(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio)
#Sample_size_OF_boundary_poisson_NI(0.2,0.4,0.15,1,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_OF_boundary_poisson_NI(0.2,0.4,0.15,2,c(1/3,2/3,1),0.025,0.1,1)

#####################################################
#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Non-inferiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_poisson_NI(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1)
#'
#' @export
power_OF_boundary_poisson_NI<-function(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given lambda_t,lambda_c,margin,direction,s,sampsz_control,rand_ratio,alpha for non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size of the control arm
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  r=rand_ratio
  T=sampsz_control/(lambda_t/r+lambda_c)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_OF_boundary_poisson_NI(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio)
#power_OF_boundary_poisson_NI(0.2,0.4,0.15,1,c(1/3,2/3,1),0.025,53,1)
#power_OF_boundary_poisson_NI(0.2,0.4,0.15,2,c(1/3,2/3,1),0.025,53,1)
######################################################


#######################
### Pocock boundary ###
#######################

#######################################################
#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Non-inferiority trial, Sample size calculation, Pocock boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_poisson_NI(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_pocock_boundary_poisson_NI<-function(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_poisson_NI(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio)
#Sample_size_pocock_boundary_poisson_NI(0.2,0.4,0.15,1,c(1/3,2/3,1),0.025,0.1,1)
#Sample_size_pocock_boundary_poisson_NI(0.2,0.4,0.15,2,c(1/3,2/3,1),0.025,0.1,1)
###################################################

#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Non-inferiority trial, Power calculation, Pocock boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_poisson_NI(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1)
#'
#' @export
power_pocock_boundary_poisson_NI<-function(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given lambda_t,lambda_c,margin,direction,s,sampsz_control,rand_ratio,alpha for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size of the control arm
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better

  if (direction==1){
    mu=abs(lambda_t-lambda_c+margin)
  }else{
    mu=abs(lambda_t-lambda_c-margin)
  }
  r=rand_ratio
  T=sampsz_control/(lambda_t/r+lambda_c)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_pocock_boundary_poisson_NI(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio)
#power_pocock_boundary_poisson_NI(0.2,0.4,0.15,1,c(1/3,2/3,1),0.025,60,1)
#power_pocock_boundary_poisson_NI(0.2,0.4,0.15,2,c(1/3,2/3,1),0.025,60,1)
##############################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################
#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Non-inferiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_poisson_NI(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
Sample_size_alsp_boundary_poisson_NI<-function(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio){
  ### This function calculates the critical boundary and saple size for given lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
   if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  if (mu<=0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(info_time*ratio)
  n_t=n_c*rand_ratio
  ######################
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_poisson_NI(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio)
#Sample_size_alsp_boundary_poisson_NI(0.2,0.4,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1)
#Sample_size_alsp_boundary_poisson_NI(0.2,0.4,0.15,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1)

#########################################################################
#' Group sequential designs / Adaptive sequential designs: Poisson distribution, Non-inferiority trial, Power calculation, Alpha-spending boundary
#'
#' @param lambda_t Mean of the test arm
#' @param lambda_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_poisson_NI(lambda_t = 0.2, lambda_c = 0.4, margin = 0.15, direction = 1,
#' s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), sampsz_control = 52, rand_ratio = 1)
#'
#' @export
power_alsp_boundary_poisson_NI<-function(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio){
  ### This function calculates the critical boundary and power for given lambda_t,lambda_c,margin,direction,s,sampsz_control,rand_ratio,alpha for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #lambda_t is the mean of the test arm
  #lambda_c is the mean of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size of the control arm
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better

  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  r=rand_ratio
  T=sampsz_control/(lambda_t/r+lambda_c)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha")
  return(out)
}
########################
#power_alsp_boundary_poisson_NI(lambda_t,lambda_c,margin,direction,s,sampsz_control,rand_ratio,alpha)
#power_alsp_boundary_poisson_NI(0.2,0.4,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),52,1)
#power_alsp_boundary_poisson_NI(0.2,0.4,0.15,2,c(1/3,2/3,1),c(0.001,0.001,0.023),52,1)
#power_alsp_boundary_poisson_NI(0.2,0.4,0,2,c(1/3,2/3,1),c(0.001,0.001,0.023),161,1)


########################################################################################
#### Part II.5. Sample size and power clculation: Negative binomial distributions   ####
########################################################################################

####################################
####     supeiority trials   #######
####################################

################################
### O'Brien_Fleming boundary ###
################################

########################################################
#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Superiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_neg_binomial(r_t = 0.2, r_c = 0.4, s = c(1/3,2/3,1), alpha_2 = 0.05,
#' beta = 0.1, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
Sample_size_OF_boundary_neg_binomial<-function(r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and sample size for given r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  alpha=alpha_2/2
  r_0=r_c
  r_1=r_t
  theta=rand_ratio
  mu=abs(log(r_1/r_0))
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(V)
  n_c=ceiling(info_time*V)
  n_t=n_c*rand_ratio
  ######################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_neg_binomial(r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa)
#Sample_size_OF_boundary_neg_binomial(0.2,0.4,c(1/3,2/3,1),0.05,0.1,1,1,1)
#Sample_size_OF_boundary_neg_binomial(0.2,0.2,c(1/3,2/3,1),0.05,0.1,1,1,1)
#Sample_size_OF_boundary_neg_binomial(1.2,0.84,1,0.05,0.1,1,1,1)
#Sample_size_OF_boundary_neg_binomial(r_t=0.85,r_c=1.25,s=c(0.5,1),alpha_2=0.05,beta=0.2,rand_ratio=1,nu_t=1,kappa=1)
#Sample_size_OF_boundary_neg_binomial(r_t=0.85,r_c=1.25,s=c(0.5,1),alpha_2=0.05,beta=0.2,rand_ratio=2,nu_t=1,kappa=1)
######################################################################

#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Superiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_neg_binomial(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, sampsz_control = 100, rand_ratio = 1, nu_t = 6, kappa = 1)
#'
#' @export
power_OF_boundary_neg_binomial<-function(r_t,r_c,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and power for given r_t,r_c,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  #sampsz_control is the sample size of the control arm
  alpha=alpha_2/2
  r_0=r_c
  r_1=r_t
  theta=rand_ratio
  mu=abs(log(r_1/r_0))

  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  T=sampsz_control/V

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
########################
#power_OF_boundary_neg_binomial(r_t,r_c,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa)
#power_OF_boundary_neg_binomial(1,0.8,1,0.05,897,1,1,1)
#power_OF_boundary_neg_binomial(1,0.8,c(1/3,2/3,1),0.05,898,1,1,1)
#power_OF_boundary_neg_binomial(0.2,0.4,c(1/3,2/3,1),0.05,212,1,1,1)
#power_OF_boundary_neg_binomial(1.2,0.84,1,0.05,333,1,1,1)
################################################################

#######################
### Pocock boundary ###
#######################

#################################################
#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Superiority trial, Sample size calculation, Pocock boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_neg_binomial(r_t = 0.2, r_c = 0.4, s = c(1/3,2/3,1),
#' alpha_2 = 0.05, beta = 0.1, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
Sample_size_pocock_boundary_neg_binomial<-function(r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and sample size for given r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa
  # the critical boundary is Pocock boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  alpha=alpha_2/2
  r_0=r_c
  r_1=r_t
  theta=rand_ratio
  mu=abs(log(r_1/r_0))
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(V)
  n_c=ceiling(info_time*V)
  n_t=n_c*rand_ratio
  ######################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_neg_binomial(r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa)
#Sample_size_pocock_boundary_neg_binomial(0.2,0.4,c(1/3,2/3,1),0.05,0.1,1,1,1)
#Sample_size_pocock_boundary_neg_binomial(0.2,0.2,c(1/3,2/3,1),0.05,0.1,1,1,1)
#Sample_size_pocock_boundary_neg_binomial(1.2,0.84,1,0.05,0.1,1,1,1)
######################################################################

#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Superiority trial, Power calculation, Pocock boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_neg_binomial(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, s = c(1/3, 2/3, 1),
#' alpha_2 = 0.05, sampsz_control = 100, rand_ratio = 1, nu_t = 6, kappa = 1)
#'
#' @export
power_pocock_boundary_neg_binomial<-function(r_t,r_c,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and power for given r_t,r_c,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  #sampsz_control is the sample size of the control arm
  alpha=alpha_2/2
  r_0=r_c
  r_1=r_t
  theta=rand_ratio
  mu=abs(log(r_1/r_0))
  # n_c=(((a+b)/mu)**2)*(lambda_t/r+lambda_c);

  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  T=sampsz_control/V

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
########################
#power_pocock_boundary_neg_binomial(r_t,r_c,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa)
#power_pocock_boundary_neg_binomial(1,0.8,1,0.05,897,1,1,1)
#power_pocock_boundary_neg_binomial(1,0.8,c(1/3,2/3,1),0.05,897,1,1,1)
#power_pocock_boundary_neg_binomial(0.2,0.4,c(1/3,2/3,1),0.05,240,1,1,1)
################################################################

###############################
### Alpha spending boundary ###
###############################

#############################################################
#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Superiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_neg_binomial(r_t = 0.2, r_c = 0.4, s = c(1/3,2/3,1),
#' alpha_2 = c(0.001,0.001,0.048), beta = 0.1, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
Sample_size_alsp_boundary_neg_binomial<-function(r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and sample size for given r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa
  # the critical boundary is alpha-spending boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  alpha=alpha_2/2
  r_0=r_c
  r_1=r_t
  theta=rand_ratio
  mu=abs(log(r_1/r_0))
  if (mu==0){return(NA)}

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(V)
  n_c=ceiling(info_time*V)
  n_t=n_c*rand_ratio
  ######################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_neg_binomial(r_t,r_c,s,alpha_2,beta,rand_ratio,nu_t,kappa)
#Sample_size_alsp_boundary_neg_binomial(0.2,0.4,c(1/3,2/3,1),c(0.001,0.001,0.048),0.1,1,1,1)
#Sample_size_alsp_boundary_neg_binomial(0.2,0.2,c(1/3,2/3,1),c(0.001,0.001,0.048),0.1,1,1,1)
#Sample_size_alsp_boundary_neg_binomial(1.2,0.84,1,0.05,0.1,1,1,1)
######################################################################

#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Superiority trial, Power calculation, Alpha-spending boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_neg_binomial(r_t = 1, r_c = 0.8, s = 1, alpha_2 = 0.05,
#' sampsz_control = 897, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
power_alsp_boundary_neg_binomial<-function(r_t,r_c,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and power for given r_t,r_c,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  #sampsz_control is the sample size of the control arm
  alpha=alpha_2/2
  r_0=r_c
  r_1=r_t
  theta=rand_ratio
  mu=abs(log(r_1/r_0))
  # n_c=(((a+b)/mu)**2)*(lambda_t/r+lambda_c);

  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  T=sampsz_control/V

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  ######################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
########################
#power_alsp_boundary_neg_binomial(r_c,r_t,s,alpha_2,sampsz_control,rand_ratio,nu_t,kappa)
#power_alsp_boundary_neg_binomial(1,0.8,1,0.05,897,1,1,1)
#power_alsp_boundary_neg_binomial(1,0.8,c(1/3,2/3,1),c(0.001,0.001,0.048),897,1,1,1)
#power_alsp_boundary_neg_binomial(0.2,0.4,c(1/3,2/3,1),c(0.001,0.001,0.048),209,1,1,1)
################################################################

#########################################
####     non-inferiority trials   #######
#########################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################
#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Non-inferiority trial, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger r_t is better, 0 = smaller r_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_OF_boundary_neg_binomial_NI(r_t = 0.2, r_c = 0.4, margin = 0.1, direction = 1,
#' s = c(1/3,2/3,1), alpha = 0.025, beta = 0.1, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
Sample_size_OF_boundary_neg_binomial_NI<-function(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and saple size for given r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger r_t is better
  #direction ne 1 if smaller r_t is better
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  theta=rand_ratio
  r=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  if (mu<=0){return(NA)}

  r_0=r_c
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(V)
  n_c=ceiling(info_time*V)
  n_t=n_c*rand_ratio
  ######################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
####################################################################
#Sample_size_OF_boundary_neg_binomial_NI(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa)
#Sample_size_OF_boundary_neg_binomial_NI(0.2,0.4,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#Sample_size_OF_boundary_neg_binomial_NI(0.2,0.4,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#Sample_size_OF_boundary_neg_binomial_NI(1.2,0.84,0,1, 1,0.025,0.1,1,1,1)
#Sample_size_OF_boundary_neg_binomial_NI(1.2,0.84,0,2, 1,0.025,0.1,1,1,1)
###################################################

#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Non-inferiority trial, Power calculation, O'Brien-Fleming boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger r_t is better, 0 = smaller r_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will display the power
#'
#' @examples
#' power_OF_boundary_neg_binomial_NI(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1,
#' nu_t = 6, kappa = 1)
#'
#' @export
power_OF_boundary_neg_binomial_NI<-function(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and power for given r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size of the control arm
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger r_t is better
  #direction ne 1 if smaller r_t is better
  #nu_t is the exposure time
  #kappa is the dispersion parameter

  r_0=r_c

  theta=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }

  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  T=sampsz_control/V

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=OF_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
########################
#power_OF_boundary_neg_binomial_NI(r_c,r_t,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa)
#power_OF_boundary_neg_binomial_NI(1,0.8,0,1,1,0.025,897,1,1,1)
#power_OF_boundary_neg_binomial_NI(1,0.8,0,2,1,0.025,897,1,1,1)
#power_OF_boundary_neg_binomial_NI(1,0.8,0,1,c(1/3,2/3,1),0.05,897,1,1,1)
#power_OF_boundary_neg_binomial_NI(1,0.8,0,2,c(1/3,2/3,1),0.05,897,1,1,1)
#power_OF_boundary_neg_binomial_NI(0.2,0.4,0.1,1,c(1/3,2/3,1),0.025,81,1,1,1)
#power_OF_boundary_neg_binomial_NI(0.2,0.4,0.1,2,c(1/3,2/3,1),0.025,81,1,1,1)
##################################################

#######################
### Pocock boundary ###
#######################

#################################################
#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Non-inferiority trial, Sample size calculation, Pocock boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger r_t is better, 0 = smaller r_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_pocock_boundary_neg_binomial_NI(r_t = 0.2, r_c = 0.4, margin = 0.1, direction = 1,
#' s = c(1/3,2/3,1), alpha = 0.025, beta = 0.1, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
Sample_size_pocock_boundary_neg_binomial_NI<-function(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and saple size for given r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger r_t is better
  #direction ne 1 if smaller r_t is better
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  theta=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  if (mu<=0){return(NA)}

  r_0=r_c
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(V)
  n_c=ceiling(info_time*V)
  n_t=n_c*rand_ratio
  ######################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
####################################################################
#Sample_size_pocock_boundary_neg_binomial_NI(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa){
#Sample_size_pocock_boundary_neg_binomial_NI(0.2,0.4,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#Sample_size_pocock_boundary_neg_binomial_NI(0.2,0.4,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#Sample_size_pocock_boundary_neg_binomial_NI(1.2,0.84,0,1, 1,0.025,0.1,1,1,1)
#Sample_size_pocock_boundary_neg_binomial_NI(1.2,0.84,0,2, 1,0.025,0.1,1,1,1)
#####################################################

#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Non-inferiority trial, Power calculation, Pocock boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger r_t is better, 0 = smaller r_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will display the power
#'
#' @examples
#' power_pocock_boundary_neg_binomial_NI(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1,
#' nu_t = 6, kappa = 1)
#'
#' @export
power_pocock_boundary_neg_binomial_NI<-function(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and power for given r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is Pocock boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size of the control arm
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger r_t is better
  #direction ne 1 if smaller r_t is better
  #nu_t is the exposure time
  #kappa is the dispersion parameter

  r_0=r_c

  theta=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }

  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  T=sampsz_control/V

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=pocock_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
########################
#power_pocock_boundary_neg_binomial_NI(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa)
#power_pocock_boundary_neg_binomial_NI(1,0.8,0,1,1,0.025,897,1,1,1)
#power_pocock_boundary_neg_binomial_NI(1,0.8,0,2,1,0.025,897,1,1,1)
#power_pocock_boundary_neg_binomial_NI(1,0.8,0,1,c(1/3,2/3,1),0.025,897,1,1,1)
#power_pocock_boundary_neg_binomial_NI(1,0.8,0,2,c(1/3,2/3,1),0.025,897,1,1,1)
#power_pocock_boundary_neg_binomial_NI(0.2,0.4,0.1,1,c(1/3,2/3,1),0.025,92,1,1,1)
#power_pocock_boundary_neg_binomial_NI(0.2,0.4,0.1,2,c(1/3,2/3,1),0.025,92,1,1,1)
######################################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################
#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Non-inferiority trial, Sample size calculation, Alpha-spending boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger r_t is better, 0 = smaller r_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' Sample_size_alsp_boundary_neg_binomial_NI(r_t = 0.2, r_c = 0.4, margin = 0.1, direction = 1,
#' s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), beta = 0.1, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
Sample_size_alsp_boundary_neg_binomial_NI<-function(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and saple size for given r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger r_t is better
  #direction ne 1 if smaller r_t is better
  #nu_t is the exposure time
  #kappa is the dispersion parameter

  theta=rand_ratio

  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  if (mu<=0){return(NA)}

  r_0=r_c
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######
  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ###################################################
  info_time=info_time(mu,s,c_bry,beta,alpha)
  #### saple size for specified distribution ####
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(V)
  n_c=ceiling(info_time*V)
  n_t=n_c*rand_ratio
  ######################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("test arm sample size","control arm sample size","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
####################################################################
#Sample_size_alsp_boundary_neg_binomial_NI(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa)
#Sample_size_alsp_boundary_neg_binomial_NI(0.2,0.4,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#Sample_size_alsp_boundary_neg_binomial_NI(0.2,0.4,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#Sample_size_alsp_boundary_neg_binomial_NI(1.2,0.84,0,1, 1,0.025,0.1,1,1,1)
#Sample_size_alsp_boundary_neg_binomial_NI(1.2,0.84,0,2, 1,0.025,0.1,1,1,1)
########################################################################

#' Group sequential designs / Adaptive sequential designs: Negative binomial distribution, Non-inferiority trial, Power calculation, Alpha-spending boundary
#'
#' @param r_t Mean of the test arm
#' @param r_c Mean of the control arm
#' @param margin Non-inferiority margin
#' @param direction Test direction: 1 = larger r_t is better, 0 = smaller r_t is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t exposure time
#' @param kappa dispersion parameter
#'
#' @return The output will display the power
#'
#' @examples
#' power_alsp_boundary_neg_binomial_NI(r_t = 1, r_c = 0.8, margin = 0, direction = 1, s = 1,
#' alpha = 0.025, sampsz_control = 897, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
power_alsp_boundary_neg_binomial_NI<-function(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa){
  ### This function calculates the critical boundary and power for given r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa for non-inferiroity trial
  # the critical boundary is alpha-spending boundary
  #r_t is the rate of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size of the control arm
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger r_t is better
  #direction ne 1 if smaller r_t is better
  #nu_t is the exposure time
  #kappa is the dispersion parameter

  r_0=r_c

  theta=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }

  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  T=sampsz_control/V

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=alsp_bundry(s,alpha)
  #######Step 1.2: Common to all boundary types ######

  c_bry=as.numeric(unlist(bry_out[1]))
  nom_alpha=as.numeric(unlist(bry_out[2]))
  alpha_s=as.numeric(unlist(bry_out[3]))
  cumu_alpha=as.numeric(unlist(bry_out[4]))
  sum_nom_alpha=as.numeric(unlist(bry_out[5]))
  ########################################
  power=time_power(mu,s,c_bry,T,alpha)
  ########################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power,c_bry,nom_alpha,alpha_s,cumu_alpha,sum_nom_alpha,2*nom_alpha,2*alpha_s,2*cumu_alpha,2*sum_nom_alpha,ref)
  names(out)=c("power","Critical boundaries","1-sided nominal alpha","1-sided alpha spending","1-sided cummulative alpha","1-sided sum_nom_alpha","2-sided nominal alpha","2-sided alpha spending","2-sided cummulative alpha","2-sided sum_nom_alpha","Reference")
  return(out)
}
########################
#power_alsp_boundary_neg_binomial_NI(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa)
#power_alsp_boundary_neg_binomial_NI(1,0.8,0,1,1,0.025,897,1,1,1)
#power_alsp_boundary_neg_binomial_NI(1,0.8,0,2,1,0.025,897,1,1,1)
#power_alsp_boundary_neg_binomial_NI(1,0.8,0,1,c(1/3,2/3,1),c(0.001,0.001,0.023),897,1,1,1)
#power_alsp_boundary_neg_binomial_NI(1,0.8,0,2,c(1/3,2/3,1),c(0.001,0.001,0.023),897,1,1,1)
#power_alsp_boundary_neg_binomial_NI(0.2,0.4,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),80,1,1,1)
#power_alsp_boundary_neg_binomial_NI(0.2,0.4,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),80,1,1,1)

