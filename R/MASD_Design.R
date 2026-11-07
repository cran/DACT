#rm(list=ls())
#library(dplyr)
#library(foreach)
#library(parallel)
#library(doParallel)
#library(bigstatsr)
unirTol = 1e-6
unirIter = 2000

###################################################
##########   Part I                    ############
###  Critical bundary calculations     ############
###################################################

#############################
MSD_alsp_bundry <- function(s, alpha, comp, dist) {
  #this function selects critical boundaries for multiple dose or endpoint comparison#
  #It provides critical boundaries for each q-level testing, 1<=q<=comp#
  #the boundary selection is based on alpha spending at eatch interm analysis#
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #alpha is the vector of alpha spending (with the same length as s). The sum of alpha should be the target one-sided type I error#
  #comp is the number of comparisons, it can be for dose, or endpoints#
  #dist is the distribution of the variable. dist==1 means that the variable is normally distributed and the comparison is for multiple dose comparison#
  #### dist ne 1 means that the distribution is not normmal #####
  ## the next section constructs the covariance matrix###
  ## cry_q calculates the critical bundary for each q level testing ###
  cry_q<-function(s,alpha,q,dist){
    cov_mtx=diag(q)
    if (dist==1){
      d_mtx=(1/2)*matrix(c(rep(1,q*q)),q,q)
      cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
      cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
    } else {
      cov_mtx=diag(q)
    }
    cov_mtx_ws=matrix(nrow=length(s)*q, ncol=length(s)*q)
    for (i in 1:length(s)){
      cov_mtx_ws[(q*(i-1)+1):(q*(i-1)+q),1:(length(s)*q)]=s[i]*cov_mtx
    }
    cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]

    #### Function for calculating critical boundaries #######################

    c_cal1<-function(e0_1,e_bry,sigma1){
      #e0_1 is the  exit boundary for the first interim analysis#
      #sigma1 is the covariance matrix at the first interim analysis#
      rep_1=length(e_bry)
      u_bry_0=NULL
      if (rep_1==0) {
        up_bry_0=NULL
      } else {
        up_bry_0=rep(0,q*rep_1)
        for (i in (1:rep_1)) {
          up_bry_0[(1+(i-1)*q):(i*q)]=e_bry[i]
        }
      }
      up_bry=c(up_bry_0,rep(e0_1,q))
      tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,length(up_bry)),sigma=sigma1,algorithm=mvtnorm::GenzBretz())[1]
      return(tmp)
    }

    ################ function for cumularive alpha ######
    c_cal2<-function(c_bry,inter_vt,sigma1){
      #e0_1 is the  exit boundary for the first interim analysis#
      #sigma1 is the covariance matrix at the first interim analysis#
      rank=inter_vt
      e_bry=c_bry*sqrt(s)
      up_bry=rep(0,q*rank)
      for (i in (1:rank)) {
        up_bry[(1+(i-1)*q):(i*q)]=e_bry[i]
      }
      sigma1=cov_mtx_ws
      sigma2=sigma1[1:(q*rank),1:(q*rank)]
      tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,(q*rank)),sigma=sigma2,algorithm=mvtnorm::GenzBretz())[1]
      return(tmp)
    }
    ##############################################################

    ## the next loop canculate the exit bpundaries bry1, bry2, and critical boundaries c_bry, when there are a to of two interim analyses##
    e_bry=NULL
    for (i in (1:length(s))){
      e_bry_i=as.numeric(stats::uniroot(function(x) c_cal1(x,e_bry,cov_mtx_ws[1:(q*i),1:(q*i)])-sum(alpha[1:i]),c(-2^8,2^128),maxiter = unirIter,extendInt="yes",tol=unirTol)[1])
      e_bry=c(e_bry,e_bry_i)
    }
    c_bry=e_bry/sqrt(s)

    ################################
    cumu_alpha=rep(0,length(s))
    for (j in 1:length(s)) {
      cumu_alpha[j]=c_cal2(c_bry,j,cov_mtx_ws)
    }

    out_q=rbind(c_bry,cumu_alpha)
    # names(out)=c("O'Brien-Fleming type boundary","alpha spent at each interim analysis")
    return(out_q)
  } ## end of q loop##

  ########################################################################################################
  ## The folowing section computes q-level critical boundary and cumulative alphs for each 1<= q <= comp##
  ########################################################################################################
  c_bry_m=matrix(nrow=comp,ncol=length(s))
  cumu_alpha=matrix(nrow=comp,ncol=length(s))
  for (q in 1:comp){
    ret=cry_q(s,alpha,q,dist)
    c_bry_m[q,]=ret[1,]
    cumu_alpha[q,]=ret[2,]
  }

  ##### clacuating alpha-spending ###############
  alpha_s=matrix(nrow=comp,ncol=length(s))
  if(length(s)==1){
    alpha_s=cumu_alpha
  } else {
    alpha_s[,1]=cumu_alpha[,1]
    for (i in 2:length(s)){
      alpha_s[,i]=cumu_alpha[,i]-cumu_alpha[,i-1]
    }
  }

  nom_p=1-stats::pnorm(c_bry_m)

  levelname=c("critical bounadry", "alpha spent at each interim analysis", "cumulative alpha", "nominal p-value")

  out = list()
  if (comp == 1) {
    level = list(c_bry_m, alpha_s, cumu_alpha, nom_p)
    names(level) = levelname
    out = append(out, list(level))
  } else {
    for (i in (comp:1)) {
      level = list(c_bry_m[i, ], alpha_s[i,], cumu_alpha[i,], nom_p[i, ])
      names(level) = levelname
      out = append(out, list(level))
    }
  }
  return(out)
}
##############################
#MSD_alsp_bundry(s,alpha,comp,dist)
#MSD_alsp_bundry(1,0.025,1,1)
#MSD_alsp_bundry(c(0.5,1),c(0.0001,0.0249),1,1)
#MSD_alsp_bundry(c(0.5,1),c(0.0001,0.0249),2,1)
#MSD_alsp_bundry(c(0.5,1),c(0.0001,0.0249),2,0)
#MSD_alsp_bundry(c(0.5,1),c(0.0001,0.0249),3,1)
#MSD_alsp_bundry(c(0.5,1),c(0.0001,0.0249),3,0)
#MSD_alsp_bundry(c(1/3,2/3,1),c(0.0001,0.0001,0.0248),2,1)
#MSD_alsp_bundry(c(1/3,2/3,1),c(0.0001,0.0001,0.0248),2,0)
#MSD_alsp_bundry(c(1/3,2/3,1),c(0.0001,0.0001,0.0248),3,1)
#MSD_alsp_bundry(c(1/3,2/3,1),c(0.0001,0.0001,0.0248),3,0)
#MSD_alsp_bundry(c(1/3,2/3,1),c(0.0001,0.0001,0.0248),4,1)
########################################

##############################
MSD_OF_bundry <- function(s, alpha, comp, dist) {
  #this function selects critical boundaries for multiple dose or endpoint comparison#
  #It provides critical boundaries for each q-level testing, 1<= q <=comp#
  #the boundary selection is O'Brien-Fleming type#
  #This function also calculates the actual alpha spending at each interim analysis########
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #alpha  is the target one-sided type I error#
  #comp is the number of comparisons, it can be for dose, or endpoints#
  #dist is the distribution of the variable. dist==1 means that the variable is normally distributed and the comparison is for multiple dose comparison#
  ## dist ne 1 means that the distributionn is not normal ###
  ## cry_q calculates the critical bundary for each q level testing ###
  cry_q<-function(s,alpha,q,dist){
    if (q==1){
      cov_mtx_ws= matrix(rep(s,length(s)),nrow=length(s), ncol=length(s))
      cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
    }else {
      cov_mtx=diag(q)
      if (dist==1){
        d_mtx=(1/2)*matrix(c(rep(1,q*q)),q,q)
        cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
        cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
      } else {
        cov_mtx=diag(q)
      }
      cov_mtx_ws = matrix(nrow=length(s)*q, ncol=length(s)*q)

      for (i in 1:length(s)){
        cov_mtx_ws[(q*(i-1)+1):(q*(i-1)+q),1:(length(s)*q)]=s[i]*cov_mtx
      }
      cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]
    }

    #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e0#
    c_cal1<-function(e0,rep_time,sigma){
      #e0 is the  exit boundary #
      #sigma is the covariance matrix, it has the rank of =rep_time x reptime#
      up_bry=rep(e0,rep_time)
      tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rep_time),sigma=sigma,algorithm=mvtnorm::GenzBretz())[1]
      return(tmp)
    }

    ################ function for cumularive alpha ######
    c_cal2<-function(c_bry,inter_vt,sigma1){
      #e0_1 is the  exit boundary for the first interim analysis#
      #sigma1 is the covariance matrix at the first interim analysis#
      rank=inter_vt
      e_bry=c_bry*sqrt(s)
      up_bry=rep(0,q*rank)
      for (i in (1:rank)) {
        up_bry[(1+(i-1)*q):(i*q)]=e_bry[i]
      }
      sigma1=cov_mtx_ws
      sigma2=sigma1[1:(q*rank),1:(q*rank)]
      tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,(q*rank)),sigma=sigma2,algorithm=mvtnorm::GenzBretz())[1]
      return(tmp)
    }

    ############## calculating the O'Brien-Fleming boundaries ############
    e_bry=as.numeric(stats::uniroot(function(x) c_cal1(x,length(s)*q,cov_mtx_ws)-alpha,c(-2^8,2^8),extendInt="yes", maxiter = unirIter,tol=unirTol)[1])
    #bry are exit boundaries. c_bry are the ctritical bundaries#
    c_bry=e_bry/sqrt(s)
    ############### calcuating the cummulative alphs spending #############
    cumu_alpha=rep(0,length(s))
    for (j in 1:length(s)) {
      cumu_alpha[j]=c_cal2(c_bry,j,cov_mtx_ws)
    }

    out_q=rbind(c_bry,cumu_alpha)
    # names(out)=c("O'Brien-Fleming type boundary","cumulative alpha at each interim analysis")
    return(out_q)
  } ## end of q loop##

  ########################################################################################################
  ## The folowing section computes q-level critical boundary and cumulative alphs for each 1<= q <= comp##
  ########################################################################################################
  c_bry_m=matrix(nrow=comp,ncol=length(s))
  cumu_alpha=matrix(nrow=comp,ncol=length(s))
  for (q in 1:comp){
    ret=cry_q(s,alpha,q,dist)
    c_bry_m[q,]=ret[1,]
    cumu_alpha[q,]=ret[2,]
  }

  ##### clacuating alpha-spending ###############
  alpha_s=matrix(nrow=comp,ncol=length(s))
  if(length(s)==1){
    alpha_s=cumu_alpha
  } else {
    alpha_s[,1]=cumu_alpha[,1]
    for (i in 2:length(s)){
      alpha_s[,i]=cumu_alpha[,i]-cumu_alpha[,i-1]
    }
  }

  nom_p=1-stats::pnorm(c_bry_m)

  levelname = c("critical bounadry", "alpha spent at each interim analysis", "cumulative alpha", "nominal p-value")
  out = list()
  if (comp == 1) {
    level = list(c_bry_m, alpha_s, cumu_alpha, nom_p)
    names(level) = levelname
    out = append(out, list(level))
  } else {
    for (i in (comp:1)) {
      level = list(c_bry_m[i, ], alpha_s[i, ], cumu_alpha[i, ], nom_p[i, ])
      names(level) = levelname
      out = append(out, list(level))
    }
  }
  return(out)
}
##############################
#MSD_OF_bundry(s,alpha,comp,dist)
#MSD_OF_bundry(1,0.025,1,1)
#MSD_OF_bundry(c(0.5,1),0.025,1,1)
#MSD_OF_bundry(c(0.5,1),0.025,2,1)
#MSD_OF_bundry(c(0.5,1),0.025,2,0)
#MSD_OF_bundry(c(1/3,2/3,1),0.025,2,1)
#MSD_OF_bundry(c(0.3,0.7,1),0.025,3,1)
#MSD_OF_bundry(c(1/3,2/3,1),0.025,2,0)
#MSD_OF_bundry(c(1/4,1/3,2/3,1),0.025,2,1)

###############################

#########################################

###################################################
####                 Part II                  #####
###  Information time and power  calculation  #####
###################################################
#######################################################################
masd_info_time<-function(mu,s,c_bry,alpha,beta,dist){
  #this function calculates the information time that provides 1-beta power#
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #alpha  is the target one-sided type I error#
  #mu is the vector of efeect sizes#
  #dist is the distribution of the efficacy variable
  #dist=1 mean normal distribution
  # dist ne 1, means non-normal distribution
  #### the next section constructs the covariance matrix ###
  ## for the vector of Brownian motions####
  #mu=c(0.1,0.2,0.3)
  #mu=c(0.2,0.2,0.2)
  #s=c(1/3,2/3,1)
  # c_bry=c(3.888374, 3.834043,2.212828)
  #alpha=0.025
  #beta=0.1
  #dist=1
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)
  ###############################
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }
  #############################################
  cov_mtx_ws = matrix(nrow=length(s)*comp, ncol=length(s)*comp)
  for (i in 1:length(s)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s)*comp)]=s[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]

  #####################################################
  e_bry = sqrt(s) * c_bry
  #############################
  ### boundary for the  vector of Brownian  motions####

  cal <- function(T, sigma) {
    #constructioj the upper exit boundary for each T#
    e_bry_m = c(rep(NA, length(s) * comp))
    for (i in 1:length(s)) {
      e_bry_m[((i - 1) * comp + 1):(i * comp)] = e_bry[i] - mu * sqrt(T) * s[i]
    }
    #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e_bry and information time T#
    tmp = 1-pmvnorm(lower=-Inf,upper=e_bry_m,mean=rep(0,comp*length(s)),sigma=sigma,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }

  time_tot=as.numeric(stats::uniroot(function(x) cal(x,cov_mtx_ws)-1+beta,c(0,2^32),extendInt="yes", maxiter = unirIter,tol=unirTol)[1])

  #################################################################
  masd_info_time_s<-function(mu,s,c_bry,beta,alpha){
    #This function calculates the infromation required for 1-beta power for a one-to-one comparison
    #s is the vector of information fractions. The length of s can be either 2 or 3#
    #c_bry is the critial boundary
    #1-beta is the desired power
    #mu is the effct size
    #alpha  is the target one-sided type I error#
    #### the next section constructs the covariance matrix ###
    #mu=0.2
    #s=c(1/3,2/3,1)
    # c_bry=c(4.115415, 2.910038, 2.376036)
    #beta=0.1
    # alpha=0.025

    ########################################
    #### the next section constructs the covariance matrix###
    cov_mtx= matrix(nrow=length(s), ncol=length(s))
    for (i in 1:length(s)){
      cov_mtx[i,1:length(s)]=s[i]
    }
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
      tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,length(s)),sigma=sigma1,algorithm=mvtnorm::GenzBretz())[1]
      return(tmp)
    }

    if (length(s)==1){
      total_time=as.numeric(stats::uniroot(function(x) cal1(x)-1+beta,c(0,2^32),extendInt="yes", maxiter = unirIter,tol=unirTol)[1])
    } else {
      total_time=as.numeric(stats::uniroot(function(x) cal2(x,cov_mtx)-1+beta,c(0,2^36),extendInt="yes", maxiter = unirIter,tol=unirTol)[1])
    }
    return(total_time)
  }

  info_time_ind=NULL
  for (i in 1:length(mu)) {
    if (mu[i]<=0){
      info_time_ind[i]=NA
    }else{
      info_time_ind[i]=masd_info_time_s(mu[i],s,c_bry,beta,alpha)
    }
  }
  out=list(time_tot,info_time_ind)
  names(out)=c("total information time","Individual information time")
  return(out)
}
########################
#masd_info_time(mu,s,c_bry,alpha,beta,dist)
#masd_info_time(c(0,0,0),c(1/3,2/3,1),c(4.115415,2.910038,2.376036),0.025,0.1,1)
#masd_info_time(c(0,0.2,0.3),c(1/3,2/3,1),c(4.115415,2.910038,2.376036),0.025,0.1,1)
#masd_info_time(c(0.1,0.2,0.3),c(1/3,2/3,1),c(4.115415,2.910038,2.376036),0.025,0.1,1)
#masd_info_time(c(0.1,0.2,0.3),c(1/3,2/3,1),c(4.115415,2.910038,2.376036),0.025,0.1,2)
#masd_info_time(c(0.2,0.2,0.2),c(1/3,2/3,1),c(4.115415,2.910038,2.376036),0.025,0.1,1)

##########################################
masd_time_power<-function(mu,s,c_bry,alpha,T,dist){
  #this function calculates the  power for a given information time T#
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #alpha  is the target one-sided type I error#
  #mu is the vector of efeect sizes#
  #dist is the distribution of the efficacy variable
  #dist=1 mean normal distribution
  # dist ne 1, means non-normal distribution
  #### the next section constructs the covariance matrix ###
  ## for the vector of Brownian motions####
  #mu=c(0.1,0.2,0.3)
  #mu=c(0.2,0.2,0.2)
  #s=c(1/3,2/3,1)
  # c_bry=c(3.888374, 3.834043,2.212828)
  #alpha=0.025
  #beta=0.1
  #dist=1
  comp=length(mu)

  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx=diag(comp)
  ###############################
  if (dist==1){
    d_mtx=(1/2)*matrix(c(rep(1,comp*comp)),comp,comp)
    cov_mtx[lower.tri(cov_mtx)]=d_mtx[lower.tri(d_mtx)]
    cov_mtx[upper.tri(cov_mtx)]=d_mtx[upper.tri(d_mtx)]
  } else {
    cov_mtx=diag(comp)
  }
  #############################################
  cov_mtx_ws <- matrix(nrow=length(s)*comp, ncol=length(s)*comp)
  for (i in 1:length(s)){
    cov_mtx_ws[(comp*(i-1)+1):(comp*(i-1)+comp),1:(length(s)*comp)]=s[i]*cov_mtx
  }
  cov_mtx_ws[lower.tri(cov_mtx_ws)]=t(cov_mtx_ws)[lower.tri(cov_mtx_ws)]

  #####################################################
  e_bry=sqrt(s)*c_bry
  #############################
  ### boundary for the  vector of Brownian  motions####

  cal<-function(time,sigma){
    #constructioj the upper exit boundary for each T#
    e_bry_m=c(rep(NA,length(s)*comp))
    for (i in 1:length(s)){
      e_bry_m[((i-1)*comp+1):(i*comp)]=e_bry[i]-mu*sqrt(time)*s[i]
    }
    #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e_bry and information time T#
    tmp = 1-pmvnorm(lower=-Inf,upper=e_bry_m,mean=rep(0,comp*length(s)),sigma=sigma,algorithm=mvtnorm::GenzBretz())[1]
    return(tmp)
  }
  power= cal(T,cov_mtx_ws)
  return(power)
}
########################
#masd_time_power(mu,s,c_bry,alpha,T,dist)
#masd_time_power(c(0.1,0.2,0.3),c(1/3,2/3,1),c(4.115415,2.910038,2.376036),0.025,138,1)
#masd_time_power(c(0.2,0.2,0.2),c(1/3,2/3,1),c(4.115415,2.910038,2.376036),0.025,209,1)

###################################################
####                 Part III                  #####
###  sample size  and power  calculation       #####
###################################################


#############################################################################
#### Part III.1 power and Sample size clculation: Normal distributions   ####
#############################################################################

#############################################################
####    supeiority trials and non-inferiority trials  #######
#############################################################


################################
### O'Brien_Fleming boundary ###
################################

###########################################################
#' Multiple comparisons Group sequential design: Normal distribution, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param mu_t Vector of means of the test arms, one entry per dose or endpoint; the function returns NA
#'   if the largest effect size is below 1e-4
#' @param mu_c Mean of the common control arm
#' @param std_t Common standard deviation of the test arms; must be > 0
#' @param std_c Standard deviation of the control arm; must be > 0
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{test arm sample size}{Sample size required for a test arm; equal to rand_ratio times the control arm sample size}
#'   \item{control arm sample size}{Sample size required for the common control arm}
#'   \item{sample size for each dose in a two arm study with same critical boudary}{Vector with one test arm sample size per dose or endpoint, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{sample size for control arm corresponding to each dose in a two arm study with same critical boudary}{Vector with the matching control arm sample size for each dose or endpoint}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_OF_boundary_normal( mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#'   margin = 0.1, direction = 1, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1,
#'   dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_OF_boundary_normal<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #mu_t is the vector of means for the test arms
  #mu_c is the mean of the control arm
  #std_t>0 is the common standard deviation of the test arms
  #std_c>0 is the standard deviation of the control arm
  #s is the vetor of infomation fraction
  #beta is type II error, power=1-beta
  #0<beta<1
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))
  # alpha_s_1_side=as.numeric(unlist(bry_out[2]))
  # cumu_alpha_1_side=as.numeric(unlist(bry_out[3]))
  # alpha_s_2_side=as.numeric(unlist(bry_out[4]))
  # cumu_alpha_2_side=as.numeric(unlist(bry_out[5]))
  # c_bry_1=as.numeric(unlist(bry_out[6]))
  # nom_p_3=as.numeric(unlist(bry_out[7]))
  # c("Critical boundary","1-sided alpha spending at each interim analysis","1-sided accumulative alpha","2-sided alpha spending at each interim analysis","2-sided accumulative alpha","Critical boundary-phase 3 only","nominal p-value for phase 3")

  #######Step 1.2: Common to all boundary types ######
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  ########## two stage sample size###
  n_c=ceiling(time_tot*ratio)
  n_t=n_c*rand_ratio
  ##############################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm corresponding to each dose in a two arm study with same critical boudary","Type I error control")
  return(out)
}
####################################################################
#MSD_Sample_size_OF_boundary_normal(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2)
#MSD_Sample_size_OF_boundary_normal(c(0,0.3,0.4),0,1,1,0,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_normal(c(0,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_normal(c(0,0.3,0.4),0,1,1,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_normal(c(0,0.3,0.4),0,1,1,0,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_normal(c(0,0.3,0.4),0,1,1,0,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_normal(c(0.1,0.3,0.4),0,1,1,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)

#####################################################################
#' Multiple comparisons Group sequential design: Normal distribution, Power calculation, O'Brien-Fleming boundary
#'
#' @param mu_t Vector of means of the test arms, one entry per dose or endpoint
#' @param mu_c Mean of the common control arm
#' @param std_t Common standard deviation of the test arms; must be > 0
#' @param std_c Standard deviation of the control arm; must be > 0
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param sampsz_control Sample size of the control arm
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given control arm sample size}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_OF_boundary_normal(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1),
#' alpha = 0.025, sampsz_control = 100, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_OF_boundary_normal<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and powerfor given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #mu_t is the vector of means for the test arms
  #mu_c is the mean of the control arm
  #std_t>0 is the common standard deviation of the test arms
  #std_c>0 is the standard deviation of the control arm
  #s is the vetor of infomation fraction
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  #sampsz_control is the sample size for the control arm
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  comp=length(mu)
  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################

  out=list(power,bry_out)
  names(out)=c("power","Type I error control")
  return(out)
}
####################################################################
#MSD_power_OF_boundary_normal(mu_t,mu_c,std_t,std_c,margin,diretion,s,alpha,sampsz_control,rand_ratio,dist,dist2)
#MSD_power_OF_boundary_normal(c(0,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),0.025,94,1,1,1)
#MSD_power_OF_boundary_normal(c(0,0.3,0.4),0,1,1,0.1,2,c(1/3,2/3,1),0.025,2661,1,1,1)
#MSD_power_OF_boundary_normal(c(0,0.3,0.4),0,1,1,0,1,c(1/3,2/3,1),0.025,151,1,1,1)
#MSD_power_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),0.025,149,1,1,1)
#MSD_power_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,2,c(1/3,2/3,1),0.025,149,1,1,1)
###########################################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################
#' Multiple comparisons Group sequential design: Normal distribution, Sample size calculation, Alpha-spending boundary
#'
#' @param mu_t Vector of means of the test arms, one entry per dose or endpoint; the function returns NA
#'   if the largest effect size is below 1e-4
#' @param mu_c Mean of the common control arm
#' @param std_t Common standard deviation of the test arms; must be > 0
#' @param std_c Standard deviation of the control arm; must be > 0
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{test arm sample size}{Sample size required for a test arm; equal to rand_ratio times the control arm sample size}
#'   \item{control arm sample size}{Sample size required for the common control arm}
#'   \item{sample size for each dose in a two arm study with same critical boudary}{Vector with one test arm sample size per dose or endpoint, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{sample size for control arm corresponding to each dose in a two arm study with same critical boudary}{Vector with the matching control arm sample size for each dose or endpoint}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_alsp_boundary_normal( mu_t = c(0.2,0.3,0.4), mu_c = 0, std_t = 1, std_c = 1,
#'   margin = 0, direction = 1, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), beta = 0.1,
#'   rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_alsp_boundary_normal<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #mu_t is the vector of means for the test arms
  #mu_c is the mean of the control arm
  #std_t>0 is the common standard deviation of the test arms
  #std_c>0 is the standard deviation of the control arm
  #s is the vetor of infomation fraction
  #beta is type II error, power=1-beta
  #0<beta<1
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  ########## two stage sample size###
  n_c=ceiling(time_tot*ratio)
  n_t=n_c*rand_ratio
  ##############################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm corresponding to each dose in a two arm study with same critical boudary","Type I error control")
  return(out)
}
####################################################################
#MSD_Sample_size_alsp_boundary_normal(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2)
#MSD_Sample_size_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_normal(c(0.3,0.3,0.3),0,1,1,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
######################################################

#' Multiple comparisons Group sequential design: Normal distribution, Power calculation, Alpha-spending boundary
#'
#' @param mu_t Vector of means of the test arms, one entry per dose or endpoint
#' @param mu_c Mean of the common control arm
#' @param std_t Common standard deviation of the test arms; must be > 0
#' @param std_c Standard deviation of the control arm; must be > 0
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param sampsz_control Sample size of the control arm
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given control arm sample size}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_alsp_boundary_normal(mu_t = c(0.2,0.3,0.4), mu_c = 0, std_t = 1, std_c = 1,
#' margin = 0, direction = 1, s = c(1/3,2/3,1),
#' alpha = c(0.001,0.001,0.023), sampsz_control = 148, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_alsp_boundary_normal<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and powerfor given mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #mu_t is the vector of means for the test arms
  #mu_c is the mean of the control arm
  #std_t>0 is the common standard deviation of the test arms
  #std_c>0 is the standard deviation of the control arm
  #s is the vetor of infomation fraction
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  #sampsz_control is the sample size for the control arm
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  comp=length(mu)
  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################

  out=list(power,bry_out)
  names(out)=c("power","Type I error control")
  return(out)
}
####################################################################
#MSD_power_alsp_boundary_normal(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2)
#MSD_power_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0,1,c(1/3,2/3,1),c(0.001,0.001,0.023),148,1,1,1)
#MSD_power_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),91,1,1,1)
#MSD_power_alsp_boundary_normal(c(0.3,0.3,0.3),0,1,1,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),91,1,1,1)
############################################################
############################################################
####                        Part III.2.                  ###
##  Multiple comparison Sample size and power clculation ###
####                     binary distributions            ###
############################################################

#############################################################
####    supeiority trials and non-inferiority trials  #######
#############################################################



################################
### O'Brien_Fleming boundary ###
################################

#########################################################
#' Multiple comparisons Group sequential design: Binary distribution, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param p_t Vector of event rates of the test arms, one entry per dose or endpoint; 0 < p_t < 1, and
#'   the function returns NA if the largest effect size is below 1e-4
#' @param p_c Event rate of the common control arm; 0 < p_c < 1
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{test arm sample size}{Sample size required for a test arm; equal to rand_ratio times the control arm sample size}
#'   \item{control arm sample size}{Sample size required for the common control arm}
#'   \item{sample size for each dose in a two arm study with same critical boudary}{Vector with one test arm sample size per dose or endpoint, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{sample size for control arm corresponding to each dose in a two arm study with same critical boudary}{Vector with the matching control arm sample size for each dose or endpoint}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_OF_boundary_binary_diff( p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 1,
#'   s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_OF_boundary_binary_diff<-function(p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  #for superiority and non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the vector of rates of the test arm
  #p_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #s is the vetor of infomation fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise

  r=rand_ratio
  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  ############################################
  n_c=ceiling(time_tot*max(ratio))
  n_t=n_c*rand_ratio
  ##############################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm corresponding to each dose in a two arm study with same critical boudary","Type I error control")
  return(out)
}
####################################################################
#MSD_Sample_size_OF_boundary_binary_diff(p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
#MSD_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_binary_diff(c(0.3,0.3,0.3),0.5,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_binary_diff(c(0.3,0.3,0.3),0.5,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_binary_diff(c(0.3,0.3,0.3),0.5,0,2,c(1/3,2/3,1),0.025,0.1,1,1,1)

##########################################################
#' Multiple comparisons Group sequential design: Binary distribution, Power calculation, O'Brien-Fleming boundary
#'
#' @param p_t Vector of event rates of the test arms, one entry per dose or endpoint; 0 < p_t < 1
#' @param p_c Event rate of the common control arm; 0 < p_c < 1
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param sampsz_control Sample size of the control arm
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given control arm sample size}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_OF_boundary_binary_diff(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100,
#' rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_OF_boundary_binary_diff<-function(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  #for superiority and non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the vector of rates of the test arm
  #p_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #s is the vetor of infomation fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size for the control arm
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise

  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  comp=length(mu)
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  T=sampsz_control/max(ratio)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################

  out=list(power,bry_out)
  names(out)=c("power","Type I error control")
  return(out)
}
####################################################################
#MSD_power_OF_boundary_binary_diff(p_t,p_c,margin,diretion,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
#MSD_power_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,1,c(1/3,2/3,1),0.025,68,1,1,1)
#MSD_power_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),0.025,68,1,1,1)
#MSD_power_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),0.025,37,1,1,1)

###########################################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################
#' Multiple comparisons Group sequential design: Binary distribution, Sample size calculation, Alpha-spending boundary
#'
#' @param p_t Vector of event rates of the test arms, one entry per dose or endpoint; 0 < p_t < 1, and
#'   the function returns NA if the largest effect size is below 1e-4
#' @param p_c Event rate of the common control arm; 0 < p_c < 1
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{test arm sample size}{Sample size required for a test arm; equal to rand_ratio times the control arm sample size}
#'   \item{control arm sample size}{Sample size required for the common control arm}
#'   \item{sample size for each dose in a two arm study with same critical boudary}{Vector with one test arm sample size per dose or endpoint, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{sample size for control arm corresponding to each dose in a two arm study with same critical boudary}{Vector with the matching control arm sample size for each dose or endpoint}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_alsp_boundary_binary_diff( p_t = c(0.2,0.3,0.4), p_c = 0.5, margin = 0.1,
#'   direction = 2, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), beta = 0.1, rand_ratio = 1,
#'   dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_alsp_boundary_binary_diff<-function(p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #p_t is the vector of rates of the test arm
  #p_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #s is the vetor of infomation fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise

  r=rand_ratio
  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  ############################################
  n_c=ceiling(time_tot*max(ratio))
  n_t=n_c*rand_ratio
  ##############################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm corresponding to each dose in a two arm study with same critical boudary","Type I error control")
  return(out)
}
####################################################################
#MSD_Sample_size_alsp_boundary_binary_diff(p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
#MSD_Sample_size_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_binary_diff(c(0.3,0.3,0.3),0.5,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_binary_diff(c(0.3,0.3,0.3),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)

######################################################

#' Multiple comparisons Group sequential design: Binary distribution, Power calculation, Alpha-spending boundary
#'
#' @param p_t Vector of event rates of the test arms, one entry per dose or endpoint; 0 < p_t < 1
#' @param p_c Event rate of the common control arm; 0 < p_c < 1
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param sampsz_control Sample size of the control arm
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given control arm sample size}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_alsp_boundary_binary_diff(p_t = c(0.2,0.3,0.4), p_c = 0.5, margin = 0.1,
#' direction = 1, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.048), sampsz_control = 68,
#' rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_alsp_boundary_binary_diff<-function(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #p_t is the vector of rates of the test arm
  #p_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #s is the vetor of infomation fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size for the control arm
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  r=rand_ratio
  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  comp=length(mu)
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  T=sampsz_control/max(ratio)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################

  out=list(power,bry_out)
  names(out)=c("power","Type I error control")
  return(out)
}
####################################################################
#MSD_power_alsp_boundary_binary_diff(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2)
#MSD_power_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.048),68,1,1,1)
#MSD_power_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.048),68,1,1,1)
#MSD_power_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),37,1,1,1)

####################################################################

#############################################################################
#### Part III.3 power and Sample size clculation: Survival analysis      ####
#############################################################################

#############################################################
####    supeiority trials and non-inferiority trials  #######
#############################################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################
#' Multiple comparisons Group sequential design: Survival analysis, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param HR Vector of hazard ratios of the test arms vs. the control arm, one entry per dose or endpoint;
#'   the function returns NA if the largest effect size is below 1e-4
#' @param margin Non-inferiority margin on the hazard ratio scale; must be > 0, with margin = 1 for a
#'   superiority trial and margin != 1 for a non-inferiority trial (0 < margin < 1 if larger HR is
#'   better, margin > 1 if smaller HR is better)
#' @param direction Test direction: 1 = larger HR is better, otherwise smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage; dist1 = 2 is the conservative choice
#' @param dist2 Distribution indicator for the second stage; dist2 = 1 is the conservative choice
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Combined number of events between a dose and control}{Total number of events, test arm and control arm combined, required for one dose vs. control comparison}
#'   \item{total number of events from all doses and control}{Total number of events required for the whole trial, i.e. all doses and the common control arm}
#'   \item{combined number of events in each dose and control in a two arm study with same critical boudary}{Vector with the combined number of events for each dose and its control, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#'   \item{reference}{Reference for the Schoenfeld event number formula}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_OF_boundary_survival_Schoenfeld( HR = 0.7, margin = 0.1, direction = 1,
#'   s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_OF_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #HR is the vector of hazards ratio of test arms vs. control arm
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
  #Margin = 1 for superiority trial
  #margin !=1 for non-inferiority trial
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  #dist1=2 (this is conservative)
  #dist2=1 (this is conservative)

  if (direction==1) {
    mu=log(HR/margin)
  } else{
    mu=-log(HR/margin)
  }
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_2=ceiling(time_tot/(p_c*p_t))
  n_tot=ceiling((n_2/2)*(length(mu)+1))
  ##############################
  n_t_individual=ceiling(info_time_ind/(p_c*p_t))
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(n_2,n_tot,n_t_individual,bry_out,ref)
  names(out)=c("Combined number of events between a dose and control","total number of events from all doses and control","combined number of events in each dose and control in a two arm study with same critical boudary","Type I error control","reference")
  return(out)
}
####################################################################
#MSD_Sample_size_OF_boundary_survival_Schoenfeld(HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
#MSD_Sample_size_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),0.8,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)

#####################################################################

#' Multiple comparisons Group sequential design: Survival analysis, Power calculation, O'Brien-Fleming boundary
#'
#' @param HR Vector of hazard ratios of the test arms vs. the control arm, one entry per dose or endpoint
#' @param margin Non-inferiority margin on the hazard ratio scale; must be > 0, with margin = 1 for a
#'   superiority trial and margin != 1 for a non-inferiority trial (0 < margin < 1 if larger HR is
#'   better, margin > 1 if smaller HR is better)
#' @param direction Test direction: 1 = larger HR is better, otherwise smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param events_tot Total number of events for the trial, used in place of the control arm sample size
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage; dist1 = 2 is the conservative choice
#' @param dist2 Distribution indicator for the second stage; dist2 = 1 is the conservative choice
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given total number of events}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#'   \item{reference}{Reference for the Schoenfeld event number formula}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_OF_boundary_survival_Schoenfeld(HR = 0.7, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, events_tot = 100,
#' rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_OF_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,events_tot,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for given HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #HR is the vector of hazards ratio of test arms vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #events_tot is total number of events
  #margin is the non-inferiority margin
  #margin>0
  #Margin = 1 for superiority trial
  #margin !=1 for non-inferiority trial
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  #dist1=2 (this is conservative)
  #dist2=1 (this is conservative)

  if (direction==1) {
    mu=log(HR/margin)
  } else{
    mu=-log(HR/margin)
  }
  comp=length(mu)
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_2=2* events_tot/(length(mu)+1)
  T=n_2*p_c*p_t
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power,bry_out,ref)
  names(out)=c("power","Type I error control","reference")
  return(out)
}
####################################################################
#MSD_power_OF_boundary_survival_Schoenfeld(HR,margin,direction,s,alpha,events_tot,rand_ratio,dist1,dist2)
#MSD_power_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),0.025,17,1,1,1)
#MSD_power_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),0.025,19,1,1,1)
#MSD_power_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),0.025,38,1,1,1)
###########################################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################
#' Multiple comparisons Group sequential design: Survival analysis, Sample size calculation, Alpha-spending boundary
#'
#' @param HR Vector of hazard ratios of the test arms vs. the control arm, one entry per dose or endpoint;
#'   the function returns NA if the largest effect size is below 1e-4
#' @param margin Non-inferiority margin on the hazard ratio scale; must be > 0, with margin = 1 for a
#'   superiority trial and margin != 1 for a non-inferiority trial (0 < margin < 1 if larger HR is
#'   better, margin > 1 if smaller HR is better)
#' @param direction Test direction: 1 = larger HR is better, otherwise smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage; dist1 = 2 is the conservative choice
#' @param dist2 Distribution indicator for the second stage; dist2 = 1 is the conservative choice
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Combined number of events between a dose and control}{Total number of events, test arm and control arm combined, required for one dose vs. control comparison}
#'   \item{total number of events from all doses and control}{Total number of events required for the whole trial, i.e. all doses and the common control arm}
#'   \item{combined number of events in each dose and control in a two arm study with same critical boudary}{Vector with the combined number of events for each dose and its control, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#'   \item{reference}{Reference for the Schoenfeld event number formula}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_alsp_boundary_survival_Schoenfeld( HR = c(0.2,0.3,0.4), margin = 1.1,
#'   direction = 2, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), beta = 0.1, rand_ratio = 1,
#'   dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_alsp_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #HR is the vector of hazards ratio of test arms vs. control arm
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
  #Margin = 1 for superiority trial
  #margin !=1 for non-inferiority trial
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  #dist1=2 (this is conservative)
  #dist2=1 (this is conservative)

  if (direction==1) {
    mu=log(HR/margin)
  } else{
    mu=-log(HR/margin)
  }
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_2=ceiling(time_tot/(p_c*p_t))
  n_tot=ceiling((n_2/2)*(length(mu)+1))
  ##############################
  n_t_individual=ceiling(info_time_ind/(p_c*p_t))
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(n_2,n_tot,n_t_individual,bry_out,ref)
  names(out)=c("Combined number of events between a dose and control","total number of events from all doses and control","combined number of events in each dose and control in a two arm study with same critical boudary","Type I error control","reference")
  return(out)
}
####################################################################
#MSD_Sample_size_alsp_boundary_survival_Schoenfeld(HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
#MSD_Sample_size_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)

########################################################################


#' Multiple comparisons Group sequential design: Survival analysis, Power calculation, Alpha-spending boundary
#'
#' @param HR Vector of hazard ratios of the test arms vs. the control arm, one entry per dose or endpoint
#' @param margin Non-inferiority margin on the hazard ratio scale; must be > 0, with margin = 1 for a
#'   superiority trial and margin != 1 for a non-inferiority trial (0 < margin < 1 if larger HR is
#'   better, margin > 1 if smaller HR is better)
#' @param direction Test direction: 1 = larger HR is better, otherwise smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param events_tot Total number of events for the trial, used in place of the control arm sample size
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage; dist1 = 2 is the conservative choice
#' @param dist2 Distribution indicator for the second stage; dist2 = 1 is the conservative choice
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given total number of events}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#'   \item{reference}{Reference for the Schoenfeld event number formula}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_alsp_boundary_survival_Schoenfeld(HR = c(0.2,0.3,0.4), margin = 1.1, direction = 2,
#' s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023),
#' events_tot = 17, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_alsp_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,events_tot,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for given HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #HR is the vector of hazards ratio of test arms vs. control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #events_tot is total number of events
  #margin is the non-inferiority margin
  #margin>0
  #Margin = 1 for superiority trial
  #margin !=1 for non-inferiority trial
  #direction=1 if larger HR is better
  #direction ne 1 if smaller HR is better
  #dist1=2 (this is conservative)
  #dist2=1 (this is conservative)
  if (direction==1) {
    mu=log(HR/margin)
  } else{
    mu=-log(HR/margin)
  }
  comp=length(mu)
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  n_2=2* events_tot/(length(mu)+1)
  T=n_2*p_c*p_t
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power,bry_out,ref)
  names(out)=c("power","Type I error control","reference")
  return(out)
}
####################################################################
#MSD_power_alsp_boundary_survival_Schoenfeld(HR,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
#MSD_power_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),17,1,1,1)
#MSD_power_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),38,1,1,1)
###########################################################


############################################################
####                        Part III.4.                  ###
##  Multiple comparison Sample size and power calculation ###
####                     Poisson distributions           ###
############################################################

############################################################
####    supeiority trials andnon-inferiority trials  #######
############################################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################

#' Multiple comparisons Group sequential design: Poisson distribution, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param lambda_t Vector of means (rates) of the test arms, one entry per dose or endpoint; the function
#'   returns NA if the largest effect size is below 1e-4
#' @param lambda_c Mean (rate) of the common control arm
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{test arm sample size}{Sample size required for a test arm; equal to rand_ratio times the control arm sample size}
#'   \item{control arm sample size}{Sample size required for the common control arm}
#'   \item{sample size for each dose in a two arm study with same critical boudary}{Vector with one test arm sample size per dose or endpoint, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{sample size for control arm corresponding to each dose in a two arm study with same critical boudary}{Vector with the matching control arm sample size for each dose or endpoint}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_OF_boundary_poisson( lambda_t = 1.1, lambda_c = 1.0, margin = 0.1,
#'   direction = 1, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1, dist1 = 1,
#'   dist2 = 1)
#'
#' @export
MSD_Sample_size_OF_boundary_poisson<-function(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #lambda_t is the vector of rates of the test arm
  #lambda_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #s is the vetor of infomation fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(time_tot*max(ratio))
  n_t=n_c*rand_ratio
  ##############################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm corresponding to each dose in a two arm study with same critical boudary","Type I error control")
  return(out)
}
####################################################################
#MSD_Sample_size_OF_boundary_poisson(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
#MSD_Sample_size_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.15,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_poisson(c(0.3,0.3,0.3),0.5,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#MSD_Sample_size_OF_boundary_poisson(c(0.3,0.3,0.3),0.5,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)

######################################################
#' Multiple comparisons Group sequential design: Poisson distribution, Power calculation, O'Brien-Fleming boundary
#'
#' @param lambda_t Vector of means (rates) of the test arms, one entry per dose or endpoint
#' @param lambda_c Mean (rate) of the common control arm
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param sampsz_control Sample size of the control arm
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given control arm sample size}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_OF_boundary_poisson(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_OF_boundary_poisson<-function(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #lambda_t is the vector of rates of the test arm
  #lambda_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #s is the vetor of infomation fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size for the control arm
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  comp=length(mu)
  ratio=lambda_c+lambda_t/rand_ratio
  T=sampsz_control/max(ratio)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################

  out=list(power,bry_out)
  names(out)=c("power","Type I error control")
  return(out)
}
####################################################################
#MSD_power_OF_boundary_poisson(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2)
#MSD_power_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.1,1,c(1/3,2/3,1),0.025,125,1,1,1)
#MSD_power_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),0.025,125,1,1,1)
#MSD_power_OF_boundary_poisson(c(0.3,0.3,0.3),0.5,0.1,2,c(1/3,2/3,1),0.025,75,1,1,1)
###########################################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################
#' Multiple comparisons Group sequential design: Poisson distribution, Sample size calculation, Alpha-spending boundary
#'
#' @param lambda_t Vector of means (rates) of the test arms, one entry per dose or endpoint; the function
#'   returns NA if the largest effect size is below 1e-4
#' @param lambda_c Mean (rate) of the common control arm
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{test arm sample size}{Sample size required for a test arm; equal to rand_ratio times the control arm sample size}
#'   \item{control arm sample size}{Sample size required for the common control arm}
#'   \item{sample size for each dose in a two arm study with same critical boudary}{Vector with one test arm sample size per dose or endpoint, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{sample size for control arm corresponding to each dose in a two arm study with same critical boudary}{Vector with the matching control arm sample size for each dose or endpoint}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_alsp_boundary_poisson( lambda_t = c(0.2,0.3,0.4), lambda_c = 0.5, margin = 0.11,
#'   direction = 2, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), beta = 0.1, rand_ratio = 1,
#'   dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_alsp_boundary_poisson<-function(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #lambda_t is the vector of rates of the test arm
  #lambda_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #s is the vetor of infomation fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  if(max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(time_tot*max(ratio))
  n_t=n_c*rand_ratio
  ##############################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm corresponding to each dose in a two arm study with same critical boudary","Type I error control")
  return(out)
}
####################################################################
#MSD_Sample_size_alsp_boundary_poisson(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
#MSD_Sample_size_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.11,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_poisson(c(0.3,0.3,0.3),0.5,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#MSD_Sample_size_alsp_boundary_poisson(c(0.3,0.3,0.3),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)

###########################################################

#' Multiple comparisons Group sequential design: Poisson distribution, Power calculation, Alpha-spending boundary
#'
#' @param lambda_t Vector of means (rates) of the test arms, one entry per dose or endpoint
#' @param lambda_c Mean (rate) of the common control arm
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param sampsz_control Sample size of the control arm
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control, with normal distributions; 2 = otherwise
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given control arm sample size}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_alsp_boundary_poisson(lambda_t = c(0.2,0.3,0.4), lambda_c = 0.5, margin = 0.1,
#' direction = 2, s = c(1/3,2/3,1),
#' alpha = c(0.001,0.001,0.023), sampsz_control = 67, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_alsp_boundary_poisson<-function(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #lambda_t is the vector of rates of the test arm
  #lambda_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #s is the vetor of infomation fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #sampsz_control is the sample size for the control arm
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control, with normal distributions
  #dist2=2 otherwise
  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  comp=length(mu)
  ratio=lambda_c+lambda_t/rand_ratio
  T=sampsz_control/max(ratio)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################

  out=list(power,bry_out)
  names(out)=c("power","Type I error control")
  return(out)
}
####################################################################
#MSD_power_alsp_boundary_poisson(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2)
#MSD_power_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),67,1,1,1)

#MSD_power_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),124,1,1,1)
#MSD_power_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),124,1,1,1)
#MSD_power_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0,2,c(1/3,2/3,1),c(0.001,0.001,0.023),124,1,1,1)
#MSD_power_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),67,1,1,1)

###################################################################
############################################################
####                  Part III.5.                        ###
##  Multiple comparison Sample size and power clculation ###
####          negative binomial  distributionn           ###
############################################################
###########################################################

#############################################################
####    supeiority trials  and non-inferiority trials #######
#############################################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################

#' Multiple comparisons Group sequential design: Negative binomial distribution, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param r_t Vector of mean rates of the test arms, one entry per dose or endpoint; the function returns NA
#'   if the largest log rate ratio is below 1e-4 or if any adjusted test arm rate r_t +/- margin is below 1e-4
#' @param r_c Mean rate of the common control arm
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param nu_t Exposure time; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control; 2 = otherwise. It has no effect here: the boundary calculation always uses dist1 = 2
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control; 2 = otherwise. It has no effect here: the information time calculation always uses dist2 = 1
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{test arm sample size}{Sample size required for a test arm; equal to rand_ratio times the control arm sample size}
#'   \item{control arm sample size}{Sample size required for the common control arm}
#'   \item{sample size for each dose in a two arm study with same critical boudary}{Vector with one test arm sample size per dose or endpoint, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{sample size for control arm corresponding to each dose in a two arm study with same critical boudary}{Vector with the matching control arm sample size for each dose or endpoint}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#'   \item{Reference}{Reference for the negative binomial sample size formula}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_OF_boundary_neg_binomial(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1),
#' alpha = 0.025, beta = 0.1, rand_ratio = 1, nu_t = 6, kappa = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_OF_boundary_neg_binomial<-function(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2 for
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the vector of the rates of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control
  #dist2=2 otherwise
  r_0=r_c
  theta=rand_ratio
  r=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  if(max(mu)<1e-4 || min(r_1)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  dist1=2
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  dist2=1
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(max(V))
  n_c=ceiling(time_tot*max(V))
  n_t=n_c*rand_ratio
  ##############################
  n_c_individual=ceiling(info_time_ind*V)
  n_t_individual=n_c_individual*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,ref)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm corresponding to each dose in a two arm study with same critical boudary","Type I error control","Reference")
  return(out)
}
####################################################################
#MSD_Sample_size_OF_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2){
#MSD_Sample_size_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1,1,1)
#MSD_Sample_size_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1,1,1)

######################################################

#' Multiple comparisons Group sequential design: Negative binomial distribution, Power calculation, O'Brien-Fleming boundary
#'
#' @param r_t Vector of mean rates of the test arms, one entry per dose or endpoint
#' @param r_c Mean rate of the common control arm
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Type I error for the one-sided test; 0 < alpha < 1
#' @param sampsz_control Sample size of the control arm
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param nu_t Exposure time; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param dist1 Distribution indicator for the first stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control; 2 = otherwise. It has no effect here: the boundary calculation always uses dist1 = 2
#' @param dist2 Distribution indicator for the second stage: 1 = the sample size calculation is for multiple
#'   dose vs. a common control; 2 = otherwise. It has no effect here: the information time calculation always uses dist2 = 1
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given control arm sample size}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#'   \item{Reference}{Reference for the negative binomial sample size formula}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_OF_boundary_neg_binomial(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1,
#' nu_t = 6, kappa = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_power_OF_boundary_neg_binomial<-function(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2){
  ### This function calculates the critical boundary and power for given r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the vector of the rates of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  #sampsz_control is the sample size for the control arm
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control
  #dist2=2 otherwise
  r_0=r_c
  theta=rand_ratio
  r=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  comp=length(mu)
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(max(V))
  T=sampsz_control/max(V)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  dist1=2
  bry_out=MSD_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  dist2=1
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power,bry_out,ref)
  names(out)=c("power","Type I error control","Reference")
  return(out)
}
####################################################################
#MSD_power_OF_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2){
#MSD_power_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,1,c(1/3,2/3,1),0.05,91,1,1,1,1,1)
#MSD_power_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),0.05,91,1,1,1,1,1)
#MSD_power_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),0.025,56,1,1,1,1,1)

###########################################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################
#' Multiple comparisons Group sequential design: Negative binomial distribution, Sample size calculation, Alpha-spending boundary
#'
#' @param r_t Vector of mean rates of the test arms, one entry per dose or endpoint; the function returns NA
#'   if the largest log rate ratio is below 1e-4 or if any adjusted test arm rate r_t +/- margin is below 1e-4
#' @param r_c Mean rate of the common control arm
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param beta Type II error rate, i.e. 1 - power; 0 < beta < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param nu_t Exposure time; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param dist1 Distribution indicator for the first stage; it has no effect here, because the boundary
#'   calculation always uses dist1 = 2
#' @param dist2 Distribution indicator for the second stage; it has no effect here, because the information
#'   time calculation always uses dist2 = 1
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{test arm sample size}{Sample size required for a test arm; equal to rand_ratio times the control arm sample size}
#'   \item{control arm sample size}{Sample size required for the common control arm}
#'   \item{sample size for each dose in a two arm study with same critical boudary}{Vector with one test arm sample size per dose or endpoint, as if a separate two arm trial were conducted for that dose with the same critical boundary}
#'   \item{sample size for control arm corresponding to each dose in a two arm study with same critical boudary}{Vector with the matching control arm sample size for each dose or endpoint}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#'   \item{Reference}{Reference for the negative binomial sample size formula}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_Sample_size_alsp_boundary_neg_binomial( r_t = c(0.2,0.3,0.4), r_c = 0.6, margin = 0.1,
#'   direction = 2, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), beta = 0.1, rand_ratio = 1,
#'   nu_t = 1, kappa = 1, dist1 = 1, dist2 = 1)
#'
#' @export
MSD_Sample_size_alsp_boundary_neg_binomial<-function(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the vector of the rates of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  #dist1=2
  #dist2=1
  r_0=r_c
  theta=rand_ratio
  r=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  if(max(mu)<1e-4 || min(r_1)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  dist1=2
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))
  #######Step 1.2: Common to all boundary types ######
  dist2=1
  info_time_out=masd_info_time(mu,s,c_bry,alpha,beta,dist2)
  ##########################################################
  time_tot=as.numeric(unlist( info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(max(V))
  n_c=ceiling(time_tot*max(V))
  n_t=n_c*rand_ratio
  ##############################
  n_c_individual=ceiling(info_time_ind*V)
  n_t_individual=n_c_individual*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,ref)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm corresponding to each dose in a two arm study with same critical boudary","Type I error control","Reference")
  return(out)
}
####################################################################
#MSD_Sample_size_alsp_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2){
#MSD_Sample_size_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1,1,1)
#MSD_Sample_size_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1,1,1)
#MSD_Sample_size_alsp_boundary_neg_binomial(c(0.3,0.3,0.3),0.6,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1,1,1)
#MSD_Sample_size_alsp_boundary_neg_binomial(c(0.3,0.3,0.3),0.6,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1,1,1)
#MSD_Sample_size_alsp_boundary_neg_binomial(c(0.3,0.3,0.3),0.6,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1,1,1)
#######################################################


#' Multiple comparisons Group sequential design: Negative binomial distribution, Power calculation, Alpha-spending boundary
#'
#' @param r_t Vector of mean rates of the test arms, one entry per dose or endpoint
#' @param r_c Mean rate of the common control arm
#' @param margin Non-inferiority margin; must be >= 0, set to 0 for a superiority trial and > 0 for a non-inferiority trial
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha Vector of alpha-spending values for the one-sided test, one per analysis and of the same
#'   length as s, summing to the overall one-sided type I error; each value must satisfy 0 < alpha < 1
#' @param sampsz_control Sample size of the control arm
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param nu_t Exposure time; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param dist Distribution indicator kept for consistency with the sibling functions; it is not used here,
#'   because the function always sets dist1 = 2 for the boundary calculation and dist2 = 1 for the information time calculation
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{power}{Power of the trial for the given control arm sample size}
#'   \item{Type I error control}{Boundary output of the closed testing procedure: a list with one element per comparison level (the M-level boundary first, down to the 1-level boundary), each element a list containing the vectors 'critical bounadry', 'alpha spent at each interim analysis', 'cumulative alpha' and 'nominal p-value' for the analyses}
#'   \item{Reference}{Reference for the negative binomial sample size formula}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' MSD_power_alsp_boundary_neg_binomial(r_t=c(0.2,0.3,0.4), r_c=0.6, margin=0.1, direction=1,
#' s=c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), sampsz_control = 90, rand_ratio = 1,
#' nu_t = 1, kappa = 1, dist = 1)
#'
#' @export
MSD_power_alsp_boundary_neg_binomial<-function(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist){
  ### This function calculates the critical boundary and power for given r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #r_t is the vector of the rates of the test arm
  #r_c is the rate of the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  # s is the information fraction
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>=0
  #Margin=0 for superiority trial
  #margin>0 for non-inferiority trial
  #direction=1 if larger lambda is better
  #direction ne 1 if smaller lambda is better
  #sampsz_control is the sample size for the control arm
  #nu_t is the exposure time
  #kappa is the dispersion parameter
  #dist1=1 if the the sample size calculation is for multiple dose vs. a common control
  #dist1=2 otherwise
  #dist2=1 if the the sample size calculation is for multiple dose vs. a common control
  #dist2=2 otherwise
  r_0=r_c
  theta=rand_ratio
  r=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  comp=length(mu)
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(max(V))
  T=sampsz_control/max(V)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  dist1=2
  bry_out=MSD_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[[1]][1]))

  #######Step 1.2: Common to all boundary types ######
  dist2=1
  power=masd_time_power(mu,s,c_bry,alpha,T,dist2)
  ##########################################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power,bry_out,ref)
  names(out)=c("power","Type I error control","Reference")
  return(out)
}
####################################################################
#MSD_power_alsp_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist)
#MSD_power_alsp_boundary_neg_binomial(r_t=c(0.2,0.3,0.4),r_c=0.6,margin=0.1,direction=1,s=c(1/3,2/3,1),c(0.001,0.001,0.023),90,1,1,1,1)
#MSD_power_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),90,1,1,1,1)
#MSD_power_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),60,1,1,1,1)
#############################################################
