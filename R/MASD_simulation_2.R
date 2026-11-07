#rm(list=ls())
#install.packages("glm")
#install.packages("glm2")
#install.packages("MASS")
#library(tidyverse)
#library(ggplot2)


#disable_code = TRUE
#disable_code = FALSE

MAX_INT = 2^32

###################################################
##########   Part I                    ############
### Component functions, to be sourced ############
###################################################
#######################################################
###########################################################
masd_c_power_com<-function(s,c_bry,theta,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist){
  #this function calculates conditional power, assuming theta, given that theta_hat_inter,theta_hat_inter_sd have been observed at inter_vt=length(s)-1#
  #s is the vector of information fractions#
  #c_bry is the critical boundary
  # theta_hat_inter is the vector of theta_hat at the inter_vt interim analysis###
  #dist=1 if the the test arms share one common contorl and the endpoint is normally distributed
  inter_vt=length(s)-1
  comp=length(theta_hat_inter_m1)
  #theta=theta_hat_inter_m1
  theta_hat=theta_hat_inter_m1

  T_inter_0=1/(max(theta_hat_inter_m1_sd))^2
  T_inter=rep(T_inter_0,length(theta_hat_inter_m1))
  score_inter=theta_hat_inter_m1*T_inter
  inter_wald=theta_hat_inter_m1*sqrt(T_inter)
  T=T_inter/s[inter_vt]
  # T=rep(T_0,length(theta_hat_inter_m1))

  # t=T*s
  ########################################
  #### the next section constructs the covariance matrix###

  cov_mtx = if (dist==1) {
    matrix(0.5, comp, comp)+diag(0.5, comp)
  } else {
    diag(comp)
  }

  #####################################################
  ### boundary for the  vector of Brownian  motions####

  u_bry_theta=((c_bry[length(s)]*sqrt(T))-score_inter)/sqrt(T-T_inter)-theta*sqrt(T-T_inter)


  #### The functions c_cal calcuate the  exit probablities for the given exit bounadries e0#
  ## the function for the vector of Brownian motions ###
  #################################################
  ## a multivariate-normal integral with non-finite boundaries must not abort
  ## the replicate: report NA (non-evaluable) instead
  c_power= tryCatch(1-pmvnorm(lower=-Inf,upper=u_bry_theta,mean=0,sigma=cov_mtx,abseps=1e-8)[1],
                    error=function(e) NA_real_)
  return(c_power)
}
########################################
#masd_c_power_com(s,c_bry,theta,theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)
#masd_c_power_com(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),theta_hat_inter_m1=c(0.3,0.4,1),theta_hat_inter_m1_sd=0.5,dist=1)
#masd_c_power_com(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),theta=c(0.3,0.4,1),theta_hat_inter_m1=c(0.3,0.4,1),theta_hat_inter_m1_sd=c(0.4,0.5,0.6),dist=1)
#masd_c_power_com(s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),theta=c(0,0,0),theta_hat_inter_m1=c(0.3,0.4,1),theta_hat_inter_m1_sd=c(0.4,0.5,0.6),dist=1)

##########
new_ebry_com<-function(T,T_inter,T_new,c_alpha,s,theta_hat_inter_m2,theta_hat_inter_m2_sd,dist){
  #this function calculates the new critical boundaries, given the conditional type I error of c_alpha#
  #only one final analysis after sample size change
  #c_alpha is tconditional type I error##

  # theta_hat_inter is the vector of the theta_hat at the interim analysis###
  #theta_hat_inter_sd is the standard error of theta_hat_inter

  #dist is the distribution of the variable. dist=1 means that the variable is normally distributed and the comparison is for multiple dose comparison#

  comp=length(theta_hat_inter_m2)
  score_inter=theta_hat_inter_m2/(theta_hat_inter_m2_sd*theta_hat_inter_m2_sd)

  ###################################
  #### the next section constructs the covariance matrix###


  cov_mtx = if (dist==1){
      matrix(0.5, comp, comp)+diag(0.5, comp)
    } else {
      diag(comp)
    }


    #############################
    #### The functions c_cal calcuates the  exit probablities for the given exit bounadries e0_1
    c_cal<-function(e_1){
      u_bry <- ifelse(T_new > T_inter, (e_1 - score_inter)/sqrt(T_new - T_inter), -Inf)
      ########################################
      tmp = 1-  pmvnorm(lower=-Inf,upper=u_bry,mean=0, sigma=cov_mtx,abseps=1e-8)[1]

      return(tmp)
    }

  ## the root search can fail when the conditional type I error is not
  ## attainable (NA statistics, degenerate boundaries); return NA instead of
  ## aborting the replicate
  e_bry=tryCatch(
    ifelse(c_cal(MAX_INT)<c_alpha, stats::uniroot(function(x) c_cal(x)-c_alpha,c(-MAX_INT,MAX_INT),tol=1e-5)$root, +Inf),
    error=function(e) NA_real_
  )

  return(e_bry)
}
########
#new_ebry_com(T,T_inter,T_new,c_alpha,s,theta_hat_inter_m2,theta_hat_inter_m2_sd,dist)
#new_e_bry_com(T_inter=c(30,31),T_new=c(30.1,55),c_alpha=0.025,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1.5,dist=1)
#new_ebry_com(T=c(35,40),T_inter=c(30,31),T_new=c(29,30),c_alpha=0.0025,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.4,0.5),dist=1)
#new_ebry_com(T=c(35,40),T_inter=c(30,31),T_new=c(29,45),c_alpha=0.0025,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.4,0.5),dist=1)
#new_ebry_com(T=c(35,40),T_inter=c(30,31),T_new=c(45,50),c_alpha=0.0025,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.4,0.5),dist=1)
#new_ebry_com(T=c(35,40),T_inter=c(30,31),T_new=c(45,50),c_alpha=0.0025,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.4,0.5),dist=2)
#new_ebry_com(T=35,T_inter=30,T_new=45,c_alpha=0.0025,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),theta_hat_inter_m2=1.2,theta_hat_inter_m2_sd=0.4,dist=2)

#########################################################
masd_c_power_new_com<-function(theta,s,new_e_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist){
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
  s_inter=s[length(s)-1]

  T_inter = 1 / (theta_hat_inter_m2_sd * theta_hat_inter_m2_sd)
  T_old=T_inter/s_inter
  T=T_old
  comp=length(theta_hat_inter_m2)
  T_new=T_old*sample_r

  score_inter = theta_hat_inter_m2 / (theta_hat_inter_m2_sd * theta_hat_inter_m2_sd)

  #################################################
  #### the next section constructs the covariance matrix###
  cov_mtx = if (dist == 1) {
    matrix(0.5, comp, comp) + diag(0.5, comp)
  } else {
    diag(comp)
  }


  #############################
  u_bry <- ifelse(T_new > T_inter,
                  (new_e_bry - score_inter)/sqrt(T_new - T_inter)-theta*sqrt(T_new-T_inter),
                  -Inf)

  ########################################

  ## the function for the vector of Brownian motions ###
  ### the function for single Brownian motion###
  c_power_new= tryCatch(
    ifelse(isTRUE(max(T_inter-T_new)>0), 0, 1-pmvnorm(lower=-Inf,upper=u_bry,mean=0,sigma=cov_mtx,abseps=1e-8)[1]),
    error=function(e) NA_real_
  )
  return(c_power_new)
}
##################################
#masd_c_power_new_com(theta,s,new_e_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)
#masd_c_power_new(theta=c(0,0.2),s=c(0.3,0.7,1),inter_vt=2,snew=1, new_c_bry=2.781890,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.4,sample_r=2,dist=1)
#masd_c_power_new(theta=c(0.3,0.2),s=c(0.3,0.7,1),inter_vt=2,snew=1,new_c_bry=c(3.310743, 2.781890),theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=1,sample_r=2,dist=1)
#masd_c_power_new(theta=c(0.3,0.2),s=c(0.3,0.7,1),inter_vt=2,snew=c(0.5,1),new_c_bry=c(3.310743, 2.781890),theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=0.4,sample_r=2,dist=1)
#masd_c_power_new(theta=c(0.01492537,0.01492537),s=c(1/3,2/3,1),inter_vt=2,snew=1,new_c_bry=c(3.310743, 2.781890),theta_hat_inter_m2=c( 0.01492537,0.01492537),theta_hat_inter_m2_sd=c(0.07917485,0.08195375),sample_r=(2/3)*1.0001,dist=1)
#masd_c_power_new_com(theta=c(0.01492537,0.01492537),s=c(1/3,2/3,1),new_e_bry=10,theta_hat_inter_m2=c( 0.01492537,0.01492537),theta_hat_inter_m2_sd=c(0.07917485,0.08195375),sample_r=2^1,dist=1)
#masd_c_power_new_com(theta=c(0.01492537,0.01492537),s=c(1/3,2/3,1),new_e_bry=Inf,theta_hat_inter_m2=c( 0.01492537,0.01492537),theta_hat_inter_m2_sd=c(0.07917485,0.08195375),sample_r=2^1,dist=1)

##################################
masd_new_info_time_com<-function(theta,s,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,dist){
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
  s_inter=s[length(s)-1]

  comp=length(theta_hat_inter_m2)

  T_inter = 1 / (theta_hat_inter_m2_sd * theta_hat_inter_m2_sd)
  T=T_inter/s_inter

  #  c_alpha=masd_c_power_com(s,c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1,theta_hat_inter_m1_sd,dist)


  power<-function(sample_r){

    T_new=T*sample_r
    new_e_bry=new_ebry_com(T,T_inter,T_new,c_alpha,s,theta_hat_inter_m2=theta_hat_inter_m2,theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=dist)
    con_power_new=masd_c_power_new_com(theta=theta_hat_inter_m2,s,new_e_bry,theta_hat_inter_m2,theta_hat_inter_m2_sd,sample_r,dist)

    return(con_power_new)
  }

  #new_sample_r=stats::uniroot(function(x) power(x)-(1-beta),c(s_inter*1.0001,s_inter * 100), extendInt="upX",tol=1e-4)$root

  new_sample_r = tryCatch(
    stats::uniroot(function(x) power(x)-(1-beta), c(s_inter*1.0001, s_inter*100), extendInt="upX", tol=1e-3)$root,
    error = function(e) 1
  )
  T_new =T*new_sample_r
  out=list(T_new,new_sample_r)

  names(out)=c("New information Time","new sample size ratio")
  return(out)
}
#################################
#masd_new_info_time_com(theta,s,c_alpha,theta_hat_inter_m2,theta_hat_inter_m2_sd,beta,dist)
#masd_new_info_time_com(theta=c(0,0.2),s=c(0.3,0.7,1),c_alpha=0.4,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.2,0.3),beta=0.1,dist=1)
#masd_new_info_time_com(theta=c(0,0.2),s=c(0.3,0.7,1),c_alpha=0.4,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.2,0.3),beta=0.1,dist=2)
##################################
#new_ebry_com(T_inter=c(25.00, 11.11),T_new=c(32.66, 14.51),c_alpha=0.04157,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.2,0.3),dist=1)
#new_ebry_com(T_inter=c(25.00, 11.11),T_new=c(32.66, 14.51),c_alpha=0.04157,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.2,0.3),dist=2)
#new_ebry_com(T_inter=c(25.00, 11.11),T_new=c(31.31, 13.92),c_alpha=0.04157,theta_hat_inter_m2=c(1.2,1.5),theta_hat_inter_m2_sd=c(0.2,0.3),dist=2)
################

######################################################
con_density <- function(theta_hat_inter,theta_hat_inter_sd,dist){

  comp=length(theta_hat_inter)
  #theta=theta_hat_inter_m1
  theta_max=max(theta_hat_inter)
  T_inter = 1 / (theta_hat_inter_sd * theta_hat_inter_sd)
  ########################################
  #### the next section constructs the covariance matrix###
  cov_mtx = if (dist == 1) {
    matrix(0.5, comp, comp) + diag(0.5, comp)
  } else {
    diag(comp)
  }

  #####################################################
  ### boundary for the  vector of Brownian  motions####

  u_bry=theta_max/sqrt(T_inter)

  #################################################
  c_density= tryCatch(1-pmvnorm(lower=-Inf,upper=u_bry,mean=0,sigma=cov_mtx,abseps=1e-8)[1],
                      error=function(e) NA_real_)
  return(c_density)


}
##################
#con_density(theta_hat_inter,theta_hat_inter_sd,dist)
#con_density(theta_hat_inter=1.2,theta_hat_inter_sd=0.2, dist=1)
#con_density(theta_hat_inter=c(1.2,1.5),theta_hat_inter_sd=c(0.2,0.3),dist=1)
#con_density(theta_hat_inter=c(1.2,1.5,1.4),theta_hat_inter_sd=c(0.2,0.3,0.3),dist=1)
#############################################################
### Part II. Simulation functions        ###################
#############################################################

#########################################################
####  II.1   Simulations for normal distributions  ######
#########################################################

####################
msd_rej_rate_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,samsz,sim_num){
  #mu is the vector of means for each dose ##
  ## the effect size for control=0=mu[length(mu)]##
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size

  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(mu_t)
  N=ceiling(samsz*s)
  t<- N/(sigma_t^2+sigma_c^2)
  T=t[length(s)]
  N1=N[1]
  NF=samsz

  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej=0; stg=0; samsz_simu_control=N[length(s)]
    samsz_simu_tot=N[length(s)]*(comp+1); arm_rej=rep(0,comp)

      #the loop specifies what is generated for each simulation
      ########################################
      ## step 1: generating random samples ###
      ########################################


      x0=stats::rnorm(NF,mu_c,sigma_c)
      x=matrix(nrow=comp,ncol=NF)
      for (i in 1:comp){
        x[i,]=stats::rnorm(NF, mu_t[i],sigma_t)
      }
      ####  analysis parameters ###

      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      var_hat_m1=matrix(nrow=comp,ncol=length(s))

      if (direction ==1){
        for (j in 1:length(s)){
          for (i in 1:comp){
            theta_hat_m1[i,j]=mean(x[i,1:N[j]])-mean(x0[1:N[j]])
            var_hat_m1[i,j]=stats::var(x[i,1:N[j]])+stats::var(x0[1:N[j]])
            t_hat_m1[i,j]=N[j]/var_hat_m1[i,j]
            theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])
           wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
          }
        }
      }else if (direction ==2){
        for (j in 1:length(s)){
          for (i in 1:comp){
            theta_hat_m1[i,j]=mean(x0[1:N[j]])-mean(x[i,1:N[j]])
            var_hat_m1[i,j]=stats::var(x[i,1:N[j]])+stats::var(x0[1:N[j]])
            t_hat_m1[i,j]=N[j]/var_hat_m1[i,j]
            theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])
            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
          }
        }
      }
       wald_max_m1=NULL

      for (j in 1:length(s) ){
        wald_max_m1[j]=max(wald_m1[,j])
      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################

      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry
  # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################

      if (isTRUE(max(test_max)>=0)){ ##null hypothesis rejected, no more random sampling

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max[test_max>0]
        stg=match(test_max_t[1],test_max)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max)<0)) {

            rej=0
            stg=0
          last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)

        }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej,stg,samsz_simu_control,samsz_simu_tot,arm_rej)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr

  rej = as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  samsz_simu_control=as.numeric(simr[,3])
  samsz_simu_tot=simr[,4]
  arm_rej=simr[,5:(4+length(mu_t))]
  # unique(case)
  ######################################

  if (direction==1){
    mu=mu_t-mu_c
  }else if (direction==2){
    mu=mu_c-mu_t
  }
  rej_rt=mean(rej)

  arm_rej_rt=NULL

  for (i in 1:length(mu_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  arm_summary=cbind(mu_t,mu_c,arm_rej_rt)
  colnames(arm_summary)=c(LABEL()$mu_t1,LABEL()$mu_c,"Rejection rate by arm")


  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  Design=cbind(c(1:length(s)),s,c_bry,rej_time)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look")

  out_1=c( rej_rt,tyI_er_rt,pos_rej_rt,samsz,tot_samsz,mean_control_samsz,mean_tot_samsz)
  names(out_1)=c("Rejection rate","typr I error","non-type I error rejection","Planned control sample size",
                 "Planned total sample size","mean control sample size","mean total sample size")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")
  return(out)

}

##############################
#msd_rej_rate_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,samsz,sim_num)
#msd_rej_rate_normal(mu_t=c(0,0,0),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.115204, 2.909888, 2.375914),samsz=100,sim_num=10000)

################################################

##############################################

#' Multiple comparisons Group sequential design: Simulations, Normal distribution, Type I error and power
#'
#' @param mu_t Vector of true means of the test arms, one per dose
#' @param mu_c Mean of the control arm
#' @param sigma_t Standard deviation of the test arm; must be > 0
#' @param sigma_c Standard deviation of the control arm; must be > 0
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param samsz Sample size of the control arm used in the simulation; must be > 0
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary and the rejection rate by look}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the test-arm mean, the control-arm mean and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate, the simulated type I error rate, the non-type I error rejection rate, the planned control and total sample sizes and the mean simulated control and total sample sizes}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' msd_power_simu_normal(mu_t = c(0,0,0), mu_c = 0, sigma_t = 1, sigma_c = 1, direction = 1,
#' s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
msd_power_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,samsz,sim_num){
   # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  mu=mu_t-mu_c
  comp=length(mu)
  if (bry_type==1) {

    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out=msd_rej_rate_normal(mu_t=mu_t,mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,
                               s=s,c_bry=c_bry,samsz=samsz,sim_num=sim_num)

  out=simu_out
  return(out)
}

##################################
#msd_power_simu_normal(mu_t,mu_c,s,bry_type,alpha_2,sigma,samsz,sim_num)
#
#msd_power_simu_normal(mu_t=c(0,0,0),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)
#msd_power_simu_normal(mu_t=c(0.1,0.2,0),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)


#######################################

###################
masd_rej_simu_adapt_com_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,sim_num,dist){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size

  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(mu_t)
  N=ceiling(samsz*s)
  t<- N/(sigma_t^2+sigma_c^2)
  T=t[length(s)]
  N1=N[1]
  NF=samsz

  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    samsz_simu_control=N[length(s)]; samsz_simu_tot=N[length(s)]*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)
      #the loop specifies what is generated for each simulation
      ########################################
      ## step 1: generating random samples ###
      ########################################
      x0=stats::rnorm(N_max,mu_c,sigma_c)
      x=matrix(nrow=comp,ncol=N_max)
      for (i in 1:comp){
        x[i,]=stats::rnorm(N_max, mu_t[i],sigma_t)
      }
      ####  analysis parameters ###

      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      var_hat_m1=matrix(nrow=comp,ncol=length(s))

      if (direction ==1){

        for (j in 1:length(s)){

          for (i in 1:comp){

            theta_hat_m1[i,j]=mean(x[i,1:N[j]])-mean(x0[1:N[j]])
            var_hat_m1[i,j]=stats::var(x[i,1:N[j]])+stats::var(x0[1:N[j]])
            t_hat_m1[i,j]=N[j]/var_hat_m1[i,j]
            #      t_hat_m1[i,j]=N[j]/2/sigma^2 #######  accurate information time
            theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])

            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])

          }
        }
      }else if (direction ==2){
        for (j in 1:length(s)){
          for (i in 1:comp){
            theta_hat_m1[i,j]=mean(x0[1:N[j]])-mean(x[i,1:N[j]])
            var_hat_m1[i,j]=stats::var(x[i,1:N[j]])+stats::var(x0[1:N[j]])
            t_hat_m1[i,j]=N[j]/var_hat_m1[i,j]
            # t_hat_m1[i,j]=N[j]/2/sigma^2 #######  accurate information time
            theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])

            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
          }
        }
      }
      T_hat=t_hat_m1[,length(s)]
      wald_max_m1=NULL

      for (j in 1:length(s) ){
        wald_max_m1[j]=max(wald_m1[,j])
      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop=rep(0,length(mu_t))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }
      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################
      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry
      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################
      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }

      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling
        futi=0
        rej=1
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        last_vt=stg
        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)
      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        #########################################################
        ### Step 4.1. If max(theta<=theta_cut),                ##
        ###     Continue to planned sample size and stop trial  ##
        #########################################################
        if (futi==1) { # If TRUE, there is not satisfactory efficacy, stop trial at planned sample size
          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            rej_adap=0
            stg=0}
          last_vt=length(s)
          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)
          mu_reduc=0
        }else if (futi==0){ #beginning of step 4.2

          ########################################################################
          ### Step 4.2. If max(theta>theta_cut), calculate conditional power    ##
          ###          and re-estiate sample size                               ##
          ########################################################################

          ########################################################################################
          ## Step 4.2.1. If futi==0. Check if there is need to drop some doses due to efficacy ###
          ### A dose is dropped for lack of efficacy if: theta_hat<theta cut                   ###
          ### or of min theta_hat is less than 1/3 of the max theta_hat                        ###
          ########################################################################################
          mu_drop=match(arm_drop[arm_drop>0],arm_drop)
          if (length(mu_drop)==0){
            mu_keep=c(1:comp)
          }else if ( length(mu_drop)>0){
            mu_keep=c(1:comp)[-mu_drop]
          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_hat_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2
          n_control_inter=N[length(s)-1]


          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1=theta_hat_inter_m1,theta_hat_inter_m1_sd=theta_hat_inter_m1_sd,dist=dist)
          if (isTRUE(c_power_m2>= 1-beta)){
            #########################################
            ### Proceed to planned sample size  #####
            ### and conduct final testing       #####
            #########################################
            last_vt=length(s)
            wald_m2=wald_m1[mu_keep,]
            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last
            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }
            samsz_simu_control=N[last_vt]
            samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)   #For q-th (collumn number is q) simulation,
          }else if (isTRUE(c_power_m2<1-beta)){
            t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                           theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)
            ## default when the re-estimation is not evaluable
            N_new_t=samsz
            if (isTRUE(as.numeric(t_new_0[2])<=1)){
              N_new_t=samsz
            }else if (isTRUE(as.numeric(t_new_0[2])>1)){
              T_new_0=as.numeric(t_new_0[[1]])
              N_new_t_0=max(ceiling(n_control_inter*(T_new_0/T_inter)))
              N_new_t_1=max(samsz,N_new_t_0)
              N_new_t= min(N_new_t_1,N_max)
            }

            if(isTRUE(N_new_t==samsz)){
              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)
              samsz_simu_control=N[last_vt]
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)
            }else  if(isTRUE(N_new_t>samsz)){
              ####################################################################################
              ###  Step 4.2.1.3.1: Generate new random samples according to the new sample size ##
              ####################################################################################
              x0_new=x0[1:N_new_t]
              x_new=x[c(mu_keep),1:N_new_t]
              if (isTRUE(length(mu_keep)==1)){
                if (direction==1){
                  theta_hat_m2_new=mean(x_new)-mean(x0_new)
                }else if (direction==2){
                  theta_hat_m2_new=mean(x0_new)-mean(x_new)
                }
                var_hat_m2_new=stats::var(x_new)+stats::var(x0_new)
                t_hat_m2_new=N_new_t/var_hat_m2_new
                theta_hat_m2_sd_new=1/sqrt(t_hat_m2_new)
                wald_m2_new=theta_hat_m2_new*sqrt(t_hat_m2_new)
                wald_max_m2_new=max(wald_m2_new)
              }else if (isTRUE(length(mu_keep)>1)){
                t_hat_m2_new=NULL
                theta_hat_m2_new=NULL
                theta_hat_m2_sd_new=NULL
                wald_m2_new=NULL
                var_hat_m2_new=NULL

                if (direction==1){
                  for (i in 1:length(mu_keep)){
                    theta_hat_m2_new[i]=mean(x_new[i,])-mean(x0_new)
                    var_hat_m2_new[i]=stats::var(x_new[i,])+stats::var(x0_new)
                    t_hat_m2_new[i]=N_new_t/var_hat_m2_new[i]
                    # t_hat_m2_new[i]=N_new_t/2/sigma^2 #accurate information time
                    theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                    wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])
                  }
                }else if (direction==2){
                  ## BUGFIX: the original code lost its 'for' loop in this branch,
                  ## so 'i' leaked from an earlier loop and could exceed nrow(x_new),
                  ## aborting the replicate with 'subscript out of bounds'.
                  for (i in 1:length(mu_keep)){
                  theta_hat_m2_new[i]=mean(x0_new)-mean(x_new[i,])
                  var_hat_m2_new[i]=stats::var(x_new[i,])+stats::var(x0_new)
                  t_hat_m2_new[i]=N_new_t/var_hat_m2_new[i]
                  theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])

                  wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])
                  }

                }

                wald_max_m2_new=max(wald_m2_new)
              }

              t_hat_m2_new_max=max(t_hat_m2_new)
              T_new=t_hat_m2_new
              T_inter=(N[length(s)-1]/N_new_t)*T_new
              T_hat_m2=T_hat[mu_keep]

              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=1)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=1)
              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=1)

              c_bry_new=NULL

              for (i in 1:length(theta_hat_inter_m2)){
                c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
              }

              test_new=wald_m2_new-max(c_bry_new)
              #################################################################################
              ## Step 4.2.1.3.3. Hypothesis testing: compare the Wald statistics with the critical bundaries ##
              #################################################################################

              samsz_simu_control=N_new_t
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N_new_t-N[length(s)-1])*(length(mu_keep)+1)

              if (isTRUE(max(test_new)>0)){

                rej=1

                stg=length(s)

              }else{
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0

              }

              arm_rej_new=NULL

              for (i in 1:length(mu_keep)){

                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }

              }

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }
        }

      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej_nossr,rej,stg,samsz_simu_control,samsz_simu_tot,futi,arm_drop,arm_rej)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr

  rej_nossr = as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  samsz_simu_control=as.numeric(simr[,4])
  samsz_simu_tot=simr[,5]
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(mu_t))]
  arm_rej=simr[,(7+length(mu_t)):(6+2*length(mu_t))]
  ######################################

  if (direction==1){
    mu=mu_t-mu_c
  }else if (direction==2){
    mu=mu_c-mu_t
  }
  rej_rt=mean(rej)
  rej_nossr_rt=mean(rej_nossr)

  arm_rej_rt=NULL

  for (i in 1:length(mu_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  arm_drop_rt=NULL

  for (i in 1:length(mu_t)){
    arm_drop_rt[i]=mean(arm_drop[,i])
  }

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
  }

  futi_rt=mean(futi)

  arm_summary=cbind(mu_t,mu_c,arm_drop_rt,arm_rej_rt)
  colnames(arm_summary)=c(LABEL()$mu_t1,LABEL()$mu_c,"drop rate by arm","Rejection rate by arm")

  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)
  max_tot_samsz=N_max*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu_control[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= samsz_simu_control[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  #########################
  Design=cbind(c(1:length(s)),s,c_bry,rej_time,samsz_icr_byLook)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look","Rate of sample size increase")

  out_1=c(rej_nossr_rt, rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,samsz,N_max, mean_control_samsz, mean_samsz_simu_control_samsz_icr)
  names(out_1)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "type I error","non-type I error rejection","futility rate","Planned control sample size",
                 "Maximum control sample size",
                 "mean control sample size","mean control sample size conditional on increase")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")
  return(out)

}

###################
#masd_rej_simu_adapt_com_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,sim_num,dist)
#masd_rej_simu_adapt_com_normal(mu_t=c(0,0,0),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.115204, 2.909888, 2.375914),theta_cut=0.02,beta=0.1,samsz=100,N_max=200,sim_num=10000,dist=1)
#masd_rej_simu_adapt_com_normal(mu_t=c(0,0),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.115204, 2.909888, 2.375914),theta_cut=0.02,beta=0.1,samsz=100,N_max=200,sim_num=10000,dist=1)
#masd_rej_simu_adapt_com_normal(mu_t=c(0,0),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.115204, 2.909888, 2.375914),theta_cut=0.02,beta=0.1,samsz=100,N_max=200,sim_num=1000000,dist=1)
#masd_rej_simu_adapt_com_normal(mu_t=c(0,0.3,0.4),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.115204, 2.909888, 2.375914),theta_cut=0.02,beta=0.1,samsz=100,N_max=200,sim_num=10000,dist=1)


#masd_rej_simu_adapt_com_normal(mu=c(0,0),s=c(1/2,1),c_bry=c(3.142524,2.222100),sigma=1,theta_cut=0.05,beta=0.1,samsz=100,N_max=450,sim_num=100000,dist=1)
################################################

##############################################

#' Multiple comparisons Group sequential design: Simulations, Normal distribution, Operating characteristics
#'
#' @param mu_t Vector of true means of the test arms, one per dose; setting mu_t = mu_c gives the simulated type I error rate
#' @param mu_c Mean of the control arm
#' @param sigma_t Standard deviation of the test arm; must be > 0
#' @param sigma_c Standard deviation of the control arm; must be > 0
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param theta_cut Lack-of-efficacy threshold on the effect-size scale: a test arm is
#'   dropped at the interim analysis when its observed effect size is below this value
#' @param samsz Sample size of the control arm used in the simulation; must be > 0
#' @param N_max Maximum sample size for the control arm allowed by the sample size
#'   re-estimation; must be >= samsz
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary, the rejection rate by look and the rate of sample size increase}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the test-arm mean, the control-arm mean, the drop rate by arm and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate without and with sample size re-estimation, the rate of sample size increase, the rejection rates conditional on meeting the sample size increase criteria, the simulated type I error rate, the non-type I error rejection rate, the futility rate, the planned and maximum control sample sizes, the mean simulated control sample size and the mean simulated control sample size conditional on an increase}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' masd_power_simu_adapt_normal(mu_t = c(0,0,0), mu_c = 0, sigma_t = 1, sigma_c = 1, direction = 1,
#' s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, beta = 0.1, theta_cut = 0.01, samsz = 100,
#' N_max = 500, sim_num = 20)
#' }
#'
#' @export
masd_power_simu_adapt_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=length(mu_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out=masd_rej_simu_adapt_com_normal(mu_t=mu_t,mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,
                                          s=s,c_bry=c_bry,theta_cut=theta_cut,beta=beta,samsz=samsz,N_max=N_max,sim_num=sim_num,dist=1)

  out=simu_out
  #names(out)=c("simulation summary","planned control","planned total")
  return(out)
}

##################################
#masd_power_simu_adapt_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sim_num)

#masd_power_simu_adapt_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=10000)

#masd_power_simu_adapt_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=10000)

#masd_power_simu_adapt_normal(mu_t=c(0,0.3,0.4),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=1000)



###################################################################################
## Simulations for Optimizing Operating characteristics - normal distributions  ###
###################################################################################
###################

masd_rej_OC_simu_adapt_com_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,sim_num,dist){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size

  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(mu_t)
  N=ceiling(samsz*s)
  t<- N/(sigma_t^2+sigma_c^2)
  T=t[length(s)]
  N1=N[1]
  NF=samsz

  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    samsz_simu_control=N[length(s)]; samsz_simu_tot=N[length(s)]*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)

      #the loop specifies what is generated for each simulation
      ########################################
      ## step 1: generating random samples ###
      ########################################

      x0=stats::rnorm(N_max,mu_c,sigma_c)
      x=matrix(nrow=comp,ncol=N_max)
      for (i in 1:comp){
        x[i,]=stats::rnorm(N_max, mu_t[i],sigma_t)
      }
      ####  analysis parameters ###

      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      var_hat_m1=matrix(nrow=comp,ncol=length(s))

      if (direction ==1){

        for (j in 1:length(s)){

          for (i in 1:comp){

            theta_hat_m1[i,j]=mean(x[i,1:N[j]])-mean(x0[1:N[j]])
            var_hat_m1[i,j]=stats::var(x[i,1:N[j]])+stats::var(x0[1:N[j]])
            t_hat_m1[i,j]=N[j]/var_hat_m1[i,j]
            #      t_hat_m1[i,j]=N[j]/2/sigma^2 #######  accurate information time
            theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])

            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])

          }
        }
      }else if (direction ==2){
        for (j in 1:length(s)){
          for (i in 1:comp){
            theta_hat_m1[i,j]=mean(x0[1:N[j]])-mean(x[i,1:N[j]])
            var_hat_m1[i,j]=stats::var(x[i,1:N[j]])+stats::var(x0[1:N[j]])
            t_hat_m1[i,j]=N[j]/var_hat_m1[i,j]
            # t_hat_m1[i,j]=N[j]/2/sigma^2 #######  accurate information time
            theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])

            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
          }
        }
      }
      T_hat=t_hat_m1[,length(s)]
      wald_max_m1=NULL

      for (j in 1:length(s) ){
        wald_max_m1[j]=max(wald_m1[,j])
      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop = ifelse(theta_hat_inter_m1<theta_cut, 1, 0)

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################

      test_all = wald_m1-c_bry
      rej_all = ifelse(test_all>=0, 1, 0)
      test_max=wald_max_m1-c_bry
      test_max_1=test_max[1:length(s)-1]

      futi = ifelse(max( theta_hat_inter_m1)<theta_cut, 1, 0)

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################
      rej_nossr = ifelse(max(test_max)>=0, 1, 0)

      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling

        futi=0
        rej=1
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        #########################################################
        ### Step 4.1. If max(theta<=theta_cut),                ##
        ###     Continue to planned sample size and stop trial  ##
        #########################################################

        if (futi==1) { # If TRUE, there is not satisfactory efficacy, stop trial at planned sample size


          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            rej_adap=0
            stg=0}
          last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)
          mu_reduc=0
        }else if (futi==0){ #beginning of step 4.2

          ########################################################################
          ### Step 4.2. If max(theta>theta_cut), calculate conditional power    ##
          ###          and re-estiate sample size                               ##
          ########################################################################

          ########################################################################################
          ## Step 4.2.1. If futi==0. Check if there is need to drop some doses due to efficacy ###
          ### A dose is dropped for lack of efficacy if: theta_hat<theta cut                   ###
          ### or of min theta_hat is less than 1/3 of the max theta_hat                        ###
          ########################################################################################

          mu_drop=match(arm_drop[arm_drop>0],arm_drop)

          if (length(mu_drop)==0){

            mu_keep=c(1:comp)

          }else if ( length(mu_drop)>0){

            mu_keep=c(1:comp)[-mu_drop]

          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_hat_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2
          n_control_inter=N[length(s)-1]


          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1=theta_hat_inter_m1,theta_hat_inter_m1_sd=theta_hat_inter_m1_sd,dist=dist)
          if (isTRUE(c_power_m2>= 1-beta)){
            #########################################
            ### Proceed to planned sample size  #####
            ### and conduct final testing       #####
            #########################################

            last_vt=length(s)

            wald_m2=wald_m1[mu_keep,]

            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last

            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }

            samsz_simu_control=N[last_vt]
            samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)   #For q-th (collumn number is q) simulation,

          }else if (isTRUE(c_power_m2<1-beta)){

            t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                           theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

            ## default when the re-estimation is not evaluable
            N_new_t=samsz
            if (isTRUE(as.numeric(t_new_0[2])<=1)){
              N_new_t=samsz
            }else if (isTRUE(as.numeric(t_new_0[2])>1)){
              T_new_0=as.numeric(t_new_0[[1]])

              N_new_t_0=max(ceiling(n_control_inter*(T_new_0/T_inter)))
              N_new_t_1=max(samsz,N_new_t_0)
              N_new_t= min(N_new_t_1,N_max)
            }

            if(isTRUE(N_new_t==samsz)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              samsz_simu_control=N[last_vt]
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)

            }else  if(isTRUE(N_new_t>samsz)){

              ####################################################################################
              ###  Step 4.2.1.3.1: Generate new random samples according to the new sample size ##
              ####################################################################################

              x0_new=x0[1:N_new_t]
              x_new=x[c(mu_keep),1:N_new_t]

              if (isTRUE(length(mu_keep)==1)){

                if (direction==1){
                  theta_hat_m2_new=mean(x_new)-mean(x0_new)
                }else if (direction==2){
                  theta_hat_m2_new=mean(x0_new)-mean(x_new)
                }
                var_hat_m2_new=stats::var(x_new)+stats::var(x0_new)

                t_hat_m2_new=N_new_t/var_hat_m2_new
                theta_hat_m2_sd_new=1/sqrt(t_hat_m2_new)

                wald_m2_new=theta_hat_m2_new*sqrt(t_hat_m2_new)

                wald_max_m2_new=max(wald_m2_new)

              }else if (isTRUE(length(mu_keep)>1)){
                t_hat_m2_new=NULL
                theta_hat_m2_new=NULL
                theta_hat_m2_sd_new=NULL
                wald_m2_new=NULL
                var_hat_m2_new=NULL

                if (direction==1){

                  for (i in 1:length(mu_keep)){

                    theta_hat_m2_new[i]=mean(x_new[i,])-mean(x0_new)
                    var_hat_m2_new[i]=stats::var(x_new[i,])+stats::var(x0_new)
                    t_hat_m2_new[i]=N_new_t/var_hat_m2_new[i]
                    # t_hat_m2_new[i]=N_new_t/2/sigma^2 #accurate information time
                    theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])

                    wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                  }

                }else if (direction==2){
                  ## BUGFIX: the original code lost its 'for' loop in this branch,
                  ## so 'i' leaked from an earlier loop and could exceed nrow(x_new),
                  ## aborting the replicate with 'subscript out of bounds'.
                  for (i in 1:length(mu_keep)){
                  theta_hat_m2_new[i]=mean(x0_new)-mean(x_new[i,])
                  var_hat_m2_new[i]=stats::var(x_new[i,])+stats::var(x0_new)
                  t_hat_m2_new[i]=N_new_t/var_hat_m2_new[i]
                  theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])

                  wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])
                  }

                }

                wald_max_m2_new=max(wald_m2_new)
              }

              t_hat_m2_new_max=max(t_hat_m2_new)
              T_new=t_hat_m2_new
              T_inter=(N[length(s)-1]/N_new_t)*T_new
              T_hat_m2=T_hat[mu_keep]

              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=1)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=1)
              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=1)

              c_bry_new=(e_bry_new-score_inter_m2)/sqrt(T_new-T_inter)

              test_new=wald_m2_new-max(c_bry_new)

              #################################################################################
              ## Step 4.2.1.3.3. Hypothesis testing: compare the Wald statistics with the critical bundaries ##
              #################################################################################

              samsz_simu_control=N_new_t
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N_new_t-N[length(s)-1])*(length(mu_keep)+1)

              if (isTRUE(max(test_new)>0)){

                rej=1

                stg=length(s)

              }else{
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0

              }

              arm_rej_new = ifelse(test_new>=0, 1, 0)

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }
        }

      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej_nossr,rej,stg,samsz_simu_control,samsz_simu_tot,futi,arm_drop,arm_rej)

  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr

  rej_nossr = as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  samsz_simu_control=as.numeric(simr[,4])
  samsz_simu_tot=simr[,5]
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(mu_t))]
  arm_rej=simr[,(7+length(mu_t)):(6+2*length(mu_t))]
  ######################################

  if (direction==1){
    mu=mu_t-mu_c
  }else if (direction==2){
    mu=mu_c-mu_t
  }
  rej_rt=mean(rej)
  rej_nossr_rt=mean(rej_nossr)

  arm_rej_rt = colMeans(arm_rej)

  futi_rt=mean(futi)

  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)
  max_tot_samsz=N_max*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  # pos_rej=NULL
  # for (i in 1: dim(simr)[1]){
  #
  #   if (sum(arm_pos_rej[i,])>0){
  #     pos_rej[i]=1
  #   }else{
  #     pos_rej[i]=0
  #   }
  # }
  pos_rej <- as.integer(rowSums(arm_pos_rej) > 0)
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL
  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr = ifelse(samsz_simu_control>samsz, 1, 0)
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= samsz_simu_control[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  #########################

  out=c(rej_nossr_rt, rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,samsz,N_max, mean_control_samsz, mean_samsz_simu_control_samsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "type I error","non-type I error rejection","futility rate","Planned control sample size",
                 "Maximum control sample size",
                 "mean control sample size","mean control sample size conditional on increase")

  return(out)

}

#masd_rej_OC_simu_adapt_com_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,sim_num,dist)
#masd_rej_OC_simu_adapt_com_normal(mu_t=c(0,0,0),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.115204, 2.909888, 2.375914),theta_cut=0.02,beta=0.1,samsz=100,N_max=200,sim_num=1000,dist=1)
#masd_rej_OC_simu_adapt_com_normal(mu_t=c(0,0),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.115204, 2.909888, 2.375914),theta_cut=0.02,beta=0.1,samsz=100,N_max=200,sim_num=10000,dist=1)
#masd_rej_OC_simu_adapt_com_normal(mu_t=c(0.2,0.3,0.4),mu_c=0, sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.115204, 2.909888, 2.375914),theta_cut=0.02,beta=0.1,samsz=100,N_max=200,sim_num=10000,dist=1)
#
################################################

masd_OC_power_simu_adapt_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=length(mu_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out=masd_rej_OC_simu_adapt_com_normal(mu_t=mu_t,mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,
                                          s=s,c_bry=c_bry,theta_cut=theta_cut,beta=beta,samsz=samsz,N_max=N_max,sim_num=sim_num,dist=1)

  mu_t_names=paste(paste(LABEL()$mu_t,"[",c(1:length(mu_t)),"]",sep=""))
  out_1= round(c(mu_t,mu_c,as.numeric(simu_out)),4)
  names(out_1)=c(mu_t_names,LABEL()$mu_c, names(simu_out))
  # out=data.frame(out_1)
  out=out_1
  return(out)
}

##################################
#masd_OC_power_simu_adapt_normal(mu_t,mu_c,s,bry_type,alpha_2,theta_cut,samsz,N_max,sigma,sim_num)
#masd_OC_power_simu_adapt_normal(mu_t=c(0.4,0.4),mu_c=c(0,0),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=100)
#masd_OC_power_simu_adapt_normal(mu_t=c(0.4,0.4),mu_c=c(0,0),s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=100)
#masd_OC_power_simu_adapt_normal(mu_t=c(0,0),mu_c=c(0,0),s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=100000)
#masd_OC_power_simu_adapt_normal(mu_t=c(0,0,0.3),mu_c=c(0,0,0),s=c(1/2,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=100000)
#masd_OC_power_simu_adapt_normal(mu_t=c(0,0,0.3),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=10000)

#masd_OC_power_simu_adapt_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=1000)

#masd_OC_power_simu_adapt_normal(mu_t=c(0.2,0.3,0.4),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=1000)
######################################################

###################
#' Multiple comparisons Group sequential design: Simulations, Normal distribution, Optimizing adaptive sequential design
#'
#' @param groups Matrix of true means of the test arms: each row is one simulation scenario and each column one test arm/dose
#' @param mu_c Mean of the control arm
#' @param sigma_t Standard deviation of the test arm; must be > 0
#' @param sigma_c Standard deviation of the control arm; must be > 0
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g. c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Two-sided type I error levels spent at each look of the schedule s1; must sum to alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error levels spent at each look of the schedule s2; must sum to alpha_2
#' @param theta_cut Lack-of-efficacy threshold on the effect-size scale: a test arm is
#'   dropped at the interim analysis when its observed effect size is below this value
#' @param samsz_1 Sample size of the control arm used in the simulation for the first sample size scenario; must be > 0
#' @param samsz_2 Sample size of the control arm used in the simulation for the second sample size scenario; must be > 0
#' @param N_max_1 Maximum sample size for the control arm allowed by the sample size re-estimation in the first N_max scenario; must be >= samsz_1
#' @param N_max_2 Maximum sample size for the control arm allowed by the sample size re-estimation in the second N_max scenario; must be >= samsz_2
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the O'Brien-Fleming group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the O'Brien-Fleming group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the alpha-spending group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the alpha-spending group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' masd_simu_OC_normal(groups = rbind(c(0.3,0.4,0.5),c(0.1,0.2,0.3),c(0.3,0.4,0.5),c(0.1,0.2,0.3)),
#' mu_c = 0.15, sigma_t = 1, sigma_c = 1, direction = 1, s1 = c(1/3,2/3,1), s2 = c(0.5,1),
#' alpha_2 = 0.05, alpha_2_alps_s1 = c(0.001,0.001,0.048), alpha_2_alps_s2 = c(0.001,0.049),
#' theta_cut = 0.01, samsz_1 = 50, samsz_2 = 75, N_max_1 = 1000, N_max_2 = 1500, beta = 0.1,
#' sim_num = 2)
#' }
#'
#' @export
masd_simu_OC_normal<-function(groups,mu_c,sigma_t,sigma_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)
  comp=length(groups[1,])

  dist=1
  ### critical boundaries for s1 ############


  s1_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=dist)[[1]][1]))
  s2_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=dist)[[1]][1]))

  s1_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry (s=s1,alpha=alpha_alsp_s1,comp=comp,dist=dist)[[1]][1]))
  s2_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry (s=s2,alpha=alpha_alsp_s2,comp=comp,dist=dist)[[1]][1]))

  ##### Two-stage GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")
  ################
  mu_t_vec = groups
  rownames( mu_t_vec)=NULL

  u_t_group_names=paste(paste(LABEL()$mu_t_group,c(1:dim(mu_t_vec)[1]),sep=" "))
  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_OF_s1_samsz_1_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s1_samsz_1_N_max_1=rbind(simu_out_OF_s1_samsz_1_N_max_1,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))

  }
  rownames(simu_out_OF_s1_samsz_1_N_max_1)=u_t_group_names


  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_1=rbind(simu_out_OF_s2_samsz_1_N_max_1,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_1)=u_t_group_names

  ## samsz_2###


  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_1=rbind(simu_out_OF_s1_samsz_2_N_max_1,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_1)=u_t_group_names


  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_1=rbind(simu_out_OF_s2_samsz_2_N_max_1,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_2_N_max_1)=u_t_group_names
  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_OF_s1_samsz_1_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s1_samsz_1_N_max_2=rbind(simu_out_OF_s1_samsz_1_N_max_2,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_2)=u_t_group_names


  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_2=rbind(simu_out_OF_s2_samsz_1_N_max_2,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_2)=u_t_group_names

  ## samsz_2###


  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_2=rbind(simu_out_OF_s1_samsz_2_N_max_2,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_2)=u_t_group_names


  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_2=rbind(simu_out_OF_s2_samsz_2_N_max_2,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_2_N_max_2)=u_t_group_names


  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s1_samsz_1_N_max_1=rbind(simu_out_alsp_s1_samsz_1_N_max_1,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                               bry_type=2,alpha_2=alpha_alsp_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_1)=u_t_group_names


  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_1=rbind(simu_out_alsp_s2_samsz_1_N_max_1,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                               bry_type=2,alpha_2=alpha_alsp_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_1)=u_t_group_names

  ## samsz_2###


  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_1=rbind(simu_out_alsp_s1_samsz_2_N_max_1,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                               bry_type=2,alpha_2=alpha_alsp_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_1)=u_t_group_names


  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_1=rbind(simu_out_alsp_s2_samsz_2_N_max_1,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                               bry_type=2,alpha_2=alpha_alsp_s2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_2_N_max_1)=u_t_group_names
  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s1_samsz_1_N_max_2=rbind(simu_out_alsp_s1_samsz_1_N_max_2,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                               bry_type=2,alpha_2=alpha_alsp_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_2)=u_t_group_names


  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_2=rbind(simu_out_alsp_s2_samsz_1_N_max_2,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                               bry_type=2,alpha_2=alpha_alsp_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_2)=u_t_group_names

  ## samsz_2###


  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_2=rbind(simu_out_alsp_s1_samsz_2_N_max_2,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                               bry_type=2,alpha_2=alpha_alsp_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_2)=u_t_group_names


  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_2=rbind(simu_out_alsp_s2_samsz_2_N_max_2,masd_OC_power_simu_adapt_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                               bry_type=2,alpha_2=alpha_alsp_s2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_2_N_max_2)=u_t_group_names

  ###############################
  ##############

  list_OF_1=list(s1_OF_design,
                 simu_out_OF_s1_samsz_1_N_max_1,
                 simu_out_OF_s1_samsz_2_N_max_1,
                 simu_out_OF_s1_samsz_1_N_max_2,
                 simu_out_OF_s1_samsz_2_N_max_2)
  names(list_OF_1)=c("Group sequential design",
                     "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                     "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                     "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                     "Simulation results: Information fractions 1, Sample size 2, N_max 2")
  list_OF_2=list(s2_OF_design,
                 simu_out_OF_s2_samsz_1_N_max_1,
                 simu_out_OF_s2_samsz_2_N_max_1,
                 simu_out_OF_s2_samsz_1_N_max_2,
                 simu_out_OF_s2_samsz_2_N_max_2)
  names(list_OF_2)=c("Group sequential design",
                     "Simulation results: Information fractions 2, Sample size 1, N_max 1",
                     "Simulation results: Information fractions 2, Sample size 2, N_max 1",
                     "Simulation results: Information fractions 2, Sample size 1, N_max 2",
                     "Simulation results: Information fractions 2, Sample size 2, N_max 2")


  list_alsp_1=list(s1_alsp_design,
                   simu_out_alsp_s1_samsz_1_N_max_1,
                   simu_out_alsp_s1_samsz_2_N_max_1,
                   simu_out_alsp_s1_samsz_1_N_max_2,
                   simu_out_alsp_s1_samsz_2_N_max_2)
  names(list_alsp_1)=c("Group sequential design",
                       "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                       "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                       "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                       "Simulation results: Information fractions 1, Sample size 2, N_max 2")

  list_alsp_2=list(s2_alsp_design,
                   simu_out_alsp_s2_samsz_1_N_max_1,
                   simu_out_alsp_s2_samsz_2_N_max_1,
                   simu_out_alsp_s2_samsz_1_N_max_2,
                   simu_out_alsp_s2_samsz_2_N_max_2)
  names(list_alsp_2)=c("Group sequential design",
                       "Simulation results: Information fractions 2, Sample size 1, N_max 1",
                       "Simulation results: Information fractions 2, Sample size 2, N_max 1",
                       "Simulation results: Information fractions 2, Sample size 1, N_max 2",
                       "Simulation results: Information fractions 2, Sample size 2, N_max 2")


  out=list(list_OF_1,list_OF_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}

#################
#masd_simu_OC_normal(mu_t_vec_1,mu_t_vec_2,mu_c,sigma_t,sigma_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num,dist)
#masd_simu_OC_normal(groups=rbind(c(0.3,0.4,0.5),c(0.1,0.2,0.3),c(0.3,0.4,0.5),c(0.1,0.2,0.3)), mu_c=0,sigma_t=1,sigma_c=1,direction=1,
#                             s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                           theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1, sim_num=100)

#########################################################
####  II.2  Simulations for survival analysis      ######
#########################################################


#######################################
msd_rej_rate_survival<-function(HR,lambda,s,c_bry,evt_num,samp_multip,surv_rate,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  comp=length(HR)
  ##########################
  evt_num_s=ceiling(evt_num*s)
  sampsz=ceiling(samp_multip*evt_num/surv_rate)
  mu=-log(HR)
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej=0; stg=0; evt_num_simu_2arm=evt_num_s[length(s)]
    evt_num_simu_tot=(evt_num_s[length(s)]/2)*(comp+1); arm_rej=rep(0,comp)

      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      ### geting random sample from exponential distributions ################
      x0=stats::rexp(sampsz,lambda) #generate exponential data for the control arm
      x=matrix(stats::rexp(sampsz*comp, HR*lambda), nrow=comp,ncol=sampsz)
      T=-log(surv_rate)/lambda  #Censoring time. Events after T will be censored
      time_x0=x0
      time_x=x
      x0_status=c(rep(1,length(x0)))
      x_status=matrix(rep(1,sampsz*comp),nrow=comp,ncol=sampsz)

      ######################## Censored data for the control arm #####################################
      for (i in (1:sampsz)){
        if (x0[i]>T){
          time_x0[i]=T  #x0[i] is censored after T
          x0_status[i]=0   #censoring indicator
        }
      }
      grp_x0=c(rep(1,length(time_x0)))#vector of zeros, length of the vector is the same as the original smaple size

      ######################## Censored data for the test  arms #####################################
      grp_x=matrix(rep(0,sampsz*comp),nrow=comp,ncol=sampsz)
      for (j in (1:comp)){
        for (i in (1:sampsz)){
          if (x[j,i]>T){
            time_x[j,i]=T #x[i] is censored after T
            x_status[j,i]=0  #censoring indicator
          }
        }
      }

      #### count the numnber of events in the control and test arms #######
      evt_count=rep(0,sampsz)
      # for (j in (1:comp)){
      #   for (i in 1:sampsz){
      #     evt_count[i]=sum(x0_status[1:i])+sum(x_status[,1:i]) #number of events after comp*i observations
      #   }
      # }
      total_evt_count = 0
      for (i in 1:sampsz){
          total_evt_count = total_evt_count + x0_status[i]+sum(x_status[,i])
          evt_count[i]=total_evt_count #number of events after comp*i observations
      }
      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      ###############################set up matrices fro parameters########################################################

      theta_hat=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      theta_hat_sd=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      wald=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      score=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      t=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      ##### Obtaningin parameters and store in the matrices ##############

      evt_short=0  ## set when a required number of events is unreachable
      for (i in (1:length(s))) {

        evt_pos=evt_count[evt_count>evt_num_s[i]]#substring of  evt_count that are larger than evt_num
        ## the required number of events may be unreachable within the available
        ## follow-up: clamp the index, and blank the statistics below so the
        ## replicate is reported as non-evaluable instead of using an NA index
        evt_short=isTRUE(evt_short) || is.na(match(evt_pos[1],evt_count))
        evt_cut=if (isTRUE(evt_short)) length(evt_count) else match(evt_pos[1],evt_count) #the first incidence that evt_count is larger than evt_num

        for (j in (1:comp)) {
          #  vector of zeros, length of the vector is the same as the original smaple size

          time=c(time_x0[1:evt_cut],time_x[j,1:evt_cut]) #combined samples before total number of events exceed evt_num
          event=c(x0_status[1:evt_cut],x_status[j,1:evt_cut]) #censoring status before total number of events exceed evt_num
          grp=c(grp_x0[1:evt_cut],grp_x[j,1:evt_cut]) #vector of zero, with length=2*evt_cut


          test= list(time=time,status=event,grp)#survival data
          coxcf=tryCatch(summary(coxph(Surv(time,status)~grp,test)), error=function(e) NULL)  ## a failed fit is reported as NA
          # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
          # test is the same as data=test
          theta_hat[j,i]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[1]
          theta_hat_sd[j,i]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[3]
          t[j,i]=(theta_hat_sd[j,i])^(-2)
          wald[j,i]=theta_hat[j,i]/ theta_hat_sd[j,i]
          score[j,i]= theta_hat[j,i]*t[j,i]
        }
      }

      ## the number of events required at a look was not reached within the
      ## available follow-up: the statistics of this replicate are not
      ## evaluable and it falls back to the non-evaluable defaults above
      if (isTRUE(evt_short)) {
        theta_hat[]=NA_real_; theta_hat_sd[]=NA_real_; t[]=NA_real_
        wald[]=NA_real_; score[]=NA_real_
      }

      wald_max=rep(0,(length(s)))
         for (i in (1:(length(s)))) {
        wald_max[i]=max(wald[,i])
      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################

      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald[i,(1:length(s))]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }

      test_max=wald_max-c_bry

      ####################################

      #any value of test>0 indicates that the wald statistics crossed the critical boundaries

      if (isTRUE(max(test_max)>=0)){ ##null hypothesis rejected, no more random sampling

        rej=1
        test_max_t=test_max[test_max>0]
        stg=match(test_max_t[1],test_max)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        evt_num_simu_2arm=evt_num_s[last_vt]
        evt_num_simu_tot= (evt_num_simu_2arm/2)*(comp+1)
      }else if (isTRUE(max(test_max)<0)) {   #beginnin of step 4
            rej=0
            stg=0
          last_vt=length(s)

          evt_num_simu_2arm=evt_num_s[last_vt]
          evt_num_simu_tot= (evt_num_simu_2arm/2)*(comp+1)
        }
          ########################################################################
          ### Step 4.2. If max(theta>theta_cut), calculate conditional power    ##
          ###          and re-estiate sample size                               ##
          ########################################################################

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }

      out=c(rej,stg,evt_num_simu_2arm,evt_num_simu_tot,arm_rej)
  }) #(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr

  rej = as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  evt_num_simu_2arm=as.numeric(simr[,3])
  evt_num_simu_tot=as.numeric(simr[,4])
  arm_rej=simr[,5:(4+length(HR))]

  rej_rt=mean(rej)

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }

  rej_rt=mean(rej)

  arm_rej_rt=NULL

  for (i in 1:length(HR)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  arm_summary=cbind(HR,arm_rej_rt)
  colnames(arm_summary)=c("HR","Rejection rate by arm")

  mean_evt_2arm=mean(evt_num_simu_2arm)
  mean_evt_tot=mean(evt_num_simu_tot)

  tot_evt=(evt_num/2)*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)


  Design=cbind(c(1:length(s)),s,c_bry,rej_time)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look")

  out_1=c( rej_rt,tyI_er_rt,pos_rej_rt,evt_num,tot_evt,mean_evt_2arm,mean_evt_tot)
  names(out_1)=c("Rejection rate","typr I error","non-type I error rejection",
                 "Planned two arm event number", "Planned total event number",
                 "mean two arm event number","mean total event number")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")

  return(out)

}

#msd_rej_rate_survival(HR,lambda,s,c_bry,evt_num,samp_multip,surv_rate,sim_num)

#msd_rej_rate_survival(HR=rep(1,2),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.931599,2.780060, 2.269910),evt_num=150,samp_multip=5,surv_rate=0.5,sim_num=100000)
#msd_rej_rate_survival(HR=c(1,0.8),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.931599,2.780060, 2.269910),evt_num=150,samp_multip=5,surv_rate=0.5,sim_num=100)

################################################
#' Multiple comparisons Group sequential design: Simulations, Survival analysis, Type I error and power
#'
#' @param HR Vector of hazard ratios of each test arm vs. the control arm; HR < 1 favours the test arm; must be > 0
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param surv_rate Survival rate at the end of follow-up; 0 < surv_rate < 1
#' @param evt_num Planned number of events between one test arm and the control arm by the final analysis; must be > 0
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary and the rejection rate by look}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the hazard ratio and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate, the simulated type I error rate, the non-type I error rejection rate, the planned two-arm and total event numbers and the mean simulated two-arm and total event numbers}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' msd_power_simu_survival(HR = c(0.8,0.8), s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05,
#' surv_rate = 0.5, evt_num = 800, sim_num = 20)
#'
#' @export
msd_power_simu_survival<-function(HR,s,bry_type,alpha_2,surv_rate,evt_num,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  comp=length(HR)
  lambda=1
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))
  }
  samp_multip=5
  sim_out=msd_rej_rate_survival(HR=HR,lambda=1,s=s,c_bry=c_bry,evt_num=evt_num,samp_multip=5,surv_rate=surv_rate,sim_num=sim_num)
  out=sim_out
  #names(out)=c("simulation summary","planned control","planned total")
  return(out)
}

##################################
##################################
#msd_power_simu_survival(HR,s,bry_type,alpha_2,surv_rate,evt_num,sim_num)
#msd_power_simu_survival(HR=c(0.8,0.8),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=800,sim_num=100)
#msd_power_simu_survival(HR=c(1,1),s=c(1/2,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=300,sim_num=10000)
#msd_power_simu_survival(HR=c(1,0.9),s=c(1/2,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=300,sim_num=10000)
#msd_power_simu_survival(HR=c(1,0.8),s=c(1/2,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=300,sim_num=10000)
#msd_power_simu_survival(HR=c(1,0.7),s=c(1/2,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=300,sim_num=1000)

##############################################

######################################################
masd_rej_simu_adapt_com_survival<-function(HR,lambda,s,c_bry,evt_num,evt_max,HR_cut,samp_multip,surv_rate,beta,dist,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(HR)
  mu=-log(HR)
  ##########################
  evt_num_s=c(ceiling(evt_num*s),evt_max)
  evt_num_m1=(comp+1)*evt_num/2
  evt_max_m1=(comp+1)*evt_max/2
  evt_num_m1_s=c(ceiling(evt_num_m1*s),evt_max_m1)

  sampsz=ceiling(samp_multip*evt_max/surv_rate)
  theta_cut=-log(HR_cut)

  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    evt_num_simu_2arm=evt_num_s[length(s)]; evt_num_simu_tot=(evt_num_s[length(s)]/2)*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)
      #the loop specifies what is generated for each simulation
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      ### geting random sample from exponential distributions ################
      x0=stats::rexp(sampsz,lambda) #generate exponential data for the control arm

      x_m1=matrix(nrow=comp,ncol=sampsz)
      for (i in 1:comp){
        x_m1[i,]=stats::rexp(sampsz, HR[i]*lambda)
      }

      T=-log(surv_rate)/lambda  #Censoring time. Events after T will be censored
      time_x0=x0
      time_x_m1=x_m1
      x0_status=c(rep(1,length(x0)))
      x_status_m1=matrix(rep(1,sampsz*comp),nrow=comp,ncol=sampsz)

      ######################## Censored data for the control arm #####################################
      for (i in (1:sampsz)){
        if (x0[i]>T){
          time_x0[i]=T  #x0[i] is censored after T
          x0_status[i]=0   #censoring indicator
        }
      }
      grp_x0=c(rep(1,length(time_x0)))#vector of zeros, length of the vector is the same as the original smaple size

      ######################## Censored data for the test  arms #####################################
      #  grp_x=c(rep(0,length(time_x))) #
      grp_x_m1=matrix(rep(0,sampsz*comp),nrow=comp,ncol=sampsz)
      for (j in (1:comp)){
        for (i in (1:sampsz)){
          if (x_m1[j,i]>T){
            time_x_m1[j,i]=T #x_m1[i] is censored after T
            x_status_m1[j,i]=0  #censoring indicator
          }
        }
      }

      #### count the numnber of events in the control and test arms #######
      evt_count_m1=rep(0,sampsz)
      # for (j in (1:comp)){
      #   for (i in 1:sampsz){
      #     evt_count_m1[i]=sum(x0_status[1:i])+sum(x_status_m1[,1:i]) #number of events after comp*i observations
      #   }
      # }

      total_evt_count_m1 = 0
      for (i in 1:sampsz){
          total_evt_count_m1 = total_evt_count_m1 + x0_status[i]+sum(x_status_m1[,i])
          evt_count_m1[i]=total_evt_count_m1 #number of events after comp*i observations
      }

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      ###############################set up matrices fro parameters########################################################

      theta_hat_m1=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      theta_hat_m1_sd=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      wald_m1=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      score_m1=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      t_m1=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      ##### Obtaningin parameters and store in the matrices ##############

      evt_short=0  ## set when a required number of events is unreachable
      for (i in (1:(length(s)+1))) {

        evt_pos=evt_count_m1[evt_count_m1>evt_num_m1_s[i]]#substring of  evt_count that are larger than evt_num
        ## the required number of events may be unreachable within the available
        ## follow-up: clamp the index, and blank the statistics below so the
        ## replicate is reported as non-evaluable instead of using an NA index
        evt_short=isTRUE(evt_short) || is.na(match(evt_pos[1],evt_count_m1))
        evt_cut_m1=if (isTRUE(evt_short)) length(evt_count_m1) else match(evt_pos[1],evt_count_m1) #the first incidence that evt_count is larger than evt_num

        for (j in (1:comp)) {
          #  vector of zeros, length of the vector is the same as the original smaple size

          time=c(time_x0[1:evt_cut_m1],time_x_m1[j,1:evt_cut_m1]) #combined samples before total number of events exceed evt_num
          event=c(x0_status[1:evt_cut_m1],x_status_m1[j,1:evt_cut_m1]) #censoring status before total number of events exceed evt_num
          grp=c(grp_x0[1:evt_cut_m1],grp_x_m1[j,1:evt_cut_m1]) #vector of zero, with length=2*evt_cut


          test= list(time=time,status=event,grp)#survival data
          coxcf=tryCatch(summary(coxph(Surv(time,status)~grp,test)), error=function(e) NULL)  ## a failed fit is reported as NA
          # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
          # test is the same as data=test
          theta_hat_m1[j,i]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[1]
          theta_hat_m1_sd[j,i]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[3]
          t_m1[j,i]=(theta_hat_m1_sd[j,i])^(-2)
          wald_m1[j,i]=theta_hat_m1[j,i]/ theta_hat_m1_sd[j,i]
          score_m1[j,i]= theta_hat_m1[j,i]*t_m1[j,i]
        }
      }

      ## the number of events required at a look was not reached within the
      ## available follow-up: the statistics of this replicate are not
      ## evaluable and it falls back to the non-evaluable defaults above
      if (isTRUE(evt_short)) {
        theta_hat_m1[]=NA_real_; theta_hat_m1_sd[]=NA_real_; t_m1[]=NA_real_
        wald_m1[]=NA_real_; score_m1[]=NA_real_
      }

      T_hat=t_m1[,length(s)]
      wald_max_m1=NULL
      for (j in 1:length(s) ){
        wald_max_m1[j]=max(wald_m1[,j])
      }
      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]
      arm_drop=rep(0,length(HR))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################

      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,(1:length(s))]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }

      test_max=wald_max_m1-c_bry
      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}

      ####################################

      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }

      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling
        futi=0
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        evt_num_simu_2arm=evt_num_s[last_vt]
        evt_num_simu_tot= (evt_num_simu_2arm/2)*(comp+1)
      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################
        #############################################################################################
        ### Step 4.1. If max(theta<=theta_cut), conditional power will not exceed conditional type I error ##
        ###     Continue to planned sampe size and stop trial                                      ##
        #############################################################################################
        if (futi==1) { # If TRUE, thiere is not satisfactory efficacy, stop trial at planned sample size
          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            stg=0
          }
          last_vt=length(s)
          evt_num_simu_2arm=evt_num_s[last_vt]
          evt_num_simu_tot= (evt_num_simu_2arm/2)*(comp+1)
        }else if (futi==0){ #beginning of step 4.2
          ########################################################################
          ### Step 4.2. If max(theta>theta_cut), calculate conditional power            ##
          ###          and re-estiate sample size                               ##
          ########################################################################
          mu_drop=match(arm_drop[arm_drop>0],arm_drop)
          if (length(mu_drop)==0){
            mu_keep=c(1:comp)
          }else if ( length(mu_drop)>0){
            mu_keep=c(1:comp)[-mu_drop]
          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          T_inter_m2=T_inter
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2

          event_inter=evt_num_s[length(s)-1]

          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=2)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1=theta_hat_inter_m1,theta_hat_inter_m1_sd=theta_hat_inter_m1_sd,dist=2)
          if (isTRUE(c_power_m2>= 1-beta)){
            #########################################
            ### Proceed to planned sample size  #####
            ### and conduct final testing       #####
            #########################################
            last_vt=length(s)
            wald_m2=wald_m1[mu_keep,]
            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last
            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }
            evt_num_simu_2arm=evt_num_s[last_vt]
            evt_num_simu_tot=(evt_num_s[length(s)-1]/2)*(comp+1)+((evt_num_s[last_vt]-evt_num_s[length(s)-1])/2)*(length(mu_keep)+1)
          }else if (isTRUE(c_power_m2<1-beta)){

            t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                           theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

            ## default when the re-estimation is not evaluable
            new_sampsz=evt_num
            if (isTRUE(as.numeric(t_new_0[2])<=1)){
              new_sampsz=evt_num
            }else if (isTRUE(as.numeric(t_new_0[2])>1)){
              T_new_0=as.numeric(t_new_0[[1]])

              event_inter=evt_num_s[length(s)-1]
              new_sampsz_0=ceiling(event_inter*T_new_0/stats::median(T_inter_m2))

              new_sampsz=min(max(new_sampsz_0,evt_num),evt_max)
            }

            if (isTRUE(new_sampsz==evt_num)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              evt_num_simu_2arm=evt_num_s[last_vt]
              evt_num_simu_tot=(evt_num_s[length(s)-1]/2)*(comp+1)+((evt_num_s[last_vt]-evt_num_s[length(s)-1])/2)*(length(mu_keep)+1)

            }else if (isTRUE(new_sampsz>evt_num)){

              theta_m2=theta_hat_inter_m1[mu_keep]
              theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
              theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
              evt_inter= evt_num_m1_s[length(s)-1]
              t_m1_inter=t_m1[,(length(s)-1)]
              t_m2_inter=t_m1_inter[mu_keep]

              inter_vt=length(s)-1
              T_inter_m2=t_m1[mu_keep,inter_vt]

              ##### Generate addtional random samples #######

              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              #### count the numnber of events in the control and test arms #######
              evt_count_m2=rep(0,sampsz)
              x_status_m2=x_status_m1[mu_keep,]
              time_x_m2=time_x_m1[mu_keep,]
              grp_x_m2=grp_x_m1[mu_keep,]

              # if (length(mu_keep)==1){
              #   for (i in 1:sampsz){
              #     evt_count_m2[i]=sum(x0_status[1:i])+sum(x_status_m2[1:i]) #number of events after comp*i observations
              #   }
              # }else if (length(mu_keep)>1) {
              #
              #   for (j in (1:length(mu_keep))){
              #     for (i in 1:sampsz){
              #       evt_count_m2[i]=sum(x0_status[1:i])+sum(x_status_m2[,1:i]) #number of events after comp*i observations
              #     }
              #   }
              # }

              if (isTRUE(length(mu_keep)==1)){
                total_evt_count_m2 = 0
                for (i in 1:sampsz){
                  total_evt_count_m2 = total_evt_count_m2 + x0_status[i]+x_status_m2[i]
                  evt_count_m2[i]=total_evt_count_m2 #number of events after comp*i observations
                }
              } else {
                total_evt_count_m2 = 0
                for (i in 1:sampsz){
                  total_evt_count_m2 = total_evt_count_m2 + x0_status[i]+sum(x_status_m2[,i])
                  evt_count_m2[i]=total_evt_count_m2 #number of events after comp*i observations
                }
              }

              theta_hat_m2_last=rep(0,length(mu_keep))
              theta_hat_m2_last_sd=rep(0,length(mu_keep))
              wald_m2_last=rep(0,length(mu_keep))
              score_m2_last=rep(0,length(mu_keep))
              t_m2_last=rep(0,length(mu_keep))
              ##### Obtaningin parameters and store in the matrices ##############

              new_evt_m2=(new_sampsz/2)*(length(mu_keep)+1)
              evt_pos=evt_count_m2[evt_count_m2>new_sampsz]#substring of  evt_count that are larger than evt_num
              ## the required number of events may be unreachable within the available
              ## follow-up: clamp the index, and blank the statistics below so the
              ## replicate is reported as non-evaluable instead of using an NA index
              evt_short=is.na(match(evt_pos[1],evt_count_m2))
              evt_cut_m2=if (isTRUE(evt_short)) length(evt_count_m2) else match(evt_pos[1],evt_count_m2) #the first incidence that evt_count is larger than evt_num

              if (isTRUE(length(mu_keep)==1)){
                time=c(time_x0[1:evt_cut_m2],time_x_m2[1:evt_cut_m2]) #combined samples before total number of events exceed evt_num
                event=c(x0_status[1:evt_cut_m2],x_status_m2[1:evt_cut_m2]) #censoring status before total number of events exceed evt_num
                grp=c(grp_x0[1:evt_cut_m2],grp_x_m2[1:evt_cut_m2]) #vector of zero, with length=2*evt_cut

                test= list(time=time,status=event,grp)#survival data
                coxcf=tryCatch(summary(coxph(Surv(time,status)~grp,test)), error=function(e) NULL)  ## a failed fit is reported as NA
                # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
                # test is the same as data=test
                theta_hat_m2_last=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[1]
                theta_hat_m2_last_sd=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[3]
                t_m2_last=(theta_hat_m2_last_sd)^(-2)
                wald_m2_last=theta_hat_m2_last/theta_hat_m2_last_sd
                score_m2_last= theta_hat_m2_last* t_m2_last
              }else if (isTRUE(length(mu_keep)>1)){

                for (j in (1:length(mu_keep))) {
                  #  vector of zeros, length of the vector is the same as the original smaple size
                  time=c(time_x0[1:evt_cut_m2],time_x_m2[j,1:evt_cut_m2]) #combined samples before total number of events exceed evt_num
                  event=c(x0_status[1:evt_cut_m2],x_status_m2[j,1:evt_cut_m2]) #censoring status before total number of events exceed evt_num
                  grp=c(grp_x0[1:evt_cut_m2],grp_x_m2[j,1:evt_cut_m2]) #vector of zero, with length=2*evt_cut

                  test= list(time=time,status=event,grp)#survival data
                  coxcf=tryCatch(summary(coxph(Surv(time,status)~grp,test)), error=function(e) NULL)  ## a failed fit is reported as NA
                  # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
                  # test is the same as data=test
                  theta_hat_m2_last[j]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[1]
                  theta_hat_m2_last_sd[j]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[3]
                  t_m2_last[j]=(theta_hat_m2_last_sd[j])^(-2)
                  wald_m2_last[j]=theta_hat_m2_last[j]/theta_hat_m2_last_sd[j]
                  score_m2_last[j]= theta_hat_m2_last[j]* t_m2_last[j]
                }
              }

              ## the number of events required after the sample size increase
              ## was not reached: the statistics after re-estimation are not
              ## evaluable and the replicate falls back to its defaults
              if (isTRUE(evt_short)) {
                theta_hat_m2_last[]=NA_real_; theta_hat_m2_last_sd[]=NA_real_
                t_m2_last[]=NA_real_; wald_m2_last[]=NA_real_; score_m2_last[]=NA_real_
              }

              T_new=t_m2_last

              T_inter= (event_inter/new_sampsz)*T_new

              T_hat_m2=T_hat[mu_keep]

              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=2)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=2)
              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=2)

              c_bry_new=NULL
              for (i in 1:length(theta_hat_inter_m2)){
                c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
              }

              test_new=wald_m2_last-max(c_bry_new)

              #####################################################################
              ## Reject the Null hypothesis if the critical boundary is crossed  ###
              #####################################################################

              evt_num_simu_2arm=new_sampsz
              # evt_num_simu_tot= (evt_num_simu_2arm/2)*(comp+1)
              evt_num_simu_tot=(evt_num_s[length(s)-1]/2)*(comp+1)+((new_sampsz-evt_num_s[length(s)-1])/2)*(length(mu_keep)+1)
              #  if ( max(test_max_new)>0){
              if (isTRUE(max(test_new)>=0)){
                rej=1
                stg=length(s)
              }else if (isTRUE(max(test_new)<0)){
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0
              }

              arm_rej_new=NULL
              for (i in 1:length(mu_keep)){
                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }
              }
              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new
              rej_all[,length(s)]=arm_rej_last
            }
          }
        }
      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }

      out=c(rej_nossr,rej,stg,evt_num_simu_2arm,evt_num_simu_tot,futi,arm_drop,arm_rej)
      #  out
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  ##############################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  evt_num_simu_2arm=as.numeric(simr[,4])
  evt_num_simu_tot=as.numeric(simr[,5])
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(HR))]
  arm_rej=simr[,(7+length(HR)):(6+2*length(HR))]
  ########################

  rej_rt=mean(rej)
  rej_nossr_rt=mean(rej_nossr)

  futi_rt=mean(futi)
  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }

  rej_rt=mean(rej)

  arm_rej_rt=NULL

  for (i in 1:length(HR)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  arm_drop_rt = NULL

  for (i in 1:length(HR)){
    arm_drop_rt[i]=mean(arm_drop[,i])
  }

  futi_rt=mean(futi)

  arm_summary=cbind(HR,arm_drop_rt,arm_rej_rt)
  colnames(arm_summary)=c("HR","drop rate by arm","Rejection rate by arm")
  mean_evt_2arm=mean(evt_num_simu_2arm)
  mean_evt_tot=mean(evt_num_simu_tot)

  tot_evt=(evt_num/2)*(comp+1)
  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){
        arm_pos_rej[j,i]=1
      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL
  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL
  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (evt_num_simu_2arm[i]>evt_num){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  evt_num_simu_2arm_samsz_icr= evt_num_simu_2arm[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_evt_num_simu_2arm_samsz_icr=mean(evt_num_simu_2arm_samsz_icr)
  ################
  Design=cbind(c(1:length(s)),s,c_bry,rej_time,samsz_icr_byLook)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look","Rate of sample size increase")

  out_1=c( rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,evt_num,evt_max,mean_evt_2arm,mean_evt_num_simu_2arm_samsz_icr)
  names(out_1)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "typr I error","non-type I error rejection","futility rate","Planned two arm event number",
                 "Maximum two arm event nmuber",
                 "mean two arm event number","mean two arm event number conditional on increase")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")

  return(out)

}

#######################################

##################################################
#masd_rej_simu_adapt_com_survival(HR,lambda,s,c_bry,evt_num,evt_max,HR_cut,samp_multip,surv_rate,beta,dist,sim_num)
#masd_rej_simu_adapt_com_survival(HR=rep(1,2),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.931599,2.780060, 2.269910),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=1000)
#masd_rej_simu_adapt_com_survival(HR=c(0.8,0.8),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.931599,2.780060, 2.269910),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=100)
#masd_rej_simu_adapt_com_survival(HR=c(0.8,0.8),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.931599,2.780060, 2.269910),evt_num=600,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=1000)
#masd_rej_simu_adapt_com_survival(HR=c(1,1,1),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),evt_num=200,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=1000)


########################################################

################################################
#' Multiple comparisons Group sequential design: Simulations, Survival analysis, Operating characteristics
#'
#' @param HR Vector of hazard ratios of each test arm vs. the control arm; HR < 1 favours the test arm; must be > 0
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param surv_rate Survival rate at the end of follow-up; 0 < surv_rate < 1
#' @param HR_cut Hazard-ratio cut-off for the lack-of-efficacy rule at the interim analysis:
#'   a test arm is dropped when its observed hazard ratio is greater than HR_cut
#' @param evt_num Planned number of events between one test arm and the control arm by the final analysis; must be > 0
#' @param evt_max Maximum number of events between one test arm and the control arm allowed by the sample size re-estimation; must be >= evt_num
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary, the rejection rate by look and the rate of sample size increase}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the hazard ratio, the drop rate by arm and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate without and with sample size re-estimation, the rate of sample size increase, the rejection rates conditional on meeting the sample size increase criteria, the simulated type I error rate, the non-type I error rejection rate, the futility rate, the planned and maximum two-arm event numbers and the mean simulated two-arm event numbers with and without conditioning on an increase}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' masd_power_simu_adapt_survival(HR = c(0.63, 0.65, 0.7), s = c(1/3,2/3,1), bry_type = 1,
#' alpha_2 = 0.05, surv_rate = 0.5, HR_cut = 0.95, evt_num = 200, evt_max = 500, beta = 0.1,
#' sim_num = 2)
#' }
#'
#' @export
masd_power_simu_adapt_survival<-function(HR,s,bry_type,alpha_2,surv_rate,HR_cut,evt_num,evt_max,beta,sim_num){

  alpha=alpha_2/2
  comp=length(HR)
  lambda=1
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))
  }
  samp_multip=5
  simu_out=masd_rej_simu_adapt_com_survival(HR=HR,lambda=1,s=s,c_bry=c_bry,evt_num=evt_num,evt_max=evt_max,HR_cut=HR_cut,samp_multip=samp_multip,surv_rate=surv_rate,beta=beta,dist=1,sim_num=sim_num)
  out=simu_out

  return(out)
}

##################################

##################################
#masd_power_simu_adapt_survival(HR,s,bry_type,alpha_2,surv_rate,HR_cut,evt_num,evt_max,beta,sim_num)
#masd_power_simu_adapt_survival(HR=c(0.8,0.8),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=800,evt_max=1200,beta=0.1,sim_num=100)
#masd_power_simu_adapt_survival(HR=c(0.8,0.8),s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),surv_rate=0.5,HR_cut=0.95,evt_num=50,evt_max=1200,beta=0.1,sim_num=100)
#masd_power_simu_adapt_survival(HR=c(1,1,1),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=800,evt_max=1200,beta=0.1,sim_num=10000)
#masd_power_simu_adapt_survival(HR=c(0.63,0.65,0.7),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=200,evt_max=500,beta=0.1,sim_num=1000)
#masd_power_simu_adapt_survival(HR=c(0.65,0.7,0.75),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=200,evt_max=500,beta=0.1,sim_num=1000)
#masd_power_simu_adapt_survival(HR=c(0.7,0.75,0.8),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=200,evt_max=500,beta=0.1,sim_num=1000)
#masd_power_simu_adapt_survival(HR=c(0.75,0.8,0.85),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=200,evt_max=500,beta=0.1,sim_num=1000)
#masd_power_simu_adapt_survival(HR=c(1,0.6,0.7),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=300,evt_max=800,beta=0.1,sim_num=1000)


######################################################

###################################################################################
## Simulations for Optimizing Operating characteristics - Survival analysis     ###
###################################################################################
masd_rej_OC_simu_adapt_com_survival<-function(HR,lambda,s,c_bry,evt_num,evt_max,HR_cut,samp_multip,surv_rate,beta,dist,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(HR)
  mu=-log(HR)
  ##########################
  evt_num_s=c(ceiling(evt_num*s),evt_max)
  evt_num_m1=(comp+1)*evt_num/2
  evt_max_m1=(comp+1)*evt_max/2
  evt_num_m1_s=c(ceiling(evt_num_m1*s),evt_max_m1)

  sampsz=ceiling(samp_multip*evt_max/surv_rate)
  theta_cut=-log(HR_cut)

  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    evt_num_simu_2arm=evt_num_s[length(s)]; evt_num_simu_tot=(evt_num_s[length(s)]/2)*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)
      #the loop specifies what is generated for each simulation

      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      ### geting random sample from exponential distributions ################
      x0=stats::rexp(sampsz,lambda) #generate exponential data for the control arm

      x_m1=matrix(nrow=comp,ncol=sampsz)
      for (i in 1:comp){
        x_m1[i,]=stats::rexp(sampsz, HR[i]*lambda)
      }

      T=-log(surv_rate)/lambda  #Censoring time. Events after T will be censored
      time_x0=x0
      time_x_m1=x_m1
      x0_status=c(rep(1,length(x0)))
      x_status_m1=matrix(rep(1,sampsz*comp),nrow=comp,ncol=sampsz)

      ######################## Censored data for the control arm #####################################
      for (i in (1:sampsz)){
        if (x0[i]>T){
          time_x0[i]=T  #x0[i] is censored after T
          x0_status[i]=0   #censoring indicator
        }
      }
      grp_x0=c(rep(1,length(time_x0)))#vector of zeros, length of the vector is the same as the original smaple size

      ######################## Censored data for the test  arms #####################################
      #  grp_x=c(rep(0,length(time_x))) #
      grp_x_m1=matrix(rep(0,sampsz*comp),nrow=comp,ncol=sampsz)
      for (j in (1:comp)){
        for (i in (1:sampsz)){
          if (x_m1[j,i]>T){
            time_x_m1[j,i]=T #x_m1[i] is censored after T
            x_status_m1[j,i]=0  #censoring indicator
          }
        }
      }

      #### count the numnber of events in the control and test arms #######
      evt_count_m1=rep(0,sampsz)

      # for (j in (1:comp)){
      #   for (i in 1:sampsz){
      #     evt_count_m1[i]=sum(x0_status[1:i])+sum(x_status_m1[,1:i]) #number of events after comp*i observations
      #   }
      # }

      total_evt_count_m1 = 0
      for (i in 1:sampsz){
          total_evt_count_m1 = total_evt_count_m1 + x0_status[i] + sum(x_status_m1[,i])
          evt_count_m1[i] = total_evt_count_m1 #number of events after comp*i observations
      }

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      ###############################set up matrices fro parameters########################################################

      theta_hat_m1=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      theta_hat_m1_sd=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      wald_m1=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      score_m1=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      t_m1=matrix(rep(0,(length(s)+1)*comp),nrow=comp,ncol=(length(s)+1))
      ##### Obtaningin parameters and store in the matrices ##############

      evt_short=0  ## set when a required number of events is unreachable
      for (i in (1:(length(s)+1))) {

        evt_pos=evt_count_m1[evt_count_m1>evt_num_m1_s[i]]#substring of  evt_count that are larger than evt_num
        ## the required number of events may be unreachable within the available
        ## follow-up: clamp the index, and blank the statistics below so the
        ## replicate is reported as non-evaluable instead of using an NA index
        evt_short=isTRUE(evt_short) || is.na(match(evt_pos[1],evt_count_m1))
        evt_cut_m1=if (isTRUE(evt_short)) length(evt_count_m1) else match(evt_pos[1],evt_count_m1) #the first incidence that evt_count is larger than evt_num

        for (j in (1:comp)) {
          #  vector of zeros, length of the vector is the same as the original smaple size

          time=c(time_x0[1:evt_cut_m1],time_x_m1[j,1:evt_cut_m1]) #combined samples before total number of events exceed evt_num
          event=c(x0_status[1:evt_cut_m1],x_status_m1[j,1:evt_cut_m1]) #censoring status before total number of events exceed evt_num
          grp=c(grp_x0[1:evt_cut_m1],grp_x_m1[j,1:evt_cut_m1]) #vector of zero, with length=2*evt_cut


          test= list(time=time,status=event,grp)#survival data
          coxcf=tryCatch(summary(coxph(Surv(time,status)~grp,test)), error=function(e) NULL)  ## a failed fit is reported as NA
          # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
          # test is the same as data=test
          theta_hat_m1[j,i]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[1]
          theta_hat_m1_sd[j,i]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[3]
          t_m1[j,i]=(theta_hat_m1_sd[j,i])^(-2)
          wald_m1[j,i]=theta_hat_m1[j,i]/ theta_hat_m1_sd[j,i]
          score_m1[j,i]= theta_hat_m1[j,i]*t_m1[j,i]
        }
      }

      ## the number of events required at a look was not reached within the
      ## available follow-up: the statistics of this replicate are not
      ## evaluable and it falls back to the non-evaluable defaults above
      if (isTRUE(evt_short)) {
        theta_hat_m1[]=NA_real_; theta_hat_m1_sd[]=NA_real_; t_m1[]=NA_real_
        wald_m1[]=NA_real_; score_m1[]=NA_real_
      }

      T_hat=t_m1[,length(s)]
      wald_max_m1=NULL

      for (j in 1:length(s) ){
        wald_max_m1[j]=max(wald_m1[,j])
      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop=rep(0,length(HR))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################

      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,(1:length(s))]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }

      test_max=wald_max_m1-c_bry
      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}

      ####################################

      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }

      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling

        futi=0
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        evt_num_simu_2arm=evt_num_s[last_vt]
        evt_num_simu_tot= (evt_num_simu_2arm/2)*(comp+1)

      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        #############################################################################################
        ### Step 4.1. If max(theta<=theta_cut), conditional power will not exceed conditional type I error ##
        ###     Continue to planned sampe size and stop trial                                      ##
        #############################################################################################

        if (futi==1) { # If TRUE, thiere is not satisfactory efficacy, stop trial at planned sample size

          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            stg=0}
          last_vt=length(s)

          evt_num_simu_2arm=evt_num_s[last_vt]
          evt_num_simu_tot= (evt_num_simu_2arm/2)*(comp+1)

        }else if (futi==0){ #beginning of step 4.2

          ########################################################################
          ### Step 4.2. If max(theta>theta_cut), calculate conditional power            ##
          ###          and re-estiate sample size                               ##
          ########################################################################

          mu_drop=match(arm_drop[arm_drop>0],arm_drop)

          if (length(mu_drop)==0){

            mu_keep=c(1:comp)

          }else if ( length(mu_drop)>0){

            mu_keep=c(1:comp)[-mu_drop]

          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          T_inter_m2=T_inter
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2

          event_inter=evt_num_s[length(s)-1]

          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=2)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1=theta_hat_inter_m1,theta_hat_inter_m1_sd=theta_hat_inter_m1_sd,dist=2)
          if (isTRUE(c_power_m2>= 1-beta)){

            #########################################
            ### Proceed to planned sample size  #####
            ### and conduct final testing       #####
            #########################################

            last_vt=length(s)

            wald_m2=wald_m1[mu_keep,]

            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last

            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }

            evt_num_simu_2arm=evt_num_s[last_vt]
            evt_num_simu_tot=(evt_num_s[length(s)-1]/2)*(comp+1)+((evt_num_s[last_vt]-evt_num_s[length(s)-1])/2)*(length(mu_keep)+1)

          }else if (isTRUE(c_power_m2<1-beta)){

            t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                           theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

            ## default when the re-estimation is not evaluable
            new_sampsz=evt_num
            if (isTRUE(as.numeric(t_new_0[2])<=1)){
              new_sampsz=evt_num
            }else if (isTRUE(as.numeric(t_new_0[2])>1)){
              T_new_0=as.numeric(t_new_0[[1]])

              event_inter=evt_num_s[length(s)-1]
              new_sampsz_0=ceiling(event_inter*T_new_0/stats::median(T_inter_m2))

              new_sampsz=min(max(new_sampsz_0,evt_num),evt_max)
            }

            if (isTRUE(new_sampsz==evt_num)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              evt_num_simu_2arm=evt_num_s[last_vt]
              evt_num_simu_tot=(evt_num_s[length(s)-1]/2)*(comp+1)+((evt_num_s[last_vt]-evt_num_s[length(s)-1])/2)*(length(mu_keep)+1)

            }else if (isTRUE(new_sampsz>evt_num)){

              theta_m2=theta_hat_inter_m1[mu_keep]
              theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
              theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
              evt_inter= evt_num_m1_s[length(s)-1]
              t_m1_inter=t_m1[,(length(s)-1)]
              t_m2_inter=t_m1_inter[mu_keep]

              inter_vt=length(s)-1
              T_inter_m2=t_m1[mu_keep,inter_vt]


              ##### Generate addtional random samples #######

              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              #### count the numnber of events in the control and test arms #######
              evt_count_m2=rep(0,sampsz)
              x_status_m2=x_status_m1[mu_keep,]
              time_x_m2=time_x_m1[mu_keep,]
              grp_x_m2=grp_x_m1[mu_keep,]

              # if (length(mu_keep)==1){
              #
              #   for (i in 1:sampsz){
              #     evt_count_m2[i]=sum(x0_status[1:i])+sum(x_status_m2[1:i]) #number of events after comp*i observations
              #   }
              # }else if (length(mu_keep)>1) {
              #
              #   for (j in (1:length(mu_keep))){
              #     for (i in 1:sampsz){
              #       evt_count_m2[i]=sum(x0_status[1:i])+sum(x_status_m2[,1:i]) #number of events after comp*i observations
              #     }
              #   }
              # }
              if (isTRUE(length(mu_keep)==1)){
                total_evt_count_m2 = 0
                for (i in 1:sampsz){
                  total_evt_count_m2 = total_evt_count_m2 + x0_status[i]+x_status_m2[i]
                  evt_count_m2[i]=total_evt_count_m2 #number of events after comp*i observations
                }
              } else {
                total_evt_count_m2 = 0
                for (i in 1:sampsz){
                  total_evt_count_m2 = total_evt_count_m2 + x0_status[i]+sum(x_status_m2[,i])
                  evt_count_m2[i]=total_evt_count_m2 #number of events after comp*i observations
                }
              }


              theta_hat_m2_last=rep(0,length(mu_keep))
              theta_hat_m2_last_sd=rep(0,length(mu_keep))
              wald_m2_last=rep(0,length(mu_keep))
              score_m2_last=rep(0,length(mu_keep))
              t_m2_last=rep(0,length(mu_keep))
              ##### Obtaningin parameters and store in the matrices ##############

              new_evt_m2=(new_sampsz/2)*(length(mu_keep)+1)
              evt_pos=evt_count_m2[evt_count_m2>new_sampsz]#substring of  evt_count that are larger than evt_num
              ## the required number of events may be unreachable within the available
              ## follow-up: clamp the index, and blank the statistics below so the
              ## replicate is reported as non-evaluable instead of using an NA index
              evt_short=is.na(match(evt_pos[1],evt_count_m2))
              evt_cut_m2=if (isTRUE(evt_short)) length(evt_count_m2) else match(evt_pos[1],evt_count_m2) #the first incidence that evt_count is larger than evt_num

              if (isTRUE(length(mu_keep)==1)){
                time=c(time_x0[1:evt_cut_m2],time_x_m2[1:evt_cut_m2]) #combined samples before total number of events exceed evt_num
                event=c(x0_status[1:evt_cut_m2],x_status_m2[1:evt_cut_m2]) #censoring status before total number of events exceed evt_num
                grp=c(grp_x0[1:evt_cut_m2],grp_x_m2[1:evt_cut_m2]) #vector of zero, with length=2*evt_cut

                test= list(time=time,status=event,grp)#survival data
                coxcf=tryCatch(summary(coxph(Surv(time,status)~grp,test)), error=function(e) NULL)  ## a failed fit is reported as NA
                # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
                # test is the same as data=test
                theta_hat_m2_last=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[1]
                theta_hat_m2_last_sd=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[3]
                t_m2_last=(theta_hat_m2_last_sd)^(-2)
                wald_m2_last=theta_hat_m2_last/theta_hat_m2_last_sd
                score_m2_last= theta_hat_m2_last* t_m2_last
              }else if (isTRUE(length(mu_keep)>1)){

                for (j in (1:length(mu_keep))) {
                  #  vector of zeros, length of the vector is the same as the original smaple size

                  time=c(time_x0[1:evt_cut_m2],time_x_m2[j,1:evt_cut_m2]) #combined samples before total number of events exceed evt_num
                  event=c(x0_status[1:evt_cut_m2],x_status_m2[j,1:evt_cut_m2]) #censoring status before total number of events exceed evt_num
                  grp=c(grp_x0[1:evt_cut_m2],grp_x_m2[j,1:evt_cut_m2]) #vector of zero, with length=2*evt_cut


                  test= list(time=time,status=event,grp)#survival data
                  coxcf=tryCatch(summary(coxph(Surv(time,status)~grp,test)), error=function(e) NULL)  ## a failed fit is reported as NA
                  # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
                  # test is the same as data=test
                  theta_hat_m2_last[j]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[1]
                  theta_hat_m2_last_sd[j]=if (is.null(coxcf)) NA_real_ else coxcf$coefficients[3]
                  t_m2_last[j]=(theta_hat_m2_last_sd[j])^(-2)
                  wald_m2_last[j]=theta_hat_m2_last[j]/theta_hat_m2_last_sd[j]
                  score_m2_last[j]= theta_hat_m2_last[j]* t_m2_last[j]
                }
              }

              ## the number of events required after the sample size increase
              ## was not reached: the statistics after re-estimation are not
              ## evaluable and the replicate falls back to its defaults
              if (isTRUE(evt_short)) {
                theta_hat_m2_last[]=NA_real_; theta_hat_m2_last_sd[]=NA_real_
                t_m2_last[]=NA_real_; wald_m2_last[]=NA_real_; score_m2_last[]=NA_real_
              }

              T_new=t_m2_last

              T_inter= (event_inter/new_sampsz)*T_new

              T_hat_m2=T_hat[mu_keep]

              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=2)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=2)
              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=2)

              c_bry_new=NULL

              for (i in 1:length(theta_hat_inter_m2)){
                c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
              }

              test_new=wald_m2_last-max(c_bry_new)

              #####################################################################
              ## Reject the Null hypothesis if the critical boundary is crossed  ###
              #####################################################################

              evt_num_simu_2arm=new_sampsz
              # evt_num_simu_tot= (evt_num_simu_2arm/2)*(comp+1)
              evt_num_simu_tot=(evt_num_s[length(s)-1]/2)*(comp+1)+((new_sampsz-evt_num_s[length(s)-1])/2)*(length(mu_keep)+1)
              #  if ( max(test_max_new)>0){
              if (isTRUE(max(test_new)>=0)){
                rej=1
                stg=length(s)

              }else if (isTRUE(max(test_new)<0)){
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0

              }

              arm_rej_new=NULL

              for (i in 1:length(mu_keep)){

                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }

              }

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }
        }

      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }

      out=c(rej_nossr,rej,stg,evt_num_simu_2arm,evt_num_simu_tot,futi,arm_drop,arm_rej)
      #  out
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  ##############################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  evt_num_simu_2arm=as.numeric(simr[,4])
  evt_num_simu_tot=as.numeric(simr[,5])
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(HR))]
  arm_rej=simr[,(7+length(HR)):(6+2*length(HR))]
  ########################

  rej_rt=mean(rej)
  rej_nossr_rt=mean(rej_nossr)

  futi_rt=mean(futi)

  rej_rt=mean(rej)

  arm_rej_rt=NULL
  for (i in 1:length(HR)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  arm_drop_rt = NULL
  for (i in 1:length(HR)){
    arm_drop_rt[i]=mean(arm_drop[,i])
  }

  futi_rt=mean(futi)

  arm_summary=cbind(HR,arm_drop_rt,arm_rej_rt)
  colnames(arm_summary)=c("HR","drop rate by arm","Rejection rate by arm")

  mean_evt_2arm=mean(evt_num_simu_2arm)
  mean_evt_tot=mean(evt_num_simu_tot)

  tot_evt=(evt_num/2)*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (evt_num_simu_2arm[i]>evt_num){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  evt_num_simu_2arm_samsz_icr= evt_num_simu_2arm[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_evt_num_simu_2arm_samsz_icr=mean(evt_num_simu_2arm_samsz_icr)
  ################

  out=c( rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,evt_num,evt_max,mean_evt_2arm,mean_evt_num_simu_2arm_samsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "typr I error","non-type I error rejection","futility rate","Planned two arm event number",
                 "Maximum two arm event nmuber",
                 "mean two arm event number","mean two arm event number conditional on increase")

  return(out)

}

#masd_rej_OC_simu_adapt_com_survival(HR,lambda,s,c_bry,evt_num,evt_max,HR_cut,samp_multip,surv_rate,beta,dist,sim_num)
#masd_rej_OC_simu_adapt_com_survival(HR=rep(1,3),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=1000)
#
#masd_rej_OC_simu_adapt_com_survival(HR=c(0.8,0.8),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.931599,2.780060, 2.269910),evt_num=600,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=100)
#masd_rej_OC_simu_adapt_com_survival(HR=c(0.8,0.8,0.8),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),evt_num=200,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=100)
#masd_rej_OC_simu_adapt_com_survival(HR=c(0.85,0.85,0.85),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),evt_num=200,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=1000)
#masd_rej_OC_simu_adapt_com_survival(HR=c(1,1,1),lambda=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),evt_num=20,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,dist=2,sim_num=10000)
################################################

################
masd_OC_power_simu_adapt_survival<-function(HR,s,c_bry,alpha_2,surv_rate,HR_cut,evt_num,evt_max,beta,sim_num){

  lambda=1
  samp_multip=5
  simu_out=masd_rej_OC_simu_adapt_com_survival(HR=HR,lambda=lambda,s=s,c_bry=c_bry,evt_num=evt_num,evt_max=evt_max,HR_cut=HR_cut,samp_multip=samp_multip,surv_rate=surv_rate,beta=beta,dist=2,sim_num=sim_num)

  HR_names=paste(paste("HR","[",c(1:length(HR)),"]",sep=""))
  out=round(c(HR,as.numeric(simu_out)),5)
  names(out)=c(HR_names, names(simu_out))

  return(out)
}

##################################
##################################
#masd_OC_power_simu_adapt_survival(HR,s,bry_type,alpha_2,surv_rate,HR_cut,evt_num,evt_max,beta,sim_num)
#masd_OC_power_simu_adapt_survival(HR=c(0.8,0.8),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=800,evt_max=1200,beta=0.1,sim_num=100)
#masd_OC_power_simu_adapt_survival(HR=c(0.8,0.8),s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),surv_rate=0.5,HR_cut=0.95,evt_num=50,evt_max=1200,beta=0.1,sim_num=100)
#masd_OC_power_simu_adapt_survival(HR=c(1,1,1),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=800,evt_max=1200,beta=0.1,sim_num=1000)
#masd_OC_power_simu_adapt_survival(HR=c(0.85,0.85,0.85),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=200,evt_max=1200,beta=0.1,sim_num=1000)
#masd_OC_power_simu_adapt_survival(HR=c(0.75,0.85,1),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,HR_cut=0.95,evt_num=400,evt_max=1200,beta=0.1,sim_num=1000)

######################################################
##################################
#' Multiple comparisons Group sequential design: Simulations, Survival analysis, Optimizing adaptive sequential design
#'
#' @param groups Matrix of true hazard ratios of each test arm vs. the control arm: each row is one simulation scenario and each column one test arm; must be > 0
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g. c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Two-sided type I error levels spent at each look of the schedule s1; must sum to alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error levels spent at each look of the schedule s2; must sum to alpha_2
#' @param surv_rate Survival rate at the end of follow-up; 0 < surv_rate < 1
#' @param evt_num_1 Planned number of events between one test arm and the control arm by the final analysis in the first event number scenario; must be > 0
#' @param evt_num_2 Planned number of events between one test arm and the control arm by the final analysis in the second event number scenario; must be > 0
#' @param evt_max_1 Maximum number of events between one test arm and the control arm in the first maximum-event-number scenario; must be >= evt_num_1
#' @param evt_max_2 Maximum number of events between one test arm and the control arm in the second maximum-event-number scenario; must be >= evt_num_2
#' @param HR_cut Hazard-ratio cut-off for the lack-of-efficacy rule at the interim analysis:
#'   a test arm is dropped when its observed hazard ratio is greater than HR_cut
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the O'Brien-Fleming group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned event numbers and the two maximum event numbers}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the O'Brien-Fleming group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned event numbers and the two maximum event numbers}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the alpha-spending group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned event numbers and the two maximum event numbers}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the alpha-spending group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned event numbers and the two maximum event numbers}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' masd_simu_OC_survival(
#' groups = rbind(c(0.85, 0.87, 0.9), c(0.8, 0.82, 0.85), c(0.75, 0.78, 0.8), c(0.7, 0.75, 0.8)),
#' s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05, alpha_2_alps_s1 = c(0.001,0.001,0.048),
#' alpha_2_alps_s2 = c(0.001,0.049), surv_rate = 0.5, evt_num_1 = 200, evt_num_2 = 300,
#' evt_max_1 = 2000, evt_max_2 = 3000, HR_cut = 0.95, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
masd_simu_OC_survival<-function(groups,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,surv_rate,evt_num_1,evt_num_2,evt_max_1,evt_max_2,HR_cut,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)
  comp=length(groups[1,])

  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=2)[[1]][1]))
  s2_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=2)[[1]][1]))

  s1_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry(s=s1,alpha=alpha_alsp_s1,comp=comp,dist=2)[[1]][1]))
  s2_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry(s=s2,alpha=alpha_alsp_s2,comp=comp,dist=2)[[1]][1]))

  ##### Two-stage GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")
  #######################
  HR_vec=groups
  rownames( HR_vec)=NULL

  HR_group_names=paste(paste("HR group",c(1:dim(HR_vec)[1]),sep=" "))

  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####

  simu_out_OF_s1_EVTnum_1_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s1_EVTnum_1_EVTmax_1=rbind(simu_out_OF_s1_EVTnum_1_EVTmax_1,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s1,c_bry=s1_OF_cbry,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                 evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_EVTnum_1_EVTmax_1)=HR_group_names

  ###############  s2 #######
  simu_out_OF_s2_EVTnum_1_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s2_EVTnum_1_EVTmax_1=rbind(simu_out_OF_s2_EVTnum_1_EVTmax_1,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s2,c_bry=s2_OF_cbry,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                 evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_EVTnum_1_EVTmax_1)=HR_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_OF_s1_EVTnum_2_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s1_EVTnum_2_EVTmax_1=rbind(simu_out_OF_s1_EVTnum_2_EVTmax_1,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s1,c_bry=s1_OF_cbry,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                 evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_EVTnum_2_EVTmax_1)=HR_group_names

  ###############  s2 #######
  simu_out_OF_s2_EVTnum_2_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s2_EVTnum_2_EVTmax_1=rbind(simu_out_OF_s2_EVTnum_2_EVTmax_1,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s2,c_bry=s2_OF_cbry,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                 evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_EVTnum_2_EVTmax_1)=HR_group_names



  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####

  simu_out_OF_s1_EVTnum_1_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s1_EVTnum_1_EVTmax_2=rbind(simu_out_OF_s1_EVTnum_1_EVTmax_2,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s1,c_bry=s1_OF_cbry,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                 evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_EVTnum_1_EVTmax_2)=HR_group_names

  ###############  s2 #######
  simu_out_OF_s2_EVTnum_1_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s2_EVTnum_1_EVTmax_2=rbind(simu_out_OF_s2_EVTnum_1_EVTmax_2,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s2,c_bry=s2_OF_cbry,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                 evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_EVTnum_1_EVTmax_2)=HR_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_OF_s1_EVTnum_2_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s1_EVTnum_2_EVTmax_2=rbind(simu_out_OF_s1_EVTnum_2_EVTmax_2,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s1,c_bry=s1_OF_cbry,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                 evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_EVTnum_2_EVTmax_2)=HR_group_names

  ###############  s2 #######
  simu_out_OF_s2_EVTnum_2_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s2_EVTnum_2_EVTmax_2=rbind(simu_out_OF_s2_EVTnum_2_EVTmax_2,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s2,c_bry=s2_OF_cbry,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                 evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_EVTnum_2_EVTmax_2)=HR_group_names


  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################


  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####

  simu_out_alsp_s1_EVTnum_1_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s1_EVTnum_1_EVTmax_1=rbind(simu_out_alsp_s1_EVTnum_1_EVTmax_1,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s1,c_bry=s1_alsp_cbry,alpha_2=alpha_2_alps_s1,surv_rate=surv_rate,
                                                                                                                     evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_EVTnum_1_EVTmax_1)=HR_group_names

  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_1_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s2_EVTnum_1_EVTmax_1=rbind(simu_out_alsp_s2_EVTnum_1_EVTmax_1,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s2,c_bry=s2_alsp_cbry,alpha_2=alpha_2_alps_s2,surv_rate=surv_rate,
                                                                                                                     evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_EVTnum_1_EVTmax_1)=HR_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_alsp_s1_EVTnum_2_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s1_EVTnum_2_EVTmax_1=rbind(simu_out_alsp_s1_EVTnum_2_EVTmax_1,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s1,c_bry=s1_alsp_cbry,alpha_2=alpha_2_alps_s1,surv_rate=surv_rate,
                                                                                                                     evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_EVTnum_2_EVTmax_1)=HR_group_names

  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_2_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s2_EVTnum_2_EVTmax_1=rbind(simu_out_alsp_s2_EVTnum_2_EVTmax_1,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s2,c_bry=s2_alsp_cbry,alpha_2=alpha_2_alps_s2,surv_rate=surv_rate,
                                                                                                                     evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_EVTnum_2_EVTmax_1)=HR_group_names



  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####

  simu_out_alsp_s1_EVTnum_1_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s1_EVTnum_1_EVTmax_2=rbind(simu_out_alsp_s1_EVTnum_1_EVTmax_2,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s1,c_bry=s1_alsp_cbry,alpha_2=alpha_2_alps_s1,surv_rate=surv_rate,
                                                                                                                     evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_EVTnum_1_EVTmax_2)=HR_group_names

  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_1_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s2_EVTnum_1_EVTmax_2=rbind(simu_out_alsp_s2_EVTnum_1_EVTmax_2,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s2,c_bry=s2_alsp_cbry,alpha_2=alpha_2_alps_s2,surv_rate=surv_rate,
                                                                                                                     evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_EVTnum_1_EVTmax_2)=HR_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_alsp_s1_EVTnum_2_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s1_EVTnum_2_EVTmax_2=rbind(simu_out_alsp_s1_EVTnum_2_EVTmax_2,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s1,c_bry=s1_alsp_cbry,alpha_2=alpha_2_alps_s1,surv_rate=surv_rate,
                                                                                                                     evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_EVTnum_2_EVTmax_2)=HR_group_names

  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_2_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s2_EVTnum_2_EVTmax_2=rbind(simu_out_alsp_s2_EVTnum_2_EVTmax_2,masd_OC_power_simu_adapt_survival(HR=HR_vec[i,],s=s2,c_bry=s2_alsp_cbry,alpha_2=alpha_2_alps_s2,surv_rate=surv_rate,
                                                                                                                     evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_EVTnum_2_EVTmax_2)=HR_group_names

  ################################
  ##############

  list_OF_1=list(s1_OF_design,
                 simu_out_OF_s1_EVTnum_1_EVTmax_1,
                 simu_out_OF_s1_EVTnum_2_EVTmax_1,
                 simu_out_OF_s1_EVTnum_1_EVTmax_2,
                 simu_out_OF_s1_EVTnum_2_EVTmax_2)

  names(list_OF_1)=c("Group sequential design",
                     "Simulation results: Information fractions 1, Sample size 1, evt_max 1",
                     "Simulation results: Information fractions 1, Sample size 2, evt_max 1",
                     "Simulation results: Information fractions 1, Sample size 1, evt_max 2",
                     "Simulation results: Information fractions 1, Sample size 2, evt_max 2")


  list_OF_2=list(s2_OF_design,
                 simu_out_OF_s2_EVTnum_1_EVTmax_1,
                 simu_out_OF_s2_EVTnum_2_EVTmax_1,
                 simu_out_OF_s2_EVTnum_1_EVTmax_2,
                 simu_out_OF_s2_EVTnum_2_EVTmax_2)

  names(list_OF_2)=c("Group sequential design",
                     "Simulation results: Information fractions 2, Sample size 1, evt_max 1",
                     "Simulation results: Information fractions 2, Sample size 2, evt_max 1",
                     "Simulation results: Information fractions 2, Sample size 1, evt_max 2",
                     "Simulation results: Information fractions 2, Sample size 2, evt_max 2")



  list_alsp_1=list(s1_alsp_design,
                   simu_out_alsp_s1_EVTnum_1_EVTmax_1,
                   simu_out_alsp_s1_EVTnum_2_EVTmax_1,
                   simu_out_alsp_s1_EVTnum_1_EVTmax_2,
                   simu_out_alsp_s1_EVTnum_2_EVTmax_2
  )
  names(list_alsp_1)=c("Group sequential design",
                       "Simulation results: Information fractions 1, Sample size 1, evt_max 1",
                       "Simulation results: Information fractions 1, Sample size 2, evt_max 1",
                       "Simulation results: Information fractions 1, Sample size 1, evt_max 2",
                       "Simulation results: Information fractions 1, Sample size 2, evt_max 2")


  list_alsp_2=list(s2_alsp_design,
                   simu_out_alsp_s2_EVTnum_1_EVTmax_1,
                   simu_out_alsp_s2_EVTnum_2_EVTmax_1,
                   simu_out_alsp_s2_EVTnum_1_EVTmax_2,
                   simu_out_alsp_s2_EVTnum_2_EVTmax_2
  )
  names(list_alsp_2)=c("Group sequential design",
                       "Simulation results: Information fractions 2, Sample size 1, evt_max 1",
                       "Simulation results: Information fractions 2, Sample size 2, evt_max 1",
                       "Simulation results: Information fractions 2, Sample size 1, evt_max 2",
                       "Simulation results: Information fractions 2, Sample size 2, evt_max 2")
  out=list(list_OF_1,list_OF_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}
##############################
#masd_simu_OC_survival(groups,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,evt_num_1,evt_num_2,evt_max_1,evt_max_2,beta,sim_num)
#masd_simu_OC_survival(groups=rbind(c(1,1,1),c(0.9,0.85,0.8),c(0.8,0.84,0.85),c(0.5,0.6,0.7)),
#                    s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                   surv_rate=0.5,evt_num_1=200,evt_num_2=300,evt_max_1=800,evt_max_2=1000,HR_cut=0.95,beta=0.1, sim_num=50)
#masd_simu_OC_survival(groups=rbind(c(0.75,0.8,0.85),c(0.7,0.75,0.8),c(0.65,0.7,0.75)),
#                    s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.01,0.04),
#                   surv_rate=0.5,evt_num_1=200,evt_num_2=300,evt_max_1=500,evt_max_2=600,HR_cut=0.95,beta=0.1, sim_num=1000)
##############################################

#########################################################
####  II.3   simulations for binary distributions  ######
#########################################################

####################################

################################################
#######################################
msd_rej_rate_binary<-function(p_t,p_c,direction,s,c_bry,samsz,sim_num){

  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size


  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(p_t)
  N=ceiling(samsz*s)
  NF=samsz
  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej=0; stg=0; samsz_simu_control=N[length(s)]
    samsz_simu_tot=N[length(s)]*(comp+1); arm_rej=rep(0,comp)


      #########################################################
      ### Step 1: Generate random samples ####
      #########################################################

      x0=stats::rbinom(NF,1,p_c)

      x=matrix(nrow=comp,ncol=NF)
      for (i in 1:comp){
        x[i,]=stats::rbinom(NF,1,p_t[i])
      }

      wald_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      score_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      t_hat_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))

      if (direction==1){

        for (j in 1:length(s) ){

          for (i in 1:length(p_t)) {

            p_t_temp=mean(x[i,1:N[j]])
            p_c_temp=mean(x0[1:N[j]])
            p_t_var_temp=p_t_temp*(1-p_t_temp)
            p_c_var_temp=p_c_temp*(1-p_c_temp)
            t_hat_m1[i,j]=N[j]/(1.15*max(p_t_var_temp,p_c_var_temp)+0.85*mean(p_t_var_temp,p_c_var_temp))
            theta_hat_m1[i,j]= p_t_temp-p_c_temp
            theta_hat_m1_sd[i,j]=sqrt(t_hat_m1[i,j])
            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
            score_m1[i,j]=theta_hat_m1[i,j]*t_hat_m1[i,j]
          }
        }

      }else if (direction==2){
        for (j in 1:length(s) ){

          for (i in 1:length(p_t)) {

            p_t_temp=mean(x[i,1:N[j]])
            p_c_temp=mean(x0[1:N[j]])
            p_t_var_temp=p_t_temp*(1-p_t_temp)
            p_c_var_temp=p_c_temp*(1-p_c_temp)
            t_hat_m1[i,j]=N[j]/(1.15*max(p_t_var_temp,p_c_var_temp)+0.85*mean(p_t_var_temp,p_c_var_temp))
            theta_hat_m1[i,j]=p_c_temp- p_t_temp
            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
            score_m1[i,j]=theta_hat_m1[i,j]*t_hat_m1[i,j]
          }
        }
      }

      wald_max_m1=NULL

      for (j in 1:length(s) ){

        wald_max_m1[j]=max(wald_m1[,j])

      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################
      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry

      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################

      if (isTRUE(max(test_max)>=0)){ ##null hypothesis rejected, no more random sampling
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max[test_max>0]
        stg=match(test_max_t[1],test_max)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

            rej=0
            stg=0
            last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)
        }
          #########################################################
          ### Step 4.1. If max(theta<=theta_cut),                ##
          ###     Continue to planned sample size and stop trial  ##
          #########################################################

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej,stg,samsz_simu_control,samsz_simu_tot,arm_rej)

  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  ##################

  rej = as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  samsz_simu_control=as.numeric(simr[,3])
  samsz_simu_tot=simr[,4]
  arm_rej=simr[,5:(4+length(p_t))]

  ######################################

  if (direction==1){
    mu=p_t-p_c
  }else if (direction==2){
    mu=p_c-p_t
  }
  rej_rt=mean(rej)

  arm_rej_rt=NULL

  for (i in 1:length(p_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }


  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
  }

  arm_summary=cbind(p_t,p_c,arm_rej_rt)
  colnames(arm_summary)=c("p<sub>t1</sub>,...","p<sub>c</sub>","Rejection rate by arm")


  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)


  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

###############################
  Design=cbind(c(1:length(s)),s,c_bry,rej_time)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look")

  out_1=c( rej_rt,tyI_er_rt,pos_rej_rt,
           samsz,tot_samsz, mean_control_samsz,mean_tot_samsz)
  names(out_1)=c("Rejection rate","typr I error","non-type I error rejection",
                 "Planned control sample size", "Planned total sample size",
                 "mean control sample size","mean total sample size")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")

  return(out)

}
#######################
#msd_rej_rate_binary(p_t,p_c,direction,s,c_bry,samsz,sim_num)
##################################
#msd_rej_rate_binary(p_t=rep(0.8,3),p_c=0.8,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),samsz=100,sim_num=10000)
#msd_rej_rate_binary(p_t=rep(0.8,3),p_c=0.6,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),samsz=100,sim_num=1000)
#msd_rej_rate_binary(p_t=rep(0.8,3),p_c=0.8,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),samsz=100,sim_num=1000)
################################################

##############################################

#' Multiple comparisons Group sequential design: Simulations, Binary distribution, Type I error and power
#'
#' @param p_t Vector of true event rates of the test arms, one per dose; 0 < p_t < 1
#' @param p_c Event rate of the control arm; 0 < p_c < 1
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param samsz Sample size of the control arm used in the simulation; must be > 0
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary and the rejection rate by look}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the test-arm event rate, the control-arm event rate and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate, the simulated type I error rate, the non-type I error rejection rate, the planned control and total sample sizes and the mean simulated control and total sample sizes}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' msd_power_simu_binary(p_t = c(0.4,0.4), p_c = 0.4, direction = 1, s = c(1/3,2/3,1),
#' bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
msd_power_simu_binary<-function(p_t,p_c,direction,s,bry_type,alpha_2,samsz,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2

  comp=length(p_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out=msd_rej_rate_binary(p_t,p_c,direction,s,c_bry,samsz,sim_num)
  #samsz_tot=samsz*(comp+1)
  out=simu_out
  #names(out)=c("simulation summary","planned control sample size","planned total sample size")
  return(out)
}

##################################
#msd_power_simu_binary(p_t,p_c,s,bry_type,alpha_2,samsz,sim_num)
#msd_power_simu_binary(p_t=c(0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)


#######################################
#############################
masd_rej_simu_adapt_com_binary<-function(p_t,p_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,sim_num,dist){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size

  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(p_t)
  mu=p_t-p_c
  N=ceiling(samsz*s)
  NF=samsz
  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    samsz_simu_control=N[length(s)]; samsz_simu_tot=N[length(s)]*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)

      #########################################################
      ### Step 1: Generate random samples ####
      #########################################################

      x0=stats::rbinom(N_max,1,p_c)

      x=matrix(nrow=comp,ncol=N_max)
      for (i in 1:comp){
        x[i,]=stats::rbinom(N_max,1,p_t[i])
      }

      wald_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      score_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      t_hat_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))

      if (direction==1){

        for (j in 1:length(s) ){

          for (i in 1:length(p_t)) {

            p_t_temp=mean(x[i,1:N[j]])
            p_c_temp=mean(x0[1:N[j]])
            p_t_var_temp=p_t_temp*(1-p_t_temp)
            p_c_var_temp=p_c_temp*(1-p_c_temp)
            t_hat_m1[i,j]=N[j]/(1.15*max(p_t_var_temp,p_c_var_temp)+0.85*mean(p_t_var_temp,p_c_var_temp))
            theta_hat_m1[i,j]= p_t_temp-p_c_temp
            theta_hat_m1_sd[i,j]=sqrt(t_hat_m1[i,j])
            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
            score_m1[i,j]=theta_hat_m1[i,j]*t_hat_m1[i,j]
          }
        }

      }else if (direction==2){
        for (j in 1:length(s) ){

          for (i in 1:length(p_t)) {

            p_t_temp=mean(x[i,1:N[j]])
            p_c_temp=mean(x0[1:N[j]])
            p_t_var_temp=p_t_temp*(1-p_t_temp)
            p_c_var_temp=p_c_temp*(1-p_c_temp)
            t_hat_m1[i,j]=N[j]/(1.15*max(p_t_var_temp,p_c_var_temp)+0.85*mean(p_t_var_temp,p_c_var_temp))
            theta_hat_m1[i,j]=p_c_temp- p_t_temp
            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
            score_m1[i,j]=theta_hat_m1[i,j]*t_hat_m1[i,j]
          }
        }
      }

      T_hat=t_hat_m1[,length(s)]
      theta_hat_m1_sd=1/sqrt(t_hat_m1)
      wald_max_m1=NULL

      for (j in 1:length(s) ){

        wald_max_m1[j]=max(wald_m1[,j])

      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop=rep(0,length(p_t))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }
      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################
      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry
      #  test=wald_m1-c_bry
      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################
      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }

      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling

        futi=0
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        if (futi==1) { # If TRUE, thiere is not satisfactory efficacy, stop trial at planned sample size

          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            stg=0}
          last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)

        }else if (futi==0){ #beginning of step 4.2

          #########################################################
          ### Step 4.1. If max(theta<=theta_cut),                ##
          ###     Continue to planned sample size and stop trial  ##
          #########################################################

          mu_drop=match(arm_drop[arm_drop>0],arm_drop)

          if (length(mu_drop)==0){

            mu_keep=c(1:comp)

          }else if ( length(mu_drop)>0){

            mu_keep=c(1:comp)[-mu_drop]

          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_hat_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2
          n_control_inter=N[length(s)-1]

          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1=theta_hat_inter_m1,theta_hat_inter_m1_sd=theta_hat_inter_m1_sd,dist=dist)

          if (isTRUE(c_power_m2>= 1-beta)){
            ########################################################################
            ### Step 4.2. If max(theta>theta_cut), calculate conditional power    ##
            ###          and re-estiate sample size                               ##
            ########################################################################
            last_vt=length(s)

            wald_m2=wald_m1[mu_keep,]

            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last

            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }
            samsz_simu_control=N[last_vt]
            samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)   #For q-th (collumn number is q) simulation,

            ########################################################################################
            ## Step 4.2.1. If futi==0. Check if there is need to drop some doses due to efficacy ###
            ### A dose is dropped for lack of efficacy if: theta_hat<theta cut                   ###
            ### or of min theta_hat is less than 1/3 of the max theta_hat                        ###
            ########################################################################################
          }else if (isTRUE(c_power_m2<1-beta)){

               t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                             theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

              ## default when the re-estimation is not evaluable
              N_new_t=samsz
              if (isTRUE(as.numeric(t_new_0[2])<=1)){
                N_new_t=samsz
              }else if (isTRUE(as.numeric(t_new_0[2])>1)){
                T_new_0=as.numeric(t_new_0[[1]])

                N_new_t_0=max(ceiling(n_control_inter*(T_new_0/T_inter)))
                N_new_t_1=max(samsz,N_new_t_0)
                N_new_t= min(N_new_t_1,N_max)
              }

            if(isTRUE(N_new_t==samsz)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              samsz_simu_control=N[last_vt]
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)

            }else  if(isTRUE(N_new_t>samsz)){

              x0_new=x0[1:N_new_t]
              x_new=x[c(mu_keep),1:N_new_t]

                p_c_new=mean(x0_new)
                p_c_var_new=p_c_new*(1-p_c_new)

                t_hat_m2_new=NULL
                theta_hat_m2_new=NULL
                theta_hat_m2_sd_new=NULL
                wald_m2_new=NULL
                p_t_new_m2=NULL
                p_t_var_new=NULL

                if (direction==1){

                  for (i in 1:length(mu_keep)){
                    p_t_new_m2[i]=mean(x_new[i,])
                    p_t_var_new[i]=p_t_new_m2[i]*(1-p_t_new_m2[i])
                    theta_hat_m2_new[i]=mean(x_new[i,])-mean(x0_new)
                    t_hat_m2_new[i]=N_new_t/( 1.15*max(p_t_var_new[i],p_c_var_new)+0.85*mean(p_t_var_new[i],p_c_var_new))

                    theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                    wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                  }

                }else if (direction==2){

                  for (i in 1:length(mu_keep)){
                    p_t_new_m2[i]=mean(x_new[i,])
                    theta_hat_m2_new[i]=mean(x0_new)-mean(x_new[i,])
                    t_hat_m2_new[i]=N_new_t/(p_t_new_m2[i]*(1-p_t_new_m2[i])+p_c_new*(1-p_c_new))
                    theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                    wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                  }

                  wald_max_m2_new=max(wald_m2_new)
                }

                t_hat_m2_new_max=max(t_hat_m2_new)

                T_new=t_hat_m2_new
                T_inter=(N[length(s)-1]/N_new_t)*T_new
                T_hat_m2=T_hat[mu_keep]

                w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=1)
                w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=1)
                c_alpha_m2=c_alpha_m1*w_m1/w_m2

                e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                        theta_hat_inter_m2=theta_hat_inter_m2,
                                        theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=1)

                c_bry_new=NULL

                for (i in 1:length(theta_hat_inter_m2)){
                  c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
                }

                test_new=wald_m2_new-max(c_bry_new)

              #################################################################################
              ## Step 4.2.1.3.3. Hypothesis testing: compare the Wald statistics with the critical boundaries ##
              #################################################################################


              samsz_simu_control=N_new_t
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N_new_t-N[length(s)-1])*(length(mu_keep)+1)

              if (isTRUE(max(test_new)>0)){

                rej=1

                stg=length(s)

              }else{
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0
              }
              arm_rej_new=NULL

              for (i in 1:length(mu_keep)){

                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }

              }

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }

        }
      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej_nossr,rej,stg,samsz_simu_control,samsz_simu_tot,futi,arm_drop,arm_rej)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr


  rej_nossr = as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  samsz_simu_control=as.numeric(simr[,4])
  samsz_simu_tot=simr[,5]
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(p_t))]
  arm_rej=simr[,(7+length(p_t)):(6+2*length(p_t))]
  # unique(case)
  ######################################

  if (direction==1){
    mu=p_t-p_c
  }else if (direction==2){
    mu=p_c-p_t
  }
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  arm_rej_rt=NULL

  for (i in 1:length(p_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  arm_drop_rt=NULL

  for (i in 1:length(p_t)){
    arm_drop_rt[i]=mean(arm_drop[,i])
  }

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
  }

  futi_rt=mean(futi)

  arm_summary=cbind(p_t,p_c,arm_drop_rt,arm_rej_rt)
  colnames(arm_summary)=c("p<sub>t1</sub>,...","p<sub>c</sub>","drop rate by arm","Rejection rate by arm")


  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)
  max_tot_samsz=N_max*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu_control[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= samsz_simu_control[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)
  #############
  Design=cbind(c(1:length(s)),s,c_bry,rej_time,samsz_icr_byLook)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look","Rate of sample size increase")

  out_1=c(rej_nossr_rt, rej_rt, samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,samsz,N_max,mean_control_samsz, mean_samsz_simu_control_samsz_icr)
  names(out_1)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "type I error","non-type I error rejection","futility rate","Planned control sample size",
                 "Maximum control sample size",
                 "mean control sample size","mean control sample size conditional on increase")
  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")

  return(out)

}

#######################################
#masd_rej_simu_adapt_com_binary(p_t,p_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,sim_num,dist)
# masd_rej_simu_adapt_com_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.05,beta=0.1,samsz=100,N_max=200,sim_num=1000,dist=1)
#masd_rej_simu_adapt_com_binary(p_t=c(0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.05,beta=0.1,samsz=100,N_max=200,sim_num=100000,dist=1)

##############################################

#' Multiple comparisons Group sequential design: Simulations, Binary distribution, Operating characteristics
#'
#' @param p_t Vector of true event rates of the test arms, one per dose; 0 < p_t < 1
#' @param p_c Event rate of the control arm; 0 < p_c < 1
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param theta_cut Lack-of-efficacy threshold on the effect-size scale: a test arm is
#'   dropped at the interim analysis when its observed effect size is below this value
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param samsz Sample size of the control arm used in the simulation; must be > 0
#' @param N_max Maximum sample size for the control arm allowed by the sample size
#'   re-estimation; must be >= samsz
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary, the rejection rate by look and the rate of sample size increase}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the test-arm event rate, the control-arm event rate, the drop rate by arm and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate without and with sample size re-estimation, the rate of sample size increase, the rejection rates conditional on meeting the sample size increase criteria, the simulated type I error rate, the non-type I error rejection rate, the futility rate, the planned and maximum control sample sizes, the mean simulated control sample size and the mean simulated control sample size conditional on an increase}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' masd_power_simu_adapt_binary(p_t = c(0.4,0.4,0.4), p_c = 0.4, direction = 1, s = c(1/3,2/3,1),
#' bry_type = 1, alpha_2 = 0.05, theta_cut = 0.01, beta = 0.1, samsz = 100,
#' N_max = 200, sim_num = 20)
#' }
#'
#' @export
masd_power_simu_adapt_binary<-function(p_t,p_c,direction,s,bry_type,alpha_2,theta_cut,beta,samsz,N_max,sim_num){


  alpha = alpha_2/2
  comp = length(p_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out = masd_rej_simu_adapt_com_binary(p_t=p_t,p_c=p_c,direction=direction,s=s,
                                            c_bry=c_bry,theta_cut=theta_cut,beta=beta,samsz=samsz,N_max=N_max,sim_num=sim_num,dist=1)

  out=simu_out
  return(out)
}

##################################
#  masd_power_simu_adapt_binary(p_t,p_c,direction,s,bry_type,alpha_2,theta_cut,beta,samsz,N_max,sim_num)
#masd_power_simu_adapt_binary(p_t=c(0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.01,beta=0.1,samsz=100,N_max=500,sim_num=100)
# masd_power_simu_adapt_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.01,beta=0.1,samsz=100,N_max=200,sim_num=10000)
#  masd_power_simu_adapt_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.01,beta=0.1,samsz=100,N_max=200,sim_num=100000)
#  masd_power_simu_adapt_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.01,beta=0.1,samsz=100,N_max=200,sim_num=10000)

#masd_power_simu_adapt_binary(p_t=c(0.6,0.7,0.8),p_c=0.6,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.01,beta=0.1,samsz=100,N_max=400,sim_num=1000)
#############################

###################################################################################
## Simulations for Optimizing Operating characteristics - Binary distributions  ###
###################################################################################

#######################################

####################
##################
masd_rej_OC_simu_adapt_com_binary<-function(p_t,p_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,sim_num,dist){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size

  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(p_t)
  mu=p_t-p_c
  N=ceiling(samsz*s)
  NF=samsz
  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    samsz_simu_control=N[length(s)]; samsz_simu_tot=N[length(s)]*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)

      #########################################################
      ### Step 1: Generate random samples ####
      #########################################################

      x0=stats::rbinom(N_max,1,p_c)

      x=matrix(nrow=comp,ncol=N_max)
      for (i in 1:comp){
        x[i,]=stats::rbinom(N_max,1,p_t[i])
      }

      wald_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      score_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      t_hat_m1=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))

      if (direction==1){

        for (j in 1:length(s) ){

          for (i in 1:length(p_t)) {

            p_t_temp=mean(x[i,1:N[j]])
            p_c_temp=mean(x0[1:N[j]])
            p_t_var_temp=p_t_temp*(1-p_t_temp)
            p_c_var_temp=p_c_temp*(1-p_c_temp)
            t_hat_m1[i,j]=N[j]/(1.15*max(p_t_var_temp,p_c_var_temp)+0.85*mean(p_t_var_temp,p_c_var_temp))
            theta_hat_m1[i,j]= p_t_temp-p_c_temp
            theta_hat_m1_sd[i,j]=sqrt(t_hat_m1[i,j])
            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
            score_m1[i,j]=theta_hat_m1[i,j]*t_hat_m1[i,j]
          }
        }

      }else if (direction==2){
        for (j in 1:length(s) ){

          for (i in 1:length(p_t)) {

            p_t_temp=mean(x[i,1:N[j]])
            p_c_temp=mean(x0[1:N[j]])
            p_t_var_temp=p_t_temp*(1-p_t_temp)
            p_c_var_temp=p_c_temp*(1-p_c_temp)
            t_hat_m1[i,j]=N[j]/(1.15*max(p_t_var_temp,p_c_var_temp)+0.85*mean(p_t_var_temp,p_c_var_temp))
            theta_hat_m1[i,j]=p_c_temp- p_t_temp
            wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
            score_m1[i,j]=theta_hat_m1[i,j]*t_hat_m1[i,j]
          }
        }
      }

      T_hat=t_hat_m1[,length(s)]
      theta_hat_m1_sd=1/sqrt(t_hat_m1)
      wald_max_m1=NULL

      for (j in 1:length(s) ){

        wald_max_m1[j]=max(wald_m1[,j])

      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop=rep(0,length(p_t))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }
      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################
      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry
      #  test=wald_m1-c_bry
      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################
      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }

      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling

        futi=0
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        if (futi==1) { # If TRUE, thiere is not satisfactory efficacy, stop trial at planned sample size

          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            stg=0}
          last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)

        }else if (futi==0){ #beginning of step 4.2

          #########################################################
          ### Step 4.1. If max(theta<=theta_cut),                ##
          ###     Continue to planned sample size and stop trial  ##
          #########################################################

          mu_drop=match(arm_drop[arm_drop>0],arm_drop)

          if (length(mu_drop)==0){

            mu_keep=c(1:comp)

          }else if ( length(mu_drop)>0){

            mu_keep=c(1:comp)[-mu_drop]

          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_hat_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2
          n_control_inter=N[length(s)-1]

          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1=theta_hat_inter_m1,theta_hat_inter_m1_sd=theta_hat_inter_m1_sd,dist=dist)

          if (isTRUE(c_power_m2>= 1-beta)){
            ########################################################################
            ### Step 4.2. If max(theta>theta_cut), calculate conditional power    ##
            ###          and re-estiate sample size                               ##
            ########################################################################
            last_vt=length(s)

            wald_m2=wald_m1[mu_keep,]

            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last

            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }
            samsz_simu_control=N[last_vt]
            samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)   #For q-th (collumn number is q) simulation,

            ########################################################################################
            ## Step 4.2.1. If futi==0. Check if there is need to drop some doses due to efficacy ###
            ### A dose is dropped for lack of efficacy if: theta_hat<theta cut                   ###
            ### or of min theta_hat is less than 1/3 of the max theta_hat                        ###
            ########################################################################################
          }else if (isTRUE(c_power_m2<1-beta)){

            t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                           theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

            ## default when the re-estimation is not evaluable
            N_new_t=samsz
            if (isTRUE(as.numeric(t_new_0[2])<=1)){
              N_new_t=samsz
            }else if (isTRUE(as.numeric(t_new_0[2])>1)){
              T_new_0=as.numeric(t_new_0[[1]])

              N_new_t_0=max(ceiling(n_control_inter*(T_new_0/T_inter)))
              N_new_t_1=max(samsz,N_new_t_0)
              N_new_t= min(N_new_t_1,N_max)
            }

            if(isTRUE(N_new_t==samsz)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              samsz_simu_control=N[last_vt]
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)

            }else  if(isTRUE(N_new_t>samsz)){

              x0_new=x0[1:N_new_t]
              x_new=x[c(mu_keep),1:N_new_t]

              p_c_new=mean(x0_new)
              p_c_var_new=p_c_new*(1-p_c_new)

              t_hat_m2_new=NULL
              theta_hat_m2_new=NULL
              theta_hat_m2_sd_new=NULL
              wald_m2_new=NULL
              p_t_new_m2=NULL
              p_t_var_new=NULL

              if (direction==1){

                for (i in 1:length(mu_keep)){
                  p_t_new_m2[i]=mean(x_new[i,])
                  p_t_var_new[i]=p_t_new_m2[i]*(1-p_t_new_m2[i])
                  theta_hat_m2_new[i]=mean(x_new[i,])-mean(x0_new)
                  t_hat_m2_new[i]=N_new_t/( 1.15*max(p_t_var_new[i],p_c_var_new)+0.85*mean(p_t_var_new[i],p_c_var_new))

                  theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                  wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                }

              }else if (direction==2){

                for (i in 1:length(mu_keep)){
                  p_t_new_m2[i]=mean(x_new[i,])
                  theta_hat_m2_new[i]=mean(x0_new)-mean(x_new[i,])
                  t_hat_m2_new[i]=N_new_t/(p_t_new_m2[i]*(1-p_t_new_m2[i])+p_c_new*(1-p_c_new))
                  theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                  wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                }

                wald_max_m2_new=max(wald_m2_new)
              }

              t_hat_m2_new_max=max(t_hat_m2_new)

              T_new=t_hat_m2_new
              T_inter=(N[length(s)-1]/N_new_t)*T_new
              T_hat_m2=T_hat[mu_keep]

              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=1)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=1)
              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=1)

              c_bry_new=NULL

              for (i in 1:length(theta_hat_inter_m2)){
                c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
              }

              test_new=wald_m2_new-max(c_bry_new)

              #################################################################################
              ## Step 4.2.1.3.3. Hypothesis testing: compare the Wald statistics with the critical boundaries ##
              #################################################################################


              samsz_simu_control=N_new_t
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N_new_t-N[length(s)-1])*(length(mu_keep)+1)

              if (isTRUE(max(test_new)>0)){

                rej=1

                stg=length(s)

              }else{
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0
              }
              arm_rej_new=NULL

              for (i in 1:length(mu_keep)){

                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }

              }

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }

        }
      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej_nossr,rej,stg,samsz_simu_control,samsz_simu_tot,futi,arm_drop,arm_rej)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr


  rej_nossr = as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  samsz_simu_control=as.numeric(simr[,4])
  samsz_simu_tot=simr[,5]
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(p_t))]
  arm_rej=simr[,(7+length(p_t)):(6+2*length(p_t))]

  ######################################

  if (direction==1){
    mu=p_t-p_c
  }else if (direction==2){
    mu=p_c-p_t
  }
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  arm_rej_rt=NULL

  for (i in 1:length(p_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  arm_drop_rt=NULL

  for (i in 1:length(p_t)){
    arm_drop_rt[i]=mean(arm_drop[,i])
  }

  futi_rt=mean(futi)

  arm_summary=cbind(p_t,p_c,arm_drop_rt,arm_rej_rt)
  colnames(arm_summary)=c("p_t","p_c","drop rate by arm","Rejection rate by arm")

  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)
  max_tot_samsz=N_max*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu_control[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= samsz_simu_control[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)
  #############

  out=c(rej_nossr_rt, rej_rt, samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,samsz,N_max,mean_control_samsz, mean_samsz_simu_control_samsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "type I error","non-type I error rejection","futility rate","Planned control sample size",
                 "Maximum control sample size",
                 "mean control sample size","mean control sample size conditional on increase")

  return(out)

}

#masd_rej_OC_simu_adapt_com_binary(p_t,p_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,sim_num,dist)
# masd_rej_OC_simu_adapt_com_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.05,beta=0.1,samsz=100,N_max=200,sim_num=1000,dist=1)
#masd_rej_OC_simu_adapt_com_binary(p_t=c(0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.05,beta=0.1,samsz=100,N_max=200,sim_num=10000,dist=1)
##############################################

masd_OC_power_simu_adapt_binary<-function(p_t,p_c,direction,s,bry_type,alpha_2,theta_cut,beta,samsz,N_max,sim_num){


  alpha = alpha_2/2
  comp = length(p_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out = masd_rej_OC_simu_adapt_com_binary(p_t=p_t,p_c=p_c,direction=direction,s=s,
                                            c_bry=c_bry,theta_cut=theta_cut,beta=beta,samsz=samsz,N_max=N_max,sim_num=sim_num,dist=1)
  p_t_names=paste(paste("p<sub>t</sub>","[",c(1:length(p_t)),"]",sep=""))
  out= round(c(p_t,p_c,as.numeric(simu_out)),4)
  names(out)=c(p_t_names,"p<sub>c</sub>",names(simu_out))
 # out=simu_out
  return(out)
}

##################################
#  masd_OC_power_simu_adapt_binary(p_t,p_c,direction,s,bry_type,alpha_2,theta_cut,beta,samsz,N_max,sim_num)
#masd_OC_power_simu_adapt_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.01,beta=0.1,samsz=100,N_max=500,sim_num=10000)
# masd_OC_power_simu_adapt_binary(p_t=c(0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.01,beta=0.1,samsz=100,N_max=200,sim_num=10000)
#  masd_OC_power_simu_adapt_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.01,beta=0.1,samsz=100,N_max=200,sim_num=10000)
#masd_OC_power_simu_adapt_binary(p_t=c(0.6,0.7,0.8),p_c=0.6,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.01,beta=0.1,samsz=100,N_max=200,sim_num=1000)


#############################

###################
#' Multiple comparisons Group sequential design: Simulations, Binary distribution, Optimizing adaptive sequential design
#'
#' @param groups Matrix of true event rates of the test arms: each row is one simulation scenario and each column one test arm/dose; 0 < groups < 1
#' @param p_c Event rate of the control arm; 0 < p_c < 1
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g. c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Two-sided type I error levels spent at each look of the schedule s1; must sum to alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error levels spent at each look of the schedule s2; must sum to alpha_2
#' @param theta_cut Lack-of-efficacy threshold on the effect-size scale: a test arm is
#'   dropped at the interim analysis when its observed effect size is below this value
#' @param samsz_1 Sample size of the control arm used in the simulation for the first sample size scenario; must be > 0
#' @param samsz_2 Sample size of the control arm used in the simulation for the second sample size scenario; must be > 0
#' @param N_max_1 Maximum sample size for the control arm allowed by the sample size re-estimation in the first N_max scenario; must be >= samsz_1
#' @param N_max_2 Maximum sample size for the control arm allowed by the sample size re-estimation in the second N_max scenario; must be >= samsz_2
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the O'Brien-Fleming group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the O'Brien-Fleming group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the alpha-spending group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the alpha-spending group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' masd_simu_OC_binary(groups = rbind(c(0.35,0.4,0.45),c(0.4,0.45,0.5),c(0.35,0.45,0.55),
#' c(0.45,0.5,0.55)), p_c = 0.3, direction = 1, s1 = c(1/3,2/3,1), s2 = c(0.5,1),
#' alpha_2 = 0.05, alpha_2_alps_s1 = c(0.001,0.001,0.048),
#' alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01, samsz_1 = 50, samsz_2 = 75,
#' N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
masd_simu_OC_binary<-function(groups,p_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)
  comp=length(groups[1,])

  dist=1
  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=dist)[[1]][1]))
  s2_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=dist)[[1]][1]))

  s1_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry (s=s1,alpha=alpha_alsp_s1,comp=comp,dist=dist)[[1]][1]))
  s2_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry (s=s2,alpha=alpha_alsp_s2,comp=comp,dist=dist)[[1]][1]))

  ##### Two-stage GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")

  ##################
  p_t_vec=groups
  rownames( p_t_vec)=NULL

  p_t_group_names=paste(paste("p<sub>t</sub> group",c(1:dim(p_t_vec)[1]),sep=" "))
  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####

  simu_out_OF_s1_samsz_1_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s1_samsz_1_N_max_1=rbind(simu_out_OF_s1_samsz_1_N_max_1,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_1)=p_t_group_names

  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_1=rbind(simu_out_OF_s2_samsz_1_N_max_1,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_1)=p_t_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_1=rbind(simu_out_OF_s1_samsz_2_N_max_1,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_1)=p_t_group_names

  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_1=rbind(simu_out_OF_s2_samsz_2_N_max_1,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_2_N_max_1)=p_t_group_names
  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####

  simu_out_OF_s1_samsz_1_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s1_samsz_1_N_max_2=rbind(simu_out_OF_s1_samsz_1_N_max_2,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_2)=p_t_group_names

  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_2=rbind(simu_out_OF_s2_samsz_1_N_max_2,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_2)=p_t_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_2=rbind(simu_out_OF_s1_samsz_2_N_max_2,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_2)=p_t_group_names

  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_2=rbind(simu_out_OF_s2_samsz_2_N_max_2,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_2_N_max_2)=p_t_group_names

  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################


  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####

  simu_out_alsp_s1_samsz_1_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s1_samsz_1_N_max_1=rbind(simu_out_alsp_s1_samsz_1_N_max_1,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                               bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_1)=p_t_group_names

  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_1=rbind(simu_out_alsp_s2_samsz_1_N_max_1,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                               bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_1)=p_t_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_1=rbind(simu_out_alsp_s1_samsz_2_N_max_1,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                               bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_1)=p_t_group_names

  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_1=rbind(simu_out_alsp_s2_samsz_2_N_max_1,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                               bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_2_N_max_1)=p_t_group_names
  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####

  simu_out_alsp_s1_samsz_1_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s1_samsz_1_N_max_2=rbind(simu_out_alsp_s1_samsz_1_N_max_2,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                               bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_2)=p_t_group_names

  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_2=rbind(simu_out_alsp_s2_samsz_1_N_max_2,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                               bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_2)=p_t_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_2=rbind(simu_out_alsp_s1_samsz_2_N_max_2,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                               bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_2)=p_t_group_names

  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_2=rbind(simu_out_alsp_s2_samsz_2_N_max_2,masd_OC_power_simu_adapt_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                               bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                               beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_2_N_max_2)=p_t_group_names

  ##################################
  ##############

  list_OF_1=list(s1_OF_design,
                 simu_out_OF_s1_samsz_1_N_max_1,
                 simu_out_OF_s1_samsz_2_N_max_1,
                 simu_out_OF_s1_samsz_1_N_max_2,
                 simu_out_OF_s1_samsz_2_N_max_2)
  names(list_OF_1)=c("Group sequential design",
                     "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                     "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                     "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                     "Simulation results: Information fractions 1, Sample size 2, N_max 2")
  list_OF_2=list(s2_OF_design,
                 simu_out_OF_s2_samsz_1_N_max_1,
                 simu_out_OF_s2_samsz_2_N_max_1,
                 simu_out_OF_s2_samsz_1_N_max_2,
                 simu_out_OF_s2_samsz_2_N_max_2)
  names(list_OF_2)=c("Group sequential design",
                     "Simulation results: Information fractions 2, Sample size 1, N_max 1",
                     "Simulation results: Information fractions 2, Sample size 2, N_max 1",
                     "Simulation results: Information fractions 2, Sample size 1, N_max 2",
                     "Simulation results: Information fractions 2, Sample size 2, N_max 2")


  list_alsp_1=list(s1_alsp_design,
                   simu_out_alsp_s1_samsz_1_N_max_1,
                   simu_out_alsp_s1_samsz_2_N_max_1,
                   simu_out_alsp_s1_samsz_1_N_max_2,
                   simu_out_alsp_s1_samsz_2_N_max_2)
  names(list_alsp_1)=c("Group sequential design",
                       "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                       "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                       "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                       "Simulation results: Information fractions 1, Sample size 2, N_max 2")

  list_alsp_2=list(s2_alsp_design,
                   simu_out_alsp_s2_samsz_1_N_max_1,
                   simu_out_alsp_s2_samsz_2_N_max_1,
                   simu_out_alsp_s2_samsz_1_N_max_2,
                   simu_out_alsp_s2_samsz_2_N_max_2)
  names(list_alsp_2)=c("Group sequential design",
                       "Simulation results: Information fractions 2, Sample size 1, N_max 1",
                       "Simulation results: Information fractions 2, Sample size 2, N_max 1",
                       "Simulation results: Information fractions 2, Sample size 1, N_max 2",
                       "Simulation results: Information fractions 2, Sample size 2, N_max 2")


  out=list(list_OF_1,list_OF_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}

#################
#masd_simu_OC_binary(groups,p_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num)
#masd_simu_OC_binary(groups=rbind(c(0.2,0.3,0.4),c(0.3,0.4,0.5),c(0.4,0.5,0.6),c(0.5,0.6,0.7)),p_c=0.2,direction=1,
#                            s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                           theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1, sim_num=100)

##############################################
##############################################
########################################################
####  II.4  Simulations for Poisson distributions  #####
########################################################

################################################
#######################################
msd_rej_rate_poisson<-function(lambda_c,lambda_t,direction,s,c_bry,samsz,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size
  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(lambda_t)
  N=ceiling(samsz*s)
  NF=samsz
  if (direction==1){
    mu=lambda_t-lambda_c
  }else if (direction==2){
    mu=lambda_c-lambda_t
  }
  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej=0; stg=0; samsz_simu_control=N[length(s)]
    samsz_simu_tot=N[length(s)]*(comp+1); arm_rej=rep(0,comp)

      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      #the loop specifies what is generated for each simulation

      #######################
      x0=stats::rpois(NF,lambda_c)
      #x0 is the random sample at the first interim analysis for the control
      x=matrix(nrow=comp,ncol=NF)
      for (i in 1:comp){
        x[i,]=stats::rpois(NF,lambda_t[i])
      }

      lambda_c_hat=NULL

      for (j in 1:length(s)){

        lambda_c_hat[j]=mean(x0[1:N[j]])
      }

      lambda_t_hat_m1=matrix(nrow=comp,ncol=length(s))

      for (i in 1:comp ){
        for (j in 1:length(s)){

          lambda_t_hat_m1[i,j]=mean(x[i,1:N[j]])
        }
      }

      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      if (direction==2){
        for (i in 1:comp ){
          for (j in 1:length(s)){

            theta_hat_m1[i,j]=lambda_c_hat[j]-lambda_t_hat_m1[i,j]
          }
        }
      } else if (direction==1){

        for (i in 1:comp ){
          for (j in 1:length(s)){
            theta_hat_m1[i,j]=lambda_t_hat_m1[i,j]-lambda_c_hat[j]
          }
        }
      }
      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))

      for (i in 1:comp ){
        for (j in 1:length(s)){
          t_hat_m1[i,j]= N[j]/(0.01+ lambda_t_hat_m1[i,j]+lambda_c_hat[j])
          wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
          theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])
        }
      }

      wald_max_m1=NULL
      for (j in 1:length(s) ){

        wald_max_m1[j]=max(wald_m1[,j])
      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################
      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################

      if (isTRUE(max(test_max)>=0)){ ##null hypothesis rejected, no more random sampling

        rej=1
        test_max_t=test_max[test_max>0]
        stg=match(test_max_t[1],test_max)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

            rej=0
            stg=0
            last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)

        }
      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej,stg,samsz_simu_control,samsz_simu_tot,arm_rej)

  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  ########################################
  rej = as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  samsz_simu_control=as.numeric(simr[,3])
  samsz_simu_tot=simr[,4]
  arm_rej=simr[,5:(4+length(lambda_t))]
  ######################################
  rej_rt=mean(rej)

  arm_rej_rt=NULL
  for (i in 1:length(lambda_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

    rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
  }

  arm_summary=cbind(lambda_t,lambda_c,arm_rej_rt)


  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL
  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  Design=cbind(c(1:length(s)),s,c_bry,rej_time)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look")

  out_1=c( rej_rt,tyI_er_rt,pos_rej_rt,samsz,tot_samsz,mean_control_samsz,mean_tot_samsz)
  names(out_1)=c("Rejection rate","typr I error","non-type I error rejection",
                 "Planned control sample size","Planned total sample size",
                 "mean control sample size","mean total sample size")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")

  return(out)

}

######################################
#msd_rej_rate_poisson(lambda_c,lambda_t,direction,s,c_bry,samsz,sim_num)

#msd_rej_rate_poisson(lambda_c=1.1,lambda_t=c(0.8,0.9,1),s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),samsz=100,sim_num=100)
#msd_rej_rate_poisson(lambda_c=1.1,lambda_t=c(1.1,1.1,1.1),s=c(1/3,2/3,1),c_bry=c(4.183799, 2.958392,2.415517),samsz=100,sim_num=100000)
#msd_rej_rate_poisson(lambda_c=1.1,lambda_t=c(1.1,1.1,1.1),direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),samsz=100,sim_num=100000)


################################################

##############################################

#' Multiple comparisons Group sequential design: Simulations, Poisson distribution, Type I error and power
#'
#' @param lambda_c Mean event count of the control arm; must be > 0
#' @param lambda_t Vector of true mean event counts of the test arms, one per dose; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param samsz Sample size of the control arm used in the simulation; must be > 0
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary and the rejection rate by look}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the test-arm mean, the control-arm mean and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate, the simulated type I error rate, the non-type I error rejection rate, the planned control and total sample sizes and the mean simulated control and total sample sizes}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' msd_power_simu_poisson(lambda_c = 0.8, lambda_t = c(0.8,0.8,0.8), direction = 1,
#' s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
msd_power_simu_poisson<-function(lambda_c,lambda_t,direction,s,bry_type,alpha_2,samsz,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2

  comp=length(lambda_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out=msd_rej_rate_poisson(lambda_c,lambda_t,direction,s=s,c_bry=c_bry,samsz=samsz,sim_num=sim_num)
  # simu_out=msd_rej_rate_normal(mu=mu,s,c_bry,sigma,samsz,sim_num)
  #samsz_tot=samsz*(comp+1)
  out=simu_out
  #names(out)=c("simulation summary","planned control","planned total")
  return(out)
}

##################################
#msd_power_simu_poisson(lambda_c,lambda_t,direction,s,bry_type,alpha_2,samsz,sim_num)
#msd_power_simu_poisson(lambda_c=0.8,lambda_t=c(0.8,0.8,0.8),direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=100)
#msd_power_simu_poisson(lambda_c=0.8,lambda_t=c(0.8,1,0.9),direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=10000)


#######################################

###############
masd_rej_simu_adapt_com_poisson<-function(lambda_t,lambda_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,dist,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size
  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(lambda_t)
  N=ceiling(samsz*s)
  # N is the vetoc of sample sizes at which interim analysis will be conducted
  N1=N[1]
  NF=samsz
  if (direction==1){
    mu=lambda_t-lambda_c
  }else if (direction==2){
    mu=lambda_c-lambda_t
  }
  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    samsz_simu_control=N[length(s)]; samsz_simu_tot=N[length(s)]*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)

      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      #the loop specifies what is generated for each simulation

      #######################
      x0=stats::rpois(N_max,lambda_c)
      #x0 is the random sample at the first interim analysis for the control
      x=matrix(nrow=comp,ncol=N_max)
      for (i in 1:comp){
        x[i,]=stats::rpois(N_max,lambda_t[i])
      }

      #x is the matrix of random samples for each dose
      #for each simulation.
      #x matirx does not include the control sample.
      #The row number is the same as the active dose
      lambda_c_hat=NULL

      for (j in 1:length(s)){

        lambda_c_hat[j]=mean(x0[1:N[j]])
      }

      lambda_t_hat_m1=matrix(nrow=comp,ncol=length(s))

      for (i in 1:comp ){
        for (j in 1:length(s)){

          lambda_t_hat_m1[i,j]=mean(x[i,1:N[j]])
        }
      }

      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      if (direction==2){
        for (i in 1:comp ){
          for (j in 1:length(s)){

            theta_hat_m1[i,j]=lambda_c_hat[j]-lambda_t_hat_m1[i,j]
          }
        }
      } else if (direction==1){

        for (i in 1:comp ){
          for (j in 1:length(s)){
            theta_hat_m1[i,j]=lambda_t_hat_m1[i,j]-lambda_c_hat[j]
          }
        }
      }
      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))

      for (i in 1:comp ){
        for (j in 1:length(s)){
          t_hat_m1[i,j]= N[j]/( 0.01+lambda_t_hat_m1[i,j]+lambda_c_hat[j])
          wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
          theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])
        }
      }

      wald_max_m1=NULL
      for (j in 1:length(s) ){

        wald_max_m1[j]=max(wald_m1[,j])
      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop=rep(0,length(lambda_t))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }

      T_hat=t_hat_m1[,length(s) ]
      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################
      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry
      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################
      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }

      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling
        futi=0
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

                samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        if (futi==1) { # If TRUE, thiere is not satisfactory efficacy, stop trial at planned sample size

          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            stg=0}
          last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)
          mu_reduc=0
        }else if (futi==0){ #beginning of step 4.2

          mu_drop=match(arm_drop[arm_drop>0],arm_drop)

          if (length(mu_drop)==0){

            mu_keep=c(1:comp)

          }else if ( length(mu_drop)>0){

            mu_keep=c(1:comp)[-mu_drop]

          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_hat_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2
          n_control_inter=N[length(s)-1]

          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1=theta_hat_inter_m1,theta_hat_inter_m1_sd=theta_hat_inter_m1_sd,dist=dist)

          if (isTRUE(c_power_m2>= 1-beta)){
            #########################################
            ### Proceed to planned sample size  #####
            ### and conduct final testing       #####
            #########################################
            last_vt=length(s)

            wald_m2=wald_m1[mu_keep,]

            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last

            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }
          samsz_simu_control=N[last_vt]
            samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)   #For q-th (collumn number is q) simulation,
          }else if (isTRUE(c_power_m2<1-beta)){

            t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                           theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

            ## default when the re-estimation is not evaluable
            N_new_t=samsz
            if (isTRUE(as.numeric(t_new_0[2])<=1)){
              N_new_t=samsz
            }else if (isTRUE(as.numeric(t_new_0[2])>1)){
              T_new_0=as.numeric(t_new_0[[1]])

              N_new_t_0=max(ceiling(n_control_inter*(T_new_0/T_inter)))
              N_new_t_1=max(samsz,N_new_t_0)
              N_new_t= min(N_new_t_1,N_max)
            }

            if(isTRUE(N_new_t==samsz)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              samsz_simu_control=N[last_vt]
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)

            }else  if(isTRUE(N_new_t>samsz)){

              ####################################################################################
              ###  Step 4.2.1.3.1: Generate new random samples according to the new sample size ##
              ####################################################################################

              ####################################################################################
              ###  Step 4.2.1.3.1: Generate new random samples according to the new sample size ##
              ####################################################################################

              x0_new=x0[1:N_new_t]
              x_new=x[c(mu_keep),1:N_new_t]
              #x0 is the random sample at the first interim analysis for the control

              if (isTRUE(length(mu_keep)==1)){

                if (direction==1){

                  theta_hat_m2_new=mean(x_new)-mean(x0_new)

                }else if (direction==2){

                  theta_hat_m2_new=mean(x0_new)-mean(x_new)

                }
                t_hat_m2_new= N_new_t/(0.01+mean(x_new)+mean(x0_new))
                theta_hat_m2_sd_new=1/sqrt(t_hat_m2_new)

                wald_m2_new=theta_hat_m2_new*sqrt(t_hat_m2_new)

                wald_max_m2_new=max(wald_m2_new)

              }else if (isTRUE(length(mu_keep)>1)){

                t_hat_m2_new=NULL
                theta_hat_m2_new=NULL
                theta_hat_m2_sd_new=NULL
                wald_m2_new=NULL


                if (direction==1){
                  for (i in 1:length(mu_keep)){

                    theta_hat_m2_new[i]=mean(x_new[i,])-mean(x0_new)
                    t_hat_m2_new[i]= N_new_t/(0.01+mean(x_new[i])+mean(x0_new))
                    theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                    wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                  }

                }else if (direction==2){

                  for (i in 1:length(mu_keep)){

                    theta_hat_m2_new[i]=mean(x0_new)-mean(x_new[i,])
                    t_hat_m2_new[i]= N_new_t/(0.01+mean(x_new[i])+mean(x0_new))
                    theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                    wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                  }

                }
                wald_max_m2_new=max(wald_m2_new)
              }

              t_hat_m2_new_max=max(t_hat_m2_new)

              t_hat_m2_new_min=min(t_hat_m2_new)

              T_new=(t_hat_m2_new_min/t_hat_m2_new_max)*t_hat_m2_new

              T_inter=(N[length(s)-1]/N_new_t)*T_new

              T_hat_m2=(N[length(s)]/N_new_t)*T_new
              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=2)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=2)

              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=2)

              c_bry_new=NULL

              for (i in 1:length(theta_hat_inter_m2)){
                c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
              }

              test_new=wald_m2_new-max(c_bry_new)
              ###############################################################
              #### Step 4.2.1.3.2. Calculate Wald statistics with the added random samples ##
              ##############################################################

              samsz_simu_control=N_new_t
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N_new_t-N[length(s)-1])*(length(mu_keep)+1)

              if (isTRUE(max(test_new)>0)){

                rej=1

                stg=length(s)

              }else{
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0
              }

              arm_rej_new=NULL

              for (i in 1:length(mu_keep)){

                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }

              }

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }
        }

      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej_nossr,rej,stg,samsz_simu_control,samsz_simu_tot,futi,arm_drop,arm_rej)

  })#(end of q loop)

  #parallel computing ends. All resuts are combined and send to simr
  ########################################
  rej_nossr = as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  samsz_simu_control=as.numeric(simr[,4])
  samsz_simu_tot=simr[,5]
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(lambda_t))]
  arm_rej=simr[,(7+length(lambda_t)):(6+2*length(lambda_t))]
  ######################################
  rej_nossr_rt=mean(rej_nossr)

  rej_rt=mean(rej)

  arm_rej_rt=NULL
  for (i in 1:length(lambda_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  arm_drop_rt=NULL
  for (i in 1:length(lambda_t)){
    arm_drop_rt[i]=mean(arm_drop[,i])
  }

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
  }

  futi_rt=mean(futi)

  arm_summary=cbind(lambda_t,lambda_c,arm_drop_rt,arm_rej_rt)
  colnames(arm_summary)=c(LABEL()$lambda_t1,LABEL()$lambda_c,"drop rate by arm","Rejection rate by arm")

  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)
  max_tot_samsz=N_max*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL
  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu_control[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))

  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= samsz_simu_control[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  Design=cbind(c(1:length(s)),s,c_bry,rej_time,samsz_icr_byLook)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look","Rate of sample size increase")

  out_1=c(rej_nossr_rt, rej_rt, samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,samsz,N_max,mean_control_samsz, mean_samsz_simu_control_samsz_icr)
  names(out_1)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "type I error","non-type I error rejection","futility rate","Planned control sample size",
                 "Maximum control sample size",
                 "mean control sample size","mean control sample size conditional on increase")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")

  return(out)

}


###########################################

#masd_rej_simu_adapt_com_poisson(lambda_c,lambda_t,direction,s,c_bry,theta_cut,beta,samsz,N_max,dist,sim_num)
# masd_rej_simu_adapt_com_poisson(lambda_t=c(0.8,0.9),lambda_c=1.1,direction=2,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.2,beta=0.1,samsz=100,N_max=200,dist=1,sim_num=1000)
# masd_rej_simu_adapt_com_poisson(lambda_t=c(0.8,0.9,1.1),lambda_c=1.1,direction=2,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=100,N_max=291,dist=1,sim_num=10000)
#masd_rej_simu_adapt_com_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=100,N_max=200,dist=1,sim_num=10000)
#masd_rej_simu_adapt_com_poisson(lambda_t=c(1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=100,N_max=200,dist=1,sim_num=10000)

#masd_rej_simu_adapt_com_normal(mu=c(0,0),s=c(1/2,1),c_bry=c(3.142524,2.222100),sigma=1,theta_cut=0.05,beta=0.1,samsz=100,N_max=450,sim_num=100000,dist=1)
################################################

##############################################

#' Multiple comparisons Group sequential design: Simulations, Poisson distribution, Operating characteristics
#'
#' @param lambda_t Vector of true mean event counts of the test arms, one per dose; must be > 0
#' @param lambda_c Mean event count of the control arm; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param theta_cut Lack-of-efficacy threshold on the effect-size scale: a test arm is
#'   dropped at the interim analysis when its observed effect size is below this value
#' @param samsz Sample size of the control arm used in the simulation; must be > 0
#' @param N_max Maximum sample size for the control arm allowed by the sample size
#'   re-estimation; must be >= samsz
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary, the rejection rate by look and the rate of sample size increase}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the test-arm mean, the control-arm mean, the drop rate by arm and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate without and with sample size re-estimation, the rate of sample size increase, the rejection rates conditional on meeting the sample size increase criteria, the simulated type I error rate, the non-type I error rejection rate, the futility rate, the planned and maximum control sample sizes, the mean simulated control sample size and the mean simulated control sample size conditional on an increase}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' masd_power_simu_adapt_poisson(lambda_t = c(1.1,1.1,1.1), lambda_c = 1.1, direction = 1,
#' s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, beta = 0.1, theta_cut = 0.01, samsz = 100,
#' N_max = 500, sim_num = 20)
#' }
#'
#' @export
masd_power_simu_adapt_poisson<-function(lambda_t,lambda_c,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=length(lambda_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out=masd_rej_simu_adapt_com_poisson(lambda_t,lambda_c,direction=direction,s=s,c_bry=c_bry,theta_cut=theta_cut,beta=beta,samsz=samsz,N_max=N_max,dist=1,sim_num=sim_num)
  #masd_rej_simu_adapt_com_normal(mu,s,c_bry,sigma,theta_cut,beta,samsz,N_max ,sim_num,dist)
  out=simu_out
  # names(out)=c("simulation summary","planned control","planned total")
  return(out)
}

##################################
#  masd_power_simu_adapt_poisson(lambda_t,lambda_c,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sigma,sim_num)
#masd_power_simu_adapt_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=10000)

#masd_power_simu_adapt_poisson(lambda_t=c(0.8,0.9,1.1),lambda_c=1.1,direction=2,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=10000)
##############################################
####################################################################################
## Simulations for Optimizing Operating characteristics - Poisson  distributions ###
####################################################################################


###############################
masd_rej_OC_simu_adapt_com_poisson<-function(lambda_t,lambda_c,direction,s,c_bry,theta_cut,beta,samsz,N_max,dist,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size
  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  comp=length(lambda_t)
  N=ceiling(samsz*s)
  N1=N[1]
  NF=samsz
  if (direction==1){
    mu=lambda_t-lambda_c
  }else if (direction==2){
    mu=lambda_c-lambda_t
  }
  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    samsz_simu_control=N[length(s)]; samsz_simu_tot=N[length(s)]*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      #the loop specifies what is generated for each simulation

      #######################
      x0=stats::rpois(N_max,lambda_c)
      #x0 is the random sample at the first interim analysis for the control
      x=matrix(nrow=comp,ncol=N_max)
      for (i in 1:comp){
        x[i,]=stats::rpois(N_max,lambda_t[i])
      }

      #x is the matrix of random samples for each dose
      #for each simulation.
      #x matirx does not include the control sample.
      #The row number is the same as the active dose
      lambda_c_hat=NULL

      for (j in 1:length(s)){

        lambda_c_hat[j]=mean(x0[1:N[j]])
      }

      lambda_t_hat_m1=matrix(nrow=comp,ncol=length(s))

      for (i in 1:comp ){
        for (j in 1:length(s)){

          lambda_t_hat_m1[i,j]=mean(x[i,1:N[j]])
        }
      }

      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      if (direction==2){
        for (i in 1:comp ){
          for (j in 1:length(s)){

            theta_hat_m1[i,j]=lambda_c_hat[j]-lambda_t_hat_m1[i,j]
          }
        }
      } else if (direction==1){

        for (i in 1:comp ){
          for (j in 1:length(s)){
            theta_hat_m1[i,j]=lambda_t_hat_m1[i,j]-lambda_c_hat[j]
          }
        }
      }
      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))

      for (i in 1:comp ){
        for (j in 1:length(s)){
          t_hat_m1[i,j]= N[j]/( 0.01+lambda_t_hat_m1[i,j]+lambda_c_hat[j])
          wald_m1[i,j]=theta_hat_m1[i,j]*sqrt(t_hat_m1[i,j])
          theta_hat_m1_sd[i,j]=1/sqrt(t_hat_m1[i,j])
        }
      }

      wald_max_m1=NULL
      for (j in 1:length(s) ){

        wald_max_m1[j]=max(wald_m1[,j])
      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop=rep(0,length(lambda_t))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }

      T_hat=t_hat_m1[,length(s) ]
      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################
      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry
      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################
      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }

      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling
        futi=0
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        if (futi==1) { # If TRUE, thiere is not satisfactory efficacy, stop trial at planned sample size

          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            stg=0}
          last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)
          mu_reduc=0
        }else if (futi==0){ #beginning of step 4.2

          mu_drop=match(arm_drop[arm_drop>0],arm_drop)

          if (length(mu_drop)==0){

            mu_keep=c(1:comp)

          }else if ( length(mu_drop)>0){

            mu_keep=c(1:comp)[-mu_drop]

          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_hat_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2
          n_control_inter=N[length(s)-1]

          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m1)),theta_hat_inter_m1=theta_hat_inter_m1,theta_hat_inter_m1_sd=theta_hat_inter_m1_sd,dist=dist)

          if (isTRUE(c_power_m2>= 1-beta)){
            #########################################
            ### Proceed to planned sample size  #####
            ### and conduct final testing       #####
            #########################################
            last_vt=length(s)

            wald_m2=wald_m1[mu_keep,]

            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last

            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }
            samsz_simu_control=N[last_vt]
            samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)   #For q-th (collumn number is q) simulation,
          }else if (isTRUE(c_power_m2<1-beta)){

            t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                           theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

            ## default when the re-estimation is not evaluable
            N_new_t=samsz
            if (isTRUE(as.numeric(t_new_0[2])<=1)){
              N_new_t=samsz
            }else if (isTRUE(as.numeric(t_new_0[2])>1)){
              T_new_0=as.numeric(t_new_0[[1]])

              N_new_t_0=max(ceiling(n_control_inter*(T_new_0/T_inter)))
              N_new_t_1=max(samsz,N_new_t_0)
              N_new_t= min(N_new_t_1,N_max)
            }

            if(isTRUE(N_new_t==samsz)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              samsz_simu_control=N[last_vt]
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)

            }else  if(isTRUE(N_new_t>samsz)){

              ####################################################################################
              ###  Step 4.2.1.3.1: Generate new random samples according to the new sample size ##
              ####################################################################################

              ####################################################################################
              ###  Step 4.2.1.3.1: Generate new random samples according to the new sample size ##
              ####################################################################################

              x0_new=x0[1:N_new_t]
              x_new=x[c(mu_keep),1:N_new_t]
              #x0 is the random sample at the first interim analysis for the control

              if (isTRUE(length(mu_keep)==1)){

                if (direction==1){

                  theta_hat_m2_new=mean(x_new)-mean(x0_new)

                }else if (direction==2){

                  theta_hat_m2_new=mean(x0_new)-mean(x_new)

                }
                t_hat_m2_new= N_new_t/(0.01+mean(x_new)+mean(x0_new))
                theta_hat_m2_sd_new=1/sqrt(t_hat_m2_new)

                wald_m2_new=theta_hat_m2_new*sqrt(t_hat_m2_new)

                wald_max_m2_new=max(wald_m2_new)

              }else if (isTRUE(length(mu_keep)>1)){

                t_hat_m2_new=NULL
                theta_hat_m2_new=NULL
                theta_hat_m2_sd_new=NULL
                wald_m2_new=NULL


                if (direction==1){
                  for (i in 1:length(mu_keep)){

                    theta_hat_m2_new[i]=mean(x_new[i,])-mean(x0_new)
                    t_hat_m2_new[i]= N_new_t/(0.01+mean(x_new[i])+mean(x0_new))
                    theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                    wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                  }

                }else if (direction==2){

                  for (i in 1:length(mu_keep)){

                    theta_hat_m2_new[i]=mean(x0_new)-mean(x_new[i,])
                    t_hat_m2_new[i]= N_new_t/(0.01+mean(x_new[i])+mean(x0_new))
                    theta_hat_m2_sd_new[i]=1/sqrt(t_hat_m2_new[i])
                    wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])

                  }

                }
                wald_max_m2_new=max(wald_m2_new)
              }

              t_hat_m2_new_max=max(t_hat_m2_new)
              t_hat_m2_new_min=min(t_hat_m2_new)

              T_new=(t_hat_m2_new_min/t_hat_m2_new_max)*t_hat_m2_new

              T_inter=(N[length(s)-1]/N_new_t)*T_new

              T_hat_m2=(N[length(s)]/N_new_t)*T_new
              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=2)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=2)

              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=2)

              c_bry_new=NULL

              for (i in 1:length(theta_hat_inter_m2)){
                c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
              }

              test_new=wald_m2_new-max(c_bry_new)
              ###############################################################
              #### Step 4.2.1.3.2. Calculate Wald statistics with the added random samples ##
              ##############################################################

              samsz_simu_control=N_new_t
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N_new_t-N[length(s)-1])*(length(mu_keep)+1)

              if (isTRUE(max(test_new)>0)){

                rej=1

                stg=length(s)

              }else{
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0
              }

              arm_rej_new=NULL

              for (i in 1:length(mu_keep)){

                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }

              }

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }
        }

      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej_nossr,rej,stg,samsz_simu_control,samsz_simu_tot,futi,arm_drop,arm_rej)

  })#(end of q loop)

  #parallel computing ends. All resuts are combined and send to simr
  ########################################
  rej_nossr = as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  samsz_simu_control=as.numeric(simr[,4])
  samsz_simu_tot=simr[,5]
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(lambda_t))]
  arm_rej=simr[,(7+length(lambda_t)):(6+2*length(lambda_t))]
  ######################################
  rej_nossr_rt=mean(rej_nossr)

  rej_rt=mean(rej)

  arm_rej_rt=NULL
  for (i in 1:length(lambda_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  futi_rt=mean(futi)

  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)
  max_tot_samsz=N_max*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL
  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu_control[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)

  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= samsz_simu_control[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  out=c(rej_nossr_rt, rej_rt, samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,samsz,N_max,mean_control_samsz, mean_samsz_simu_control_samsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "type I error","non-type I error rejection","futility rate","Planned control sample size",
                 "Maximum control sample size",
                 "mean control sample size","mean control sample size conditional on increase")

  return(out)

}

###########################################
#masd_rej_OC_simu_adapt_com_poisson(lambda_c,lambda_t,direction,s,c_bry,theta_cut,beta,samsz,N_max,dist,sim_num)
#masd_rej_OC_simu_adapt_com_poisson(lambda_t=c(0.8,0.9),lambda_c=1.1,direction=2,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.2,beta=0.1,samsz=100,N_max=200,dist=1,sim_num=1000)
#masd_rej_OC_simu_adapt_com_poisson(lambda_t=c(0.8,0.9,1),lambda_c=1.1,direction=2,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=291,N_max=291,dist=1,sim_num=1000)
#masd_rej_OC_simu_adapt_com_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=100,N_max=500,dist=1,sim_num=10000)
#  masd_rej_OC_simu_adapt_com_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.114805,2.909606,2.375684),theta_cut=0.002,beta=0.1,samsz=100,N_max=500,dist=1,sim_num=10000)

################################################
##############################################

masd_OC_power_simu_adapt_poisson<-function(lambda_t,lambda_c,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sigma,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=length(lambda_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)[[1]][1]))
  }
  simu_out=masd_rej_OC_simu_adapt_com_poisson(lambda_t,lambda_c,direction=direction,s=s,c_bry=c_bry,theta_cut=theta_cut,beta=beta,samsz=samsz,N_max=N_max,dist=1,sim_num=sim_num)
  #masd_rej_simu_adapt_com_normal(mu,s,c_bry,sigma,theta_cut,beta,samsz,N_max ,sim_num,dist)
  lambda_t_names=paste(paste(LABEL()$lambda_t,"[",c(1:length(lambda_t)),"]",sep=""))
  out= round(c(lambda_t,lambda_c,as.numeric(simu_out)),4)
  names(out)=c(lambda_t_names,LABEL()$lambda_c, names(simu_out))
  return(out)
}


#  masd_OC_power_simu_adapt_poisson(lambda_t,lambda_c,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sigma,sim_num)
#masd_OC_power_simu_adapt_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=10000)
#masd_OC_power_simu_adapt_poisson(lambda_t=c(0.8,0.9,1.1),lambda_c=1.1,direction=2,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=10000)
#masd_OC_power_simu_adapt_poisson(lambda_t=c(0.8,0.9,1.1),lambda_c=1.1,direction=2,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=10000)

#masd_OC_power_simu_adapt_poisson(lambda_t=c(1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sigma=1,sim_num=10000)

##############################################
###################
#' Multiple comparisons Group sequential design: Simulations, Poisson distribution, Optimizing adaptive sequential design
#'
#' @param groups Matrix of true mean event counts of the test arms: each row is one simulation scenario and each column one test arm/dose; must be > 0
#' @param lambda_c Mean event count of the control arm; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g. c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Two-sided type I error levels spent at each look of the schedule s1; must sum to alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error levels spent at each look of the schedule s2; must sum to alpha_2
#' @param theta_cut Lack-of-efficacy threshold on the effect-size scale: a test arm is
#'   dropped at the interim analysis when its observed effect size is below this value
#' @param samsz_1 Sample size of the control arm used in the simulation for the first sample size scenario; must be > 0
#' @param samsz_2 Sample size of the control arm used in the simulation for the second sample size scenario; must be > 0
#' @param N_max_1 Maximum sample size for the control arm allowed by the sample size re-estimation in the first N_max scenario; must be >= samsz_1
#' @param N_max_2 Maximum sample size for the control arm allowed by the sample size re-estimation in the second N_max scenario; must be >= samsz_2
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the O'Brien-Fleming group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the O'Brien-Fleming group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the alpha-spending group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the alpha-spending group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' masd_simu_OC_poisson(
#' groups = rbind(c(0.2,0.3,0.4),c(0.3,0.4,0.5),c(0.4,0.5,0.6),c(0.5,0.6,0.7)), lambda_c = 0.3,
#' direction = 1, s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05,
#' alpha_2_alps_s1 = c(0.001,0.001,0.048), alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01,
#' samsz_1 = 50, samsz_2 = 75, N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
masd_simu_OC_poisson<-function(groups,lambda_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)
  comp=length(groups[1,])


  dist=1
  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=dist)[[1]][1]))
  s2_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=dist)[[1]][1]))

  s1_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry (s=s1,alpha=alpha_alsp_s1,comp=comp,dist=dist)[[1]][1]))
  s2_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry (s=s2,alpha=alpha_alsp_s2,comp=comp,dist=dist)[[1]][1]))

  ##### Two-stage GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")

  ##################
  lambda_t_vec=groups
  rownames( lambda_t_vec)=NULL

  lambda_t_group_names=paste(paste(LABEL()$lambda_t_group,c(1:dim(lambda_t_vec)[1]),sep=" "))
  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####


  simu_out_OF_s1_samsz_1_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s1_samsz_1_N_max_1=rbind(simu_out_OF_s1_samsz_1_N_max_1,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_1)=lambda_t_group_names



  ###############  s2 #######


  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_1=rbind(simu_out_OF_s2_samsz_1_N_max_1,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_1)=lambda_t_group_names


  ## samsz_2###

  ##### s1 ####


  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_1=rbind(simu_out_OF_s1_samsz_2_N_max_1,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_1)=lambda_t_group_names



  ###############  s2 #######


  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_1=rbind(simu_out_OF_s2_samsz_2_N_max_1,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_2_N_max_1)=lambda_t_group_names
  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####


  simu_out_OF_s1_samsz_1_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s1_samsz_1_N_max_2=rbind(simu_out_OF_s1_samsz_1_N_max_2,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_2)=lambda_t_group_names



  ###############  s2 #######


  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_2=rbind(simu_out_OF_s2_samsz_1_N_max_2,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_2)=lambda_t_group_names


  ## samsz_2###

  ##### s1 ####


  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_2=rbind(simu_out_OF_s1_samsz_2_N_max_2,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_2)=lambda_t_group_names



  ###############  s2 #######


  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_2=rbind(simu_out_OF_s2_samsz_2_N_max_2,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_2_N_max_2)=lambda_t_group_names



  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################


  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####


  simu_out_alsp_s1_samsz_1_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s1_samsz_1_N_max_1=rbind(simu_out_alsp_s1_samsz_1_N_max_1,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                                bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_1)=lambda_t_group_names



  ###############  s2 #######


  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_1=rbind(simu_out_alsp_s2_samsz_1_N_max_1,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                                bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_1)=lambda_t_group_names


  ## samsz_2###

  ##### s1 ####


  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_1=rbind(simu_out_alsp_s1_samsz_2_N_max_1,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                                bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                                beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_1)=lambda_t_group_names



  ###############  s2 #######


  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_1=rbind(simu_out_alsp_s2_samsz_2_N_max_1,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                                bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                                beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_2_N_max_1)=lambda_t_group_names
  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####


  simu_out_alsp_s1_samsz_1_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s1_samsz_1_N_max_2=rbind(simu_out_alsp_s1_samsz_1_N_max_2,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                                bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_2)=lambda_t_group_names



  ###############  s2 #######


  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_2=rbind(simu_out_alsp_s2_samsz_1_N_max_2,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                                bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_2)=lambda_t_group_names


  ## samsz_2###

  ##### s1 ####


  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_2=rbind(simu_out_alsp_s1_samsz_2_N_max_2,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                                bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                                beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_2)=lambda_t_group_names



  ###############  s2 #######


  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_2=rbind(simu_out_alsp_s2_samsz_2_N_max_2,masd_OC_power_simu_adapt_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                                bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                                beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_2_N_max_2)=lambda_t_group_names




  ##############

  list_OF_1=list(s1_OF_design,
                 simu_out_OF_s1_samsz_1_N_max_1,
                 simu_out_OF_s1_samsz_2_N_max_1,
                 simu_out_OF_s1_samsz_1_N_max_2,
                 simu_out_OF_s1_samsz_2_N_max_2)
  names(list_OF_1)=c("Group sequential design",
                     "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                     "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                     "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                     "Simulation results: Information fractions 1, Sample size 2, N_max 2")
  list_OF_2=list(s2_OF_design,
                 simu_out_OF_s2_samsz_1_N_max_1,
                 simu_out_OF_s2_samsz_2_N_max_1,
                 simu_out_OF_s2_samsz_1_N_max_2,
                 simu_out_OF_s2_samsz_2_N_max_2)
  names(list_OF_2)=c("Group sequential design",
                     "Simulation results: Information fractions 2, Sample size 1, N_max 1",
                     "Simulation results: Information fractions 2, Sample size 2, N_max 1",
                     "Simulation results: Information fractions 2, Sample size 1, N_max 2",
                     "Simulation results: Information fractions 2, Sample size 2, N_max 2")


  list_alsp_1=list(s1_alsp_design,
                   simu_out_alsp_s1_samsz_1_N_max_1,
                   simu_out_alsp_s1_samsz_2_N_max_1,
                   simu_out_alsp_s1_samsz_1_N_max_2,
                   simu_out_alsp_s1_samsz_2_N_max_2)
  names(list_alsp_1)=c("Group sequential design",
                       "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                       "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                       "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                       "Simulation results: Information fractions 1, Sample size 2, N_max 2")

  list_alsp_2=list(s2_alsp_design,
                   simu_out_alsp_s2_samsz_1_N_max_1,
                   simu_out_alsp_s2_samsz_2_N_max_1,
                   simu_out_alsp_s2_samsz_1_N_max_2,
                   simu_out_alsp_s2_samsz_2_N_max_2)
  names(list_alsp_2)=c("Group sequential design",
                       "Simulation results: Information fractions 2, Sample size 1, N_max 1",
                       "Simulation results: Information fractions 2, Sample size 2, N_max 1",
                       "Simulation results: Information fractions 2, Sample size 1, N_max 2",
                       "Simulation results: Information fractions 2, Sample size 2, N_max 2")


  out=list(list_OF_1,list_OF_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}

#################
#masd_simu_OC_poisson(groups,lambda_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num)
#masd_simu_OC_poisson(groups=rbind(c(0.2,0.3,0.4),c(0.3,0.4,0.5),c(0.4,0.5,0.6),c(0.5,0.6,0.7)),lambda_c=0,direction=1,
#                           s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                            theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1, sim_num=100)

###################################################################
####  III.5  Simulations for Negative binomial distributions  #####
###################################################################

################################################
#######################################
msd_rej_rate_neg_binomial<-function(r_t,r_c,kappa,direction,s,c_bry,samsz,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size

  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################

  N=ceiling(samsz*s)

  N1=N[1]
  NF=samsz

  if (direction==1){
    mu=log(r_t/r_c)
  }else if (direction==2){
    mu=-log(r_t/r_c)
  }

  comp=length(r_t)

  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej=0; stg=0; samsz_simu_control=N[length(s)]
    samsz_simu_tot=N[length(s)]*(comp+1); arm_rej=rep(0,comp)
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      #the loop specifies what is generated for each simulation

      x_c=cbind(rep(1,NF),stats::rnbinom(n=NF, size=kappa, mu=r_c))

      x=matrix(nrow=NF,ncol=comp)
      for (i in 1:comp){
        x[,i]=stats::rnbinom(n=NF, size=kappa, mu=r_t[i])
      }


      x_t=cbind(rep(2,NF),x)

      ####  analysis parameters ###

      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      for (i in 1:length(s)){

        for (j in 1:comp){

          x_c_temp=x_c[1:N[i],]
          x_t_temp=x_t[1:N[i],c(1,j+1)]

          x_temp=rbind(x_c_temp,x_t_temp)

          negglm_temp=tryCatch(glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"), error=function(e) NULL)  ## a failed fit is reported as NA
          sum_temp=if (is.null(negglm_temp)) NULL else summary(negglm_temp)
          se_temp=if (is.null(sum_temp)) NA_real_ else sum_temp$coefficients[2,2]
          r_c_temp=mean(x_c_temp[,2])
          r_t_temp=mean(x_t_temp[,2])

          if (direction ==1){ #higher rate r is better
            theta_hat_m1[j,i]=log(r_t_temp/r_c_temp)
          } else if (direction ==2){ #lower rate r is better
            theta_hat_m1[j,i]=-log(r_t_temp/r_c_temp)
          }

          t_hat_m1[j,i]=1/se_temp^2
          wald_m1[j,i]=theta_hat_m1[j,i]*sqrt(t_hat_m1[j,i])
          theta_hat_m1_sd[j,i]=1/sqrt(t_hat_m1[j,i])
        }
      }

      wald_max_m1=NULL

      for (j in 1:length(s) ){
        wald_max_m1[j]=max(wald_m1[,j])
      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################

      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################
      if (isTRUE(max(test_max)>=0)){ ##null hypothesis rejected, no more random sampling

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max[test_max>0]
        stg=match(test_max_t[1],test_max)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else  if (isTRUE(max(test_max)<0)){   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################
            rej=0
            stg=0
            last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)

        }
      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej,stg,samsz_simu_control,samsz_simu_tot,arm_rej)

  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr

  rej = as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  samsz_simu_control=as.numeric(simr[,3])
  samsz_simu_tot=simr[,4]
  arm_rej=simr[,5:(4+length(r_t))]

  ######################################

  rej_rt=mean(rej)

  arm_rej_rt=NULL

  for (i in 1:length(r_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/dim(simr)[1]
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }

  arm_summary=cbind(r_t,r_c,arm_rej_rt)
  colnames(arm_summary)=c("r<sub>t1</sub>,...","r<sub>c</sub>","Rejection rate by arm")


  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)


###############
  Design=cbind(c(1:length(s)),s,c_bry,rej_time)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look")

  out_1=c( rej_rt,tyI_er_rt,pos_rej_rt,samsz,tot_samsz,mean_control_samsz,mean_tot_samsz)
  names(out_1)=c("Rejection rate","typr I error","non-type I error rejection",
                 "Planned control sample size", "Planned total sample size",
                 "mean control sample size","mean total sample size")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")

  return(out)

}

###################
#msd_rej_rate_neg_binomial(r_t,r_c,kappa,direction,s,c_bry,samsz,sim_num)
#msd_rej_rate_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),samsz=100,sim_num=10000)
#msd_rej_rate_neg_binomial(r_t=c(1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298, 3.961146, 2.349682),samsz=100,sim_num=10000)
#msd_rej_rate_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.117,2.911,2.377),samsz=100,sim_num=10000)

################################################

##############################################

#' Multiple comparisons Group sequential design: Simulations, Negative binomial distribution, Type I error and power
#'
#' @param r_t Vector of true mean rates of the test arms, one per dose; must be > 0
#' @param r_c Mean rate of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param samsz Sample size of the control arm used in the simulation; must be > 0
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary and the rejection rate by look}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the test-arm mean rate, the control-arm mean rate and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate, the simulated type I error rate, the non-type I error rejection rate, the planned control and total sample sizes and the mean simulated control and total sample sizes}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' msd_power_simu_neg_binomial(r_t = c(1.1,1.1,1.1), r_c = 1.1, kappa = 1, direction = 1,
#' s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
msd_power_simu_neg_binomial<-function(r_t,r_c,kappa,direction,s,bry_type,alpha_2,samsz,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=length(r_t)
  if (bry_type==1) {
    #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
    c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))

  } else if (bry_type==2) {
    c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))
  }
  simu_out=msd_rej_rate_neg_binomial(r_t=r_t,r_c=r_c,kappa=kappa,direction=direction,s=s,c_bry=c_bry,samsz=samsz,sim_num=sim_num)
  # simu_out=msd_rej_rate_normal(mu=mu,s,c_bry,sigma,samsz,sim_num)
  #samsz_tot=samsz*(comp+1)
  out=simu_out
#  names(out)=c("simulation summary","planned control","planned total")
  return(out)
}

##################################
#msd_power_simu_neg_binomial(r_t,r_c,kappa,direction,s,bry_type,alpha_2,samsz,sim_num)
#msd_power_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=100)
#msd_power_simu_neg_binomial(r_t=c(1.2,1.3,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)

#msd_power_simu_neg_binomial(r_t=c(1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=10000)
##########################


  #############################
masd_rej_simu_adapt_com_neg_binomial<-function(r_t,r_c,kappa,direction,s,c_bry,theta_cut,beta,samsz,N_max,dist,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size

  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  N=ceiling(samsz*s)

  N1=N[1]
  NF=samsz

  if (direction==1){
    mu=log(r_t/r_c)
  }else if (direction==2){
    mu=-log(r_t/r_c)
  }
  comp=length(r_t)
  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################
  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    samsz_simu_control=N[length(s)]; samsz_simu_tot=N[length(s)]*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      #the loop specifies what is generated for each simulation

      ####################
      x_c=cbind(rep(1,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_c))

      x=matrix(nrow=N_max,ncol=comp)
      for (i in 1:comp){
        x[,i]=stats::rnbinom(n=N_max, size=kappa, mu=r_t[i])
      }

      x_t=cbind(rep(2,N_max),x)

      ####  analysis parameters ###

      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      for (i in 1:length(s)){

        for (j in 1:comp){

          x_c_temp=x_c[1:N[i],]
          x_t_temp=x_t[1:N[i],c(1,j+1)]

          x_temp=rbind(x_c_temp,x_t_temp)

          negglm_temp=tryCatch(glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"), error=function(e) NULL)  ## a failed fit is reported as NA
          sum_temp=if (is.null(negglm_temp)) NULL else summary(negglm_temp)
          se_temp=if (is.null(sum_temp)) NA_real_ else sum_temp$coefficients[2,2]
          r_c_temp=mean(x_c_temp[,2])
          r_t_temp=mean(x_t_temp[,2])

          if (direction ==1){ #higher rate r is better
            theta_hat_m1[j,i]=log(r_t_temp/r_c_temp)
          } else if (direction ==2){ #lower rate r is better
            theta_hat_m1[j,i]=-log(r_t_temp/r_c_temp)
          }

          t_hat_m1[j,i]=1/se_temp^2
          wald_m1[j,i]=theta_hat_m1[j,i]*sqrt(t_hat_m1[j,i])
          theta_hat_m1_sd[j,i]=1/sqrt(t_hat_m1[j,i])
        }
      }

      wald_max_m1=NULL

      for (j in 1:length(s) ){
        wald_max_m1[j]=max(wald_m1[,j])
      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop=rep(0,length(r_t))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################

      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry

      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################

      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1  ## rejection rate withoput SSR
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }


      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling

        futi=0
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        #########################################################
        ### Step 4.1. If max(theta<=theta_cut),                ##
        ###     Continue to planned sample size and stop trial  ##
        #########################################################

        if (futi==1) { # If TRUE, thiere is not satisfactory efficacy, stop trial at planned sample size

          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            stg=0}
          last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)

        }else if (futi==0){ #beginning of step 4.2

          ########################################################################
          ### Step 4.2. If max(theta>theta_cut), calculate conditional power    ##
          ###          and re-estiate sample size                               ##
          ########################################################################

          ########################################################################################
          ## Step 4.2.1. If futi==0. Check if there is need to drop some doses due to efficacy ###
          ### A dose is dropped for lack of efficacy if: theta_hat<theta cut                   ###
          ### or of min theta_hat is less than 1/3 of the max theta_hat                        ###
          ########################################################################################

          mu_drop=match(arm_drop[arm_drop>0],arm_drop)

          if (length(mu_drop)==0){

            mu_keep=c(1:comp)

          }else if ( length(mu_drop)>0){

            mu_keep=c(1:comp)[-mu_drop]

          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_hat_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2
          n_control_inter=N[length(s)-1]

          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m2)),theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          if (isTRUE(c_power_m2>= 1-beta)){
            #########################################
            ### Proceed to planned sample size  #####
            ### and conduct final testing       #####
            #########################################
            last_vt=length(s)

            wald_m2=wald_m1[mu_keep,]

            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last

            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }

            samsz_simu_control=N[last_vt]
            samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)   #For q-th (collumn number is q) simulation,

          }else if (isTRUE(c_power_m2<1-beta)){

              t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                             theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

              ## default when the re-estimation is not evaluable
              N_new_t=samsz
              if (isTRUE(as.numeric(t_new_0[2])<=1)){
                N_new_t=samsz
              }else if (isTRUE(as.numeric(t_new_0[2])>1)){
                T_new_0=as.numeric(t_new_0[[1]])

                N_new_t_0=max(ceiling(n_control_inter*(T_new_0/T_inter)))
                N_new_t_1=max(samsz,N_new_t_0)
                N_new_t= min(N_new_t_1,N_max)
              }

            if(isTRUE(N_new_t==samsz)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              samsz_simu_control=N[last_vt]
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)

            }else  if(isTRUE(N_new_t>samsz)){
              ####################################################################################
              ###  Step 4.2.1.3.1: Generate new random samples according to the new sample size ##
              ####################################################################################

              x_c_new=x_c[1:N_new_t,]
              x_t_new=x_t[1:N_new_t,c(1,(mu_keep+1))]

              r_c_new=mean(x_c_new[,2])
              #x0 is the random sample at the first interim analysis for the control
              t_hat_m2_new=NULL
              theta_hat_m2_new=NULL
              theta_hat_m2_sd_new=NULL
              wald_m2_new=NULL

              for (i in 1:length(mu_keep)){

                x_t_temp=x_t_new[,c(1,i+1)]

                x_temp=rbind(x_c_new,x_t_temp)

                negglm_temp=tryCatch(glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"), error=function(e) NULL)  ## a failed fit is reported as NA
                sum_temp=if (is.null(negglm_temp)) NULL else summary(negglm_temp)
                s_e_temp=if (is.null(sum_temp)) NA_real_ else as.numeric(sum_temp$coefficients[2,2])  ## NA se, not a crash

                r_t_temp=mean(x_t_temp[,2])

                if (direction==1){
                  theta_hat_m2_new[i]=log(r_t_temp/r_c_new)
                }else if (direction==2){
                  theta_hat_m2_new[i]=-log(r_t_temp/r_c_new)
                }

                t_hat_m2_new[i]=s_e_temp^(-2)
                wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])
                theta_hat_m2_sd_new[i]=s_e_temp
              }

              wald_max_m2_new=max(wald_m2_new)

              t_hat_m2_new_max=max(t_hat_m2_new)

              t_hat_m2_new_min=min(t_hat_m2_new)

              T_new=t_hat_m2_new

              T_inter=(N[length(s)-1]/N_new_t)*T_new
              T_hat_m2=(N[length(s)]/N_new_t)*T_new

              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=2)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=2)

              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=2)

              c_bry_new=NULL

              for (i in 1:length(theta_hat_inter_m2)){
                c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
              }

              mean_min_c_bry=mean(mean(c_bry_new),min(c_bry_new))

              test_new=wald_m2_new-max(mean_min_c_bry,c_bry_new)

              ###############################################################
              #### Step 4.2.1.3.2. Calculate Wald statistics with the added random samples ##
              ###############################################################

              samsz_simu_control=N_new_t
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N_new_t-N[length(s)-1])*(length(mu_keep)+1)

              if (isTRUE(max(test_new)>0)){
                rej=1
                stg=length(s)
              }else{
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0
              }

              arm_rej_new=NULL

              for (i in 1:length(mu_keep)){

                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }

              }

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }
        }

      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej_nossr,rej,stg,samsz_simu_control,samsz_simu_tot,futi,arm_drop,arm_rej)

  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr

   rej_nossr = as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  samsz_simu_control=as.numeric(simr[,4])
  samsz_simu_tot=simr[,5]
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(r_t))]
  arm_rej=simr[,(7+length(r_t)):(6+2*length(r_t))]

  ######################################

  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  arm_rej_rt=NULL

  for (i in 1:length(r_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  arm_drop_rt=NULL

  for (i in 1:length(r_t)){
    arm_drop_rt[i]=mean(arm_drop[,i])
  }

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/dim(simr)[1]
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }

  futi_rt=mean(futi)

  arm_summary=cbind(r_t,r_c,arm_drop_rt,arm_rej_rt)
  colnames(arm_summary)=c("r<sub>t1</sub>,...","r<sub>c</sub>","drop rate by arm","Rejection rate by arm")

  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)
  max_tot_samsz=N_max*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu_control[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= samsz_simu_control[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  Design=cbind(c(1:length(s)),s,c_bry,rej_time,samsz_icr_byLook)
  colnames(Design)=c("Look","Information fraction","Critical boundaries","rejection rate by look","Rate of sample size increase")

  out_1=c(rej_nossr_rt, rej_rt, samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,samsz,N_max,mean_control_samsz, mean_samsz_simu_control_samsz_icr)
  names(out_1)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "type I error","non-type I error rejection","futility rate","Planned control sample size",
                 "Maximum control sample size",
                 "mean control sample size","mean control sample size conditional on increase")

  out=list(Design,arm_summary,out_1)
  names(out)=c("Design","Summary by arm","simulation results")

  return(out)

}

  #######################################
  # masd_rej_simu_adapt_com_neg_binomial(r_t,r_c,kappa,direction,s,c_bry,theta_cut,beta,samsz,N_max,dist,sim_num)
 #masd_rej_simu_adapt_com_neg_binomial(r_t=c(0.8,0.9,1.1),r_c=1.1,kappa=1,direction=2,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=291,N_max=500,dist=1,sim_num=10000)
 # masd_rej_simu_adapt_com_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.117,2.911,2.377),theta_cut=0.01,beta=0.1,samsz=100,N_max=500,dist=1,sim_num=100000)
 # masd_rej_simu_adapt_com_neg_binomial(r_t=c(1.3,1.3,1.3),r_c=1.3,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.117,2.911,2.377),theta_cut=0.01,beta=0.1,samsz=100,N_max=500,dist=1,sim_num=100000)

 # masd_rej_simu_adapt_com_neg_binomial(r_t=c(0.6,0.6,0.6),r_c=0.6,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(4.117,2.911,2.377),theta_cut=0.01,beta=0.1,samsz=100,N_max=500,dist=1,sim_num=100000)
 # masd_rej_simu_adapt_com_neg_binomial(r_t=c(1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=291,N_max=400,dist=1,sim_num=100000)

 # masd_rej_simu_adapt_com_neg_binomial(r_t=c(1.3,1.3),r_c=1.3,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=291,N_max=400,dist=1,sim_num=100000)
 # masd_rej_simu_adapt_com_neg_binomial(r_t=c(0.6,0.6),r_c=0.6,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=291,N_max=400,dist=1,sim_num=100000)

  ################################################

  ##############################################

#' Multiple comparisons Group sequential design: Simulations, Negative binomial distribution, Operating characteristics
#'
#' @param r_t Vector of true mean rates of the test arms, one per dose; must be > 0
#' @param r_c Mean rate of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; when bry_type = 2 this is instead the
#'   vector of two-sided type I error levels spent at each look, which must sum to the total
#'   two-sided type I error
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param theta_cut Lack-of-efficacy threshold on the effect-size scale: a test arm is
#'   dropped at the interim analysis when its observed effect size is below this value
#' @param samsz Sample size of the control arm used in the simulation; must be > 0
#' @param N_max Maximum sample size for the control arm allowed by the sample size
#'   re-estimation; must be >= samsz
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Design}{Matrix with one row per look giving the look number, the information fraction, the critical boundary, the rejection rate by look and the rate of sample size increase}
#'   \item{Summary by arm}{Matrix with one row per test arm giving the test-arm mean rate, the control-arm mean rate, the drop rate by arm and the rejection rate by arm}
#'   \item{simulation results}{Named numeric vector holding the rejection rate without and with sample size re-estimation, the rate of sample size increase, the rejection rates conditional on meeting the sample size increase criteria, the simulated type I error rate, the non-type I error rejection rate, the futility rate, the planned and maximum control sample sizes, the mean simulated control sample size and the mean simulated control sample size conditional on an increase}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' masd_power_simu_adapt_neg_binomial(r_t = c(1.1,1.1,1.1), r_c = 1.1, kappa = 1, direction = 1,
#' s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, beta = 0.1, theta_cut = 0.01, samsz = 100,
#' N_max = 500, sim_num = 20)
#' }
#'
#' @export
  masd_power_simu_adapt_neg_binomial<-function(r_t,r_c,kappa,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sim_num){

    alpha=alpha_2/2
    comp=length(r_t)
    if (bry_type==1) {
      #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
      c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))

    } else if (bry_type==2) {
      c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))
    }
    simu_out=masd_rej_simu_adapt_com_neg_binomial(r_t,r_c,kappa,direction,s=s,c_bry=c_bry,theta_cut=theta_cut,beta=beta,samsz=samsz,N_max=N_max,dist=1,sim_num=sim_num)
    out=simu_out

    return(out)
  }

  ##################################
  # masd_power_simu_adapt_neg_binomial(r_t,r_c,kappa,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sim_num)
#masd_power_simu_adapt_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=10000)
#masd_power_simu_adapt_neg_binomial(r_t=c(1.2,1.3,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=10000)

 # masd_power_simu_adapt_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=100000)

  ##############################################
  ##############################################
  #############################################################################################
  ## Simulations for Optimizing Operating characteristics - Negative binomial distributions ###
  #############################################################################################
  #####################################

masd_rej_OC_simu_adapt_com_neg_binomial<-function(r_t,r_c,kappa,direction,s,c_bry,theta_cut,beta,samsz,N_max,dist,sim_num){
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size

  ################################################
  ##                                         #####
  ##              Preparations               #####
  ################################################
  N=ceiling(samsz*s)

  N1=N[1]
  NF=samsz

  if (direction==1){
    mu=log(r_t/r_c)
  }else if (direction==2){
    mu=-log(r_t/r_c)
  }
  comp=length(r_t)
  ############################################################
  ##                                                        ##
  ##   Simulation loops using parallel computing            ##
  ############################################################

  simr <- simuloop(sim_num, "MASD", function() {
    ## ---- non-evaluable default (added 2026-10-31) ----------------------
    ## A replicate whose statistics cannot be evaluated is recorded as
    ## 'no rejection, no sample size increase, finish at the planned size'.
    ## Without these defaults an NA decision left the output variables
    ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0; stg=0
    samsz_simu_control=N[length(s)]; samsz_simu_tot=N[length(s)]*(comp+1); futi=0
    arm_drop=rep(0,comp); arm_rej=rep(0,comp)
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      #the loop specifies what is generated for each simulation

      ####################
      x_c=cbind(rep(1,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_c))

      x=matrix(nrow=N_max,ncol=comp)
      for (i in 1:comp){
        x[,i]=stats::rnbinom(n=N_max, size=kappa, mu=r_t[i])
      }

      x_t=cbind(rep(2,N_max),x)

      ####  analysis parameters ###

      t_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1=matrix(nrow=comp,ncol=length(s))
      theta_hat_m1_sd=matrix(nrow=comp,ncol=length(s))
      wald_m1=matrix(nrow=comp,ncol=length(s))
      for (i in 1:length(s)){

        for (j in 1:comp){

          x_c_temp=x_c[1:N[i],]
          x_t_temp=x_t[1:N[i],c(1,j+1)]

          x_temp=rbind(x_c_temp,x_t_temp)

          negglm_temp=tryCatch(glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"), error=function(e) NULL)  ## a failed fit is reported as NA
          sum_temp=if (is.null(negglm_temp)) NULL else summary(negglm_temp)
          se_temp=if (is.null(sum_temp)) NA_real_ else sum_temp$coefficients[2,2]
          r_c_temp=mean(x_c_temp[,2])
          r_t_temp=mean(x_t_temp[,2])

          if (direction ==1){ #higher rate r is better
            theta_hat_m1[j,i]=log(r_t_temp/r_c_temp)
          } else if (direction ==2){ #lower rate r is better
            theta_hat_m1[j,i]=-log(r_t_temp/r_c_temp)
          }

          t_hat_m1[j,i]=1/se_temp^2
          wald_m1[j,i]=theta_hat_m1[j,i]*sqrt(t_hat_m1[j,i])
          theta_hat_m1_sd[j,i]=1/sqrt(t_hat_m1[j,i])
        }
      }
      wald_max_m1=NULL

      for (j in 1:length(s) ){
        wald_max_m1[j]=max(wald_m1[,j])
      }

      inter_vt=length(s)-1
      wald_inter_m1=wald_m1[,inter_vt]
      theta_hat_inter_m1=theta_hat_m1[,inter_vt]
      theta_hat_inter_m1_sd=theta_hat_m1_sd[,inter_vt]

      arm_drop=rep(0,length(r_t))
      for (i in 1:comp){
        if (isTRUE(theta_hat_inter_m1[i]<theta_cut)){
          arm_drop[i]=1
        }else{
          arm_drop[i]=0
        }
      }

      ######################################################################
      ## Step 2: hypothesis testing at interim analyses up to length(s)-1 ##
      ######################################################################

      test_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        test_all[i,]=wald_m1[i,]-c_bry
      }
      rej_all=matrix(nrow=comp,ncol=length(s))
      for (i in 1:comp){
        for (j in 1:length(s)){
          if (isTRUE(test_all[i,j]>=0)){
            rej_all[i,j]=1
          }else if (isTRUE(test_all[i,j]<0)){
            rej_all[i,j]=0
          }
        }
      }
      test_max=wald_max_m1-c_bry

      test_max_1=test_max[1:length(s)-1]

      if (isTRUE(max( theta_hat_inter_m1)<theta_cut)) {futi=1}else {futi=0}

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      ############################################################################################
      ### step 3. If Wald statistics cross the boundary at up to length(s)-1, terminte trial ##
      ############################################################################################

      if (isTRUE(max(test_max)>=0)){
        rej_nossr=1  ## rejection rate withoput SSR
      }else if (isTRUE(max(test_max)<0)){
        rej_nossr=0
      }
      #############################
      if (isTRUE(max(test_max_1)>=0)){ ##null hypothesis rejected, no more random sampling

        futi=0
        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_max_t=test_max_1[test_max_1>0]
        stg=match(test_max_t[1],test_max_1)
        #stg[q] identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        samsz_simu_control=N[last_vt]
        samsz_simu_tot=N[last_vt]*(comp+1)

      }else if (isTRUE(max(test_max_1)<0)) {   #beginnin of step 4
        ############################################################################################
        ### Step 4: If Wald statistics does not cross the boundary at up to length(s)-1         ##
        ## Calculate conditional power. If conditional power > 1-beta, no sample size change      ##
        ## if CP < 1-beta, re-estimate sample size to provide CP of 1-beta                        ##
        ############################################################################################

        #########################################################
        ### Step 4.1. If max(theta<=theta_cut),                ##
        ###     Continue to planned sample size and stop trial  ##
        #########################################################

        if (futi==1) { # If TRUE, thiere is not satisfactory efficacy, stop trial at planned sample size

          if (isTRUE(test_max[length(s)]>0)){
            rej=1
            stg=length(s)
          } else {
            rej=0
            stg=0}
          last_vt=length(s)

          samsz_simu_control=N[last_vt]
          samsz_simu_tot=N[last_vt]*(comp+1)

        }else if (futi==0){ #beginning of step 4.2

          ########################################################################
          ### Step 4.2. If max(theta>theta_cut), calculate conditional power    ##
          ###          and re-estiate sample size                               ##
          ########################################################################

          ########################################################################################
          ## Step 4.2.1. If futi==0. Check if there is need to drop some doses due to efficacy ###
          ### A dose is dropped for lack of efficacy if: theta_hat<theta cut                   ###
          ### or of min theta_hat is less than 1/3 of the max theta_hat                        ###
          ########################################################################################

          mu_drop=match(arm_drop[arm_drop>0],arm_drop)

          if (length(mu_drop)==0){

            mu_keep=c(1:comp)

          }else if ( length(mu_drop)>0){

            mu_keep=c(1:comp)[-mu_drop]

          }

          theta=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2=theta_hat_inter_m1[mu_keep]
          theta_hat_inter_m2_sd=theta_hat_inter_m1_sd[mu_keep]
          t_hat_inter_m1=t_hat_m1[,inter_vt]
          t_hat_inter_m2=t_hat_inter_m1[mu_keep]
          T_inter=t_hat_inter_m2
          wald_inter_m2=theta_hat_inter_m2*sqrt(t_hat_inter_m2)
          score_inter_m2=theta_hat_inter_m2*t_hat_inter_m2
          n_control_inter=N[length(s)-1]

          c_power_m2=masd_c_power_com(s=s,c_bry=c_bry,theta=theta_hat_inter_m2,theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          c_alpha_m1=masd_c_power_com(s=s,c_bry=c_bry,theta=rep(0,length(theta_hat_inter_m2)),theta_hat_inter_m1=theta_hat_inter_m2,theta_hat_inter_m1_sd=theta_hat_inter_m2_sd,dist=dist)
          if (isTRUE(c_power_m2>= 1-beta)){
            #########################################
            ### Proceed to planned sample size  #####
            ### and conduct final testing       #####
            #########################################
            last_vt=length(s)

            wald_m2=wald_m1[mu_keep,]

            if (isTRUE(length(mu_keep)==1)){
              wald_m2_last=wald_m2[last_vt]
              wald_max_m2=wald_m2_last

            }else if (isTRUE(length(mu_keep)>1)){
              wald_m2_last=wald_m2[,last_vt]
              wald_max_m2=max(wald_m2_last)
            }

            test_m2=wald_max_m2-c_bry[length(s)]
            arm_test_m2=wald_m2_last-c_bry[length(s)]
            if (isTRUE(test_m2>0)) {
              rej=1
              stg=length(s)
              last_vt=stg
            }else{
              rej=0
              stg=0
              last_vt=length(s)
            }

            samsz_simu_control=N[last_vt]
            samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)   #For q-th (collumn number is q) simulation,

          }else if (isTRUE(c_power_m2<1-beta)){

            t_new_0=masd_new_info_time_com(theta=theta_hat_inter_m2,s=s,c_alpha=c_alpha_m1,theta_hat_inter_m2=theta_hat_inter_m2,
                                           theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,beta=beta,dist=1)

            ## default when the re-estimation is not evaluable
            N_new_t=samsz
            if (isTRUE(as.numeric(t_new_0[2])<=1)){
              N_new_t=samsz
            }else if (isTRUE(as.numeric(t_new_0[2])>1)){
              T_new_0=as.numeric(t_new_0[[1]])

              N_new_t_0=max(ceiling(n_control_inter*(T_new_0/T_inter)))
              N_new_t_1=max(samsz,N_new_t_0)
              N_new_t= min(N_new_t_1,N_max)
            }

            if(isTRUE(N_new_t==samsz)){

              if (isTRUE(test_max[length(s)]>0)){
                rej=1
                stg=length(s)
              } else {
                rej=0
                stg=0}
              last_vt=length(s)

              samsz_simu_control=N[last_vt]
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N[last_vt]-N[length(s)-1])*(length(mu_keep)+1)

            }else  if(isTRUE(N_new_t>samsz)){
              ####################################################################################
              ###  Step 4.2.1.3.1: Generate new random samples according to the new sample size ##
              ####################################################################################

              x_c_new=x_c[1:N_new_t,]
              x_t_new=x_t[1:N_new_t,c(1,(mu_keep+1))]

              r_c_new=mean(x_c_new[,2])
              #x0 is the random sample at the first interim analysis for the control
              t_hat_m2_new=NULL
              theta_hat_m2_new=NULL
              theta_hat_m2_sd_new=NULL
              wald_m2_new=NULL

              for (i in 1:length(mu_keep)){

                x_t_temp=x_t_new[,c(1,i+1)]

                x_temp=rbind(x_c_new,x_t_temp)

                negglm_temp=tryCatch(glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"), error=function(e) NULL)  ## a failed fit is reported as NA
                sum_temp=if (is.null(negglm_temp)) NULL else summary(negglm_temp)
                s_e_temp=if (is.null(sum_temp)) NA_real_ else as.numeric(sum_temp$coefficients[2,2])  ## NA se, not a crash

                r_t_temp=mean(x_t_temp[,2])

                if (direction==1){
                  theta_hat_m2_new[i]=log(r_t_temp/r_c_new)
                }else if (direction==2){
                  theta_hat_m2_new[i]=-log(r_t_temp/r_c_new)
                }

                t_hat_m2_new[i]=s_e_temp^(-2)
                wald_m2_new[i]=theta_hat_m2_new[i]*sqrt(t_hat_m2_new[i])
                theta_hat_m2_sd_new[i]=s_e_temp
              }

              wald_max_m2_new=max(wald_m2_new)

              t_hat_m2_new_max=max(t_hat_m2_new)

              t_hat_m2_new_min=min(t_hat_m2_new)

              T_new=t_hat_m2_new

              T_inter=(N[length(s)-1]/N_new_t)*T_new
              T_hat_m2=(N[length(s)]/N_new_t)*T_new

              w_m1=comp*con_density(theta_hat_inter=theta_hat_inter_m1,theta_hat_inter_sd=theta_hat_inter_m1_sd,dist=2)
              w_m2=length(mu_keep)*con_density(theta_hat_inter=theta_hat_inter_m2,theta_hat_inter_sd=theta_hat_inter_m2_sd,dist=2)

              c_alpha_m2=c_alpha_m1*w_m1/w_m2

              e_bry_new= new_ebry_com(T=T_hat_m2,T_inter=T_inter,T_new=T_new,c_alpha=c_alpha_m2,s=s,
                                      theta_hat_inter_m2=theta_hat_inter_m2,
                                      theta_hat_inter_m2_sd=theta_hat_inter_m2_sd,dist=2)

              c_bry_new=NULL

              for (i in 1:length(theta_hat_inter_m2)){
                c_bry_new[i]=(e_bry_new-score_inter_m2[i])/sqrt(T_new[i]-T_inter[i])
              }

              mean_min_c_bry=mean(mean(c_bry_new),min(c_bry_new))

              test_new=wald_m2_new-max(mean_min_c_bry,c_bry_new)

              ###############################################################
              #### Step 4.2.1.3.2. Calculate Wald statistics with the added random samples ##
              ###############################################################

              samsz_simu_control=N_new_t
              samsz_simu_tot=N[length(s)-1]*(comp+1)+(N_new_t-N[length(s)-1])*(length(mu_keep)+1)

              if (isTRUE(max(test_new)>0)){
                rej=1
                stg=length(s)
              }else{
                #################################################################################
                ## Fail to reject the Null hypothesis if the critical boundary is not crossed  ###
                #################################################################################
                rej=0
                stg=0
              }

              arm_rej_new=NULL

              for (i in 1:length(mu_keep)){

                if (isTRUE(test_new[i]>=0)){
                  arm_rej_new[i]=1
                }else {
                  arm_rej_new[i]=0
                }

              }

              arm_rej_last=rep(0,comp)
              arm_rej_last[mu_keep]=arm_rej_new

              rej_all[,length(s)]=arm_rej_last
            }
          }
        }

      }

      if (isTRUE(stg==0)) {
        arm_rej=rep(0,comp)
      }else if (isTRUE(stg>0)) {
        arm_rej=rej_all[,stg]
      }
      out=c(rej_nossr,rej,stg,samsz_simu_control,samsz_simu_tot,futi,arm_drop,arm_rej)

  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr

  rej_nossr = as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  samsz_simu_control=as.numeric(simr[,4])
  samsz_simu_tot=simr[,5]
  futi=as.numeric(simr[,6])
  arm_drop=simr[,7:(6+length(r_t))]
  arm_rej=simr[,(7+length(r_t)):(6+2*length(r_t))]
  ######################################

  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  arm_rej_rt=NULL

  for (i in 1:length(r_t)){
    arm_rej_rt[i]=mean(arm_rej[,i])
  }

  futi_rt=mean(futi)

  mean_control_samsz=mean(samsz_simu_control)
  mean_tot_samsz=mean(samsz_simu_tot)
  tot_samsz=samsz*(comp+1)
  max_tot_samsz=N_max*(comp+1)

  arm_pos_rej=matrix(ncol=comp,nrow=dim(simr)[1])

  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){

      if (isTRUE((arm_rej[j,i]==1)&(mu[i]>0))){

        arm_pos_rej[j,i]=1

      }else {
        arm_pos_rej[j,i]=0
      }
    }
  }

  pos_rej=NULL

  for (i in 1: dim(simr)[1]){

    if (sum(arm_pos_rej[i,])>0){
      pos_rej[i]=1
    }else{
      pos_rej[i]=0
    }
  }
  pos_rej_rt=mean(pos_rej)

  arm_tyI_er=matrix(ncol=comp,nrow=dim(simr)[1])
  for (i in 1:comp){
    for (j in 1: dim(simr)[1]){
      if (isTRUE(arm_rej[j,i]==1&mu[i]<=0)){
        arm_tyI_er[j,i]=1
      }else {
        arm_tyI_er[j,i]=0
      }
    }
  }

  tyI_er=NULL

  for (i in 1: dim(simr)[1]){
    if ((sum(arm_tyI_er[i,])>0)&(pos_rej[i]==0)){
      tyI_er[i]=1
    }else{
      tyI_er[i]=0
    }
  }
  tyI_er_rt=mean(tyI_er)

  ######### sample size increase ########
  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu_control[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= samsz_simu_control[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  out=c(rej_nossr_rt, rej_rt, samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er_rt,pos_rej_rt,futi_rt,samsz,N_max,mean_control_samsz, mean_samsz_simu_control_samsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)",
                 "Rejection rate with SSR/(sample size incerase)",
                 "type I error","non-type I error rejection","futility rate","Planned control sample size",
                 "Maximum control sample size",
                 "mean control sample size","mean control sample size conditional on increase")
  return(out)

}

####################

# masd_rej_OC_simu_adapt_com_neg_binomial(r_t,r_c,kappa,direction,s,c_bry,theta_cut,beta,samsz,N_max,dist,sim_num)
# masd_rej_OC_simu_adapt_com_neg_binomial(r_t=c(0.8,0.9,1.1),r_c=1.1,kappa=1,direction=2,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=291,N_max=500,dist=1,sim_num=1000)
# masd_rej_OC_simu_adapt_com_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=291,N_max=400,dist=1,sim_num=10000)
# masd_rej_OC_simu_adapt_com_neg_binomial(r_t=c(1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.984298,3.961146,2.349682),theta_cut=0.002,beta=0.1,samsz=291,N_max=400,dist=1,sim_num=10000)

  ################################################
  ##############################################

  masd_OC_power_simu_adapt_neg_binomial<-function(r_t,r_c,kappa,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sim_num){

    alpha=alpha_2/2
    comp=length(r_t)
    if (bry_type==1) {
      #  design_out=MSD_OF_bundry(s,alpha,comp,dist=1)[1]
      c_bry=as.numeric(unlist(MSD_OF_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))

    } else if (bry_type==2) {
      c_bry=as.numeric(unlist(MSD_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=2)[[1]][1]))
    }
    simu_out=masd_rej_OC_simu_adapt_com_neg_binomial(r_t,r_c,kappa,direction,s=s,c_bry=c_bry,theta_cut=theta_cut,beta=beta,samsz=samsz,N_max=N_max,dist=1,sim_num=sim_num)
    r_t_names=paste(paste("r<sub>t</sub>","[",c(1:length(r_t)),"]",sep=""))
    out= round(c(r_t,r_c,as.numeric(simu_out)),4)
    names(out)=c(r_t_names,"r<sub>c</sub>", names(simu_out))

  #   out=simu_out

    return(out)
  }

  ##################################
  # masd_OC_power_simu_adapt_neg_binomial(r_t,r_c,kappa,direction,s,bry_type,alpha_2,beta,theta_cut,samsz,N_max,sim_num)
 # masd_OC_power_simu_adapt_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=10000)

#masd_OC_power_simu_adapt_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=10000)
#masd_OC_power_simu_adapt_neg_binomial(r_t=c(1.3,1.4,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),beta=0.1,theta_cut=0.01,samsz=100,N_max=500,sim_num=10000)

  ##############################################

  ###################
#' Multiple comparisons Group sequential design: Simulations, Negative binomial distribution, Optimizing adaptive sequential design
#'
#' @param groups Matrix of true mean rates of the test arms: each row is one simulation scenario and each column one test arm/dose; must be > 0
#' @param r_c Mean rate of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g. c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Two-sided type I error levels spent at each look of the schedule s1; must sum to alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error levels spent at each look of the schedule s2; must sum to alpha_2
#' @param theta_cut Lack-of-efficacy threshold on the effect-size scale: a test arm is
#'   dropped at the interim analysis when its observed effect size is below this value
#' @param samsz_1 Sample size of the control arm used in the simulation for the first sample size scenario; must be > 0
#' @param samsz_2 Sample size of the control arm used in the simulation for the second sample size scenario; must be > 0
#' @param N_max_1 Maximum sample size for the control arm allowed by the sample size re-estimation in the first N_max scenario; must be >= samsz_1
#' @param N_max_2 Maximum sample size for the control arm allowed by the sample size re-estimation in the second N_max scenario; must be >= samsz_2
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the
#'   sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the O'Brien-Fleming group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the O'Brien-Fleming group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the alpha-spending group sequential design for the schedule s1 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the alpha-spending group sequential design for the schedule s2 and the simulation results for the four combinations of the two planned control sample sizes and the two maximum control sample sizes}
#' }
#'
#' @references P. Gao & Y. Li (2024) Adaptive Multiple Comparison Sequential Design (AMCSD) for clinical trials, Journal of Biopharmaceutical Statistics, 34:3, 424-440, DOI:10.1080/10543406.2023.2233590. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2023.2233590}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' masd_simu_OC_neg_binomial(
#' groups = rbind(c(0.2,0.3,0.4),c(0.3,0.4,0.5),c(0.4,0.5,0.6),c(0.5,0.6,0.7)), r_c = 0.35,
#' kappa = 1, direction = 1, s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05,
#' alpha_2_alps_s1 = c(0.001,0.001,0.048), alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01,
#' samsz_1 = 50, samsz_2 = 75, N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
  masd_simu_OC_neg_binomial<-function(groups,r_c,kappa,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
    # This simulation function checks the actua power and sample size for different scenerios##
    # alpha_s is two sided type I error rate

    alpha=alpha_2/2
    alpha_alsp_s1=alpha_2_alps_s1/2
    alpha_alsp_s2=alpha_2_alps_s2/2

    a=stats::qnorm(1-alpha)
    b=stats::qnorm(1-beta)
    comp=length(groups[1,])


    dist=1
    ### critical boundaries for s1 ############

    s1_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=dist)[[1]][1]))
    s2_OF_cbry=as.numeric(unlist(MSD_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=dist)[[1]][1]))

    s1_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry (s=s1,alpha=alpha_alsp_s1,comp=comp,dist=dist)[[1]][1]))
    s2_alsp_cbry=as.numeric(unlist(MSD_alsp_bundry (s=s2,alpha=alpha_alsp_s2,comp=comp,dist=dist)[[1]][1]))

    ##### Two-stage GSD designs #####
    s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
    names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

    s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
    names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

    s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
    names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

    s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
    names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")

    ###################
    r_t_vec=groups
    rownames( r_t_vec)=NULL

    r_t_group_names=paste(paste("r<sub>t</sub> group",c(1:dim(r_t_vec)[1]),sep=" "))
    #########################################
    #### OC simulations for OF boundaries ###
    #########################################

    ##############
    ## N_max_1 ###
    ##############

    ## samsz_1###

    ##### s1 ####

    simu_out_OF_s1_samsz_1_N_max_1=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_OF_s1_samsz_1_N_max_1=rbind(simu_out_OF_s1_samsz_1_N_max_1,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                   bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                   beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_OF_s1_samsz_1_N_max_1)=r_t_group_names



    ###############  s2 #######

    simu_out_OF_s2_samsz_1_N_max_1=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_OF_s2_samsz_1_N_max_1=rbind(simu_out_OF_s2_samsz_1_N_max_1,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                   bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                   beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_OF_s2_samsz_1_N_max_1)=r_t_group_names


    ## samsz_2###


    ##### s1 ####

    simu_out_OF_s1_samsz_2_N_max_1=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_OF_s1_samsz_2_N_max_1=rbind(simu_out_OF_s1_samsz_2_N_max_1,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                   bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                                   beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_OF_s1_samsz_2_N_max_1)=r_t_group_names



    ###############  s2 #######

    simu_out_OF_s2_samsz_2_N_max_1=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_OF_s2_samsz_2_N_max_1=rbind(simu_out_OF_s2_samsz_2_N_max_1,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                   bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                                   beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_OF_s2_samsz_2_N_max_1)=r_t_group_names

    ##############
    ## N_max_2 ###
    ##############


    ## samsz_1###

    ##### s1 ####

    simu_out_OF_s1_samsz_1_N_max_2=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_OF_s1_samsz_1_N_max_2=rbind(simu_out_OF_s1_samsz_1_N_max_2,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                   bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                   beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_OF_s1_samsz_1_N_max_2)=r_t_group_names



    ###############  s2 #######

    simu_out_OF_s2_samsz_1_N_max_2=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_OF_s2_samsz_1_N_max_2=rbind(simu_out_OF_s2_samsz_1_N_max_2,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                   bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                   beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_OF_s2_samsz_1_N_max_2)=r_t_group_names


    ## samsz_2###


    ##### s1 ####

    simu_out_OF_s1_samsz_2_N_max_2=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_OF_s1_samsz_2_N_max_2=rbind(simu_out_OF_s1_samsz_2_N_max_2,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                   bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                                   beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_OF_s1_samsz_2_N_max_2)=r_t_group_names



    ###############  s2 #######

    simu_out_OF_s2_samsz_2_N_max_2=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_OF_s2_samsz_2_N_max_2=rbind(simu_out_OF_s2_samsz_2_N_max_2,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                   bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                                   beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_OF_s2_samsz_2_N_max_2)=r_t_group_names


    ###########################################
    #### OC simulations for alsp boundaries ###
    ###########################################


    ##############
    ## N_max_1 ###
    ##############

    ## samsz_1###

    ##### s1 ####

    simu_out_alsp_s1_samsz_1_N_max_1=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_alsp_s1_samsz_1_N_max_1=rbind(simu_out_alsp_s1_samsz_1_N_max_1,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                       bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                       beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_alsp_s1_samsz_1_N_max_1)=r_t_group_names



    ###############  s2 #######

    simu_out_alsp_s2_samsz_1_N_max_1=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_alsp_s2_samsz_1_N_max_1=rbind(simu_out_alsp_s2_samsz_1_N_max_1,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                       bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                       beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_alsp_s2_samsz_1_N_max_1)=r_t_group_names


    ## samsz_2###


    ##### s1 ####

    simu_out_alsp_s1_samsz_2_N_max_1=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_alsp_s1_samsz_2_N_max_1=rbind(simu_out_alsp_s1_samsz_2_N_max_1,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                       bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                                       beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_alsp_s1_samsz_2_N_max_1)=r_t_group_names



    ###############  s2 #######

    simu_out_alsp_s2_samsz_2_N_max_1=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_alsp_s2_samsz_2_N_max_1=rbind(simu_out_alsp_s2_samsz_2_N_max_1,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                       bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                                       beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_alsp_s2_samsz_2_N_max_1)=r_t_group_names

    ##############
    ## N_max_2 ###
    ##############


    ## samsz_1###

    ##### s1 ####

    simu_out_alsp_s1_samsz_1_N_max_2=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_alsp_s1_samsz_1_N_max_2=rbind(simu_out_alsp_s1_samsz_1_N_max_2,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                       bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                       beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_alsp_s1_samsz_1_N_max_2)=r_t_group_names



    ###############  s2 #######

    simu_out_alsp_s2_samsz_1_N_max_2=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_alsp_s2_samsz_1_N_max_2=rbind(simu_out_alsp_s2_samsz_1_N_max_2,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                       bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                       beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_alsp_s2_samsz_1_N_max_2)=r_t_group_names


    ## samsz_2###


    ##### s1 ####

    simu_out_alsp_s1_samsz_2_N_max_2=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_alsp_s1_samsz_2_N_max_2=rbind(simu_out_alsp_s1_samsz_2_N_max_2,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                       bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                                       beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_alsp_s1_samsz_2_N_max_2)=r_t_group_names



    ###############  s2 #######

    simu_out_alsp_s2_samsz_2_N_max_2=NULL
    for (i in 1:dim(r_t_vec)[1]){
      simu_out_alsp_s2_samsz_2_N_max_2=rbind(simu_out_alsp_s2_samsz_2_N_max_2,masd_OC_power_simu_adapt_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                       bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                                       beta=beta,sim_num=sim_num))
    }
    rownames(simu_out_alsp_s2_samsz_2_N_max_2)=r_t_group_names



    ##############

    list_OF_1=list(s1_OF_design,
                   simu_out_OF_s1_samsz_1_N_max_1,
                   simu_out_OF_s1_samsz_2_N_max_1,
                   simu_out_OF_s1_samsz_1_N_max_2,
                   simu_out_OF_s1_samsz_2_N_max_2)
    names(list_OF_1)=c("Group sequential design",
                       "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                       "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                       "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                       "Simulation results: Information fractions 1, Sample size 2, N_max 2")
    list_OF_2=list(s2_OF_design,
                   simu_out_OF_s2_samsz_1_N_max_1,
                   simu_out_OF_s2_samsz_2_N_max_1,
                   simu_out_OF_s2_samsz_1_N_max_2,
                   simu_out_OF_s2_samsz_2_N_max_2)
    names(list_OF_2)=c("Group sequential design",
                       "Simulation results: Information fractions 2, Sample size 1, N_max 1",
                       "Simulation results: Information fractions 2, Sample size 2, N_max 1",
                       "Simulation results: Information fractions 2, Sample size 1, N_max 2",
                       "Simulation results: Information fractions 2, Sample size 2, N_max 2")


    list_alsp_1=list(s1_alsp_design,
                     simu_out_alsp_s1_samsz_1_N_max_1,
                     simu_out_alsp_s1_samsz_2_N_max_1,
                     simu_out_alsp_s1_samsz_1_N_max_2,
                     simu_out_alsp_s1_samsz_2_N_max_2)
    names(list_alsp_1)=c("Group sequential design",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 2")

    list_alsp_2=list(s2_alsp_design,
                     simu_out_alsp_s2_samsz_1_N_max_1,
                     simu_out_alsp_s2_samsz_2_N_max_1,
                     simu_out_alsp_s2_samsz_1_N_max_2,
                     simu_out_alsp_s2_samsz_2_N_max_2)
    names(list_alsp_2)=c("Group sequential design",
                         "Simulation results: Information fractions 2, Sample size 1, N_max 1",
                         "Simulation results: Information fractions 2, Sample size 2, N_max 1",
                         "Simulation results: Information fractions 2, Sample size 1, N_max 2",
                         "Simulation results: Information fractions 2, Sample size 2, N_max 2")


    out=list(list_OF_1,list_OF_2,
             list_alsp_1,list_alsp_2)



    names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
                 "O'Brien-Fleming boundaries, interim analysis schedule=s2",
                 "Alpha-spending boundaries, interim analysis schedule=s1",
                 "Alpha-spending boundaries, interim analysis schedule=s2"
    )
    return(out)
  }

  #################
  # masd_simu_OC_neg_binomial(r_t_vec_1,r_t_vec_2,r_c,kappa,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num)
 # masd_simu_OC_neg_binomial(groups=rbind(c(0.2,0.3,0.4),c(0.3,0.4,0.5),c(0.4,0.5,0.6),c(0.5,0.6,0.7)),r_c=0.1,kappa=1,direction=1,
 #                                   s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
 #                                   theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1,sim_num=10)


  