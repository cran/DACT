# rm(list=ls())

#library(dplyr)
#library(bigstatsr)
#library(tidyverse)
#library("coin")

################################################

###########################################################
####  1  power  simulations for normal distributions  ########
###########################################################

#' Fixed sample designs: Normal distribution, Simulations
#'
#' @param mu_t Mean of the test arm; the null hypothesis is mu_t <= mu_c and the alternative
#'   hypothesis is mu_t > mu_c
#' @param mu_c Mean of the control arm
#' @param sigma_t Standard deviation of the test arm
#' @param sigma_c Standard deviation of the control arm
#' @param alpha_2 Type I error for the two-sided test
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param samsz Sample size of the control arm used in the simulation
#' @param sim_num Number of simulation replicates
#'
#' @return A named numeric value with the following component:
#' \describe{
#'   \item{rejection rate}{Simulated rejection rate of the two-sided fixed-sample Wald test of mu_t = mu_c; this is the simulated power, and the type I error rate when mu_t = mu_c}
#' }
#'
#' @examples
#' fixed_rej_rate_normal(mu_t = 0, mu_c = 0, sigma_t = 1, sigma_c = 1, alpha_2 = 0.05,
#' rand_ratio = 1, samsz = 82, sim_num = 20)
#'
#' @export
fixed_rej_rate_normal<-function(mu_t,mu_c,sigma_t,sigma_c,alpha_2,rand_ratio,samsz,sim_num){
  #This function simulates the rejection rate for normal distributions
  #mu is the effect size
  #sim_num=simulation number
  #samsaz is the sample size
  #sigma is the common stdard deviation
  ## When mu=0, the rejection rate is the type I error
  ## When mu>0, the rejection rate is the power

  r=rand_ratio
  ######################################
  ###  Preparing the parameters ########
  ######################################
  t<- samsz*(sigma_t^2/r+sigma_c^2)^(-1)   # t is the vector of information times at each interim anbalysis

  a=stats::qnorm(1-alpha_2/2)
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "FIXED", function() {

    ########################################
    ### Step 1: Generate random samples ####
    ########################################
    x0=stats::rnorm(samsz,mu_c,sigma_c)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rnorm(r*samsz,mu_t,sigma_t)

    #####################################################################
    ### Step 2: calculate the Wald statistics for hypothesis testing ####
    #####################################################################

    Wald=(mean(x)-mean(x0))*sqrt(t)
    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    test=Wald-a
    #at each simulation, test compares the Wald statistics with the critical boundaries
    #any value of test>0 indicates that the wald statistics crossed the critical boundaries
    if ( test>0){
      rej=1
    }else{
      rej=0
    }

    #######################################################
    ### Step 4: output and finish the simjulation loop ####
    #######################################################

    out=rej
  })#(end of q loop)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr)

  rej_rt=mean(rej)
  #rej_rt is the rejection rate

  #############################
  ### Step 6: Final output ####
  #############################
  out=rej_rt
  names(out)="rejection rate"
  return(out)
}
#######################################
#fixed_rej_rate_normal(mu_t,mu_c,sigma_t,sigma_c,alpha_2,rand_ratio,samsz,sim_num)
#fixed_rej_rate_normal(mu_t=0,mu_c=0,sigma_t=1,sigma_c=1,alpha_2=0.05,rand_ratio=1,samsz=82,sim_num=100000)
#fixed_rej_rate_normal(mu_t=0.5,mu_c=0,sigma_t=1,sigma_c=1,alpha_2=0.05,rand_ratio=1,samsz=82,sim_num=100000)


################################################################################
###########################################################
#### 2  Power simulations for survival analysis    ########
###########################################################

#####################################################
#' Fixed sample designs: Survival analysis, Simulations
#'
#' @param HR Hazard ratio of the test arm vs. the control arm; the null hypothesis is HR >= 1 and the
#'   alternative hypothesis is HR < 1; must be > 0
#' @param evt_num Total number of events at the final analysis; the simulated data are truncated at the
#'   first observation at which this number of events is reached
#' @param alpha_2 Type I error for the two-sided test
#' @param sim_num Number of simulation replicates
#'
#' @return A named numeric value with the following component:
#' \describe{
#'   \item{rejection CPH}{Simulated rejection rate of the two-sided Wald test from the Cox proportional hazards model; this is the simulated power, and the type I error rate when HR = 1}
#' }
#'
#' @examples
#' fixed_rej_rate_surv_cox(HR = 0.8, evt_num = 845, alpha_2 = 0.05, sim_num = 20)
#'
#' @export
fixed_rej_rate_surv_cox<-function(HR,evt_num,alpha_2,sim_num){
  lambda=1
  surv_rate=0.5
  samp_multip=5

  a=stats::qnorm(1-alpha_2/2)

  ###############################################################
  ###############################################
  ############## generate the parameters ########
  ###############################################

  sampsz=ceiling(samp_multip*evt_num/surv_rate)

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "FIXED", function() {
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
    total = 0
    for (i in 1:sampsz){
      # CBQ modify 2024/4/2, performance optimization
      #evt_count[i]=sum(x0_status[1:i])+sum(x_status[1:i]) #number of events after 2*i observations
      evt_count[i] = total + x0_status[i] + x_status[i]
      total = evt_count[i]
    }

    #################################################################################################

    evt_pos=evt_count[evt_count>evt_num]#substring of  evt_count that are larger than evt_num
    evt_cut=match(evt_pos[1],evt_count) #the first incidence that evt_count is larger than evt_num

    #  vector of zeros, length of the vector is the same as the original smaple size
    time=c(time_x0[1:evt_cut],time_x[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
    event=c(x0_status[1:evt_cut],x_status[1:evt_cut]) #censoring status before total number of events exceed evt_num
    grp=c(grp_x0[1:evt_cut],grp_x[1:evt_cut]) #vector of zero, with length=2*evt_cut


    test= list(time=time,status=event,grp)#survival data
    ########################################################################################################
    ####################### genrate score function , Wald testr stratistics from Cox partial likelihood regression model #####
    ###################################################################################################################

    coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
    # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
    # test is the same as data=test
    theta_hat= coxcf$coefficients[1] #estimated log(HR)
    theta_hat_sd= coxcf$coefficients[3] #estimated standard error of log(HR)
    wald_cox=theta_hat/theta_hat_sd
    ########################################################################################################
    #######################   genrate score function , Wald testr stratistics from log-rank test ###########
    ########################################################################################################

    # logrk=survdiff(Surv(time,status) ~ grp, data = test, rho = 0)
    # logrk=survdiff(Surv(time,status) ~ grp, data = test)
    # logrank=logrank.test(time,event,grp)
    # kafang_sta=logrank$chisq#kafang test
    # kafang_sta=as.numeric((logrank$obs-logrank$exp)^2/logrank$var)[1]
    # t=(theta_hat_sd)^(-2)
    # wald_log_rank=sqrt(kafang_sta)

    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    test_cox=wald_cox-a

    # test compares the Wald statistics with the critical boundaries
    #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
    if (test_cox>0){
      rej_cox=1
    }else {
      rej_cox=0
    }

    ######################################################
    #   test_log_rank=wald_log_rank-a
    # test compares the Wald statistics with the critical boundaries
    #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
    #if ( test_log_rank>0){
    #   rej_log_rank=1
    #}else {
    #   rej_log_rank=0
    #}

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    #  out=c(rej_cox,rej_log_rank)
    # out=rej_cox
  })

  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_cox= as.numeric(simr)
  # rej_log_rank= as.numeric(simr)

  power_cox=mean(rej_cox)
  # power_log_rank=mean(rej_log_rank)

  ###################
  #############################
  ### Step 7: Final output ####
  #############################

  #out=list(power_cox,power_log_rank)
  out=power_cox
  #names(out)=c("rejection rate_cox","rejection log-rank")
  names(out)="rejection CPH"

  # out=list(c_bry,power_cox,power_log_rank,mean_bias, stats)
  # names(out)=c("c_bry","rejection rate_cox","rejection log-rank","mean_bias","stats")
  return(out)
}
#######################################
#fixed_rej_rate_surv_cox(HR,evt_num,alpha_2,sim_num)
#fixed_rej_rate_surv_cox(HR=0.8,evt_num=845,alpha_2=0.05,sim_num=1000)
#fixed_rej_rate_surv_cox(HR=1,evt_num=845,alpha_2=0.05,sim_num=10000)

##############################################################

##############################################################
#### 3  Power simulations for binary distributions    ########
##############################################################

################################################
#' Fixed sample designs: Binary distribution, Simulations
#'
#' @param p_t Event rate (probability of an event) of the test arm; the null hypothesis is p_t <= p_c
#'   and the alternative hypothesis is p_t > p_c; 0 < p_t < 1
#' @param p_c Event rate (probability of an event) of the control arm; 0 < p_c < 1
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param alpha_2 Type I error for the two-sided test
#' @param samsz Sample size of the control arm used in the simulation
#' @param sim_num Number of simulation replicates
#'
#' @return A named numeric value with the following component:
#' \describe{
#'   \item{rejection rate}{Simulated rejection rate of the two-sided fixed-sample Wald test of p_t = p_c; this is the simulated power, and the type I error rate when p_t = p_c}
#' }
#'
#' @examples
#' fixed_rej_rate_binary(p_t = 0.4, p_c = 0.2, rand_ratio = 1, alpha_2 = 0.05,
#' samsz = 106, sim_num = 20)
#'
#' @export
fixed_rej_rate_binary<-function(p_t,p_c,rand_ratio,alpha_2,samsz,sim_num){
  #This function simulates the rejection rate, given p_t, p_c,s,sim_num
  #sim_num=simulation number
  #samsz is the sample size
  r=rand_ratio
  a=stats::qnorm(1-alpha_2/2)
  ######################################
  ###  Preparing the parameters ########
  ######################################

  t=samsz/(p_t*(1-p_t)/r+p_c*(1-p_c))  # t is the vector of information times at each interim anbalysis

  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "FIXED", function() {
    ########################################
    ### Step 1: Generate random samples ####
    ########################################
    x0=stats::rbinom(samsz,1,p_c)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rbinom(r*samsz,1,p_t)

    #####################################################################
    ### Step 2: calculate the Wald statistics for hypothesis testing ####
    #####################################################################

    Wald=(mean(x)-mean(x0))*sqrt(t)

    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    test=Wald-a
    #at each simulation, test compares the Wald statistics with the critical boundaries
    #any value of test>0 indicates that the wald statistics crossed the critical boundaries
    if ( test>0){
      rej=1
    }else{
      rej=0
    }

    #######################################################
    ### Step 4: output and finish the simjulation loop ####
    #######################################################

    out=rej
  })#(end of q loop)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr)

  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  #############################
  ### Step 6: Final output ####
  #############################
  out=rej_rt
  names(out)="rejection rate"
  return(out)
}
#######################################
#fixed_rej_rate_binary(p_t,p_c,rand_ratio,alpha_2,samsz,sim_num)
#fixed_rej_rate_binary(p_t=0.4,p_c=0.2,rand_ratio=1,alpha_2=0.05,samsz=106,sim_num=10000)
#fixed_rej_rate_binary(p_t=0.4,p_c=0.2,rand_ratio=2,alpha_2=0.05,samsz=74,sim_num=10000)
#fixed_rej_rate_binary(p_t=0.45,p_c=0.25,rand_ratio=1,alpha_2=0.05,samsz=115,sim_num=10000)
#fixed_rej_rate_binary(p_t=0.5,p_c=0.5,rand_ratio=1,alpha_2=0.05,samsz=100,sim_num=10000)
#fixed_rej_rate_binary(p_t=0.7,p_c=0.2,rand_ratio=1,alpha_2=0.05,samsz=100,sim_num=100)
##############################################


##############################################################
####  4  Power simulations for Poisson  distributions    ########
##############################################################
################################################
#' Fixed sample designs: Poisson distribution, Simulations
#'
#' @param lambda_t Event rate of the test arm; the null hypothesis is lambda_t <= lambda_c and the
#'   alternative hypothesis is lambda_t > lambda_c; must be > 0
#' @param lambda_c Event rate of the control arm; must be > 0
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param alpha_2 Type I error for the two-sided test
#' @param samsz Sample size of the control arm used in the simulation
#' @param sim_num Number of simulation replicates
#'
#' @return A named numeric value with the following component:
#' \describe{
#'   \item{rejection rate}{Simulated rejection rate of the two-sided fixed-sample Wald test of lambda_t = lambda_c; this is the simulated power, and the type I error rate when lambda_t = lambda_c}
#' }
#'
#' @examples
#' fixed_rej_rate_poisson(lambda_t = 11, lambda_c = 10, rand_ratio = 2, alpha_2 = 0.05,
#' samsz = 163, sim_num = 20)
#'
#' @export
fixed_rej_rate_poisson<-function(lambda_t,lambda_c,rand_ratio,alpha_2,samsz,sim_num){
  #This function simulates the rejection rate, given p_t, p_c,s,sim_num
  #sim_num=simulation number
  #samsz is the sample size

  a=stats::qnorm(1-alpha_2/2)
  ######################################
  ###  Preparing the parameters ########
  ######################################
  r=rand_ratio
  t=samsz/(lambda_t/r+lambda_c)  # t is the vector of information times at each interim anbalysis

  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "FIXED", function() {
    ########################################
    ### Step 1: Generate random samples ####
    ########################################
    x0=stats::rpois(samsz,lambda_c)
    #x0 is the random sample at the first interim analysis for the control
    x=stats::rpois(r*samsz,lambda_t)

    #####################################################################
    ### Step 2: calculate the Wald statistics for hypothesis testing ####
    #####################################################################

    Wald=(mean(x)-mean(x0))*sqrt(t)

    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    test=Wald-a
    #at each simulation, test compares the Wald statistics with the critical boundaries
    #any value of test>0 indicates that the wald statistics crossed the critical boundaries
    if ( test>0){
      rej=1
    }else{
      rej=0
    }

    #######################################################
    ### Step 4: output and finish the simjulation loop ####
    #######################################################

    out=rej
  })#(end of q loop)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr)

  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  #############################
  ### Step 6: Final output ####
  #############################
  out=rej_rt
  names(out)="rejection rate"
  return(out)
}
##############################################
#fixed_rej_rate_poisson(lambda_t,lambda_c,rand_ratio,alpha_2,samsz,sim_num)
#fixed_rej_rate_poisson(lambda_t=11,lambda_c=10,rand_ratio=2,alpha_2=0.05,samsz=163,sim_num=10000)
#fixed_rej_rate_poisson(lambda_t=10,lambda_c=10,rand_ratio=1,alpha_2=0.05,samsz=221,sim_num=100000)
##############################################


####################################################################
###  5  Power simulations for negative binomial distributions    ###
####################################################################
#' Fixed sample designs: Negative binomial distribution, Simulations
#'
#' @param alpha_2 Type I error for the two-sided test
#' @param r_c Event rate of the control arm; must be > 0
#' @param r_t Event rate of the test arm; the null hypothesis is r_t >= r_c and the alternative
#'   hypothesis is r_t < r_c; must be > 0
#' @param rand_ratio Randomization ratio: number of subjects in the test arm vs. the number of subjects in the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial distribution; the negative binomial size
#'   parameter used in the simulation is 1/kappa, so kappa must be > 0 and a larger kappa means more overdispersion
#' @param samsz Sample size of the control arm used in the simulation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{rejection rate}{Simulated rejection rate of the two-sided Wald test from the negative binomial log-link GLM; this is the simulated power, and the type I error rate when r_t = r_c}
#'   \item{mean size}{Mean over the simulation replicates of the estimated negative binomial size (dispersion) parameter theta of the fitted GLM}
#' }
#'
#' @examples
#' fixed_rej_nbglm_simu(alpha_2 = 0.05, r_c = 1, r_t = 1*0.6, rand_ratio = 1, kappa = 1,
#' samsz = 188, sim_num = 20)
#'
#' @export
fixed_rej_nbglm_simu<-function(alpha_2,r_c,r_t,rand_ratio,kappa,samsz,sim_num){
  ### simulation with dispersion parameter ######
  ###############################################
  ############## generate the parameters ########
  ###############################################
  dispersion=1/kappa
  #sim_num=simulation number
  z_alpha=stats::qnorm(1-alpha_2/2)

  N=samsz

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "FIXED", function() {
    ########################################
    ### Step 1: Generate random samples ####
    ########################################

    x_c=cbind(rep(1,samsz),stats::rnbinom(n=samsz, size=dispersion, mu=r_c))
    x_t=cbind(rep(2,samsz),stats::rnbinom(n=rand_ratio*samsz, size=dispersion, mu=r_t))
    x=rbind(x_c,x_t)

    ##### negative binomial glm ####
    nbglm_f=glm.nb(as.numeric(x[,2])~as.factor(x[,1]),link="log")
    sum_f=summary(nbglm_f)
    s_e_f=as.numeric(sum_f$coefficients[2,2])
    mu_nbglm=-as.numeric(sum_f$coefficients[2,1])
    theta_hat_f=sum_f$theta
    Wald_nbglm=-as.numeric(sum_f$coefficients[2,3])
    if (Wald_nbglm>=z_alpha){rej_nbglm=1}else{rej_nbglm=0}
    ##### Final analysis ######

    out=c(rej_nbglm,theta_hat_f)
  }) #(end of simulation loops)
  #############termination timepoint####
  ##################################################
  ### Step 6: Generate summary statistics --- PZ####
  ##################################################

  rej_nbglm = as.numeric(simr[,1])
  theta_hat_f =as.numeric(simr[,2])
  ####### power ######

  rej_rt_nbglm=mean(rej_nbglm)
  mean_size=mean(theta_hat_f)
  #################################################################################################

  out=list(rej_rt_nbglm ,mean_size)
  names(out)=c("rejection rate","mean size")
  return(out)
}
#######################################
#fixed_rej_nbglm_simu(alpha_2,r_c,r_t,rand_ratio,kappa,samsz,sim_num)
#fixed_rej_nbglm_simu(alpha_2=0.05,r_c=1,r_t=1*0.6,rand_ratio=1,kappa=1,samsz=188,sim_num=10000)
#fixed_rej_nbglm_simu(alpha_2=0.05,r_c=1,r_t=1*0.6,rand_ratio=1,kappa=1.3,samsz=213,sim_num=10000)
#fixed_rej_nbglm_simu(alpha_2=0.05,r_c=1,r_t=1*0.6,rand_ratio=2,kappa=1.3,samsz=166,sim_num=10000)
#fixed_rej_nbglm_simu(alpha_2=0.05,r_c=1,r_t=1,rand_ratio=2,kappa=1.3,samsz=166,sim_num=10000)
