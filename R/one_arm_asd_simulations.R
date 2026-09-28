#rm(list = ls())

#library(clinfun)

#library(dplyr)
#library(mvtnorm)
#library(foreach)
#library(parallel)
#library(doParallel)

##########################################

##########################################################################
####              Two stage design  simulations                       ####
##########################################################################

#' Simulations: two stage fixed Simon's design
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param p_1   The assumed response rate.
#' @param n_1   Number of patients at first look.
#' @param r_1   If no more than r_1 responses are observed at the first look, then the trial would be stopped for futility.
#' @param n_2   Number of patients at second look.
#' @param r_2   If no more than r_2 responses are observed at the second look, then the null hypothesis will be rejected.
#' @param sim_num Purpose of simulation: To verify type I error control by setting p_1=p_0, To verify the parameters n_1,r_1,n_2,r_2
#' are correctly chosen, such that the trial has 1-beta power if the true response rate is p_1.
#'
#' @return Average sample size is EN(p_1).
#' @return Rate of termination is PET(p_1)
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references P. Gao, J.H. Ware, C. Mehta, (2008), Sample size re-estimation for adaptive sequential designs.
#' Journal of Biopharmaceutical Statistics, 18: 1184–1196, 2008.
#'
#' @examples
#' \donttest{
#' one_arm_rej_simon_binary(p_1=0.4,p_0=0.2,n_1=37,r_1=9,n_2=53,r_2=16,sim_num=20)
#' }
#'
#' @export
one_arm_rej_simon_binary<-function(p_1,p_0,n_1,r_1,n_2,r_2,sim_num){
  # p_1 is the event rate for the test arm
  # p_0 is the undesirable evnet rate
  #sim_num=simulation number
  #samsz is sample size
  #N_max is themaximum sample size

  ######################################
  ###  Preparing the parameters ########
  ######################################

  N=c(n_1,n_2)  # N is the vetoc of sample sizes at which interim analysis will be conducted

  # t=N/(p_1*(1-p_1))  # t is the vector of information times at each interim anbalysis
  # c_2=(r_2/n_tot-p_0)*sqrt(n_tot^3/(r_2*(n_tot-r_2)))
  # b=stats::qnorm(1-beta)
  #t is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "SINGLE", function() {
    ########################################
    ### Step 1: Generate random samples ####
    ########################################

    x=stats::rbinom(n_2,1,p_1)

    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    # test compares the Wald statistics with the critical boundaries
    #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
    if (sum(x[1:N[1]])<=r_1) { #futility termination
      futi=1
      rej=0
      samsz_simu=n_1
    }else if (sum(x[1:N[1]])>r_1){
      futi=0
      #################################################################################
      ###                   Step 4.2: Compare theta_hat with theta_cut             ####
      ###              if theta_hat>theta_cut, re-estimate sampe size              ####
      #################################################################################
      #sample size re-estimation procedure##
      samsz_simu=N[2]
      if (sum(x[1:N[2]])>r_2){ # the NULL hypothesis is rejected at the final analysis
        rej=1
      }else { # the NULL hypothesis is NOT rejected at the final analysis
        rej=0
      }
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out=c(rej,samsz_simu,futi)
  }) #(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  sampsz_simu=as.numeric(simr[,2])
  futi=as.numeric(simr[,3])
  ###########################

  power=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)

  ###################
  #############################
  ### Step 7: Final output ####
  #############################

  out=list(power,av_samsz,futi_rt)
  names(out)=c("Rate of rejection","Average sample size","Rate of early termination for futility")
  return(out)
}
#######################################
#one_arm_rej_simon_binary(p_1,p_0,n_1,r_1,n_2,r_2,sim_num)
#one_arm_rej_simon_binary(p_1=0.4,p_0=0.2,n_1=37,r_1=9,n_2=53,r_2=16,sim_num=100000)
#one_arm_rej_simon_binary(p_1=0.2,p_0=0.2,n_1=37,r_1=9,n_2=53,r_2=16,sim_num=100000)
#one_arm_rej_simon_binary(p_1=0.33,p_0=0.2,n_1=100,r_1=25,n_2=119,r_2=32,sim_num=100000)
#one_arm_rej_simon_binary(p_1=0.2,p_0=0.2,n_1=100,r_1=25,n_2=119,r_2=32,sim_num=100000)
############################################
#######################################

#' Simulations: Two stage fixed expanded Simon's design: either hybrid or mid-point design
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param p_1   The assumed response rate.
#' @param n_1   Number of patients at first look.
#' @param r_1_f   If no more than r_1 responses are observed at the first look, then the trial would be stopped for futility.
#' @param r_1_e   The trial would be stopped for superiority if at least r_1_e responses are observed at the first look.
#' If r_1_e is not used in the design, enter NA.
#' @param n_2   Number of patients at second look.
#' @param r_2   If no more than r_2 responses are observed at the second look, then the null hypothesis will be rejected.
#' @param sim_num Purpose of simulation: To verify type I error control by setting p_1=p_0,
#' To verify the parameters n_1,r_1_f,r_1_e,n_2,r_2 are correctly chosen, such that the rejection rate matches that
#' from the design table from either the hybrid design or the mid-point design.
#'
#' @return Average sample size is EN(p_1).
#' @return Rate of termination is PET(p_1)
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references P. Gao, J.H. Ware, C. Mehta, (2008), Sample size re-estimation for adaptive sequential designs.
#' Journal of Biopharmaceutical Statistics, 18: 1184–1196, 2008.
#'
#' @examples
#' \donttest{
#' one_arm_rej_2_stg_binary(p_1=0.2,p_0=0.2,n_1=100,r_1_f=25,r_1_e=33,
#' n_2=105,r_2=29,sim_num=20)
#' }
#'
#' @export
one_arm_rej_2_stg_binary<-function(p_1,p_0,n_1,r_1_f,r_1_e=NA,n_2,r_2,sim_num){
  # p_1 is the event rate for the test arm
  # p_0 is the undesirable evnet rate
  #sim_num=simulation number
  #samsz is sample size
  ######################################
  ###  Preparing the parameters ########
  ######################################
  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "SINGLE", function() {
    ########################################
    ### Step 1: Generate random samples ####
    ########################################

    x=stats::rbinom(n_2,1,p_1)

    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    # test compares the Wald statistics with the critical boundaries
    #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
    x_1=sum(x[1:n_1])
    x_2=sum(x[1:n_2])

    if (!is.na(r_1_e)){
      if (x_1<=r_1_f) { #futility termination
        futi=1
        rej=0
        samsz_simu=n_1
      }else if (x_1>r_1_e){
        futi=0
        rej=1
        samsz_simu=n_1
      }else if ((x_1>r_1_f)&(x_1<=r_1_e)){
        samsz_simu=n_2
        if(x_2>r_2){
          futi=0
          rej=1
        }else {
          futi=0
          rej=0
        }
      }
    }else if (is.na(r_1_e)){
      if (x_1<=r_1_f) { #futility termination
        futi=1
        rej=0
        samsz_simu=n_1
      }else if (x_1>r_1_f){
        samsz_simu=n_2
        if(x_2>r_2){
          futi=0
          rej=1
        }else {
          futi=0
          rej=0
        }
      }
    }

    #################################################################################
    ###                   Step 4.2: Compare theta_hat with theta_cut             ####
    ###              if theta_hat>theta_cut, re-estimate sampe size              ####
    #################################################################################
    #sample size re-estimation procedure##
    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    #    out=cbind(x_1,r_1,n_1,x_2,r_2_f,r_2_e,n_2,x_3,r_3,n_3,rej,futi,samsz_simu)
    #    out
    out=c(rej,samsz_simu,futi)
  }) #(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  sampsz_simu=as.numeric(simr[,2])
  futi=as.numeric(simr[,3])
  ###########################

  power=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)

  ###################
  #############################
  ### Step 7: Final output ####
  #############################

  out=list(power,av_samsz,futi_rt)
  names(out)=c("Rate of rejection","Average sample size","Rate of early termination for futility")
  return(out)
}

############################
#one_arm_rej_2_stg_binary(p_1,p_0,n_1,r_1_f,r_1_e=NA,n_2,r_2,sim_num)
#########  Minimax ###
#one_arm_rej_2_stg_binary(p_1=0.2,p_0=0.2,n_1=100,r_1_f=25,r_1_e=33,n_2=105,r_2=29,sim_num=100000)
#one_arm_rej_2_stg_binary(p_1=0.4,p_0=0.2,n_1=100,r_1_f=25,r_1_e=33,n_2=105,r_2=29,sim_num=100000)
#one_arm_rej_2_stg_binary(p_1=0.33,p_0=0.2,n_1=100,r_1_f=25,r_1_e=33,n_2=105,r_2=29,sim_num=100000)

##### Optimal ####
#one_arm_rej_2_stg_binary(p_1=0.2,p_0=0.2,n_1=50,r_1_f=11,r_1_e=NA,n_2=56,r_2=17,sim_num=100000)
#one_arm_rej_2_stg_binary(p_1=0.4,p_0=0.2,n_1=50,r_1_f=11,r_1_e=NA,n_2=56,r_2=17,sim_num=100000)
#one_arm_rej_2_stg_binary(p_1=0.33,p_0=0.2,n_1=50,r_1_f=11,r_1_e=NA,n_2=56,r_2=17,sim_num=100000)

##### Average ####
#one_arm_rej_2_stg_binary(p_1=0.2,p_0=0.2,n_1=75,r_1_f=18,r_1_e=26,n_2=80,r_2=23,sim_num=100000)
#one_arm_rej_2_stg_binary(p_1=0.4,p_0=0.2,n_1=75,r_1_f=18,r_1_e=26,n_2=80,r_2=23,sim_num=100000)
#one_arm_rej_2_stg_binary(p_1=0.33,p_0=0.2,n_1=75,r_1_f=18,r_1_e=26,n_2=80,r_2=23,sim_num=100000)


############################

#' Simulations: adaptive two stage design
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param p_1   The assumed response rate.
#' @param n_1   Number of patients at first look.
#' @param r_1   If no more than r_1 responses are observed at the first look, then the trial would be stopped for futility.
#' @param r_1_e   The trial would be stopped for superiority if at least r_1_e responses are observed at the first look.
#' If r_1_e is not used in the design, enter NA.
#' @param n_2   Number of patients at second look.
#' @param r_2   If no more than r_2 responses are observed at the second look, then the null hypothesis will be rejected.
#' @param N_max Maximum sample size for the trial.
#' @param beta  Type II error. The target power is 1-beta.
#' @param sim_num Purpose of simulation: To verify type I error control by setting p_1=p_0,
#' To evaluate the operating characteristics of the adaptive design and to determine if continuity correction is necessary.
#'
#' @return Average sample size is EN(p_1).
#' @return Rate of termination is PET(p_1)
#' @return If rejection rate without continuity correction does not exceed alpha, the continuity correction is not needed for this group of
#' parameters p_0,n_1,r_1,r_1_e,n_2,r_2. Otherwise, continuity correction should be applied.
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references P. Gao, J.H. Ware, C. Mehta, (2008), Sample size re-estimation for adaptive sequential designs.
#' Journal of Biopharmaceutical Statistics, 18: 1184–1196, 2008.
#'
#' @examples
#' \donttest{
#' one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=37,r_1=9,r_1_e=NA,n_2=53,r_2=16,
#' N_max=140,beta=0.1,sim_num=20)
#' }
#'
#' @export
one_arm_ad_two_stg_binary<-function(p_1,p_0,n_1,r_1,r_1_e=NA,n_2,r_2,N_max,beta,sim_num){
  # p_1 is the event rate for the test arm
  # p_0 is the undesirable evnet rate
  #sim_num=simulation number
  #samsz is sample size
  #N_max is themaximum sample size

  ######################################
  ###  Preparing the parameters ########
  ######################################

  N=c(n_1,n_2)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  s=N/n_2
  t=N/(p_1*(1-p_1))
  t_1=n_1/(p_1*(1-p_1))
  t_2=n_2/(p_1*(1-p_1))
  # t is the vector of information times at each interim anbalysis
  c_2=(r_2/n_2-p_0)*sqrt(n_2^3/(r_2*(n_2-r_2)))
  r_2_0=r_2+0.5   #continuity correction
  c_2_0=(r_2_0/n_2-p_0)*sqrt(n_2^3/(r_2_0*(n_2-r_2_0)))  # critical boundary with continuity correction
  # c_2=c_2_0
  b=stats::qnorm(1-beta)
  # r_1=r_1_f
  #t is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "SINGLE", function() {
    ########################################
    ### Step 1: Generate random samples ####
    ########################################

    x=stats::rbinom(n_2,1,p_1)
    #############################################################################################
    ### Step 2: Obtain the Wald statistics and theta_hat for hypothesis testing and estimate ####
    #############################################################################################

    Wald=NULL
    score=NULL
    theta_hat=NULL
    for (j in 1:length(s)){
      Wald[j]= (mean(x[1:N[j]])-p_0)*sqrt(t[j])
      score[j]= (mean(x[1:N[j]])-p_0)*t[j]
      theta_hat[j]=mean(x[1:N[j]])-p_0
    }

    theta_hat_sd=1/sqrt(t)
    #Wald containes the Wald statistics
    p_hat=mean(x[1:N[1]])
    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    # test compares the Wald statistics with the critical boundaries
    #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1

    if ((sum(x[1:N[1]])>=r_1_e)&&(!is.na(r_1_e))){ #futility termination
      cond_power=1
      futi=0
      rej=1
      rej_0=1
      rej_2=1
      samsz_simu=n_1
      case=4
    }else if (sum(x[1:N[1]])<=r_1) { #futility termination
      cond_power=0
      futi=1
      rej=0
      rej_0=0
      rej_2=0
      samsz_simu=n_1
      case=1
    }else if (sum(x[1:N[1]])>r_1){
      futi=0

      #################################################################################
      ###                   Step 4.2: Compare theta_hat with theta_cut             ####
      ###              if theta_hat>theta_cut, re-estimate sampe size              ####
      #################################################################################
      #sample size re-estimation procedure##

      cond_power=stats::pnorm((p_hat-p_0)*sqrt(t[2]-t[1])-(c_2*sqrt(t[2])-score[1])/sqrt(t[2]-t[1]))

      if (cond_power>=1-beta){
        case=2
        samsz_simu=n_2
        if (sum(x[1:N[2]])>r_2){ # the NULL hypothesis is rejected at the final analysis
          rej=1
          rej_0=1
          rej_2=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej=0
          rej_0=0
          rej_2=0
        }
      }else{

        #################################################################################
        ###                   Step 4.2: Compare theta_hat with theta_cut             ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
        #sample size re-estimation procedure##
        #         T_new_0=(p_hat-p_0)^(-2)*(b+(c_2*sqrt(t[2])-score[1])/sqrt(t[2]-t[1]))+t[1]
        #       new_sampsz=ceiling(T_new_0*2*sigma^2)
        T_new_0=(p_hat-p_0)^(-2)*(b+(c_2_0*sqrt(t[2])-score[1])/sqrt(t[2]-t[1]))+t[1]
        ### new information time with continuity correction ###
        n_new_0=ceiling(T_new_0*p_hat*(1-p_hat))
        #        T_new=n_new/2/sigma^2
        n_new=min(N_max,n_new_0)
        T_new=n_new/(p_hat*(1-p_hat))
        #####  new critical boundary without continuoty correction ####
        e_bry_new=(sqrt(T_new-t[1])/sqrt(t[2]-t[1]))*(c_2*sqrt(t[2])-score[1])+score[1]

        c_bry_new=e_bry_new/sqrt(T_new)

        ####################################

        wald_2_new_sq=NULL  # Wald statistics
        for (i in 1:n_new){
          wald_2_new_sq[i]=(i/n_new-p_0)*sqrt(n_new^3/(i*(n_new-i)))
        }  # sq=sequence ###

        r_2_new_sq=cbind(1:n_new,wald_2_new_sq)
        r_2_new_cut_sq=r_2_new_sq[r_2_new_sq[,2]>=c_bry_new,]
        r_2_new=r_2_new_cut_sq[1,1]  ### The number of response that corresponds to c_bry_new. This will be used as the rejection criterion ###

        ###################### c new critical boundary with continuity correction ####
        e_bry_new_0=(sqrt(T_new-t[1])/sqrt(t[2]-t[1]))*(c_2_0*sqrt(t[2])-score[1])+score[1]

        c_bry_new_0=e_bry_new_0/sqrt(T_new)

        ##############################################

        r_2_new_cut_sq_0=r_2_new_sq[r_2_new_sq[,2]>=c_bry_new_0,]
        r_2_new_0=r_2_new_cut_sq_0[1,1]### The number of response that corresponds to c_bry_new_0. This will be used as the rejection criterion ###

        ###########################################################
        #################################################################################
        ###                   Step 4.4: Compare theta_hat with theta_cut             ####
        ###                     if new_sampsz<=N_max, proceed to new sample size     ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        case=3
        samsz_simu=n_new
        ##### Generate addtional random samples #######

        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################

        x_add=stats::rbinom((n_new-N[1]),1,p_1)
        x_new=c(x[1:N[1]],x_add)


        score_last=(mean(x_new)-p_0)*T_new
        wald_last=(mean(x_new)-p_0)*sqrt(T_new)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        test_asd=wald_last-c_bry_new

        if (sum(x_new)>r_2_new){
          rej=1 ## rejection without continuity correction
        }else{
          rej=0
        }
        if (sum(x_new)>r_2_new_0){
          rej_0=1  ## rejection with continuity correction
        }else{
          rej_0=0
        }

        if (test_asd>=0){
          rej_2=1 ### rejection with critical boundary , instead of number of responses , without continuity correction ######

        }else{
          rej_2=0
        }
      }
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out=c(rej,samsz_simu,futi,case,cond_power,p_hat,rej_2,rej_0)
  }) #(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  sampsz_simu=as.numeric(simr[,2])
  futi=as.numeric(simr[,3])
  case=as.numeric(simr[,4])
  cond_power=as.numeric(simr[,5])
  p_hat=as.numeric(simr[,6])
  rej_2=as.numeric(simr[,7])
  rej_0=as.numeric(simr[,8])
  ###########################

  power=mean(rej)
  power_2=mean(rej_2)
  power_0=mean(rej_0)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)

  case_rt=NULL
  sz_case=NULL
  power_case=NULL
  power_2_case=NULL
  power_0_case=NULL
  cond_power_case=NULL
  p_hat_gt_ptrue_case=NULL
  for (i in 1:4){
    sz_case[i]=mean(sampsz_simu[case==i])
    case_rt[i]=length(case[case==i])/sim_num
    power_case[i]=mean(rej[case==i])
    power_2_case[i]=mean(rej_2[case==i])
    power_0_case[i]=mean(rej_0[case==i])
    cond_power_case[i]=mean(cond_power[case==i])
    p_hat_gt_ptrue_case[i]=length((p_hat[case==i])&(p_hat[case==i]>p_1))/length(p_hat[case==i])
  }

  case_summary=round(rbind(1:4,case_rt,sz_case,power_case,cond_power_case,p_hat_gt_ptrue_case,power_2_case,power_0_case),4)
  rownames(case_summary)=c("case","case rate","mean sample size by case", "rejection rate by case", "conditional power by case", "proportion of p_hat>p_true, by case","rejection 2 rate by case","corrected rejection rate by case")

  #############################
  ### Step 7: Final output ####
  #############################

  out=list(power,power_0, av_samsz,futi_rt, case_summary)
  names(out)=c("rate of rejection-- without continity correction","rejection rate with continuity correction","Average sample size","Rate of early termination for futility","summary by case")
  return(out)
}
#######################################
#one_arm_ad_two_stg_binary(p_1,p_0,n_1,r_1,r_1_e=NA,n_2,r_2,N_max,beta,sim_num)
#Simon's design'#######
# Minimax
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=37,r_1=9,r_1_e=NA,n_2=53,r_2=16,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=37,r_1=9,r_1_e=NA,n_2=53,r_2=16,N_max=140,beta=0.1,sim_num=100)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=37,r_1=9,r_1_e=23,n_2=53,r_2=16,N_max=140,beta=0.1,sim_num=100000)
# optimal
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=22,r_1=5,r_1_e=NA,n_2=66,r_2=19,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=22,r_1=5,r_1_e=NA,n_2=66,r_2=19,N_max=140,beta=0.1,sim_num=100000)
# average
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=30,r_1=7,r_1_e=NA,n_2=57,r_2=17,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=30,r_1=7,r_1_e=NA,n_2=57,r_2=17,N_max=140,beta=0.1,sim_num=100000)

### the hybrid  design ###

##### p_0=0.2, p_low=0.33, p_target=0.4#####
#Minimax
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=100,r_1=25,r_1_e=33,n_2=105,r_2=29,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=100,r_1=25,r_1_e=33,n_2=105,r_2=29,N_max=140,beta=0.1,sim_num=100000)
#Optimal
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=50,r_1=11,r_1_e=NA,n_2=56,r_2=17,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=50,r_1=11,r_1_e=NA,n_2=56,r_2=17,N_max=140,beta=0.1,sim_num=100000)

#average
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=75,r_1=18,r_1_e=26,n_2=80,r_2=23,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=75,r_1=18,r_1_e=26,n_2=80,r_2=23,N_max=140,beta=0.1,sim_num=100000)

### the midpoint design  design , q=0.5###

##### p_0=0.2, p_low=0.33, p_target=0.4#####
#Minimax
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=69,r_1=17,r_1_e=25,n_2=74,r_2=22,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=69,r_1=17,r_1_e=25,n_2=74,r_2=22,N_max=140,beta=0.1,sim_num=100000)

#Optimal
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=36,r_1=7,r_1_e=NA,n_2=56,r_2=17,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=36,r_1=7,r_1_e=NA,n_2=56,r_2=17,N_max=140,beta=0.1,sim_num=100000)

#average
#one_arm_ad_two_stg_binary(p_1=0.2,p_0=0.2,n_1=53,r_1=12,r_1_e=NA,n_2=59,r_2=18,N_max=140,beta=0.1,sim_num=100000)
#one_arm_ad_two_stg_binary(p_1=0.33,p_0=0.2,n_1=53,r_1=12,r_1_e=NA,n_2=59,r_2=18,N_max=140,beta=0.1,sim_num=100000)



####################################

##########################################################################
####              Three stage design  simulations                       ####
##########################################################################

#######################################

#' Simulations: Expanded Simon's design: the three stage design, no sample size change
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param p_1   The assumed response rate.
#' @param n_1   Number of patients at first look.
#' @param r_1   If no more than r_1 responses are observed at the first look, then the trial would be stopped for futility.
#' @param n_2   Number of patients at second look.
#' @param r_2_f   If no more than responses are observed at the second look, then the trial would be stopped for futility.
#' @param r_2_e The trial would be stopped for superiority if at least r_2_e responses are observed at the first look.
#' If r_2_e is not used in the design, enter NA.
#' @param n_3   Number of patients at third look.
#' @param r_3   If at least r_2 responses are observed at the third look, then the null hypothesis will be rejected.
#' @param sim_num Purpose of simulation: To verify type I error control by setting p_1=p_0,
#' To verify the parameters n_1,r_1,n_2,r_2_f,r_2_e,n_3,r_3 are correctly chosen, such that the rejection rate matches that
#' from the design table of the three-stage design.
#'
#' @return Average sample size is EN(p_1).
#' @return Rate of termination is PET(p_1)
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references P. Gao, J.H. Ware, C. Mehta, (2008), Sample size re-estimation for adaptive sequential designs.
#' Journal of Biopharmaceutical Statistics, 18: 1184–1196, 2008.
#'
#' @examples
#' \donttest{
#' one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=37,r_1=8,n_2=100,r_2_f=22,r_2_e=33,
#' n_3=105,r_3=29,sim_num=20)
#' }
#'
#' @export
one_arm_rej_3_stg_binary<-function(p_1,p_0,n_1,r_1,n_2,r_2_f,r_2_e=NA,n_3,r_3,sim_num){
  # p_1 is the event rate for the test arm
  # p_0 is the undesirable evnet rate
  #sim_num=simulation number
  #samsz is sample size
  ######################################
  ###  Preparing the parameters ########
  ######################################
  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "SINGLE", function() {
    ########################################
    ### Step 1: Generate random samples ####
    ########################################

    x=stats::rbinom(n_3,1,p_1)

    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    # test compares the Wald statistics with the critical boundaries
    #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
    x_1=sum(x[1:n_1])
    x_2=sum(x[1:n_2])
    x_3=sum(x[1:n_3])

    if (!is.na(r_2_e)){
      if (x_1<=r_1) { #futility termination
        futi=1
        rej=0
        samsz_simu=n_1
      }else if (x_1>r_1){
        if (x_2<=r_2_f){ #futility termination
          futi=1
          rej=0
          samsz_simu=n_2
        }else if  (x_2>r_2_e){
          futi=0
          rej=1
          samsz_simu=n_2
        }else if ((x_2>r_2_f)&(x_2<=r_2_e)){
          samsz_simu=n_3
          if(x_3>r_3){
            futi=0
            rej=1
          }else {
            futi=0
            rej=0
          }
        }
      }
    }else if (is.na(r_2_e)){
      if (x_1<=r_1) { #futility termination
        futi=1
        rej=0
        samsz_simu=n_1
      }else if (x_1>r_1){
        if (x_2<=r_2_f){ #futility termination
          futi=1
          rej=0
          samsz_simu=n_2
        }else if (x_2>r_2_f){
          samsz_simu=n_3
          if(x_3>r_3){
            futi=0
            rej=1
          }else {
            futi=0
            rej=0
          }
        }
      }
    }

    #################################################################################
    ###                   Step 4.2: Compare theta_hat with theta_cut             ####
    ###              if theta_hat>theta_cut, re-estimate sampe size              ####
    #################################################################################
    #sample size re-estimation procedure##
    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out=c(rej,samsz_simu,futi)
  }) #(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  sampsz_simu=as.numeric(simr[,2])
  futi=as.numeric(simr[,3])
  ###########################

  power=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)

  ###################
  #############################
  ### Step 7: Final output ####
  #############################

  out=list(power,av_samsz,futi_rt)
  names(out)=c("Rate of rejection","Average sample size","Rate of early termination for futility")
  return(out)
}
############################
#one_arm_rej_3_stg_binary(p_1,p_0,n_1,r_1,n_2,r_2_f,r_2_e=NA,n_3,r_3,sim_num)
######## Minimax ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=37,r_1=8,n_2=100,r_2_f=22,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.4,p_0=0.2,n_1=37,r_1=8,n_2=100,r_2_f=22,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=37,r_1=8,n_2=100,r_2_f=22,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
######## optimal ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=22,r_1=3,n_2=50,r_2_f=11,r_2_e=NA,n_3=56,r_3=17,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.4,p_0=0.2,n_1=22,r_1=3,n_2=50,r_2_f=11,r_2_e=NA,n_3=56,r_3=17,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=22,r_1=3,n_2=50,r_2_f=11,r_2_e=NA,n_3=56,r_3=17,sim_num=100000)
##########average  ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=30,r_1=6,n_2=75,r_2_f=16,r_2_e=26,n_3=80,r_3=23,sim_num=1000000)
#one_arm_rej_3_stg_binary(p_1=0.4,p_0=0.2,n_1=30,r_1=6,n_2=75,r_2_f=16,r_2_e=26,n_3=80,r_3=23,sim_num=1000000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=30,r_1=6,n_2=75,r_2_f=16,r_2_e=26,n_3=80,r_3=23,sim_num=1000000)
######## Minimax ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=25,r_1=4,n_2=100,r_2_f=25,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.41,p_0=0.2,n_1=25,r_1=4,n_2=100,r_2_f=25,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=25,r_1=4,n_2=100,r_2_f=25,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
######## optimal ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=24,r_1=4,n_2=50,r_2_f=10,r_2_e=NA,n_3=55,r_3=17,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.41,p_0=0.2,n_1=24,r_1=4,n_2=50,r_2_f=10,r_2_e=NA,n_3=55,r_3=17,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=24,r_1=4,n_2=50,r_2_f=10,r_2_e=NA,n_3=55,r_3=17,sim_num=100000)
##########average  ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=25,r_1=4,n_2=75,r_2_f=18,r_2_e=26,n_3=80,r_3=23,sim_num=1000000)
#one_arm_rej_3_stg_binary(p_1=0.41,p_0=0.2,n_1=25,r_1=4,n_2=75,r_2_f=18,r_2_e=26,n_3=80,r_3=23,sim_num=1000000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=25,r_1=4,n_2=75,r_2_f=18,r_2_e=26,n_3=80,r_3=23,sim_num=1000000)
############################
######## Minimax ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=48,r_1=11,n_2=100,r_2_f=23,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.38,p_0=0.2,n_1=48,r_1=11,n_2=100,r_2_f=23,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=48,r_1=11,n_2=100,r_2_f=23,r_2_e=33,n_3=105,r_3=29,sim_num=100000)
######## optimal ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=27,r_1=4,n_2=50,r_2_f=11,r_2_e=NA,n_3=68,r_3=20,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.38,p_0=0.2,n_1=27,r_1=4,n_2=50,r_2_f=11,r_2_e=NA,n_3=68,r_3=20,sim_num=100000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=27,r_1=4,n_2=50,r_2_f=11,r_2_e=NA,n_3=68,r_3=20,sim_num=100000)
##########average  ###
#one_arm_rej_3_stg_binary(p_1=0.2,p_0=0.2,n_1=38,r_1=8,n_2=75,r_2_f=17,r_2_e=NA,n_3=80,r_3=23,sim_num=1000000)
#one_arm_rej_3_stg_binary(p_1=0.38,p_0=0.2,n_1=38,r_1=8,n_2=75,r_2_f=17,r_2_e=NA,n_3=80,r_3=23,sim_num=1000000)
#one_arm_rej_3_stg_binary(p_1=0.33,p_0=0.2,n_1=38,r_1=8,n_2=75,r_2_f=17,r_2_e=NA,n_3=80,r_3=23,sim_num=1000000)
############################

##########################################

#' Simulations: adaptive three stage design
#'
#' @param p_0   Response rate indicating low activity with insufficient clinical benefit.
#' @param p_1   The assumed response rate.
#' @param n_1   Number of patients at first look.
#' @param r_1   If no more than r_1 responses are observed at the first look, then the trial would be stopped for futility.
#' @param n_2   Number of patients at second look.
#' @param r_2   If no more than r_2 responses are observed at the second look, then the trial would be stopped for futility.
#' @param r_2_e   The trial would be stopped for superiority if at least r_2_e responses are observed at the second look.
#' If r_2_e is not used in the design, enter NA.
#' @param n_3   Number of patients at third look.
#' @param r_3   If no more than r_3 responses are observed at the third look, then the null hypothesis will not be rejected.
#' @param N_max Maximum sample size for the trial.
#' @param beta  Type II error. The target power is 1-beta.
#' @param sim_num Purpose of simulation: To verify type I error control by setting p_1=p_0,
#' To evaluate the operating characteristics of the adaptive design and to determine if continuity correction is necessary.
#'
#' @return Average sample size is EN(p_1).
#' @return Rate of termination is PET(p_1)
#' @return If rejection rate without continuity correction does not exceed alpha, the continuity correction is not needed for this group of
#' parameters p_0,n_1,r_1,n_2,r_2,r_2_e,n_3,r_3. Otherwise, continuity correction should be applied.
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials: an expansion of Simon’s design,
#' Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#' @references P. Gao, J.H. Ware, C. Mehta, (2008), Sample size re-estimation for adaptive sequential designs.
#' Journal of Biopharmaceutical Statistics, 18: 1184–1196, 2008.
#'
#' @examples
#' \donttest{
#' one_arm_ad_3_stg_binary(p_1=0.2,p_0=0.2,n_1=37,r_1=8,n_2=100,r_2=22,r_2_e=33,
#' n_3=105,r_3=29,N_max=150,beta=0.1,sim_num=20)
#' }
#'
#' @export
one_arm_ad_3_stg_binary<-function(p_1,p_0,n_1,r_1,n_2,r_2,r_2_e=NA,n_3,r_3,N_max,beta,sim_num){
  # p_1 is the event rate for the test arm
  # p_0 is the undesirable evnet rate
  #sim_num=simulation number
  #N_max is the maximum sample size

  ######################################
  ###  Preparing the parameters ########
  ######################################

  N=c(n_1,n_2,n_3)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  s=N/n_3
  t=N/(p_1*(1-p_1))
  t_1=n_1/(p_1*(1-p_1))
  t_2=n_2/(p_1*(1-p_1))
  t_3=n_3/(p_1*(1-p_1))
  # t is the vector of information times at each interim anbalysis
  c_3=(r_3/n_3-p_0)*sqrt(n_3^3/(r_3*(n_3-r_3)))
  r_3_0=r_3+0.5   #continuity correction
  c_3_0=(r_3_0/n_3-p_0)*sqrt(n_3^3/(r_3_0*(n_3-r_3_0)))  # critical boundary with continuity correction
  # c_2=c_2_0
  b=stats::qnorm(1-beta)
  # r_1=r_1_f
  #t is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "SINGLE", function() {
    ########################################
    ### Step 1: Generate random samples ####
    ########################################

    x=stats::rbinom(n_3,1,p_1)
    #############################################################################################
    ### Step 2: Obtain the Wald statistics and theta_hat for hypothesis testing and estimate ####
    #############################################################################################

    Wald=NULL
    score=NULL
    theta_hat=NULL
    for (j in 1:length(s)){
      Wald[j]= (mean(x[1:N[j]])-p_0)*sqrt(t[j])
      score[j]= (mean(x[1:N[j]])-p_0)*t[j]
      theta_hat[j]=mean(x[1:N[j]])-p_0
    }

    theta_hat_sd=1/sqrt(t)
    #Wald containes the Wald statistics
    p_hat=mean(x[1:N[2]])
    ###################################################################################################
    ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
    ###################################################################################################

    # test compares the Wald statistics with the critical boundaries
    #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
    if (sum(x[1:N[1]])<=r_1) { #futility termination
      cond_power=0
      futi=1
      rej=0
      rej_0=0
      rej_2=0
      samsz_simu=n_1
      case=1
    }else if ((sum(x[1:N[2]])>=r_2_e)&&(!is.na(r_2_e))){ #futility termination
      cond_power=1
      futi=0
      rej=1
      rej_0=1
      rej_2=1
      samsz_simu=n_2
      case=5
    }else if (sum(x[1:N[2]])<=r_2) { #futility termination
      cond_power=0
      futi=1
      rej=0
      rej_0=0
      rej_2=0
      samsz_simu=n_2
      case=2
    }else if (sum(x[1:N[2]])>r_2){
      futi=0

      #################################################################################
      ###                   Step 4.2: Compare theta_hat with theta_cut             ####
      ###              if theta_hat>theta_cut, re-estimate sampe size              ####
      #################################################################################
      #sample size re-estimation procedure##

      cond_power=stats::pnorm((p_hat-p_0)*sqrt(t[3]-t[2])-(c_3*sqrt(t[3])-score[2])/sqrt(t[3]-t[2]))

      if (cond_power>=1-beta){
        case=3
        samsz_simu=n_3
        if (sum(x[1:N[3]])>r_3){ # the NULL hypothesis is rejected at the final analysis
          rej=1
          rej_0=1
          rej_2=1
        }else { # the NULL hypothesis is NOT rejected at the final analysis
          rej=0
          rej_0=0
          rej_2=0
        }
      }else{
        #################################################################################
        ###                   Step 4.2: Compare theta_hat with theta_cut             ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
        #sample size re-estimation procedure##
        # T_new_0=(p_hat-p_0)^(-2)*(b+(c_2*sqrt(t[2])-score[1])/sqrt(t[2]-t[1]))+t[1]
        # new_sampsz=ceiling(T_new_0*2*sigma^2)
        T_new_0=(p_hat-p_0)^(-2)*(b+(c_3_0*sqrt(t[3])-score[2])/sqrt(t[3]-t[2]))+t[2]
        ### new information time with continuity correction ###
        n_new_0=ceiling(T_new_0*p_hat*(1-p_hat))
        # T_new=n_new/2/sigma^2
        n_new=min(N_max,n_new_0)
        T_new=n_new/(p_hat*(1-p_hat))
        #####  new critical boundary without continuoty correction ####
        e_bry_new=(sqrt(T_new-t[2])/sqrt(t[3]-t[2]))*(c_3*sqrt(t[3])-score[2])+score[2]

        c_bry_new=e_bry_new/sqrt(T_new)

        ####################################

        wald_3_new_sq=NULL  # Wald statistics
        for (i in 1:n_new){
          wald_3_new_sq[i]=(i/n_new-p_0)*sqrt(n_new^3/(i*(n_new-i)))
        }  # sq=sequence ###

        r_3_new_sq=cbind(1:n_new,wald_3_new_sq)
        r_3_new_cut_sq=r_3_new_sq[r_3_new_sq[,2]>=c_bry_new,]
        r_3_new=r_3_new_cut_sq[1,1]  ### The number of response that corresponds to c_bry_new. This will be used as the rejection criterion ###

        ###################### c new critical boundary with continuity correction ####
        e_bry_new_0=(sqrt(T_new-t[2])/sqrt(t[3]-t[2]))*(c_3_0*sqrt(t[3])-score[2])+score[2]

        c_bry_new_0=e_bry_new_0/sqrt(T_new)

        ##############################################

        r_3_new_cut_sq_0=r_3_new_sq[r_3_new_sq[,2]>=c_bry_new_0,]
        r_3_new_0=r_3_new_cut_sq_0[1,1]### The number of response that corresponds to c_bry_new_0. This will be used as the rejection criterion ###

        ###########################################################
        #################################################################################
        ###                   Step 4.4: Compare theta_hat with theta_cut             ####
        ###                     if new_sampsz<=N_max, proceed to new sample size     ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        case=4
        samsz_simu=n_new
        ##### Generate addtional random samples #######

        #################################################################################
        ###                   Step 4.5: generate new random samples                  ####
        ###                         and calculate statistics                         ####
        #################################################################################

        x_add=stats::rbinom((n_new-N[2]),1,p_1)
        x_new=c(x[1:N[2]],x_add)


        score_last=(mean(x_new)-p_0)*T_new
        wald_last=(mean(x_new)-p_0)*sqrt(T_new)

        #################################################################################
        ###                   Step 4.6: Hypothesis testing                           ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        test_asd=wald_last-c_bry_new

        if (sum(x_new)>r_3_new){
          rej=1 ## rejection without continuity correction
        }else{
          rej=0
        }
        if (sum(x_new)>r_3_new_0){
          rej_0=1  ## rejection with continuity correction
        }else{
          rej_0=0
        }

        if (test_asd>=0){
          rej_2=1 ### rejection with critical boundary , instead of number of responses , without continuity correction ######
        }else{
          rej_2=0
        }
      }
    }

    #######################################################
    ### Step 5: output and finish the simjulation loop ####
    #######################################################
    out=c(rej,samsz_simu,futi,case,cond_power,p_hat,rej_2,rej_0)
  }) #(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  sampsz_simu=as.numeric(simr[,2])
  futi=as.numeric(simr[,3])
  case=as.numeric(simr[,4])
  cond_power=as.numeric(simr[,5])
  p_hat=as.numeric(simr[,6])
  rej_2=as.numeric(simr[,7])
  rej_0=as.numeric(simr[,8])
  ###########################

  power=mean(rej)
  power_2=mean(rej_2)
  power_0=mean(rej_0)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)

  case_rt=NULL
  sz_case=NULL
  power_case=NULL
  power_2_case=NULL
  power_0_case=NULL
  cond_power_case=NULL
  p_hat_gt_ptrue_case=NULL
  for (i in 1:5){
    sz_case[i]=mean(sampsz_simu[case==i])
    case_rt[i]=length(case[case==i])/sim_num
    power_case[i]=mean(rej[case==i])
    power_2_case[i]=mean(rej_2[case==i])
    power_0_case[i]=mean(rej_0[case==i])
    cond_power_case[i]=mean(cond_power[case==i])
    p_hat_gt_ptrue_case[i]=length((p_hat[case==i])&(p_hat[case==i]>p_1))/length(p_hat[case==i])
  }

  case_summary=round(rbind(1:5,case_rt,sz_case,power_case,cond_power_case,p_hat_gt_ptrue_case,power_2_case,power_0_case),4)
  rownames(case_summary)=c("case","case rate","mean sample size by case", "rejection rate by case", "conditional power by case", "proportion of p_hat>p_true, by case","rejection 2 rate by case","corrected rejection rate by case")

  #############################
  ### Step 7: Final output ####
  #############################

  out=list(power,power_0,av_samsz,futi_rt)
  # names(out)=c("rate of rejection-- without continity correction","rejection rate with continuity correction","rate of rejection using critical boundary --without continuity correction ","mean sample size","futility stop rate",)
  names(out)=c("rate of rejection-- without continity correction","rejection rate with continuity correction","Average sample size","Rate of early termination for futility")
  return(out)
}
#######################################
#one_arm_ad_3_stg_binary(p_1,p_0,n_1,r_1,n_2,r_2,r_2_e=NA,n_3,r_3,N_max,beta,sim_num)

## Minimax
#one_arm_ad_3_stg_binary(p_1=0.2,p_0=0.2,n_1=37,r_1=8,n_2=100,r_2=22,r_2_e=33,n_3=105,r_3=29,N_max=150,beta=0.1,sim_num=100000)
#one_arm_ad_3_stg_binary(p_1=0.33,p_0=0.2,n_1=37,r_1=8,n_2=100,r_2=22,r_2_e=33,n_3=105,r_3=29,N_max=150,beta=0.1,sim_num=100000)

# Optimal
#one_arm_ad_3_stg_binary(p_1=0.2,p_0=0.2,n_1=22,r_1=3,n_2=50,r_2=11,r_2_e=NA,n_3=56,r_3=17,N_max=150,beta=0.1,sim_num=100000)
#one_arm_ad_3_stg_binary(p_1=0.33,p_0=0.2,n_1=22,r_1=3,n_2=50,r_2=11,r_2_e=NA,n_3=56,r_3=17,N_max=150,beta=0.1,sim_num=100000)

## average
#one_arm_ad_3_stg_binary(p_1=0.2,p_0=0.2,n_1=30,r_1=6,n_2=75,r_2=16,r_2_e=26,n_3=80,r_3=23,N_max=150,beta=0.1,sim_num=100000)
#one_arm_ad_3_stg_binary(p_1=0.33,p_0=0.2,n_1=30,r_1=6,n_2=75,r_2=16,r_2_e=26,n_3=80,r_3=23,N_max=150,beta=0.1,sim_num=100000)
