
#install.packages("glm")
#install.packages("glm2")
#install.packages("MASS")
#library(tidyverse)
#library(ggplot2)


#############################################

###################################################
##########   Part II                    ############
###       Simulations                  ############
###################################################

############################ Simulations #########################

##########################################################################
### II.1 .Simulations  -normal distributions                         #####
##########################################################################

#########################################
two_stg_rej_rate_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry_two_stg,samsz,sim_num){

  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better
  N=ceiling(samsz*s)
  comp=length(mu_t)
  mu=mu_t-mu_c
  # N is the vetoc of sample sizes at which interim analysis will be conducted
  T<- samsz/(sigma_t^2+sigma_c^2)
  t=T*s
  #t is the vector of information times at the interim analysis
  N1=N[1]
  #N1 is the sample at the end of phase 2 stage
  NF=N[length(s)]
  ######################## exit boundaries
  e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined

  ############ the simulation loops #
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    rej=0; stg=0
    s_slct=0; control_samsz_simu=N[length(s)]
    tot_samsz_simu=N[1]*(comp-1)+2*control_samsz_simu
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################

      x0=stats::rnorm(samsz,mu_c,sigma_c)
      x=matrix(stats::rnorm(samsz*comp, mu_t,sigma_t), nrow=comp,ncol=samsz)

      #x matirx does not include the control sample.
      #The row number is the same as the active dose

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################
      S_ph2=NULL
      for (j in 1:length(mu_t)){
        S_ph2[j]=(mean(x[j,1:N1])-mean(x0[1:N1]))*t[1]
      }
      wald_ph2=S_ph2/sqrt(t[1])
      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(max(S_ph2),S_ph2)
      #s_slct identifies which dose was selected for the i-th simulation
      ###########################################################
      ### Step 4: generate  random samples for phase 3 stage ####
      ###########################################################
      # x0_ph3=x0[(N1+1):NF]
      # x_slct_ph3=x[s_slct,(N1+1):NF]
      #y are the samples for the selected dose for the phase 3 stage
      x0_two_stg=x0[1:NF]
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_slct_two_stg=x[s_slct,1:NF]
      #x_two_stg contains all the random samples(including samples from the phase 2 stage) for the seected dose

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      theta_hat=NULL

      if(direction==1){ #larger mu is better
      for (j in 1:length(s) ){
        theta_hat[j]=mean(x_slct_two_stg[1:N[j]])-mean(x0_two_stg[1:N[j]])  # theta_hat at each interim analysis
      }
      }else if (direction==2){ #smaller mu is better
        for (j in 1:length(s) ){
          theta_hat[j]=mean(x0_two_stg[1:N[j]])-mean(x_slct_two_stg[1:N[j]]) # theta_hat at each interim analysis
        }
      }

      theta_hat_sd=1/sqrt(t)
      wald_two_stg=theta_hat*sqrt(t)
      ######  calculating the statistics for phase 3 only #######

      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test=wald_two_stg-c_bry_two_stg
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if (isTRUE(max(test)>0)){

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        last_vt=stg
        control_samsz_simu=N[last_vt]
        tot_samsz_simu=N[1]*(comp-1)+2*control_samsz_simu

        #stg identifies the interim analysis at which the null hypothesis is rejected

      }else{

        rej=0
        stg=0
        last_vt=length(s)
        control_samsz_simu=N[last_vt]
        tot_samsz_simu=N[1]*(comp+1)+2*(control_samsz_simu-N[1])
      }
      #######################################################
      ### Step 7: output and finish the simjulation loop ####
      #######################################################
      out=c(rej,stg,s_slct,control_samsz_simu,tot_samsz_simu)

  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr

  ############################################
  ### Step 8: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  s_slct= as.numeric(simr[,3])
  control_samsz_simu= as.numeric(simr[,4])
  tot_samsz_simu= as.numeric(simr[,5])
  rej_rt=mean(rej)

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=sum(rej[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  select_rt=rep(0,length(mu))
  for (j in 1:length(mu)){
    select_rt[j]=length( s_slct[ s_slct==j])/sim_num
  }

  rej_rt_indiv=rep(0,length(mu))
  for (j in 1:length(mu)){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/sim_num
  }

  con_rej_rt_indiv=rep(0,length(mu))
  for (j in 1:length(mu)){
    if (isTRUE(select_rt[j]==0)){
      con_rej_rt_indiv[j]=NA
    }else {
      con_rej_rt_indiv[j]=rej_rt_indiv[j]/select_rt[j]
    }
  }

  results_by_arm=cbind(mu_t,mu_c,select_rt, rej_rt_indiv,con_rej_rt_indiv)
  colnames(results_by_arm)=c(LABEL()$mu_t1,LABEL()$mu_c,"selection rate"," rejection rate by each mu_t","conditional rejection rate")

  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look")

  tyI_er=sum(rej[mu[s_slct]==0])/sim_num
  #  tyI_er=sum(rej_0)/sim_num

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/sim_num
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(mu_t)-1)+2*samsz
  sim_results=c(rej_rt,tyI_er, non_tyI_er_rej_rt,samsz,mean_control_samsz,tot_samsz,mean_tot_samsz)
  names(sim_results)=c("rejection rate","Type I error","Non-Type I error rejection rate",
                       "Planned control sample size","Mean control sample size","Planned total sample size","Mean total sample size")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(seq_design,results_by_arm,sim_results)
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)

}


#######################################
#two_stg_rej_rate_simu_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry_two_stg,samsz,sim_num)
#two_stg_rej_rate_simu_normal(mu_t=c(0.3,0.3),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.746142,2.162836),samsz=100,sim_num=10000)
#two_stg_rej_rate_simu_normal(mu_t=c(0,0.3),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.746142,2.162836),samsz=100,sim_num=1000)

####################################

##############################################
##############################################

#' Phase 2/3 seamless combination Group sequential design: Simulations, Normal distribution, Type I error and power
#'
#' @param mu_t Vector of means of the test arms, one entry per dose
#' @param mu_c Mean of the control arm
#' @param sigma_t Standard deviation of the test arms
#' @param sigma_c Standard deviation of the control arm
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param samsz Planned sample size of the control arm used in the simulation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries" and "rejection by look"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "mu_t1,...", "mu_c", "selection rate", " rejection rate by each mu_t" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the overall "rejection rate", the "Type I error", the "Non-Type I error rejection rate", the "Planned control sample size", the "Mean control sample size", the "Planned total sample size" and the "Mean total sample size"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_power_simu_normal(mu_t = c(0,0,0), mu_c = 0, sigma_t = 1, sigma_c = 1, direction = 1,
#' s = c(1/3,1), bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
two_stg_power_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,samsz,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better

  alpha=alpha_2/2
  mu=mu_t-mu_c
  comp=length(mu)
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))
  simu_out=two_stg_rej_rate_simu_normal(mu_t=mu_t,mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s,c_bry_two_stg=c_bry_two_stg,samsz=samsz,sim_num=sim_num)
 # planned_tot_samsz=ceiling(samsz*s[1]*(comp+1)+samsz*(s[length(s)]-s[1])*2)
 # out=list(simu_out,samsz, planned_tot_samsz)
#  names(out)=c("simulation summary","planned control sample size","planned total sample size")
  out=simu_out
  return(out)
}

##################################
#two_stg_power_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=100000)
#two_stg_power_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=100000)
#two_stg_power_simu_normal(mu_t=c(0,0.1,0.2),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)
#two_stg_power_simu_normal(mu_t=c(0,0.1,0.2),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),samsz=100,sim_num=1000)
#two_stg_power_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),samsz=100,sim_num=100000)


###############Estimates#####################
#two_stg_ci_simu_normal(mu=c(rep(0,3)),s=c(1/3,2/3,1),c_bry_two_stg=c(3.938325, 2.784816,2.273793),c_bry_p3=c(3.954862,1.977431),sim_num=100000,samsz=100,sigma=1,alpha_2=0.05,dist=1)

##########################################################################
### II.1.2 .Simulations with planned TSSSD -normal distributions     #####
#### This will not be included in the software                       #####
##########################################################################

##########################################################################
### II.1.2 .Simulations with planned TSSSD -normal distributions     #####
##          type I error and power                                   #####
####    This will  be included in the software                       #####
##########################################################################

#############################################
Two_stage_adapt_rej_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist){

  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better

  N=ceiling(samsz*s)
  comp=length(mu_t)
  mu=mu_t-mu_c
   # N is the vetoc of sample sizes at which interim analysis will be conducted
  T<- samsz/(sigma_t^2+sigma_c^2)
  t=T*s
   #t is the vector of information times at the interim analysis
  N1=N[1]
  #N1 is the sample at the end of phase 2 stage
  NF=N[length(s)]
  ######################## exit boundaries
 # e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined

  ############ the simulation loops #
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=N[length(s)]
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    control_sampsz_simu=sampsz_simu; tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
    futi=0
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################

      x0=stats::rnorm(N_max,mu_c,sigma_c)
      x=matrix(stats::rnorm(N_max*comp, mu_t,sigma_t), nrow=comp,ncol=N_max)

      #x matirx does not include the control sample.
      #The row number is the same as the active dose

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################
      S_ph2=NULL

      if (direction==1){ #larger mu is better
      for (j in 1:length(mu_t)){
        S_ph2[j]=(mean(x[j,1:N1])-mean(x0[1:N1]))*t[1]
      }
      }else if (direction==2){ #smaller mu is better
        for (j in 1:length(mu_t)){
          S_ph2[j]=(mean(x0[1:N1])-mean(x[j,1:N1]))*t[1]
        }
      }

      wald_ph2=S_ph2/sqrt(t[1])
      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(max(S_ph2),S_ph2)
      #s_slct identifies which dose was selected for the i-th simulation
      ###########################################################
      ### Step 4: generate  random samples for phase 3 stage ####
      ###########################################################
    #  x0_ph3=x0[(N1+1):NF]
     # x_slct_ph3=x[s_slct,(N1+1):NF]
      #y are the samples for the selected dose for the phase 3 stage
      x0_two_stg=x0[1:NF]
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_slct_two_stg=x[s_slct,1:NF]
      #x_two_stg contains all the random samples(including samples from the phase 2 stage) for the seected dose

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
         theta_hat=NULL
      for (j in 1:length(s) ){
         theta_hat[j]=mean(x_slct_two_stg[1:N[j]])-mean(x0_two_stg[1:N[j]])  # theta_hat at each interim analysis
      }
      theta_hat_sd=1/sqrt(t)
      wald_two_stg=theta_hat*sqrt(t)
      ######  calculating the statistics for phase 3 only #######

      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test_1=test[-length(s)]  ## test statistics excluding the final analysis

      ##########################################
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries
      ############################################

      if (isTRUE(theta_hat[length(s)-1]<theta_cut)) {futi=1} else {futi=0}

      if (isTRUE(max(test)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
      }

      #If the NUll hypothesis is rejected at up to length(s)-1, terminate trial and estimate CI ####
      if (isTRUE(max(test_1)>=0)){
        #  if ( max(tes_1)>0){ #null hypothesis rejected at the next to last interim analysis. no more random sampling
         rej=1

        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
      #  theta_hat_last=theta_hat[last_vt]
      #  theta_hat_last_sd=1/sqrt(t[last_vt])

        ###########################
        sampsz_simu=N[last_vt]
        control_sampsz_simu=sampsz_simu
        tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu


      }else if (isTRUE(max(test_1)<0)){ ## (theta_p3_hat[length(s_1)-1]<theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation


        if (futi==1){ #terminate trial ath the final analysis due to futility. No sample size change
          if (isTRUE(test[length(s)]>0)){
            rej=1
            #The null hypothesis is rejected at the final analysis
            stg=length(s)
          }else{
            rej=0
            stg=0
          }
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=1/sqrt(t[last_vt])

          ###########################
          sampsz_simu=N[last_vt]
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu


        }else if (futi==0){ # no futility ter,mination. proceed to sample size re-estimation
          inter_vt=length(s)-1

          theta=theta_hat[inter_vt]

          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=1/sqrt(t[inter_vt])
          score_inter=theta_hat_inter*t[inter_vt]

         #c_bry=c_bry_two_stg[length(s)]

      #    cond_power=as.numeric(asd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd))
         cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))

          if (isTRUE(cond_power>=1-beta)){

            #stg identifies the interim analysis at which the null hypothesis is rejected
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=1/sqrt(t[last_vt])

            sampsz_simu=N[last_vt]
            control_sampsz_simu=sampsz_simu
            tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }
          }else  if (isTRUE(cond_power<1-beta)){

            ## proceed to sample size re-estimation #####

            ############ Sample size re-estimation with respect toe the combined phase 2/3 test ##

            T_new_0=theta_hat_inter^(-2)*((c_bry_two_stg[length(s)]*sqrt(T)-score_inter)/sqrt(T-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_sampsz=ceiling(T_new_0*(sigma_t^2+sigma_c^2))
            T_new=new_sampsz/(sigma_t^2+sigma_c^2)
            e_bry_new_two_stg=
              (sqrt(T_new-t[inter_vt])/sqrt(T-t[inter_vt]))*(c_bry_two_stg[length(s)]*sqrt(T)-score_inter)+ score_inter
            #new exit boundary for the combined phase 2/3 test. the NULL is rejected if the score crossed the exit boundary

            c_bry_new_two_stg=e_bry_new_two_stg/sqrt(T_new)  #new critical  boundary for the combined phase 2/3 test. the NULL is rejected if the Wald statistics crossed the critical boundary

            ###############################################

            if (isTRUE(new_sampsz>N_max)){ ## new sample size is too large, the sample size will not change. Trial ends at the planned sample size
              if (isTRUE(test[length(test)]>0)){ ### the NULL hypothesis is rejected at the final anbalysis, with the planned sample size
                rej=1
                stg=length(s)
              }else{### the NULL hypothesis is NOT rejected at the final anbalysis, with the planned sample size
                rej=0
                stg=0
              }
              ###################### Statistics for estimates ####
              ### combined phase 2/3 data ####
              last_vt=length(s)

              sampsz_simu=N[last_vt]
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              theta_hat_last=theta_hat[last_vt] ## theta_hat for the combined phase 2/3 data. Used for estimating the lower bound of the CI
              theta_hat_last_sd=1/sqrt(T)


              #################
              ####################### this completes the estimation when there is no saple size change ########

            } else if (isTRUE(new_sampsz<= N_max)) {## new sample size does not exceed the pre-selected maximum. Will change the sample  size to this new sample size
              #here

              sampsz_simu=new_sampsz
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              x0_two_stg_new=x0[1:new_sampsz] # control arm [combined phase 2/3] data
              x_slct_two_stg_new=x[s_slct,1:new_sampsz]  # data of the selected dose arm

              ### statistics from combined data ####
              if (direction==1){ #larger mu is better
              theta_hat_last_ad=(mean(x_slct_two_stg_new)-mean(x0_two_stg_new)) # theta_hat at the adjusted final analysis
              } else if (direction==2){ #smaller  mu is better
                theta_hat_last_ad=mean(x0_two_stg_new)-mean(x_slct_two_stg_new) # theta_hat at the adjusted final analysis
              }
              wald_last_ad=theta_hat_last_ad*sqrt(T_new)  # Wald statistics at the adjusted final analysis
              score_last_ad=theta_hat_last_ad*T_new  # score statistics at the adjusted final analysis

              theta_hat_last_ad_sd=1/sqrt(T_new)

              ##########################################

              test_add=wald_last_ad-c_bry_new_two_stg  # compare the final Wald statistics with the new adjusted critical boundary

              if (isTRUE(test_add>0)){ # the new critical boundary is crossed, the NULL is rejected

                rej=1
                stg=length(s)

              }else{  # the new critical boundary is NOT crossed, the NULL is rejected

                rej=0
                stg=0

              }

              ###################### Parameter estimate with sample size change ##################

            }
          }
        }
      }
      out=c(rej_nossr,rej,stg,s_slct,control_sampsz_simu,tot_sampsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  control_samsz_simu = as.numeric(simr[,5])
  tot_samsz_simu = as.numeric(simr[,6])
  futi=as.numeric(simr[,7])

  #############termination timepoint
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  rej_num=sum(rej>0)
  # rej_num1=length(rej[rej>0])
  # rej_num2=sum(stg>0)

  rej_time=rep(0,(length(s)))
  for (j in 1:(length(s))){

    count_j=c(1:sim_num)[stg==j]
    rej_time[j]=sum(rej[count_j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  select_rt=rep(0,length(mu))
  for (j in 1:length(mu)){
    select_rt[j]=length( s_slct[ s_slct==j])/sim_num
  }


  rej_rt_indiv=rep(0,length(mu))
  for (j in 1:length(mu)){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/sim_num
  }

  con_rej_rt_indiv=rep(0,length(mu))
  for (j in 1:length(mu)){
    if (isTRUE(select_rt[j]==0)){
      con_rej_rt_indiv[j]=NA
    }else {
      con_rej_rt_indiv[j]=rej_rt_indiv[j]/select_rt[j]
    }
  }

  samsz_icr=NULL
  for (i in 1:length(control_samsz_simu)){
    if (control_samsz_simu[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }

  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_rt_byLook=c(rep(0,(length(s)-1)),samsz_icr_rt)

  futi_rt=mean(futi)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= control_samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  #############################
  ### Step 9: Final output ####
  #############################

  results_by_arm=cbind(mu_t,mu_c,select_rt, rej_rt_indiv,con_rej_rt_indiv)
  colnames(results_by_arm)=c(LABEL()$mu_t1,LABEL()$mu_c,"selection rate","rejection rate by each mu_t","conditional rejection rate")
  #  rej_0=rej[mu[s_slct]==0]

  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time,samsz_icr_rt_byLook)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look","Sample size increase")

  tyI_er=sum(rej[mu[s_slct]==0])/sim_num
  #  tyI_er=sum(rej_0)/sim_num

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/sim_num
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(mu_t)-1)+2*samsz
  sim_results=c(rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er, non_tyI_er_rej_rt,futi_rt,samsz,
                N_max,mean_control_samsz,tot_samsz,mean_tot_samsz, mean_samsz_simu_control_samsz_icr)
  names(sim_results)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned control sample size","Maximum control sample size","Mean control sample size",
                       "Planned total sample size","Mean total sample size","mean control sample size conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(seq_design,results_by_arm,sim_results)
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)

}

#################################################
#Two_stage_adapt_rej_simu_normal(mu_t,mu_c,sigma_t,sigma_c,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist)

#Two_stage_adapt_rej_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000,dist=1)
#Two_stage_adapt_rej_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.605339, 3.505080, 2.240072),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=10000,dist=1)
#Two_stage_adapt_rej_simu_normal(mu_t=c(0,0.1,0.2),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.605339, 3.505080, 2.240072),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=500,sim_num=10000,dist=1)


#########################

#################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Normal distribution, Operating characteristics
#'
#' @param mu_t Vector of means of the test arms, one entry per dose
#' @param mu_c Mean of the control arm
#' @param sigma_t Standard deviation of the test arms
#' @param sigma_c Standard deviation of the control arm
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param theta_cut Lack-of-efficacy threshold: an effect size below which the test drug is not considered to provide clinically meaningful efficacy
#' @param samsz Planned sample size of the control arm used in the simulation
#' @param N_max Maximum sample size (control arm) allowed after sample size re-estimation
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named numeric vector with the following components:
#' \describe{
#'   \item{Rejection rate without SSR}{Rejection rate of the trial if the sample size is not re-estimated}
#'   \item{Rejection rate with SSR}{Rejection rate of the trial under the adaptive procedure with sample size re-estimation}
#'   \item{rate of sample size increase}{Rate of simulations in which the sample size is increased}
#'   \item{Rejection rate without SSR/(meeting sample size incerase criteria)}{Rejection rate without sample size re-estimation among the simulations that meet the sample size increase criteria}
#'   \item{Rejection rate with SSR/(sample size incerase)}{Rejection rate with sample size re-estimation among the simulations in which the sample size is increased}
#'   \item{Type I error}{Simulated one-sided type I error rate}
#'   \item{Non-Type I error rejection rate}{Simulated rejection rate when the true effect is non-zero}
#'   \item{Futility rate}{Rate of simulations stopped for lack of efficacy, i.e. for which the observed interim effect size is below theta_cut}
#'   \item{Planned control sample size}{Planned sample size of the control arm}
#'   \item{Maximum control sample size}{Maximum permitted sample size of the control arm}
#'   \item{Mean control sample size}{Mean simulated sample size of the control arm}
#'   \item{Planned total sample size}{Planned total sample size of the trial}
#'   \item{Mean total sample size}{Mean simulated total sample size of the trial}
#'   \item{mean control sample size conditional on increase}{Mean simulated control sample size among the simulations in which the sample size is increased}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_adapt_power_simu_normal(mu_t = c(0,0,0), mu_c = 0, sigma_t = 1, sigma_c = 1,
#' direction = 1, s = c(1/3,2/3,1), bry_type = 2, alpha_2 = c(0.001,0.001,0.048),
#' theta_cut = 0.05, samsz = 100, N_max = 500, beta = 0.1, sim_num = 20)
#'
#' @export
two_stg_adapt_power_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,theta_cut,samsz,N_max,beta,sim_num){
  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better

  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  mu=mu_t-mu_c
  comp=length(mu)

  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }


  c_bry_two_stg=as.numeric(unlist(design_out[1]))

  simu_out=Two_stage_adapt_rej_simu_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist=1)

  out= simu_out
  return(out)
}
##################################
#two_stg_adapt_power_simu_normal(mu_t=c(0,0.4,0.5),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=1000,beta=0.1,sim_num=10000)
#two_stg_adapt_power_simu_normal(mu_t=c(0.4,0.5),mu_c=0,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=1000,sigma=1,beta=0.1,sim_num=100)
#two_stg_adapt_power_simu_normal(mu_t=c(0.4,0.5),mu_c=0.4,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=1000,sigma=1,beta=0.1,sim_num=100)
#two_stg_adapt_power_simu_normal(mu_t=c(0.4,0.5),mu_c=0.3,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=1000,sigma=1,beta=0.1,sim_num=100)
#two_stg_adapt_power_simu_normal(mu_t=c(0.4,0.5),mu_c=0.3,s=c(1/3,1),bry_type=2,alpha_2=c(0.001,0.049),theta_cut=0.05,samsz=100,N_max=1000,sigma=1,beta=0.1,sim_num=100)
#two_stg_adapt_power_simu_normal(mu_t=c(0,0.4,0.5),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=500,beta=0.1,sim_num=1000000)
#two_stg_adapt_power_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=500,beta=0.1,sim_num=1000000)
#two_stg_adapt_power_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=500,beta=0.1,sim_num=1000000)

#two_stg_adapt_power_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,1),bry_type=2,alpha_2=c(0.001,0.049),theta_cut=0.05,samsz=100,N_max=500,beta=0.1,sim_num=1000000)
#two_stg_adapt_power_simu_normal(mu_t=c(0,0.1,0.2),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=100,beta=0.1,sim_num=10000)
#two_stg_adapt_power_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=500,beta=0.1,sim_num=100000)

####################################

############################################
## Operating characteristics simulations ###
############################################
Two_stage_adapt_rej_OC_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist){

  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better

  N=ceiling(samsz*s)
  comp=length(mu_t)
  mu=mu_t-mu_c
  # N is the vetoc of sample sizes at which interim analysis will be conducted
  T<- samsz/(sigma_t^2+sigma_c^2)
  t=T*s
  #t is the vector of information times at the interim analysis
  N1=N[1]
  #N1 is the sample at the end of phase 2 stage
  NF=N[length(s)]
  ######################## exit boundaries
  # e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined

  ############ the simulation loops #
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=N[length(s)]
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    control_sampsz_simu=sampsz_simu; tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
    futi=0
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################

      x0=stats::rnorm(N_max,mu_c,sigma_c)
      x=matrix(stats::rnorm(N_max*comp, mu_t,sigma_t), nrow=comp,ncol=N_max)

      #x matirx does not include the control sample.
      #The row number is the same as the active dose

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################
      S_ph2=NULL

      if (direction==1){ #larger mu is better
        for (j in 1:length(mu_t)){
          S_ph2[j]=(mean(x[j,1:N1])-mean(x0[1:N1]))*t[1]
        }
      }else if (direction==2){ #smaller mu is better
        for (j in 1:length(mu_t)){
          S_ph2[j]=(mean(x0[1:N1])-mean(x[j,1:N1]))*t[1]
        }
      }

      wald_ph2=S_ph2/sqrt(t[1])
      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(max(S_ph2),S_ph2)
      #s_slct identifies which dose was selected for the i-th simulation
      ###########################################################
      ### Step 4: generate  random samples for phase 3 stage ####
      ###########################################################
      #  x0_ph3=x0[(N1+1):NF]
      # x_slct_ph3=x[s_slct,(N1+1):NF]
      #y are the samples for the selected dose for the phase 3 stage
      x0_two_stg=x0[1:NF]
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_slct_two_stg=x[s_slct,1:NF]
      #x_two_stg contains all the random samples(including samples from the phase 2 stage) for the seected dose

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      theta_hat=NULL
      for (j in 1:length(s) ){
        theta_hat[j]=mean(x_slct_two_stg[1:N[j]])-mean(x0_two_stg[1:N[j]])  # theta_hat at each interim analysis
      }
      theta_hat_sd=1/sqrt(t)
      wald_two_stg=theta_hat*sqrt(t)
      ######  calculating the statistics for phase 3 only #######

      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test_1=test[-length(s)]  ## test statistics excluding the final analysis

      ##########################################
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries
      ############################################

      if (isTRUE(theta_hat[length(s)-1]<theta_cut)) {futi=1} else {futi=0}

      if (isTRUE(max(test)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
      }

      #If the NUll hypothesis is rejected at up to length(s)-1, terminate trial and estimate CI ####
      if (isTRUE(max(test_1)>=0)){
        #  if ( max(tes_1)>0){ #null hypothesis rejected at the next to last interim analysis. no more random sampling
        rej=1

        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        #  theta_hat_last=theta_hat[last_vt]
        #  theta_hat_last_sd=1/sqrt(t[last_vt])

        ###########################
        sampsz_simu=N[last_vt]
        control_sampsz_simu=sampsz_simu
        tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu


      }else if (isTRUE(max(test_1)<0)){ ## (theta_p3_hat[length(s_1)-1]<theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation


        if (futi==1){ #terminate trial ath the final analysis due to futility. No sample size change
          if (isTRUE(test[length(s)]>0)){
            rej=1
            #The null hypothesis is rejected at the final analysis
            stg=length(s)
          }else{
            rej=0
            stg=0
          }
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=1/sqrt(t[last_vt])

          ###########################
          sampsz_simu=N[last_vt]
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu


        }else if (futi==0){ # no futility ter,mination. proceed to sample size re-estimation
          inter_vt=length(s)-1

          theta=theta_hat[inter_vt]

          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=1/sqrt(t[inter_vt])
          score_inter=theta_hat_inter*t[inter_vt]

          #c_bry=c_bry_two_stg[length(s)]

          #    cond_power=as.numeric(asd_c_power_theta(theta,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd))
          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))

          if (isTRUE(cond_power>=1-beta)){

            #stg identifies the interim analysis at which the null hypothesis is rejected
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=1/sqrt(t[last_vt])

            sampsz_simu=N[last_vt]
            control_sampsz_simu=sampsz_simu
            tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }
          }else  if (isTRUE(cond_power<1-beta)){

            ## proceed to sample size re-estimation #####

            ############ Sample size re-estimation with respect toe the combined phase 2/3 test ##

            T_new_0=theta_hat_inter^(-2)*((c_bry_two_stg[length(s)]*sqrt(T)-score_inter)/sqrt(T-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_sampsz=ceiling(T_new_0*(sigma_t^2+sigma_c^2))
            T_new=new_sampsz/(sigma_t^2+sigma_c^2)
            e_bry_new_two_stg=
              (sqrt(T_new-t[inter_vt])/sqrt(T-t[inter_vt]))*(c_bry_two_stg[length(s)]*sqrt(T)-score_inter)+ score_inter
            #new exit boundary for the combined phase 2/3 test. the NULL is rejected if the score crossed the exit boundary

            c_bry_new_two_stg=e_bry_new_two_stg/sqrt(T_new)  #new critical  boundary for the combined phase 2/3 test. the NULL is rejected if the Wald statistics crossed the critical boundary

            ###############################################

            if (isTRUE(new_sampsz>N_max)){ ## new sample size is too large, the sample size will not change. Trial ends at the planned sample size
              if (isTRUE(test[length(test)]>0)){ ### the NULL hypothesis is rejected at the final anbalysis, with the planned sample size
                rej=1
                stg=length(s)
              }else{### the NULL hypothesis is NOT rejected at the final anbalysis, with the planned sample size
                rej=0
                stg=0
              }
              ###################### Statistics for estimates ####
              ### combined phase 2/3 data ####
              last_vt=length(s)

              sampsz_simu=N[last_vt]
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              theta_hat_last=theta_hat[last_vt] ## theta_hat for the combined phase 2/3 data. Used for estimating the lower bound of the CI
              theta_hat_last_sd=1/sqrt(T)


              #################
              ####################### this completes the estimation when there is no saple size change ########

            } else if (isTRUE(new_sampsz<= N_max)) {## new sample size does not exceed the pre-selected maximum. Will change the sample  size to this new sample size
              #here

              sampsz_simu=new_sampsz
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              x0_two_stg_new=x0[1:new_sampsz] # control arm [combined phase 2/3] data
              x_slct_two_stg_new=x[s_slct,1:new_sampsz]  # data of the selected dose arm

              ### statistics from combined data ####
              if (direction==1){ #larger mu is better
                theta_hat_last_ad=(mean(x_slct_two_stg_new)-mean(x0_two_stg_new)) # theta_hat at the adjusted final analysis
              } else if (direction==2){ #smaller  mu is better
                theta_hat_last_ad=mean(x0_two_stg_new)-mean(x_slct_two_stg_new) # theta_hat at the adjusted final analysis
              }
              wald_last_ad=theta_hat_last_ad*sqrt(T_new)  # Wald statistics at the adjusted final analysis
              score_last_ad=theta_hat_last_ad*T_new  # score statistics at the adjusted final analysis

              theta_hat_last_ad_sd=1/sqrt(T_new)

              ##########################################

              test_add=wald_last_ad-c_bry_new_two_stg  # compare the final Wald statistics with the new adjusted critical boundary

              if (isTRUE(test_add>0)){ # the new critical boundary is crossed, the NULL is rejected

                rej=1
                stg=length(s)

              }else{  # the new critical boundary is NOT crossed, the NULL is rejected

                rej=0
                stg=0

              }

              ###################### Parameter estimate with sample size change ##################

            }
          }
        }
      }
      out=c(rej_nossr,rej,stg,s_slct,control_sampsz_simu,tot_sampsz_simu,futi)

  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  control_samsz_simu = as.numeric(simr[,5])
  tot_samsz_simu = as.numeric(simr[,6])
  futi=as.numeric(simr[,7])

  #############termination timepoint
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  #rej_rt is the rejection rate

  rej_rt_indiv=rep(0,length(mu))
  for (j in 1:length(mu)){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/sim_num
  }

  samsz_icr=NULL
  for (i in 1:length(control_samsz_simu)){
    if (control_samsz_simu[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }

  samsz_icr_rt=mean(samsz_icr)

  futi_rt=mean(futi)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= control_samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  #############################
  ### Step 9: Final output ####
  #############################


  tyI_er=sum(rej[mu[s_slct]==0])/sim_num
  #  tyI_er=sum(rej_0)/sim_num

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/sim_num
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(mu_t)-1)+2*samsz
  sim_results=c(rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er, non_tyI_er_rej_rt,futi_rt,samsz,
                N_max,mean_control_samsz,tot_samsz,mean_tot_samsz, mean_samsz_simu_control_samsz_icr)
  names(sim_results)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned control sample size","Maximum control sample size","Mean control sample size",
                       "Planned total sample size","Mean total sample size","mean control sample size conditional on increase")  #############################
  ### Step 9: Final output ####
  #############################

  out=sim_results
  return(out)

}
#################################################

#Two_stage_adapt_rej_simu_normal(mu_t,mu_c,sigma_t,sigma_c,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist)

#Two_stage_adapt_rej_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000,dist=1)
#Two_stage_adapt_rej_OC_simu_normal(mu_t=c(0,0,0),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.605339, 3.505080, 2.240072),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=10000,dist=1)
#Two_stage_adapt_rej_OC_simu_normal(mu_t=c(0,0.1,0.2),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.605339, 3.505080, 2.240072),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000,dist=1)

#################

two_stg_adapt_power_OC_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,theta_cut,samsz,N_max,beta,sim_num){
  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better

  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
 # mu=mu_t-mu_c
  comp=length(mu_t)
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))

  simu_out=Two_stage_adapt_rej_OC_simu_normal(mu_t=mu_t,mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,
                                           s=s,c_bry_two_stg=c_bry_two_stg,alpha=alpha,beta=beta,theta_cut=theta_cut,samsz=samsz,N_max=N_max,sim_num=sim_num,dist=1)

  mu_t_names=paste(paste(LABEL()$mu_t,"[",c(1:length(mu_t)),"]",sep=""))
  out_1= round(c(mu_t,mu_c,as.numeric(simu_out)),4)
   names(out_1)=c(mu_t_names,LABEL()$mu_c,names(simu_out))
   out=out_1
  return(out)
}
#############
#two_stg_adapt_power_OC_simu_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,theta_cut,samsz,N_max,beta,sim_num)
#two_stg_adapt_power_OC_simu_normal(mu_t=c(0.4,0.5),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=1000,beta=0.1,sim_num=100)
#two_stg_adapt_power_OC_simu_normal(mu_t=c(0,0.2,0.3),mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=1000,beta=0.1,sim_num=10000)

####################


###################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Normal distribution, Optimizing adaptive sequential design
#'
#' @param groups Matrix of assumed effect sizes for the test arms under study: one row per candidate test-arm scenario and one column per dose, so that the number of doses is the number of columns
#' @param mu_c Mean of the control arm
#' @param sigma_t Standard deviation of the test arms
#' @param sigma_c Standard deviation of the control arm
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s1 Vector of information fractions of the first interim analysis schedule; the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule; the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s1; each element must be > 0 and the elements must sum to alpha_2
#' @param alpha_2_alps_s2 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s2; each element must be > 0 and the elements must sum to alpha_2
#' @param theta_cut Lack-of-efficacy threshold: an effect size below which the test drug is not considered to provide clinically meaningful efficacy
#' @param samsz_1 Planned sample size of the control arm for the first sample size scenario
#' @param samsz_2 Planned sample size of the control arm for the second sample size scenario
#' @param N_max_1 Maximum sample size (control arm) allowed after sample size re-estimation for the first maximum sample size scenario
#' @param N_max_2 Maximum sample size (control arm) allowed after sample size re-estimation for the second maximum sample size scenario
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list of four elements, one per boundary type and interim analysis schedule, with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with O'Brien-Fleming boundaries}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with O'Brien-Fleming boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with alpha-spending boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with alpha-spending boundaries}
#' }
#' Each of these four lists has the components "Group sequential design" (data frame with the columns "Looks", "Information fractions" and the critical boundaries, plus "Alpha-spending" and "Alpha-spending boundaries" for the alpha-spending schedules) and "Simulation results: Information fractions ..., Sample size ..., N_max ..." (the operating characteristics of the corresponding adaptive design).
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' Two_stage_asd_simu_OC_normal(groups = rbind(c(0.2, 0.3, 0.4), c(0.3, 0.4, 0.5)), mu_c = 0,
#' sigma_t = 1, sigma_c = 1, direction = 1, s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05,
#' alpha_2_alps_s1 = c(0.001,0.001,0.048), alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01,
#' samsz_1 = 100, samsz_2 = 150, N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
Two_stage_asd_simu_OC_normal<-function(groups,mu_c,sigma_t,sigma_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
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

  s1_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=dist)[1]))
  s2_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=dist)[1]))

  s1_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s1,alpha=alpha_alsp_s1,comp=comp,dist=dist)[1]))
  s2_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s2,alpha=alpha_alsp_s2,comp=comp,dist=dist)[1]))

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
    simu_out_OF_s1_samsz_1_N_max_1=rbind(simu_out_OF_s1_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))

  }
  rownames(simu_out_OF_s1_samsz_1_N_max_1)=u_t_group_names


  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_1=rbind(simu_out_OF_s2_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_1)=u_t_group_names

    ## samsz_2###


  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_1=rbind(simu_out_OF_s1_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_1)=u_t_group_names


  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_1=rbind(simu_out_OF_s2_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
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
    simu_out_OF_s1_samsz_1_N_max_2=rbind(simu_out_OF_s1_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_2)=u_t_group_names


  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_2=rbind(simu_out_OF_s2_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_2)=u_t_group_names

  ## samsz_2###


  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_2=rbind(simu_out_OF_s1_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_2)=u_t_group_names


  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_2=rbind(simu_out_OF_s2_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
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
    simu_out_alsp_s1_samsz_1_N_max_1=rbind(simu_out_alsp_s1_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=2,alpha_2=alpha_alsp_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_1)=u_t_group_names


  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_1=rbind(simu_out_alsp_s2_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                           bry_type=2,alpha_2=alpha_alsp_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_1)=u_t_group_names

  ## samsz_2###


  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_1=rbind(simu_out_alsp_s1_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=2,alpha_2=alpha_alsp_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_1)=u_t_group_names


  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_1=rbind(simu_out_alsp_s2_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
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
    simu_out_alsp_s1_samsz_1_N_max_2=rbind(simu_out_alsp_s1_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=2,alpha_2=alpha_alsp_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_2)=u_t_group_names


  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_2=rbind(simu_out_alsp_s2_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                                                           bry_type=2,alpha_2=alpha_alsp_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_2)=u_t_group_names

  ## samsz_2###


  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_2=rbind(simu_out_alsp_s1_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                                                           bry_type=2,alpha_2=alpha_alsp_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_2)=u_t_group_names


  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(mu_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_2=rbind(simu_out_alsp_s2_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_normal(mu_t=mu_t_vec[i,],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
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
#Two_stage_asd_simu_OC_normal(mu_t_vec_1,mu_t_vec_2,mu_c,sigma_t,sigma_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num,dist)
#Two_stage_asd_simu_OC_normal(mu_t_vec_1=c(0.3,0.4,0.5),mu_t_vec_2=c(0.1,0.2,0.3),mu_t_vec_3=c(0.3,0.4,0.5),mu_t_vec_4=c(0.1,0.2,0.3),mu_t_vec_5=NULL,mu_c=0,sigma_t=1,sigma_c=1,direction=1,
#                             s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                             theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1, sim_num=100)
#x=Two_stage_asd_simu_OC_normal(groups=rbind(c(0.2,0.3,0.4),c(0.3,0.4,0.5),c(0.4,0.5,0.6),c(0.5,0.6,0.7)),mu_c=0,sigma_t=1,sigma_c=1,direction=1,
#                                                          s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                                                        theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1, sim_num=1000)
#print(x)

#######################################################################
### II.2 .Simulations  -survival analysis                         #####
#######################################################################


##########################################################################
### II.2.1 .Simulations with planned TSSSD -survival analysis        #####
##          Type I error, power simulaitons             #####
#### This will be included in the software                           #####
##########################################################################

#######################################

Two_stage_rej_simu_survival<-function(HR,lambda,s,c_bry_two_stg,evt_num,samp_multip,surv_rate,beta,sim_num){

  #######################################


  comp=length(HR)
  ##########################
  evt_num_s=ceiling(evt_num*s)
  sampsz=ceiling(samp_multip*evt_num/surv_rate)

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=evt_num_s[length(s)]
    rej=0; stg=0
    s_slct=0; sampsz_simu=evt_num_s[length(s)]
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################


      ###########################################
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
      #  grp_x=c(rep(0,length(time_x))) #
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
      evt_count=matrix(rep(0,sampsz*comp),nrow=comp,ncol=sampsz)
      # for (j in (1:comp)){
      #   for (i in 1:sampsz){
      #     evt_count[j,i]=sum(x0_status[1:i])+sum(x_status[j,1:i]) #number of events after 2*i observations
      #   }
      # }
      evt_count_x0 = rep(0, sampsz)
      total_count_x0 = 0
      for (i in 1:sampsz){
        total_count_x0 = total_count_x0 + x0_status[i]
        evt_count_x0[i] = total_count_x0
      }
      for (j in (1:comp)){
        total_count_x = 0
        for (i in 1:sampsz){
          total_count_x = total_count_x + x_status[j,i]
          evt_count[j,i] = evt_count_x0[i] + total_count_x
        }
      }

      #  evt_count[,1:100]
      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      ###############################set up matrices fro parameters########################################################

      theta_hat_0=rep(0,comp)
      theta_hat_0_sd=rep(0,comp)
      wald_0=rep(0,comp)
      score_0=rep(0,comp)
      t_0=rep(0,comp)
      ##### Obtaningin parameters and store in the matrices ##############
      evt_cut_0=rep(0,comp)

      for (j in (1:comp)) {
        evt_pos=evt_count[j,evt_count[j,]>evt_num_s[1]]#substring of  evt_count that are larger than evt_num
        evt_cut_0[j]=match(evt_pos[1],evt_count[j,]) #the first incidence that evt_count is larger than evt_num
      }


      evt_cut=min(evt_cut_0)
      ## the required number of events may be unreachable within the available
      ## follow-up: the replicate is marked non-evaluable instead of using an NA index
      if (!is.finite(evt_cut)) evt_cut=length(evt_count[1,])

      for (j in (1:comp)) {
        #  vector of zeros, length of the vector is the same as the original smaple size
        time=c(time_x0[1:evt_cut],time_x[j,1:evt_cut]) #combined samples befroe total number of events exceed evt_num
        event=c(x0_status[1:evt_cut],x_status[j,1:evt_cut]) #censoring status before total number of events exceed evt_num
        grp=c(grp_x0[1:evt_cut],grp_x[j,1:evt_cut]) #vector of zero, with length=2*evt_cut


        test= list(time=time,status=event,grp)#survival data
        coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
        # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
        # test is the same as data=test
        theta_hat_0[j]= coxcf$coefficients[1] #estimated log(HR)
        theta_hat_0_sd[j]= coxcf$coefficients[3] #estimated standard error of log(HR)
        t_0[j]=(theta_hat_0_sd[j])^(-2)
        wald_0[j]=theta_hat_0[j]/ theta_hat_0_sd[j]
        score_0[j]= theta_hat_0[j]*t_0[j]
      }

      ## the required events for the phase-2 look were unreachable: blank the
      ## statistics so the replicate is reported as non-evaluable
      if (anyNA(evt_cut_0)) { theta_hat_0[]=NA_real_; theta_hat_0_sd[]=NA_real_; t_0[]=NA_real_; wald_0[]=NA_real_; score_0[]=NA_real_ }
      wald_p2_max=max(wald_0)
      test_1=wald_p2_max-c_bry_two_stg[1]
      ####################################

      #############################################################################################
      ### Step 3: the dose with the largest wald_0 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(wald_p2_max,wald_0)
      #s_slct identifies which dose was selected for the i-th simulation
      #############################################################
      ### Step 4: Identify  test statistics  for phase 3 stage ####
      #############################################################
      theta_hat=rep(0,length(s))
      theta_hat_sd=rep(0,length(s))
      wald=rep(0,length(s))
      score=rep(0,length(s))
      t=rep(0,length(s))
      #################################
      evt_count_slct=evt_count[s_slct,]
      time_x_slct=time_x[s_slct,]
      x_status_slct=x_status[s_slct,]
      grp_x_slct=grp_x[s_slct,]
      #################################
      for (i in (1: length(s))) {
        #  vector of zeros, length of the vector is the same as the original smaple size
        evt_pos=evt_count[s_slct,evt_count_slct>evt_num_s[i]]#substring of  evt_count that are larger than evt_num
        evt_cut=match(evt_pos[1],evt_count_slct) #the first incidence that evt_count is larger than evt_num
        ## events required at this look are unreachable within follow-up:/n        ## the replicate is left non-evaluable instead of indexing with NA
        if (is.na(evt_cut)) { evt_cut=length(evt_count_slct); theta_hat[i]=NA_real_; theta_hat_sd[i]=NA_real_; t[i]=NA_real_; wald[i]=NA_real_; score[i]=NA_real_; next }
        time=c(time_x0[1:evt_cut],time_x_slct[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
        event=c(x0_status[1:evt_cut],x_status_slct[1:evt_cut]) #censoring status before total number of events exceed evt_num
        grp=c(grp_x0[1:evt_cut],grp_x_slct[1:evt_cut]) #vector of zero, with length=2*evt_cut


        test= list(time=time,status=event,grp)#survival data
        coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
        # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
        # test is the same as data=test
        theta_hat[i]= coxcf$coefficients[1] #estimated log(HR)
        theta_hat_sd[i]= coxcf$coefficients[3] #estimated standard error of log(HR)
        t[i]=(theta_hat_sd[i])^(-2)
        wald[i]=theta_hat[i]/ theta_hat_sd[i]
        score[i]= theta_hat[i]*t[i]
      }

      #######################################

      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald-c_bry_two_stg
      # test compares the Wald statistics with the critical boundaries

      ##########################################################################################
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries
      ##########################################################################################

      if (isTRUE(max(test)>0)) {

        rej=1
        #If max(tes_1)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        sampsz_simu=evt_num_s[last_vt]



      }else{

        rej=0
        stg=0
        last_vt=length(s)

        sampsz_simu=evt_num_s[last_vt]

      }



      #######################################################
      ### Step 5: output and finish the simjulation loop ####
      #######################################################
      out=c(rej,stg,s_slct,sampsz_simu)
      out
      #      out=c(rej,stg,sampsz_simu,theta_ci_dif)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################

  rej = as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  s_slct = as.numeric(simr[,3])
  sampsz_simu = as.numeric(simr[,4])

  ###################################
  sampsz_simu_two_arm = as.numeric(simr[,4])
  sampsz_simu_tot=(evt_num_s[1]/2)*(comp-1)+sampsz_simu_two_arm
  mean_two_arm_evt=mean(sampsz_simu_two_arm)
  mean_tot_evt_num=mean(sampsz_simu_tot)

  #############termination timepoint####
  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  ################# Overall ci results######

  rej_rt=mean(rej)
  ##############################

  #####################################
  slct_HR=rep(0,length(HR))
  for (j in 1:length(HR)){
    slct_HR[j]=length( s_slct[ s_slct==j])/(dim(simr)[1])
  }


  rej_HR=rep(0,length(HR))
  for (j in 1:length(HR)){
    rej_HR[j]=sum( rej[ s_slct==j])/(dim(simr)[1])
  }

  cond_rej_HR=rep(0,length(HR))
  for (j in 1:length(HR)){
    if (slct_HR[j]==0){
      cond_rej_HR[j]=NA
    }else {
      cond_rej_HR[j]=rej_HR[j]/slct_HR[j]
    }
  }

  ##################################

  mean_two_arm_evt=mean(sampsz_simu_two_arm)
  mean_tot_evt_num=mean(sampsz_simu_tot)

  results_by_HR=cbind(HR,slct_HR,rej_HR,cond_rej_HR)

  colnames( results_by_HR)=c("HR","selection rate","rejection rate by each HR","conditional rejection rate")
#  out=list( rej_rt,rej_by_HR,tyI_er,power, rej_time, mean_two_arm_evt,mean_tot_evt_num,c_bry_two_stg)
#  names(out)=c( "rejection rate","rej_by_HR","tyI_er","pos_power", "rej_time","mean_two_arm_evt","mean_tot_evt_num","c_bry")


seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time)
colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look")

tyI_er=sum(rej_HR[HR==1])

non_tyI_er_rej_rt=sum(rej[HR[s_slct]<1])/(dim(simr)[1])
mean_two_arm_evt=mean(sampsz_simu_two_arm)
mean_tot_evt_num=mean(sampsz_simu_tot)
evt_tot=(evt_num_s[1]/2)*(comp-1)+evt_num
sim_results=c(rej_rt,tyI_er, non_tyI_er_rej_rt,evt_num,mean_two_arm_evt,evt_tot,mean_tot_evt_num)
names(sim_results)=c("rejection rate","Type I error","Non-Type I error rejection rate",
                     "Planned event number for two arms","Mean event number for two arms","Planned total number of events","total number of events")
#c( "rejection rate","rej_by_HR","tyI_er","pos_power", "rej_time","mean_two_arm_evt","mean_tot_evt_num","c_bry")

#############################
### Step 9: Final output ####
#############################

out=list(seq_design,results_by_HR,sim_results)
names(out)=c("Sequential design","Results by arm","Simulation results")
return(out)

}

#################################################
#Two_stage_rej_simu_survival(HR,lambda,s,c_bry_two_stg,evt_num,samp_multip,surv_rate,sim_num)
#Two_stage_rej_simu_survival(HR=rep(1,2),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=150,samp_multip=5,surv_rate=0.5,sim_num=1000)
#Two_stage_rej_simu_survival(HR=c(0.8,0.8),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=700,samp_multip=5,surv_rate=0.5,sim_num=100)
#Two_stage_rej_simu_survival(HR=c(0.3,0.8),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=500,samp_multip=5,surv_rate=0.5,sim_num=100)
#Two_stage_rej_simu_survival(HR=c(1,0.95,0.8),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=100,samp_multip=5,surv_rate=0.5,sim_num=1000)

#######################################


#' Phase 2/3 seamless combination Group sequential design: Simulations, Survival analysis, Type I error and power
#'
#' @param HR Vector of hazard ratios of each test arm vs. the control arm, one entry per dose; HR < 1 indicates a beneficial effect
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param surv_rate Survival rate at the end of follow-up, used to translate the number of events into a sample size
#' @param evt_num Number of events between a test arm and the control arm at the final analysis, used as the planned number of events
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries" and "rejection by look"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "HR", "selection rate", "rejection rate by each HR" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the overall "rejection rate", the "Type I error", the "Non-Type I error rejection rate", the "Planned event number for two arms", the "Mean event number for two arms", the "Planned total number of events" and the "total number of events"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_power_simu_survival(HR = c(0.8,0.8), s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05,
#' surv_rate = 0.5, evt_num = 800, sim_num = 20)
#'
#' @export
two_stg_power_simu_survival<-function(HR,s,bry_type,alpha_2,surv_rate,evt_num,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  comp=length(HR)
  lambda=1
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s,alpha,comp,dist=2)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s,alpha,comp,dist=2)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))
  samp_multip=5
  simu_out=Two_stage_rej_simu_survival(HR=HR,lambda=1,s=s,c_bry_two_stg=c_bry_two_stg,evt_num=evt_num,samp_multip=samp_multip,surv_rate=surv_rate,sim_num=sim_num)

#  evt_tot_planned=ceiling((evt_num*s[1]/2)*(comp+1)+evt_num*(s[length(s)]-s[1]))
  out=simu_out
 # names(out)=c("simulation summary","planned number of events","planned total number of events")
  return(out)
}

##################################
#two_stg_power_simu_survival<-function(HR,s,bry_type,alpha_2,surv_rate,evt_num,sim_num)
#two_stg_power_simu_survival(HR=c(0.8,0.8),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=800,sim_num=100)
#two_stg_power_simu_survival(HR=c(1,1,1),s=c(1/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=800,sim_num=100)
#two_stg_power_simu_survival(HR=c(1,1,1),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=800,sim_num=100)
#two_stg_power_simu_survival(HR=c(1,1,0.8),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=800,sim_num=1000)

##############################################

########################################################################

##############################################


##########################################################################
### II.2.2 .Simulations with adaptive  ATSSSD -survival analysis      #####
##          Type I error, power, simulaitons                         #####
#### This will be included in the software                           #####
##########################################################################

#######################################

Two_stage_adapt_rej_simu_survival<-function(HR,lambda,s,c_bry_two_stg,evt_num,evt_max,HR_cut,samp_multip,surv_rate,beta,sim_num){
  #######################################


  comp=length(HR)
  ##########################
  evt_num_s=c(ceiling(evt_num*s),evt_max)
  samsz=ceiling(samp_multip*evt_max/surv_rate)
  theta_cut=-log(HR_cut)
  ######################## exit boundaries
  #e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    samsz_simu=evt_num_s[length(s)]; futi=0
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################

      ###########################################
      ### geting random sample from exponential distributions ################
      x0=stats::rexp(samsz,lambda) #generate exponential data for the control arm
      x=matrix(stats::rexp(samsz*comp, HR*lambda), nrow=comp,ncol=samsz)
      T=-log(surv_rate)/lambda  #Censoring time. Events after T will be censored
      time_x0=x0
      time_x=x
      x0_status=c(rep(1,length(x0)))
      x_status=matrix(rep(1,samsz*comp),nrow=comp,ncol=samsz)

      ######################## Censored data for the control arm #####################################
      for (i in (1:samsz)){
        if (x0[i]>T){
          time_x0[i]=T  #x0[i] is censored after T
          x0_status[i]=0   #censoring indicator
        }
      }
      grp_x0=c(rep(1,length(time_x0)))#vector of zeros, length of the vector is the same as the original smaple size

      ######################## Censored data for the test  arms #####################################
      #  grp_x=c(rep(0,length(time_x))) #
      grp_x=matrix(rep(0,samsz*comp),nrow=comp,ncol=samsz)
      for (j in (1:comp)){
        for (i in (1:samsz)){
          if (x[j,i]>T){
            time_x[j,i]=T #x[i] is censored after T
            x_status[j,i]=0  #censoring indicator
          }
        }
      }

      #### count the numnber of events in the control and test arms #######
      evt_count=matrix(rep(0,samsz*comp),nrow=comp,ncol=samsz)
      # for (j in (1:comp)){
      #   for (i in 1:samsz){
      #     evt_count[j,i]=sum(x0_status[1:i])+sum(x_status[j,1:i]) #number of events after 2*i observations
      #   }
      # }
      evt_count_x0 = rep(0, samsz)
      total_count_x0 = 0
      for (i in 1:samsz){
        total_count_x0 = total_count_x0 + x0_status[i]
        evt_count_x0[i] = total_count_x0
      }
      for (j in (1:comp)){
        total_count_x = 0
        for (i in 1:samsz){
          total_count_x = total_count_x + x_status[j,i]
          evt_count[j,i] = evt_count_x0[i] + total_count_x
        }
      }

      #  evt_count[,1:100]
      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      ###############################set up matrices fro parameters########################################################

      theta_hat_0=rep(0,comp)
      theta_hat_0_sd=rep(0,comp)
      wald_0=rep(0,comp)
      score_0=rep(0,comp)
      t_0=rep(0,comp)
      ##### Obtaningin parameters and store in the matrices ##############
      evt_cut_0=rep(0,comp)

      for (j in (1:comp)) {
        evt_pos=evt_count[j,evt_count[j,]>evt_num_s[1]]#substring of  evt_count that are larger than evt_num_s[1]
        evt_cut_0[j]=match(evt_pos[1],evt_count[j,]) #the first incidence that evt_count is larger than evt_num
      }


      evt_cut=min(evt_cut_0)
      ## the required number of events may be unreachable within the available
      ## follow-up: the replicate is marked non-evaluable instead of using an NA index
      if (!is.finite(evt_cut)) evt_cut=length(evt_count[1,])

      for (j in (1:comp)) {
        #  vector of zeros, length of the vector is the same as the original smaple size
        time=c(time_x0[1:evt_cut],time_x[j,1:evt_cut]) #combined samples befroe total number of events exceed evt_num
        event=c(x0_status[1:evt_cut],x_status[j,1:evt_cut]) #censoring status before total number of events exceed evt_num
        grp=c(grp_x0[1:evt_cut],grp_x[j,1:evt_cut]) #vector of zero, with length=2*evt_cut


        test= list(time=time,status=event,grp)#survival data
        coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
        # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
        # test is the same as data=test
        theta_hat_0[j]= coxcf$coefficients[1] #estimated log(HR)
        theta_hat_0_sd[j]= coxcf$coefficients[3] #estimated standard error of log(HR)
        t_0[j]=(theta_hat_0_sd[j])^(-2)
        wald_0[j]=theta_hat_0[j]/ theta_hat_0_sd[j]
        score_0[j]= theta_hat_0[j]*t_0[j]
      }

      ## the required events for the phase-2 look were unreachable: blank the
      ## statistics so the replicate is reported as non-evaluable
      if (anyNA(evt_cut_0)) { theta_hat_0[]=NA_real_; theta_hat_0_sd[]=NA_real_; t_0[]=NA_real_; wald_0[]=NA_real_; score_0[]=NA_real_ }
      wald_p2_max=max(wald_0)
      test_1=wald_p2_max-c_bry_two_stg[1]
      ####################################

      #############################################################################################
      ### Step 3: the dose with the largest wald_0 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(wald_p2_max,wald_0)
      #s_slct identifies which dose was selected for the i-th simulation
      #############################################################
      ### Step 4: Identify  test statistics  for phase 3 stage ####
      #############################################################
      theta_hat_two_stg=rep(0,length(s))
      theta_hat_two_stg_sd=rep(0,length(s))
      wald=rep(0,length(s))
      score=rep(0,length(s))
      t=rep(0,length(s))
      #################################
      evt_count_slct=evt_count[s_slct,]
      time_x_slct=time_x[s_slct,]
      x_status_slct=x_status[s_slct,]
      grp_x_slct=grp_x[s_slct,]
      #################################
      for (i in (1: length(s))) {
        #  vector of zeros, length of the vector is the same as the original smaple size
        evt_pos=evt_count[s_slct,evt_count_slct>evt_num_s[i]]#substring of  evt_count that are larger than evt_num
        # evt_pos=evt_count[s_slct,evt_count_slct>evt_num_s[i+1]]#substring of  evt_count that are larger than evt_num

        evt_cut=match(evt_pos[1],evt_count_slct) #the first incidence that evt_count is larger than evt_num
        ## events required at this look are unreachable within follow-up:/n        ## the replicate is left non-evaluable instead of indexing with NA
        if (is.na(evt_cut)) { evt_cut=length(evt_count_slct); theta_hat_two_stg[i]=NA_real_; theta_hat_two_stg_sd[i]=NA_real_; t[i]=NA_real_; wald[i]=NA_real_; score[i]=NA_real_; next }
        time=c(time_x0[1:evt_cut],time_x_slct[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
        event=c(x0_status[1:evt_cut],x_status_slct[1:evt_cut]) #censoring status before total number of events exceed evt_num
        grp=c(grp_x0[1:evt_cut],grp_x_slct[1:evt_cut]) #vector of zero, with length=2*evt_cut


        test= list(time=time,status=event,grp)#survival data
        coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
        # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
        # test is the same as data=test
        theta_hat_two_stg[i]= coxcf$coefficients[1] #estimated log(HR)
        theta_hat_two_stg_sd[i]= coxcf$coefficients[3] #estimated standard error of log(HR)
        t[i]=(theta_hat_two_stg_sd[i])^(-2)
        wald[i]=theta_hat_two_stg[i]/theta_hat_two_stg_sd[i]
        score[i]= theta_hat_two_stg[i]*t[i]
      }

      #######################################

      evt_count_two_stg=evt_count[s_slct,]

      t_two_stg=t
      x_status_two_stg=x_status[s_slct,]
      wald_two_stg=wald
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg[1:length(s)]-c_bry_two_stg
      test_2=test[1:(length(s)-1)]
      # test compares the Wald statistics with the critical boundaries

      ##########################################################################################
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries
      ##########################################################################################

      inter_vt=length(s)-1
      if (isTRUE(theta_hat_two_stg[inter_vt]<theta_cut)) {futi=1} else {futi=0}
##########################################

      if (isTRUE(max(test)>=0)){
        rej_nossr=1
        test_p=test[test>0]
        stg_1=match(test_p[1],test)
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
        stg_1=0
      }

      ## (theta_hat[inter_vt]<theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation
      if (isTRUE(max(test_2)>0)) {

        rej=1
        #If max(tes_1)>0, te nulll hypothesis is rejected
        test_t=test_2[test_2>0]
        stg=match(test_t[1],test_2)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        theta_hat_last=theta_hat_two_stg[last_vt]
        theta_hat_last_sd=1/sqrt(t_two_stg[last_vt])

        samsz_simu=evt_num_s[last_vt]

      }else if (isTRUE(max(test_2)<=0)) {


        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation

          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat_two_stg[last_vt]
          theta_hat_last_sd=1/sqrt(t_two_stg[last_vt])


          samsz_simu=evt_num_s[length(s)]
          if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat_two_stg[inter_vt]
          theta_hat_inter_sd=theta_hat_two_stg_sd[inter_vt]
          t_hat_inter=theta_hat_inter_sd^(-2)
          score_inter=score[length(s)-1]
          event_inter=evt_num_s[inter_vt]
          t_final_0=t_hat_inter*(evt_num_s[length(s)]/evt_num_s[inter_vt])

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))


          if (isTRUE(cond_power>=1-beta)){

            last_vt=length(s)
            theta_hat_last=theta_hat_two_stg[last_vt]
            theta_hat_last_sd=theta_hat_two_stg_sd[last_vt]


            samsz_simu=evt_num_s[length(s)]
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (isTRUE(cond_power<1-beta)){
            #  c_bry=c_bry_two_stg
            #sample size re-estimation procedure##
            inter_vt=length(s)-1
            event_inter=evt_num_s[inter_vt]
            T_new_0=theta_hat_inter^(-2)*((c_bry_two_stg[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_samsz=ceiling(event_inter*T_new_0/t[inter_vt])

              if  (new_samsz>evt_max){

              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (isTRUE(test[length(s)]>0)){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              samsz_simu=evt_num_s[length(s)]
            } else  if  (new_samsz<=evt_max){
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              samsz_simu=new_samsz
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              evt_pos=evt_count_two_stg[evt_count_two_stg>new_samsz]#substring of  evt_count that are larger than evt_num
              ## events required after the sample size increase may be unreachable:/n              ## clamp the index so the subset below is valid, then blank the
              ## statistics so the replicate is reported as non-evaluable
              evt_short=is.na(match(evt_pos[1],evt_count_slct))
              evt_cut=if (evt_short) length(evt_count_slct) else match(evt_pos[1],evt_count_slct)

              time=c(time_x0[1:evt_cut],time_x_slct[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
              event=c(x0_status[1:evt_cut],x_status_slct[1:evt_cut]) #censoring status before total number of events exceed evt_num
              grp=c(grp_x0[1:evt_cut],grp_x_slct[1:evt_cut]) #vector of zero, with length=2*evt_cut


              test= list(time=time,status=event,grp)#survival data
              coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
              # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
              theta_hat_last= coxcf$coefficients[1] #estimated -log(HR)
              theta_hat_last_sd= coxcf$coefficients[3] #estimated standard error of log(HR)
              t_last=(theta_hat_last_sd)^(-2)
              T_new=t_last
              wald_last=theta_hat_last/ theta_hat_last_sd
              if (evt_short) { theta_hat_last=NA_real_; theta_hat_last_sd=NA_real_; t_last=NA_real_; wald_last=NA_real_ }
              e_bry_new_two_stg=
                (sqrt(T_new-t[inter_vt])/sqrt(t[length(s)]-t[inter_vt]))*(c_bry_two_stg[length(s)]*sqrt(t[length(s)])-score_inter)+ score_inter

              new_c_bry=e_bry_new_two_stg/sqrt(T_new)  #new critical  boundary for the combined phase 2/3 test. the NULL is rejected if the Wald statistics crossed the critical boundary

              #################################################################################
              ###                   Step 4.6: Hypothesis testing                           ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              test_add=wald_last-new_c_bry

              if (isTRUE(test_add>0)){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }

            }
          }
        }
      }
      #######################################################
      ### Step 5: output and finish the simjulation loop ####
      #######################################################
      out=c(rej_nossr,rej,stg,s_slct,samsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  samsz_simu = as.numeric(simr[,5])
  futi = as.numeric(simr[,6])
  ###################################
  samsz_simu_two_arm = samsz_simu
  samsz_simu_tot=(evt_num_s[1]/2)*(comp+1)+(samsz_simu_two_arm-evt_num_s[1])
  mean_two_arm_evt=mean(samsz_simu_two_arm)
  mean_tot_evt_num=mean(samsz_simu_tot)

  #############termination timepoint####
  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  #############termination timepoint####

  #################  results######
  av_samsz=mean(samsz_simu)
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  futi_rt=mean(futi)

  #####################################
  slct_HR=rep(0,length(HR))
  for (j in 1:length(HR)){
    slct_HR[j]=length( s_slct[ s_slct==j])/(dim(simr)[1])
  }

  rej_HR=rep(0,length(HR))
  for (j in 1:length(HR)){
    rej_HR[j]=sum( rej[ s_slct==j])/(dim(simr)[1])
  }

  cond_rej_HR=rep(0,length(HR))
  for (j in 1:length(HR)){
    if (slct_HR[j]==0){
      cond_rej_HR[j]=NA
    }else {
      cond_rej_HR[j]=rej_HR[j]/slct_HR[j]
    }
  }

  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu[i]>evt_num){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }

  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_rt_byLook=c(rep(0,(length(s)-1)),samsz_icr_rt)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  two_arm_evt_samsz_icr= samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_two_arm_evt_samsz_icr=mean(two_arm_evt_samsz_icr)

  ################################

  mean_two_arm_evt=mean(samsz_simu_two_arm)
  mean_tot_evt_num=mean(samsz_simu_tot)

  results_by_HR=cbind(HR,slct_HR,rej_HR,cond_rej_HR)

  colnames( results_by_HR)=c("HR","selection rate","rejection rate by each HR","conditional rejection rate")
  #  out=list( rej_rt,rej_by_HR,tyI_er,power, rej_time, mean_two_arm_evt,mean_tot_evt_num,c_bry_two_stg)
  #  names(out)=c( "rejection rate","rej_by_HR","tyI_er","pos_power", "rej_time","mean_two_arm_evt","mean_tot_evt_num","c_bry")


  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg, rej_time,samsz_icr_rt_byLook)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look","Sample size increase")

  # tyI_er=sum(rej_HR[HR==1])/(dim(simr)[1])
  tyI_er=sum(rej_HR[HR==1])

  non_tyI_er_rej_rt=sum(rej[HR[s_slct]<1])/(dim(simr)[1])
  mean_two_arm_evt=mean(samsz_simu_two_arm)
  mean_tot_evt_num=mean(samsz_simu_tot)
  evt_tot=(evt_num_s[1]/2)*(comp-1)+evt_num
  sim_results=round(c( rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,
                       tyI_er, non_tyI_er_rej_rt,futi_rt,evt_num,evt_max,
                       mean_two_arm_evt,evt_tot,mean_tot_evt_num, mean_two_arm_evt_samsz_icr),5)
  names(sim_results)=c("rejection rate without SSR","rejection rate with SSR", "Rate of sample size increase",
                       "Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned event number for two arms", "Maximum number of events for two arms",
                       "Mean event number for two arms","Planned total number of events",
                       "total number of events","Mean event number for two arms conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(round(seq_design,4),round(results_by_HR,4),round(sim_results,4))
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)
}


#################################################
#Two_stage_adapt_rej_simu_survival(HR,lambda,s,c_bry_two_stg,evt_num,evt_max,HR_cut,samp_multip,surv_rate,beta,sim_num)
#######################################
#Two_stage_adapt_rej_simu_survival(HR=rep(1,3),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=10000)
#Two_stage_adapt_rej_simu_survival(HR=rep(1,3),lambda=1,s=c(1/3,1),c_bry_two_stg=c(3.852050,  2.223982),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=10000)

 #$Two_stage_adapt_rej_simu_survival(HR=c(0.3,0.8),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=500,evt_max=1500,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)
#Two_stage_adapt_rej_simu_survival(HR=c(0.3,0.8),lambda=1,s=c(1/3,1),c_bry_two_stg=c(3.852050, 2.223982),evt_num=500,evt_max=1500,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)
#Two_stage_adapt_rej_simu_survival(HR=c(0.6,0.6,0.8),lambda=1,s=c(0.3,0.7,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=500,evt_max=1500,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)
#Two_stage_adapt_rej_simu_survival(HR=rep(1,2),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=200,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=1000)
#Two_stage_adapt_rej_simu_survival(HR=c(0.9,1),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=200,evt_max=800,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)
#Two_stage_adapt_rej_simu_survival(HR=c(0.8,1),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=200,evt_max=500,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=10000)
#Two_stage_adapt_rej_simu_survival(HR=c(0.7,1),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=200,evt_max=500,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=1000)
#Two_stage_adapt_rej_simu_survival(HR=rep(1,2),lambda=1,s=c(1/3,1),c_bry_two_stg=c(3.852050,  2.223982),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)
#Two_stage_adapt_rej_simu_survival(HR=rep(1,2),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=150, evt_max=1500,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)

#Two_stage_adapt_rej_simu_survival(HR=c(0.7,1),lambda=1,s=c(1/3,1),c_bry_two_stg=c(3.852050, 2.223982),evt_num=500,evt_max=1500,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)

##############################################################

#####################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Survival analysis, Operating characteristics
#'
#' @param HR Vector of hazard ratios of each test arm vs. the control arm, one entry per dose; must be > 0, with HR < 1 indicating a beneficial effect and HR = 1 the null hypothesis
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param evt_num Planned number of events between a test arm and the control arm at the final analysis
#' @param evt_max Maximum number of events between a test arm and the control arm allowed after sample size re-estimation
#' @param HR_cut Hazard-ratio cut-off for lack of efficacy at the interim analysis: the test drug is not considered to provide clinically meaningful efficacy when the observed hazard ratio is greater than HR_cut
#' @param surv_rate Survival rate at the end of follow-up, used to translate the number of events into a sample size
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned adaptive group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries", "rejection by look" and "Sample size increase"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "HR", "selection rate", "rejection rate by each HR" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the "rejection rate without SSR", the "rejection rate with SSR", the "Rate of sample size increase", the "Rejection rate without SSR/(meeting sample size incerase criteria)", the "Rejection rate with SSR/(sample size incerase)", the "Type I error", the "Non-Type I error rejection rate", the "Futility rate", the "Planned event number for two arms", the "Maximum number of events for two arms", the "Mean event number for two arms", the "Planned total number of events", the "total number of events" and the "Mean event number for two arms conditional on increase"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_adapt_power_simu_survival(HR = rep(1,2), s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05,
#' evt_num = 150, evt_max = 1500, HR_cut = 0.95, surv_rate = 0.5, beta = 0.1, sim_num = 20)
#'
#' @export
two_stg_adapt_power_simu_survival<-function(HR,s,bry_type,alpha_2,evt_num,evt_max,HR_cut,surv_rate,beta,sim_num){
  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2

  comp=length(HR)

  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=2)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=2)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))

  simu_out=Two_stage_adapt_rej_simu_survival(HR=HR,lambda=1,s,c_bry_two_stg,evt_num=evt_num,evt_max=evt_max,HR_cut=HR_cut,samp_multip=5,surv_rate=surv_rate,beta=beta,sim_num=sim_num)

  out=simu_out


  return(out)
}
##################
#two_stg_adapt_power_simu_survival(HR,s,bry_type,alpha_2,evt_num,evt_max,HR_cut,surv_rate,beta,sim_num)
#two_stg_adapt_power_simu_survival(HR=rep(1,2),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,evt_num=150, evt_max=1500,HR_cut=0.95,surv_rate=0.5,beta=0.1,sim_num=100)
#two_stg_adapt_power_simu_survival(HR=c(0.7,0.8,1),s=c(1/3,1), bry_type=2,alpha_2=c(0.001,0.049),evt_num=550,evt_max=2000,HR_cut=0.95, surv_rate=0.5,beta=0.1,sim_num=100)

#######################################

###################################
Two_stage_adapt_rej_OC_simu_survival<-function(HR,lambda,s,c_bry_two_stg,evt_num,evt_max,HR_cut,samp_multip,surv_rate,beta,sim_num){
  #######################################
  comp=length(HR)
  ##########################
  evt_num_s=c(ceiling(evt_num*s),evt_max)
  samsz=ceiling(samp_multip*evt_max/surv_rate)
  theta_cut=-log(HR_cut)
  ######################## exit boundaries
  #e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    samsz_simu=evt_num_s[length(s)]; futi=0
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################

      ###########################################
      ### geting random sample from exponential distributions ################
      x0=stats::rexp(samsz,lambda) #generate exponential data for the control arm
      x=matrix(stats::rexp(samsz*comp, HR*lambda), nrow=comp,ncol=samsz)
      T=-log(surv_rate)/lambda  #Censoring time. Events after T will be censored
      time_x0=x0
      time_x=x
      x0_status=c(rep(1,length(x0)))
      x_status=matrix(rep(1,samsz*comp),nrow=comp,ncol=samsz)

      ######################## Censored data for the control arm #####################################
      for (i in (1:samsz)){
        if (x0[i]>T){
          time_x0[i]=T  #x0[i] is censored after T
          x0_status[i]=0   #censoring indicator
        }
      }
      grp_x0=c(rep(1,length(time_x0)))#vector of zeros, length of the vector is the same as the original smaple size

      ######################## Censored data for the test  arms #####################################
      #  grp_x=c(rep(0,length(time_x))) #
      grp_x=matrix(rep(0,samsz*comp),nrow=comp,ncol=samsz)
      for (j in (1:comp)){
        for (i in (1:samsz)){
          if (x[j,i]>T){
            time_x[j,i]=T #x[i] is censored after T
            x_status[j,i]=0  #censoring indicator
          }
        }
      }

      #### count the numnber of events in the control and test arms #######
      evt_count=matrix(rep(0,samsz*comp),nrow=comp,ncol=samsz)
      # for (j in (1:comp)){
      #   for (i in 1:samsz){
      #     evt_count[j,i]=sum(x0_status[1:i])+sum(x_status[j,1:i]) #number of events after 2*i observations
      #   }
      # }
      evt_count_x0 = rep(0, samsz)
      total_count_x0 = 0
      for (i in 1:samsz){
        total_count_x0 = total_count_x0 + x0_status[i]
        evt_count_x0[i] = total_count_x0
      }
      for (j in (1:comp)){
        total_count_x = 0
        for (i in 1:samsz){
          total_count_x = total_count_x + x_status[j,i]
          evt_count[j,i] = evt_count_x0[i] + total_count_x
        }
      }

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      ###############################set up matrices fro parameters########################################################

      theta_hat_0=rep(0,comp)
      theta_hat_0_sd=rep(0,comp)
      wald_0=rep(0,comp)
      score_0=rep(0,comp)
      t_0=rep(0,comp)
      ##### Obtaningin parameters and store in the matrices ##############
      evt_cut_0=rep(0,comp)

      for (j in (1:comp)) {
        evt_pos=evt_count[j,evt_count[j,]>evt_num_s[1]]#substring of  evt_count that are larger than evt_num_s[1]
        evt_cut_0[j]=match(evt_pos[1],evt_count[j,]) #the first incidence that evt_count is larger than evt_num
      }

      evt_cut=min(evt_cut_0)
      ## the required number of events may be unreachable within the available
      ## follow-up: the replicate is marked non-evaluable instead of using an NA index
      if (!is.finite(evt_cut)) evt_cut=length(evt_count[1,])

      for (j in (1:comp)) {
        #  vector of zeros, length of the vector is the same as the original smaple size
        time=c(time_x0[1:evt_cut],time_x[j,1:evt_cut]) #combined samples befroe total number of events exceed evt_num
        event=c(x0_status[1:evt_cut],x_status[j,1:evt_cut]) #censoring status before total number of events exceed evt_num
        grp=c(grp_x0[1:evt_cut],grp_x[j,1:evt_cut]) #vector of zero, with length=2*evt_cut


        test= list(time=time,status=event,grp)#survival data
        coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
        # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
        # test is the same as data=test
        theta_hat_0[j]= coxcf$coefficients[1] #estimated log(HR)
        theta_hat_0_sd[j]= coxcf$coefficients[3] #estimated standard error of log(HR)
        t_0[j]=(theta_hat_0_sd[j])^(-2)
        wald_0[j]=theta_hat_0[j]/ theta_hat_0_sd[j]
        score_0[j]= theta_hat_0[j]*t_0[j]
      }

      ## the required events for the phase-2 look were unreachable: blank the
      ## statistics so the replicate is reported as non-evaluable
      if (anyNA(evt_cut_0)) { theta_hat_0[]=NA_real_; theta_hat_0_sd[]=NA_real_; t_0[]=NA_real_; wald_0[]=NA_real_; score_0[]=NA_real_ }
      wald_p2_max=max(wald_0)
      test_1=wald_p2_max-c_bry_two_stg[1]
      ####################################

      #############################################################################################
      ### Step 3: the dose with the largest wald_0 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(wald_p2_max,wald_0)
      #s_slct identifies which dose was selected for the i-th simulation
      #############################################################
      ### Step 4: Identify  test statistics  for phase 3 stage ####
      #############################################################
      theta_hat_two_stg=rep(0,length(s))
      theta_hat_two_stg_sd=rep(0,length(s))
      wald=rep(0,length(s))
      score=rep(0,length(s))
      t=rep(0,length(s))
      #################################
      evt_count_slct=evt_count[s_slct,]
      time_x_slct=time_x[s_slct,]
      x_status_slct=x_status[s_slct,]
      grp_x_slct=grp_x[s_slct,]
      #################################
      for (i in (1: length(s))) {
        #  vector of zeros, length of the vector is the same as the original smaple size
        evt_pos=evt_count[s_slct,evt_count_slct>evt_num_s[i]]#substring of  evt_count that are larger than evt_num
        # evt_pos=evt_count[s_slct,evt_count_slct>evt_num_s[i+1]]#substring of  evt_count that are larger than evt_num

        evt_cut=match(evt_pos[1],evt_count_slct) #the first incidence that evt_count is larger than evt_num
        ## events required at this look are unreachable within follow-up:/n        ## the replicate is left non-evaluable instead of indexing with NA
        if (is.na(evt_cut)) { evt_cut=length(evt_count_slct); theta_hat_two_stg[i]=NA_real_; theta_hat_two_stg_sd[i]=NA_real_; t[i]=NA_real_; wald[i]=NA_real_; score[i]=NA_real_; next }
        time=c(time_x0[1:evt_cut],time_x_slct[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
        event=c(x0_status[1:evt_cut],x_status_slct[1:evt_cut]) #censoring status before total number of events exceed evt_num
        grp=c(grp_x0[1:evt_cut],grp_x_slct[1:evt_cut]) #vector of zero, with length=2*evt_cut


        test= list(time=time,status=event,grp)#survival data
        coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
        # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
        # test is the same as data=test
        theta_hat_two_stg[i]= coxcf$coefficients[1] #estimated log(HR)
        theta_hat_two_stg_sd[i]= coxcf$coefficients[3] #estimated standard error of log(HR)
        t[i]=(theta_hat_two_stg_sd[i])^(-2)
        wald[i]=theta_hat_two_stg[i]/theta_hat_two_stg_sd[i]
        score[i]= theta_hat_two_stg[i]*t[i]
      }

      #######################################

      evt_count_two_stg=evt_count[s_slct,]

      t_two_stg=t
      x_status_two_stg=x_status[s_slct,]
      wald_two_stg=wald
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg[1:length(s)]-c_bry_two_stg
      test_2=test[1:(length(s)-1)]
      # test compares the Wald statistics with the critical boundaries

      ##########################################################################################
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries
      ##########################################################################################

      inter_vt=length(s)-1
      if (isTRUE(theta_hat_two_stg[inter_vt]<theta_cut)) {futi=1} else {futi=0}
      ##########################################

      if (isTRUE(max(test)>=0)){
        rej_nossr=1
        test_p=test[test>0]
        stg_1=match(test_p[1],test)
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
        stg_1=0
      }

      ## (theta_hat[inter_vt]<theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation
      if (isTRUE(max(test_2)>0)) {

        rej=1
        #If max(tes_1)>0, te nulll hypothesis is rejected
        test_t=test_2[test_2>0]
        stg=match(test_t[1],test_2)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        theta_hat_last=theta_hat_two_stg[last_vt]
        theta_hat_last_sd=1/sqrt(t_two_stg[last_vt])

        samsz_simu=evt_num_s[last_vt]

      }else if (isTRUE(max(test_2)<=0)) {


        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation

          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat_two_stg[last_vt]
          theta_hat_last_sd=1/sqrt(t_two_stg[last_vt])


          samsz_simu=evt_num_s[length(s)]
          if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat_two_stg[inter_vt]
          theta_hat_inter_sd=theta_hat_two_stg_sd[inter_vt]
          t_hat_inter=theta_hat_inter_sd^(-2)
          score_inter=score[length(s)-1]
          event_inter=evt_num_s[inter_vt]
          t_final_0=t_hat_inter*(evt_num_s[length(s)]/evt_num_s[inter_vt])

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))


          if (isTRUE(cond_power>=1-beta)){

            last_vt=length(s)
            theta_hat_last=theta_hat_two_stg[last_vt]
            theta_hat_last_sd=theta_hat_two_stg_sd[last_vt]


            samsz_simu=evt_num_s[length(s)]
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (isTRUE(cond_power<1-beta)){
            #  c_bry=c_bry_two_stg
            #sample size re-estimation procedure##
            inter_vt=length(s)-1
            event_inter=evt_num_s[inter_vt]
            T_new_0=theta_hat_inter^(-2)*((c_bry_two_stg[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_samsz=ceiling(event_inter*T_new_0/t[inter_vt])

            if  (new_samsz>evt_max){

              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (isTRUE(test[length(s)]>0)){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              samsz_simu=evt_num_s[length(s)]
            } else  if  (new_samsz<=evt_max){
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              samsz_simu=new_samsz
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              evt_pos=evt_count_two_stg[evt_count_two_stg>new_samsz]#substring of  evt_count that are larger than evt_num
              ## events required after the sample size increase may be unreachable:/n              ## clamp the index so the subset below is valid, then blank the
              ## statistics so the replicate is reported as non-evaluable
              evt_short=is.na(match(evt_pos[1],evt_count_slct))
              evt_cut=if (evt_short) length(evt_count_slct) else match(evt_pos[1],evt_count_slct)

              time=c(time_x0[1:evt_cut],time_x_slct[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
              event=c(x0_status[1:evt_cut],x_status_slct[1:evt_cut]) #censoring status before total number of events exceed evt_num
              grp=c(grp_x0[1:evt_cut],grp_x_slct[1:evt_cut]) #vector of zero, with length=2*evt_cut


              test= list(time=time,status=event,grp)#survival data
              coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
              # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
              theta_hat_last= coxcf$coefficients[1] #estimated -log(HR)
              theta_hat_last_sd= coxcf$coefficients[3] #estimated standard error of log(HR)
              t_last=(theta_hat_last_sd)^(-2)
              T_new=t_last
              wald_last=theta_hat_last/ theta_hat_last_sd
              if (evt_short) { theta_hat_last=NA_real_; theta_hat_last_sd=NA_real_; t_last=NA_real_; wald_last=NA_real_ }
              e_bry_new_two_stg=
                (sqrt(T_new-t[inter_vt])/sqrt(t[length(s)]-t[inter_vt]))*(c_bry_two_stg[length(s)]*sqrt(t[length(s)])-score_inter)+ score_inter

              new_c_bry=e_bry_new_two_stg/sqrt(T_new)  #new critical  boundary for the combined phase 2/3 test. the NULL is rejected if the Wald statistics crossed the critical boundary

              #################################################################################
              ###                   Step 4.6: Hypothesis testing                           ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              test_add=wald_last-new_c_bry

              if (isTRUE(test_add>0)){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }

            }
          }
        }
      }
      #######################################################
      ### Step 5: output and finish the simjulation loop ####
      #######################################################
      out=c(rej_nossr,rej,stg,s_slct,samsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  samsz_simu = as.numeric(simr[,5])
  futi = as.numeric(simr[,6])
  ###################################
  samsz_simu_two_arm = samsz_simu
  samsz_simu_tot=(evt_num_s[1]/2)*(comp+1)+(samsz_simu_two_arm-evt_num_s[1])
  mean_two_arm_evt=mean(samsz_simu_two_arm)
  mean_tot_evt_num=mean(samsz_simu_tot)

  #################  results######
  av_samsz=mean(samsz_simu)
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  futi_rt=mean(futi)

  #####################################
  rej_HR=rep(0,length(HR))
  for (j in 1:length(HR)){
    rej_HR[j]=sum( rej[ s_slct==j])/(dim(simr)[1])
  }

  samsz_icr=NULL
  for (i in 1:dim(simr)[1]){
    if (samsz_simu[i]>evt_num){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }

  samsz_icr_rt=mean(samsz_icr)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  two_arm_evt_samsz_icr= samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_two_arm_evt_samsz_icr=mean(two_arm_evt_samsz_icr)

  ################################

  mean_two_arm_evt=mean(samsz_simu_two_arm)
  mean_tot_evt_num=mean(samsz_simu_tot)

  # tyI_er=sum(rej_HR[HR==1])/(dim(simr)[1])
  tyI_er=sum(rej_HR[HR==1])

  non_tyI_er_rej_rt=sum(rej[HR[s_slct]<1])/(dim(simr)[1])
  mean_two_arm_evt=mean(samsz_simu_two_arm)
  mean_tot_evt_num=mean(samsz_simu_tot)
  evt_tot=(evt_num_s[1]/2)*(comp-1)+evt_num
  sim_results=round(c( rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,
                       tyI_er, non_tyI_er_rej_rt,futi_rt,evt_num,evt_max,
                       mean_two_arm_evt,evt_tot,mean_tot_evt_num, mean_two_arm_evt_samsz_icr),5)
  names(sim_results)=c("rejection rate without SSR","rejection rate with SSR", "Rate of sample size increase",
                       "Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned event number for two arms", "Maximum number of events for two arms",
                       "Mean event number for two arms","Planned total number of events",
                       "total number of events","Mean event number for two arms conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=sim_results
  return(out)
}

#################################################
#Two_stage_adapt_rej_OC_simu_survival(HR,lambda,s,c_bry_two_stg,evt_num,evt_max,HR_cut,samp_multip,surv_rate,beta,sim_num) #######################################
#Two_stage_adapt_rej_OC_simu_survival(HR=rep(1,2),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050, 2.723811, 2.223982),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=10000)
#Two_stage_adapt_rej_OC_simu_survival(HR=rep(1,2),lambda=1,s=c(1/3,1),c_bry_two_stg=c(3.852050,  2.223982),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=10000)
#Two_stage_adapt_rej_OC_simu_survival(HR=c(0.8,0.8),lambda=1,s=c(1/3,1),c_bry_two_stg=c(3.852050,  2.223982),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)

#Two_stage_adapt_rej_OC_simu_survival(HR=c(0.8,0.8),lambda=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.852050 ,2.723811, 2.223982),evt_num=150,evt_max=1000,HR_cut=0.95,samp_multip=5,surv_rate=0.5,beta=0.1,sim_num=100)

##############################################################

#######################

two_stg_adapt_power_OC_simu_survival<-function(HR,s,bry_type,alpha_2,surv_rate,evt_num,evt_max,HR_cut,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  comp=length(HR)
  lambda=1
  dist=2
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=2)

  } else if (bry_type!=1) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=2)

  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))
  samp_multip=5
  simu_out=Two_stage_adapt_rej_OC_simu_survival(HR=HR,lambda=1,s=s,c_bry_two_stg=c_bry_two_stg,evt_num=evt_num,evt_max=evt_max,HR_cut=HR_cut,samp_multip=5,surv_rate=surv_rate,beta=beta,sim_num=sim_num)

  HR_names=paste(paste("HR","[",c(1:length(HR)),"]",sep=""))
  out=round(c(HR,as.numeric(simu_out)),5)
  names(out)=c(HR_names, names(simu_out))

  return(out)
}

##################################
#two_stg_adapt_power_OC_simu_survival(HR,s,bry_type,alpha_2,surv_rate,evt_num,evt_max,HR_cut,beta,sim_num)
#two_stg_adapt_power_OC_simu_survival(HR=c(0.8,0.8),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=200,evt_max=1500,HR_cut=0.95,beta=0.1,sim_num=1000)
#two_stg_adapt_power_OC_simu_survival(HR=c(1,1,1),s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=800,evt_max=2000,HR_cut=0.95,beta=0.1,sim_num=100)
#two_stg_adapt_power_OC_simu_survival(HR=c(1,1,0.8),s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),surv_rate=0.5,evt_num=800,evt_max=2000,HR_cut=0.95,beta=0.1,sim_num=100)
#two_stg_adapt_power_OC_simu_survival(HR=rep(1,2),s=c(1/3,1),bry_type=1,alpha_2=0.05,evt_num=150,evt_max=1000,HR_cut=0.95, surv_rate=0.5,beta=0.1,sim_num=100)
#two_stg_adapt_power_OC_simu_survival(HR=rep(1,2),s=c(1/3,1), bry_type=2,alpha_2=c(0.001,0.049),evt_num=150,evt_max=1000,HR_cut=0.95, surv_rate=0.5,beta=0.1,sim_num=100)
#two_stg_adapt_power_OC_simu_survival(HR=c(0.7,0.8,0.9),s=c(1/3,1), bry_type=2,alpha_2=c(0.001,0.049),evt_num=550,evt_max=2000,HR_cut=0.95, surv_rate=0.5,beta=0.1,sim_num=100)
#two_stg_adapt_power_OC_simu_survival(HR=c(0.7,0.8,1),s=c(1/3,1), bry_type=2,alpha_2=c(0.001,0.049),evt_num=550,evt_max=2000,HR_cut=0.95, surv_rate=0.5,beta=0.1,sim_num=100)


##################################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Survival analysis, Optimizing adaptive sequential design
#'
#' @param groups Matrix of assumed hazard ratios for the test arms under study: one row per candidate test-arm scenario and one column per dose, so that the number of doses is the number of columns; each hazard ratio is relative to the control arm, must be > 0, with HR < 1 indicating a beneficial effect and HR = 1 the null hypothesis
#' @param s1 Vector of information fractions of the first interim analysis schedule; the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule; the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s1; each element must be > 0 and the elements must sum to alpha_2
#' @param alpha_2_alps_s2 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s2; each element must be > 0 and the elements must sum to alpha_2
#' @param surv_rate Survival rate at the end of follow-up, used to translate the number of events into a sample size
#' @param evt_num_1 Planned number of events between a test arm and the control arm for the first event number scenario
#' @param evt_num_2 Planned number of events between a test arm and the control arm for the second event number scenario
#' @param evt_max_1 Maximum number of events between a test arm and the control arm allowed after sample size re-estimation for the first maximum event number scenario
#' @param evt_max_2 Maximum number of events between a test arm and the control arm allowed after sample size re-estimation for the second maximum event number scenario
#' @param HR_cut Hazard-ratio cut-off for lack of efficacy at the interim analysis: the test drug is not considered to provide clinically meaningful efficacy when the observed hazard ratio is greater than HR_cut
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with O'Brien-Fleming boundaries}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with O'Brien-Fleming boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with alpha-spending boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with alpha-spending boundaries}
#' }
#' Each of these four lists has the components "Group sequential design" (data frame with the columns "Looks", "Information fractions" and the critical boundaries, plus "Alpha-spending" and "Alpha-spending boundaries" for the alpha-spending schedules) and the "Simulation results: Information fractions ..., Sample size ..., evt_max ..." operating characteristics of the corresponding adaptive design, one data frame per combination of interim analysis schedule, planned event number and maximum number of events, with one row per assumed hazard-ratio scenario.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' Two_stage_asd_simu_OC_survival(
#' groups = rbind(c(0.8, 0.83, 0.85), c(0.78, 0.8, 0.82), c(0.75, 0.78, 0.8)), s1 = c(1/3,2/3,1),
#' s2 = c(0.5,1), alpha_2 = 0.05, alpha_2_alps_s1 = c(0.001,0.001,0.048),
#' alpha_2_alps_s2 = c(0.01,0.04), surv_rate = 0.5, evt_num_1 = 200, evt_num_2 = 300,
#' evt_max_1 = 2000, evt_max_2 = 3000, HR_cut = 0.95, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
Two_stage_asd_simu_OC_survival<-function(groups,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,surv_rate,evt_num_1,evt_num_2,evt_max_1,evt_max_2,HR_cut,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)
  comp=length(groups[1,])

  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=2)[1]))
  s2_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=2)[1]))

  s1_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s1,alpha=alpha_alsp_s1,comp=comp,dist=2)[1]))
  s2_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s2,alpha=alpha_alsp_s2,comp=comp,dist=2)[1]))

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
    simu_out_OF_s1_EVTnum_1_EVTmax_1=rbind(simu_out_OF_s1_EVTnum_1_EVTmax_1,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s1,bry_type=1,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                                     evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_EVTnum_1_EVTmax_1)=HR_group_names

  ###############  s2 #######
  simu_out_OF_s2_EVTnum_1_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s2_EVTnum_1_EVTmax_1=rbind(simu_out_OF_s2_EVTnum_1_EVTmax_1,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s2,bry_type=1,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_EVTnum_1_EVTmax_1)=HR_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_OF_s1_EVTnum_2_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s1_EVTnum_2_EVTmax_1=rbind(simu_out_OF_s1_EVTnum_2_EVTmax_1,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s1,bry_type=1,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_EVTnum_2_EVTmax_1)=HR_group_names

  ###############  s2 #######
  simu_out_OF_s2_EVTnum_2_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s2_EVTnum_2_EVTmax_1=rbind(simu_out_OF_s2_EVTnum_2_EVTmax_1,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s2,bry_type=1,alpha_2=alpha_2,surv_rate=surv_rate,
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
    simu_out_OF_s1_EVTnum_1_EVTmax_2=rbind(simu_out_OF_s1_EVTnum_1_EVTmax_2,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s1,bry_type=1,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_EVTnum_1_EVTmax_2)=HR_group_names

  ###############  s2 #######
  simu_out_OF_s2_EVTnum_1_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s2_EVTnum_1_EVTmax_2=rbind(simu_out_OF_s2_EVTnum_1_EVTmax_2,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s2,bry_type=1,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_EVTnum_1_EVTmax_2)=HR_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_OF_s1_EVTnum_2_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s1_EVTnum_2_EVTmax_2=rbind(simu_out_OF_s1_EVTnum_2_EVTmax_2,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s1,bry_type=1,alpha_2=alpha_2,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_EVTnum_2_EVTmax_2)=HR_group_names

  ###############  s2 #######
  simu_out_OF_s2_EVTnum_2_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_OF_s2_EVTnum_2_EVTmax_2=rbind(simu_out_OF_s2_EVTnum_2_EVTmax_2,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s2,bry_type=1,alpha_2=alpha_2,surv_rate=surv_rate,
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
    simu_out_alsp_s1_EVTnum_1_EVTmax_1=rbind(simu_out_alsp_s1_EVTnum_1_EVTmax_1,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s1,bry_type=2,alpha_2=alpha_2_alps_s1,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_EVTnum_1_EVTmax_1)=HR_group_names

  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_1_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s2_EVTnum_1_EVTmax_1=rbind(simu_out_alsp_s2_EVTnum_1_EVTmax_1,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s2,bry_type=2,alpha_2=alpha_2_alps_s2,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_EVTnum_1_EVTmax_1)=HR_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_alsp_s1_EVTnum_2_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s1_EVTnum_2_EVTmax_1=rbind(simu_out_alsp_s1_EVTnum_2_EVTmax_1,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s1,bry_type=2,alpha_2=alpha_2_alps_s1,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_EVTnum_2_EVTmax_1)=HR_group_names

  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_2_EVTmax_1=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s2_EVTnum_2_EVTmax_1=rbind(simu_out_alsp_s2_EVTnum_2_EVTmax_1,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s2,bry_type=2,alpha_2=alpha_2_alps_s2,surv_rate=surv_rate,
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
    simu_out_alsp_s1_EVTnum_1_EVTmax_2=rbind(simu_out_alsp_s1_EVTnum_1_EVTmax_2,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s1,bry_type=2,alpha_2=alpha_2_alps_s1,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_EVTnum_1_EVTmax_2)=HR_group_names

  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_1_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s2_EVTnum_1_EVTmax_2=rbind(simu_out_alsp_s2_EVTnum_1_EVTmax_2,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s2,bry_type=2,alpha_2=alpha_2_alps_s2,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_EVTnum_1_EVTmax_2)=HR_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_alsp_s1_EVTnum_2_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s1_EVTnum_2_EVTmax_2=rbind(simu_out_alsp_s1_EVTnum_2_EVTmax_2,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s1,bry_type=2,alpha_2=alpha_2_alps_s1,surv_rate=surv_rate,
                                                                                                                   evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_EVTnum_2_EVTmax_2)=HR_group_names

  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_2_EVTmax_2=NULL
  for (i in 1:dim(HR_vec)[1]){
    simu_out_alsp_s2_EVTnum_2_EVTmax_2=rbind(simu_out_alsp_s2_EVTnum_2_EVTmax_2,two_stg_adapt_power_OC_simu_survival(HR=HR_vec[i,],s=s2,bry_type=2,alpha_2=alpha_2_alps_s2,surv_rate=surv_rate,
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
#Two_stage_asd_simu_OC_survival(HR_vec_1,HR_vec_2,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,evt_num_1,evt_num_2,evt_max_1,evt_max_2,beta,sim_num)
#Two_stage_asd_simu_OC_survival(groups=rbind(c(0.75,0.8,0.85),c(0.7,0.75,0.8),c(0.65,0.7,0.75)),
#                      s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.01,0.04),
#                      surv_rate=0.5,evt_num_1=200,evt_num_2=300,evt_max_1=500,evt_max_2=600,HR_cut=0.95,beta=0.1, sim_num=100)

##########################################################################
### II.3 .Simulations  -binary distributions                         #####
##########################################################################


##########################################################################
### II.3.1 .Simulations with planned TSSSD -binary distributions     #####
##          Type I error, power, performance simulaitons             #####
#### This will be included in the software                           #####
##########################################################################


#########################################
two_stg_rej_rate_simu_binary<-function(p_t,p_c,direction,s,c_bry_two_stg,samsz,sim_num){

  # direction==1 means that larger rate p is better
  # direction==2 means that smaller rate p is better
  #mu is the vector of means for each dose and control
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsaz is sample size

  N=ceiling(samsz*s)
  comp=length(p_t)
  # N is the vetoc of sample sizes at which interim analysis will be conducted
  #s_slct is the vector indicating which comparison (e.g. which dose) was selected at the end of the


  t=matrix(rep(0,comp*length(s)),nrow=comp,ncol=length(s))
  for (i in 1:length(p_t)) {
    t[i,]<- N/(p_t[i]*(1-p_t[i])+p_c*(1-p_c))
  }
  ##########################
  N1=N[1]
  NF=N[length(s)]
  #t_two_stg is the vector of information times at the interim analysis
  if (direction ==1){ #higher rate is better
    mu=p_t-p_c
  }else if (direction ==2){ #lower rate is better
    mu=p_c-p_t
  }

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    rej=0; stg=0
    s_slct=0; control_samsz_simu=N[length(s)]
    tot_samsz_simu=N[1]*comp+control_samsz_simu
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################

      x0=stats::rbinom(NF,1,p_c)
      x0_p2=x0[1:N1]
      #x0_p2 is the random sample at from phase 2 stagefor the control
      x=matrix(stats::rbinom(NF*comp,1,p_t), nrow=comp,ncol=NF)
      x_p2=x[,1:N1]
      #x is the matrix of random samples for each dose from the phase 2 stage
      #for each simulation.
      #x matirx does not include the control sample.
      #The row number is the same as the active dose

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################
      S_p2=NULL
      wald_p2=NULL
      if (direction==1){# higher event rate is better
      for (j in 1:length(p_t)){
        S_p2[j]=(mean(x_p2[j,])-mean(x0_p2))*t[j,1]
        wald_p2[j]=(mean(x_p2[j,])-mean(x0_p2))*sqrt(t[j,1])
      }
      }else if (direction==2){# lower event rate is better

        for (j in 1:length(p_t)){
          S_p2[j]=(mean(x0_p2)-mean(x_p2[j,]))*t[j,1]
          wald_p2[j]=(mean(x0_p2)-mean(x_p2[j,]))*sqrt(t[j,1])
        }
      }

      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(max(S_p2),S_p2)

      t_two_stg=N/(p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c))
      #s_slct identifies which dose was selected for the i-th simulation
      ###########################################################
      ### Step 4: generate  random samples for phase 3 stage ####
      ###########################################################
   #   x0_p3=stats::rbinom((NF-N1),1,p_c)
      x_p3=stats::rbinom((NF-N1),1,p_t[s_slct])
      #y are the samples for the selected dose for the phase 3 stage
      x0_two_stg=x0
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_two_stg=x[s_slct,]
      #x_two_stg contains all the random samples(including samples from the phase 2 stage) for the seected dose

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      wald_two_stg=NULL
      for (j in 1:length(s) ){
        if (direction==1){# higher event rate is better
        wald_two_stg[j]= (mean(x_two_stg[1:N[j]])-mean(x0_two_stg[1:N[j]]))*sqrt(t_two_stg[j])
        } else  if (direction==2){# lower event rate is better
          wald_two_stg[j]= (mean(x0_two_stg[1:N[j]])-mean(x_two_stg[1:N[j]]))*sqrt(t_two_stg[j])
        }
      }
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if (isTRUE(max(test)>0)){

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        last_vt=stg
        control_samsz_simu=N[last_vt]
        tot_samsz_simu=N[1]*comp+control_samsz_simu

        #stg identifies the interim analysis at which the null hypothesis is rejected

      }else{

        rej=0
        stg=0
        last_vt=length(s)
        control_samsz_simu=N[last_vt]
        tot_samsz_simu=N[1]*(comp+1)+2*(control_samsz_simu-N[1])
      }
      #######################################################
      ### Step 7: output and finish the simjulation loop ####
      #######################################################
      out=c(rej,stg,s_slct,control_samsz_simu,tot_samsz_simu)
  })#(end of q loop)

  ############################################
  ### Step 8: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  s_slct= as.numeric(simr[,3])
  control_samsz_simu= as.numeric(simr[,4])
  tot_samsz_simu= as.numeric(simr[,5])
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  rej_num=sum(rej>0)
  # rej_num1=length(rej[rej>0])
  # rej_num2=sum(stg>0)

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=sum(rej[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  select_rt=rep(0,length(p_t))
  for (j in 1:length(p_t)){
    select_rt[j]=length( s_slct[ s_slct==j])/sim_num
  }


  rej_rt_indiv=rep(0,length(p_t))
  for (j in 1:length(p_t)){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/sim_num
  }

  con_rej_rt_indiv=rep(0,length(p_t))
  for (j in 1:length(p_t)){
    if (isTRUE(select_rt[j]==0)){
      con_rej_rt_indiv[j]=NA
    }else {
      con_rej_rt_indiv[j]=rej_rt_indiv[j]/select_rt[j]
    }
  }


  results_by_arm=cbind(p_t,p_c,select_rt, rej_rt_indiv,con_rej_rt_indiv)
  colnames(results_by_arm)=c("p<sub>t1</sub>,...","p<sub>c</sub>","selection rate"," rejection rate by each p_t","conditional rejection rate")
  #  rej_0=rej[mu[s_slct]==0]

  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look")

  tyI_er=sum(rej[mu[s_slct]==0])/sim_num
  #  tyI_er=sum(rej_0)/sim_num

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/sim_num
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(p_t)-1)+2*samsz
  sim_results=c(rej_rt,tyI_er, non_tyI_er_rej_rt,samsz,mean_control_samsz,tot_samsz,mean_tot_samsz)
  names(sim_results)=c("rejection rate","Type I error","Non-Type I error rejection rate",
                       "Planned control sample size","Mean control sample size","Planned total sample size","Mean total sample size")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(seq_design,results_by_arm,sim_results)
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)
}
#######################################
#two_stg_rej_rate_simu_binary(p_t,p_c,s,c_bry_two_stg,samsz,sim_num)
#two_stg_rej_rate_simu_binary(p_t=c(0.75,0.6,0.8),p_c=0.6,direction=1,s=c(0.3,0.7,1),c_bry_two_stg=c( 3.938123,2.784673,2.273676),samsz=100,sim_num=10000)
#two_stg_rej_rate_simu_binary(p_t=rep(0.75,3),p_c=0.6,s=c(0.3,0.7,1),c_bry_two_stg=c( 3.938123,2.784673,2.273676),samsz=100,sim_num=10000)
#two_stg_rej_rate_simu_binary(p_t=rep(0.6,3),p_c=0.6,direction=1,s=c(0.3,0.7,1),c_bry_two_stg=c( 3.938123, 2.784673,2.273676),samsz=100,sim_num=100000)
###############Estimates#####################


##############################################

#' Phase 2/3 seamless combination Group sequential design: Simulations, Binary distribution, Type I error and power
#'
#' @param p_t Event rate of the test arm
#' @param p_c Event rate of the control arm
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param samsz Planned sample size of the control arm used in the simulation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries" and "rejection by look"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "p<sub>t1</sub>,...", "p<sub>c</sub>", "selection rate", " rejection rate by each p_t" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the overall "rejection rate", the "Type I error", the "Non-Type I error rejection rate", the "Planned control sample size", the "Mean control sample size", the "Planned total sample size" and the "Mean total sample size"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_power_simu_binary(p_t = c(0.4,0.4,0.4), p_c = 0.4, direction = 1, s = c(1/3,2/3,1),
#' bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
two_stg_power_simu_binary<-function(p_t,p_c,direction,s,bry_type,alpha_2,samsz,sim_num){
  # direction==1 means that larger rate p is better
  # direction==2 means that smaller rate p is better

  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=length(p_t)
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))
  simu_out=two_stg_rej_rate_simu_binary(p_t,p_c,direction,s,c_bry_two_stg,samsz,sim_num)
  out=simu_out
  return(out)
}

##################################
#two_stg_power_simu_binary(p_t,p_c,s,bry_type,alpha_2,samsz,sim_num)
#two_stg_power_simu_binary(p_t=c(0.4,0.4),p_c=0.2,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)
#two_stg_power_simu_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000000)
#two_stg_power_simu_binary(p_t=c(0.4,0.5,0.6),p_c=0.4,direction=1,s=c(1/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)


#######################################

##############################################


##########################################################################
### II.3 .Simulations with adaptive TSSSD -binary distributions     #####
##          Type I error, power, simulations                        #####
##########################################################################
#############################################
Two_stage_adapt_rej_simu_binary<-function(p_t,p_c,direction,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist){

  # direction==1 means that larger rate p is better
  # direction==2 means that smaller rate p is better
  ### parameters for two-stage #########

  N=ceiling(samsz*s)
  comp=length(p_t)
  if (direction ==1){ #higher rate is better
    mu=p_t-p_c
  }else if (direction ==2){ #lower rate is better
    mu=p_c-p_t
  }
  s_p3=(s[-1]-s[1])/(s[length(s)]-s[1])
  # N is the vetoc of sample sizes at which interim analysis will be conducted
  T=NULL
  for (i in 1:comp){
    T[i]=samsz/(p_t[i]*(1-p_t[i])+p_c*(1-p_c))
  }
  t=matrix(ncol=length(s),nrow=comp)
  for (i in 1:comp){
    t[i,]=T[i]*s
  }
  #t_p3=(T-t[1])*s_p3
  #t is the vector of information times at the interim analysis
  N1=N[1]
  #N1 is the sample at the end of phase 2 stage
  NF=N[length(s)]
  #s_0=s[1] #information fraction at the end of phase 2 stage
  #information fraction conditonal on phase 2 sthage
  s_1=(s[-1]-s[1])
  ################### Parameters for phase 3 only ##################
  N_p3=N[-1]-N1

  ## s_1 is the scaled information fraction for selection the ##
  ## critical boundary for phase 3 stage, this boundary will be only used for point estiante and upper limit of confiedence intgerval##
  ######################## exit boundaries##########
  e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined
  #  e_bry_p3=c_bry_p3*sqrt(s_p3) # exit boundary for phase 3  only

  ############ the simulation loops ##########
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=N[length(s)]
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    control_sampsz_simu=sampsz_simu; tot_sampsz_simu=(comp+1)*sampsz_simu
    futi=0
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################
      futi=0


      x0=stats::rbinom(N_max,1,p_c)
      #x0_p2 is the random sample at from phase 2 stagefor the control
      x=matrix(stats::rbinom(N_max*comp,1,p_t), nrow=comp,ncol=N_max)
      #x is the matrix of random samples for each dose from the phase 2 stage
      #for each simulation.
      #x matirx does not include the control sample.
      #The row number is the same as the active dose

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################


      S_ph2=NULL
      wald_ph2=NULL
      if (direction ==1){ #higher rate p is better
        for (j in 1:length(p_t)){
          S_ph2[j]=(mean(x[j,1:N1])-mean(x0[1:N1]))*t[j,1]
          wald_ph2[j]=(mean(x[j,1:N1])-mean(x0[1:N1]))*sqrt(t[j,1])
        }
      }else if (direction ==2){ #lower rate p is better
        for (j in 1:length(p_t)){
          S_ph2[j]=(mean(x0[1:N1])-mean(x[j,1:N1]))*t[j,1]
          wald_ph2[j]=(mean(x0[1:N1])-mean(x[j,1:N1]))*sqrt(t[j,1])
        }
      }

      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      # s_slct=match(max(S_ph2),S_ph2)
      s_slct=match(max(wald_ph2),wald_ph2)
      #    s_slct=match(wald_ph2,max(wald_ph2))

      T_slct<- samsz/(p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c))
      t_slct=T_slct*s
      t_p3=t_slct[-1]-t_slct[1]
      #t_p3=(T_slct-t_slct[1])*s_p3
      #s_slct identifies which dose was selected for the i-th simulation
      ###########################################################
      ### Step 4: generate  random samples for phase 3 stage ####
      ###########################################################

      x0_two_stg=x0[1:NF]
      x0_p3=x0[(N1+1):NF]
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_slct_two_stg=x[s_slct,1:NF]
      x_slct_p3=x[s_slct,(N1+1):NF]

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################

      theta_hat=NULL
      wald_two_stg=NULL
      for (j in 1:length(s) ){

        theta_hat[j]=mean(x_slct_two_stg[1:N[j]])-mean(x0_two_stg[1:N[j]])  # theta_hat at each interim analysis
        wald_two_stg[j]=theta_hat[j]*sqrt(t_slct[j])
      }
      theta_hat_sd=1/sqrt(t_slct)
      theta_hat_two_stg=theta_hat
      theta_hat_two_stg_sd=theta_hat_sd
      ######  calculating the statistics for phase 3 only #######
      theta_hat=NULL  #theta_hat for phase 3 only
      for (j in 1:length(s)){
        theta_hat[j]= mean(x[s_slct,1:N[j]])-mean(x0[1:N[j]])  # theta_hat for phase 3 only, at each interim analysis
      }
      #theta_hat_p3_sd=1/sqrt(t_p3)
      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test_1=test[-length(s)]  ## test statistics excluding the final analysis

      ##########################################
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries
      ############################################
      if (isTRUE(theta_hat[length(s)-1]<theta_cut)) {futi=1} else {futi=0}
      #If the NUll hypothesis is rejected at up to length(s)-1, terminate trial and estimate CI ####

      #############################

      if (isTRUE(max(test)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
      }
      ##############

      if  (isTRUE(max(test_1)>=0)){ #null hypothesis rejected, no more random sampling ot theta_hat at length(s_p3)-1 is almost negative


        #  if ( max(tes_1)>0){ #null hypothesis rejected, no more random sampling
        rej=1
        #If max(tes_1)>0, te nulll hypothesis is rejected
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=1/sqrt(t_slct[last_vt])
        ###########################
        sampsz_simu=N[last_vt]
        if (last_vt==1){
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*sampsz_simu
        }else {
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*N[1]+2*(sampsz_simu-N[1])
        }

      }else if  (isTRUE(max(test_1)<0)){ ## (theta_p3_hat[length(s_1)-1]<theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation

        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat_two_stg[last_vt]
          theta_hat_last_sd=1/sqrt(t_slct[last_vt])

          sampsz_simu=N[length(s)]
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*N[1]+2*(sampsz_simu-N[1])
          if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }

          ############################

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat_two_stg[inter_vt]
          theta_hat_inter_sd=theta_hat_two_stg_sd[inter_vt]
          t_inter=t_slct[inter_vt]
          score_inter=theta_hat_inter*t_inter
          c_bry=c_bry_two_stg

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t_slct[length(s)]-t_slct[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t_slct[length(s)])-score_inter)/sqrt(t_slct[length(s)]-t_slct[inter_vt]))

          if (isTRUE(cond_power>=1-beta)){
            last_vt=length(s)
            theta_hat_last=theta_hat_two_stg[last_vt]
            theta_hat_last_sd=1/sqrt(t_slct[last_vt])


            sampsz_simu=N[length(s)]
            control_sampsz_simu=sampsz_simu
            tot_sampsz_simu=(comp+1)*N[1]+2*(sampsz_simu-N[1])
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (isTRUE(cond_power<1-beta)){

            ## proceed to sample size re-estimation #####

            ######### if the NULL hypothesis is not rejected at up to length(s)-1, and theta_p3_hat[length(s_1)-1]>=theta_cut, indicating potential efficacy ###
            ## re-estimate sample size#####
            #sample size re-estimation procedure##
            ## statistic sfor combined phae 2/3 data#####
            inter_vt=length(s)-1
            theta_hat_inter=theta_hat[inter_vt]
            theta_hat_inter_sd=1/sqrt(t_slct[inter_vt])

            score_inter=theta_hat_inter*t_slct[inter_vt]
            #########################
            T_new_0=theta_hat_inter^(-2)*((e_bry_two_stg[length(s)]*sqrt(T_slct)-score_inter)/sqrt(T_slct-t_slct[inter_vt])+stats::qnorm(1-beta))^2+t_slct[inter_vt]
            new_sampsz=ceiling(T_new_0*((p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c))))
            T_new=new_sampsz/((p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c)))
            e_bry_new_two_stg=
              (sqrt(T_new-t_slct[inter_vt])/sqrt(T_slct-t_slct[inter_vt]))*(e_bry_two_stg[length(s)]*sqrt(T_slct)-score_inter)+ score_inter
            c_bry_new_two_stg=e_bry_new_two_stg/sqrt(T_new)
            ###############################################
            if (isTRUE(new_sampsz>N_max)){ ## new sample size is too large, the sample size will not change. Trial ends at the planned sample size
              if (isTRUE(test[length(test)]>0)){ ### the NULL hypothesis is rejected at the final anbalysis, with the planned sample size
                rej=1
                stg=length(s)
              }else{### the NULL hypothesis is NOT rejected at the final anbalysis, with the planned sample size
                rej=0
                stg=0
              }
              ###################### Statistics for estimates ####
              ### combined phase 2/3 data ####
              last_vt=length(s)
              sampsz_simu=N[last_vt]
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              ####################### this completes the estimation when there is no saple size change ########

            } else if (isTRUE(new_sampsz<= N_max)) {## new sample size does not exceed the pre-selected maximum. Will change the sample  size to this new sample size
              #here
              sampsz_simu=new_sampsz
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
              x0_two_stg_new=x0[1:new_sampsz] # combine the previous [combined phase 2/3] data of the control arm (up to inter_vt), and the newly generated additional random sample

              x_two_stg_new=x[s_slct,1:new_sampsz] # combine the previous [combined phase 2/3] data of the selectd dose arm (up to inter_vt), and the newly generated additional random sample

              T_new=new_sampsz/(p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c))   # the new information time for the combined phase 2/3 test


              ### statistics from combined data ####

              if (direction ==1){ #higher rate p is better
                wald_last_ad=(mean(x_two_stg_new)-mean(x0_two_stg_new))*sqrt(T_new)  # Wald statistics at the adjusted final analysis
              } else   if (direction ==2){ #lower rate p is better
                wald_last_ad=(mean(x0_two_stg_new)-mean(x_two_stg_new))*sqrt(T_new)  # Wald statistics at the adjusted final analysis
              }
              ##########################################

              test_add=wald_last_ad-c_bry_new_two_stg  # compare the final Wald statistics with the new adjusted critical boundary

              if (isTRUE(test_add>0)){ # the new critical boundary is crossed, the NULL is rejected

                rej=1
                stg=length(s)

              }else{  # the new critical boundary is NOT crossed, the NULL is rejected

                rej=0
                stg=0

              }

              ###################### Parameter estimate with sample size change ##################

            }
          }
        }
      }
      out=c(rej_nossr,rej,stg,s_slct,control_sampsz_simu,tot_sampsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  control_samsz_simu = as.numeric(simr[,5])
  tot_samsz_simu = as.numeric(simr[,6])
  futi=as.numeric(simr[,7])
  #############termination timepoint####
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  rej_time=rep(0,(length(s)))
  for (j in 1:(length(s))){

    count_j=c(1:sim_num)[stg==j]
    rej_time[j]=sum(rej[count_j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  select_rt=rep(0,length(p_t))
  for (j in 1:length(p_t)){
    select_rt[j]=length( s_slct[ s_slct==j])/sim_num
  }

  rej_rt_indiv=rep(0,length(p_t))
  for (j in 1:length(p_t)){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/sim_num
  }

  con_rej_rt_indiv=rep(0,length(p_t))
  for (j in 1:length(p_t)){
    if (isTRUE(select_rt[j]==0)){
      con_rej_rt_indiv[j]=NA
    }else {
      con_rej_rt_indiv[j]=rej_rt_indiv[j]/select_rt[j]
    }
  }

  samsz_icr=NULL
  for (i in 1:length(control_samsz_simu)){
    if (control_samsz_simu[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_rt_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  futi_rt=mean(futi)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= control_samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  #############################
  ### Step 9: Final output ####
  #############################

  results_by_arm=cbind(p_t,p_c,select_rt, rej_rt_indiv,con_rej_rt_indiv)
  colnames(results_by_arm)=c("p<sub>t1</sub>,...","p<sub>c</sub>","selection rate","rejection rate by each p_t","conditional rejection rate")
  #  rej_0=rej[mu[s_slct]==0]

  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time, samsz_icr_rt_byLook)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look","Sample size increase")

  tyI_er=sum(rej[mu[s_slct]==0])/sim_num
  #  tyI_er=sum(rej_0)/sim_num

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/sim_num
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(p_t)-1)+2*samsz
  sim_results=c(rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er, non_tyI_er_rej_rt,futi_rt,samsz,
                N_max,mean_control_samsz,tot_samsz,mean_tot_samsz, mean_samsz_simu_control_samsz_icr)
  names(sim_results)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned control sample size","Maximum control sample size","Mean control sample size",
                       "Planned total sample size","Mean total sample size","mean control sample size conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(seq_design,results_by_arm,sim_results)
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)

}


#################################################
#Two_stage_adapt_rej_simu_binary(p_t,p_c,direction,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist)

#Two_stage_adapt_rej_simu_binary(p_t=c(0.6,0.6,0.6),p_c=0.6,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000,dist=1)
#wo_stage_adapt_rej_simu_binary(p_t=c(0.6,0.6,0.6),p_c=0.6,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000,dist=1)

#Two_stage_adapt_rej_simu_binary(p_t=c(0.7,0.6,0.6),p_c=0.6,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000,dist=1)
#Two_stage_adapt_rej_simu_binary(p_t=c(0.7,0.6,0.6),p_c=0.6,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000,dist=1)

###############################################################

#' Phase 2/3 seamless combination Group sequential design: Simulations, Binary distribution, Operating characteristics
#'
#' @param p_t Event rate of the test arm; each element must be > 0 and < 1
#' @param p_c Event rate of the control arm; must be > 0 and < 1
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param theta_cut Lack-of-efficacy threshold: an effect size below which the test drug is not considered to provide clinically meaningful efficacy
#' @param N_max Maximum sample size (control arm) allowed after sample size re-estimation
#' @param samsz Planned sample size of the control arm used in the simulation
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned adaptive group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries", "rejection by look" and "Sample size increase"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "p<sub>t1</sub>,...", "p<sub>c</sub>", "selection rate", "rejection rate by each p_t" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the "Rejection rate without SSR", the "Rejection rate with SSR", the "rate of sample size increase", the "Rejection rate without SSR/(meeting sample size incerase criteria)", the "Rejection rate with SSR/(sample size incerase)", the "Type I error", the "Non-Type I error rejection rate", the "Futility rate", the "Planned control sample size", the "Maximum control sample size", the "Mean control sample size", the "Planned total sample size", the "Mean total sample size" and the "mean control sample size conditional on increase"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_adapt_power_simu_binary(p_t = c(0.4,0.5), p_c = 0.3, direction = 1, s = c(1/3,2/3,1),
#' bry_type = 1, alpha_2 = 0.05, theta_cut = 0.05, N_max = 1000, samsz = 100,
#' beta = 0.1, sim_num = 20)
#'
#' @export
two_stg_adapt_power_simu_binary<-function(p_t,p_c,direction,s,bry_type,alpha_2,theta_cut,N_max,samsz,beta,sim_num){
  # direction==1 means that larger rate p is better
  # direction==2 means that smaller rate p is better

  # This simulation function checks the actua power and sample size for different scenerios##
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2

  comp=length(p_t)

  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))

  simu_out=Two_stage_adapt_rej_simu_binary(p_t=p_t,p_c=p_c,direction=direction,s=s,c_bry_two_stg=c_bry_two_stg,alpha=alpha,beta=beta,theta_cut=theta_cut,samsz=samsz,N_max=N_max,sim_num=sim_num,dist=1)
  out=simu_out
  return(out)
}

##################################
#two_stg_adapt_power_simu_binary(p_t,p_c,s,bry_type,alpha_2,theta_cut,N_max,samsz,beta,sim_num)
#two_stg_adapt_power_simu_binary(p_t=c(0.4,0.5),p_c=0.3,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,N_max=1000,samsz=100,beta=0.1,sim_num=100)
#two_stg_adapt_power_simu_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,N_max=1000,samsz=100,beta=0.1,sim_num=100)
#two_stg_adapt_power_simu_binary(p_t=c(0.4,0.4,0.4),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,N_max=1000,samsz=100,beta=0.1,sim_num=100)

#two_stg_adapt_power_simu_binary(p_t=c(0.4,0.5,0.6),p_c=0.4,direction=1,s=c(1/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,N_max=1000,samsz=100,beta=0.1,sim_num=10000)


#############################################
Two_stage_adapt_rej_OC_simu_binary<-function(p_t,p_c,direction,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist){

  # direction==1 means that larger rate p is better
  # direction==2 means that smaller rate p is better
  ### parameters for two-stage #########

  N=ceiling(samsz*s)
  comp=length(p_t)
  if (direction ==1){ #higher rate is better
    mu=p_t-p_c
  }else if (direction ==2){ #lower rate is better
    mu=p_c-p_t
  }
  s_p3=(s[-1]-s[1])/(s[length(s)]-s[1])
  # N is the vetoc of sample sizes at which interim analysis will be conducted
  T=NULL
  for (i in 1:comp){
    T[i]=samsz/(p_t[i]*(1-p_t[i])+p_c*(1-p_c))
  }
  t=matrix(ncol=length(s),nrow=comp)
  for (i in 1:comp){
    t[i,]=T[i]*s
  }
  #t_p3=(T-t[1])*s_p3
  #t is the vector of information times at the interim analysis
  N1=N[1]
  #N1 is the sample at the end of phase 2 stage
  NF=N[length(s)]
  #s_0=s[1] #information fraction at the end of phase 2 stage
  #information fraction conditonal on phase 2 sthage
  s_1=(s[-1]-s[1])
  ################### Parameters for phase 3 only ##################
  N_p3=N[-1]-N1

  ## s_1 is the scaled information fraction for selection the ##
  ## critical boundary for phase 3 stage, this boundary will be only used for point estiante and upper limit of confiedence intgerval##
  ######################## exit boundaries##########
  e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined
  #  e_bry_p3=c_bry_p3*sqrt(s_p3) # exit boundary for phase 3  only

  ############ the simulation loops ##########
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=N[length(s)]
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    control_sampsz_simu=sampsz_simu; tot_sampsz_simu=(comp+1)*sampsz_simu
    futi=0
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################
      futi=0

      x0=stats::rbinom(N_max,1,p_c)
      #x0_p2 is the random sample at from phase 2 stagefor the control
      x=matrix(stats::rbinom(N_max*comp,1,p_t), nrow=comp,ncol=N_max)
      #x is the matrix of random samples for each dose from the phase 2 stage
      #for each simulation.
      #x matirx does not include the control sample.
      #The row number is the same as the active dose

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      S_ph2=NULL
      wald_ph2=NULL
      if (direction ==1){ #higher rate p is better
        for (j in 1:length(p_t)){
          S_ph2[j]=(mean(x[j,1:N1])-mean(x0[1:N1]))*t[j,1]
          wald_ph2[j]=(mean(x[j,1:N1])-mean(x0[1:N1]))*sqrt(t[j,1])
        }
      }else if (direction ==2){ #lower rate p is better
        for (j in 1:length(p_t)){
          S_ph2[j]=(mean(x0[1:N1])-mean(x[j,1:N1]))*t[j,1]
          wald_ph2[j]=(mean(x0[1:N1])-mean(x[j,1:N1]))*sqrt(t[j,1])
        }
      }
      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      # s_slct=match(max(S_ph2),S_ph2)
      s_slct=match(max(wald_ph2),wald_ph2)
      #    s_slct=match(wald_ph2,max(wald_ph2))

      T_slct<- samsz/(p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c))
      t_slct=T_slct*s
      t_p3=t_slct[-1]-t_slct[1]
      #t_p3=(T_slct-t_slct[1])*s_p3
      #s_slct identifies which dose was selected for the i-th simulation
      ###########################################################
      ### Step 4: generate  random samples for phase 3 stage ####
      ###########################################################

      x0_two_stg=x0[1:NF]
      x0_p3=x0[(N1+1):NF]
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_slct_two_stg=x[s_slct,1:NF]
      x_slct_p3=x[s_slct,(N1+1):NF]

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################

      theta_hat=NULL
      wald_two_stg=NULL
      for (j in 1:length(s) ){

        theta_hat[j]=mean(x_slct_two_stg[1:N[j]])-mean(x0_two_stg[1:N[j]])  # theta_hat at each interim analysis
        wald_two_stg[j]=theta_hat[j]*sqrt(t_slct[j])
      }
      theta_hat_sd=1/sqrt(t_slct)
      theta_hat_two_stg=theta_hat
      theta_hat_two_stg_sd=theta_hat_sd
      ######  calculating the statistics for phase 3 only #######
      theta_hat=NULL  #theta_hat for phase 3 only
      for (j in 1:length(s)){
        theta_hat[j]= mean(x[s_slct,1:N[j]])-mean(x0[1:N[j]])  # theta_hat for phase 3 only, at each interim analysis
      }
      #theta_hat_p3_sd=1/sqrt(t_p3)
      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test_1=test[-length(s)]  ## test statistics excluding the final analysis

      ##########################################
      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the Wald statistics crossed the critical boundaries
      ############################################
      if (isTRUE(theta_hat[length(s)-1]<theta_cut)) {futi=1} else {futi=0}
      #If the NUll hypothesis is rejected at up to length(s)-1, terminate trial and estimate CI ####

      #############################

      if (isTRUE(max(test)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
      }
      ##############

      if  (isTRUE(max(test_1)>=0)){ #null hypothesis rejected, no more random sampling ot theta_hat at length(s_p3)-1 is almost negative


        #  if ( max(tes_1)>0){ #null hypothesis rejected, no more random sampling
        rej=1
        #If max(tes_1)>0, te nulll hypothesis is rejected
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=1/sqrt(t_slct[last_vt])
        ###########################
        sampsz_simu=N[last_vt]
        if (last_vt==1){
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*sampsz_simu
        }else {
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*N[1]+2*(sampsz_simu-N[1])
        }

      }else if  (isTRUE(max(test_1)<0)){ ## (theta_p3_hat[length(s_1)-1]<theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation

        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat_two_stg[last_vt]
          theta_hat_last_sd=1/sqrt(t_slct[last_vt])

          sampsz_simu=N[length(s)]
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*N[1]+2*(sampsz_simu-N[1])
          if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }

          ############################

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat_two_stg[inter_vt]
          theta_hat_inter_sd=theta_hat_two_stg_sd[inter_vt]
          t_inter=t_slct[inter_vt]
          score_inter=theta_hat_inter*t_inter
          c_bry=c_bry_two_stg

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t_slct[length(s)]-t_slct[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t_slct[length(s)])-score_inter)/sqrt(t_slct[length(s)]-t_slct[inter_vt]))

          if (isTRUE(cond_power>=1-beta)){
            last_vt=length(s)
            theta_hat_last=theta_hat_two_stg[last_vt]
            theta_hat_last_sd=1/sqrt(t_slct[last_vt])


            sampsz_simu=N[length(s)]
            control_sampsz_simu=sampsz_simu
            tot_sampsz_simu=(comp+1)*N[1]+2*(sampsz_simu-N[1])
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (isTRUE(cond_power<1-beta)){

            ## proceed to sample size re-estimation #####

            ######### if the NULL hypothesis is not rejected at up to length(s)-1, and theta_p3_hat[length(s_1)-1]>=theta_cut, indicating potential efficacy ###
            ## re-estimate sample size#####
            #sample size re-estimation procedure##
            ## statistic sfor combined phae 2/3 data#####
            inter_vt=length(s)-1
            theta_hat_inter=theta_hat[inter_vt]
            theta_hat_inter_sd=1/sqrt(t_slct[inter_vt])

            score_inter=theta_hat_inter*t_slct[inter_vt]
            #########################
            T_new_0=theta_hat_inter^(-2)*((e_bry_two_stg[length(s)]*sqrt(T_slct)-score_inter)/sqrt(T_slct-t_slct[inter_vt])+stats::qnorm(1-beta))^2+t_slct[inter_vt]
            new_sampsz=ceiling(T_new_0*((p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c))))
            T_new=new_sampsz/((p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c)))
            e_bry_new_two_stg=
              (sqrt(T_new-t_slct[inter_vt])/sqrt(T_slct-t_slct[inter_vt]))*(e_bry_two_stg[length(s)]*sqrt(T_slct)-score_inter)+ score_inter
            c_bry_new_two_stg=e_bry_new_two_stg/sqrt(T_new)
            ###############################################
            if (isTRUE(new_sampsz>N_max)){ ## new sample size is too large, the sample size will not change. Trial ends at the planned sample size
              if (isTRUE(test[length(test)]>0)){ ### the NULL hypothesis is rejected at the final anbalysis, with the planned sample size
                rej=1
                stg=length(s)
              }else{### the NULL hypothesis is NOT rejected at the final anbalysis, with the planned sample size
                rej=0
                stg=0
              }
              ###################### Statistics for estimates ####
              ### combined phase 2/3 data ####
              last_vt=length(s)
              sampsz_simu=N[last_vt]
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              ####################### this completes the estimation when there is no saple size change ########

            } else if (isTRUE(new_sampsz<= N_max)) {## new sample size does not exceed the pre-selected maximum. Will change the sample  size to this new sample size
              #here
              sampsz_simu=new_sampsz
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
              x0_two_stg_new=x0[1:new_sampsz] # combine the previous [combined phase 2/3] data of the control arm (up to inter_vt), and the newly generated additional random sample

              x_two_stg_new=x[s_slct,1:new_sampsz] # combine the previous [combined phase 2/3] data of the selectd dose arm (up to inter_vt), and the newly generated additional random sample

              T_new=new_sampsz/(p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c))   # the new information time for the combined phase 2/3 test


              ### statistics from combined data ####

              if (direction ==1){ #higher rate p is better
                wald_last_ad=(mean(x_two_stg_new)-mean(x0_two_stg_new))*sqrt(T_new)  # Wald statistics at the adjusted final analysis
              } else   if (direction ==2){ #lower rate p is better
                wald_last_ad=(mean(x0_two_stg_new)-mean(x_two_stg_new))*sqrt(T_new)  # Wald statistics at the adjusted final analysis
              }
              ##########################################

              test_add=wald_last_ad-c_bry_new_two_stg  # compare the final Wald statistics with the new adjusted critical boundary

              if (isTRUE(test_add>0)){ # the new critical boundary is crossed, the NULL is rejected

                rej=1
                stg=length(s)

              }else{  # the new critical boundary is NOT crossed, the NULL is rejected

                rej=0
                stg=0

              }

              ###################### Parameter estimate with sample size change ##################

            }
          }
        }
      }
      out=c(rej_nossr,rej,stg,s_slct,control_sampsz_simu,tot_sampsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  control_samsz_simu = as.numeric(simr[,5])
  tot_samsz_simu = as.numeric(simr[,6])
  futi=as.numeric(simr[,7])
  #############termination timepoint####
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  samsz_icr=NULL
  for (i in 1:length(control_samsz_simu)){
    if (control_samsz_simu[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  futi_rt=mean(futi)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= control_samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  #############################
  ### Step 9: Final output ####
  #############################
  tyI_er=sum(rej[mu[s_slct]==0])/sim_num
  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/sim_num
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(p_t)-1)+2*samsz
  sim_results=c(rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er, non_tyI_er_rej_rt,futi_rt,samsz,
                N_max,mean_control_samsz,tot_samsz,mean_tot_samsz, mean_samsz_simu_control_samsz_icr)
  names(sim_results)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned control sample size","Maximum control sample size","Mean control sample size",
                       "Planned total sample size","Mean total sample size","mean control sample size conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=sim_results
  return(out)

}

#################################################
#Two_stage_adapt_rej_OC_simu_binary(p_t,p_c,direction,s,c_bry_two_stg,alpha,beta,theta_cut,samsz,N_max,sim_num,dist)
#Two_stage_adapt_rej_simu_binary(p_t=c(0.9,0.8,0.7),p_c=0.6,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000,dist=1)
#Two_stage_adapt_rej_simu_binary(p_t=c(0.6,0.6,0.6),p_c=0.6,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000,dist=1)
#Two_stage_adapt_rej_OC_simu_binary(p_t=c(0.6,0.6,0.6),p_c=0.6,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000,dist=1)

#Two_stage_adapt_rej_simu_binary(p_t=c(0.7,0.6,0.6),p_c=0.6,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000,dist=1)
#Two_stage_adapt_rej_simu_binary(p_t=c(0.7,0.6,0.6),p_c=0.6,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),alpha=0.025,beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000,dist=1)

#############################################
two_stg_adapt_power_OC_simu_binary<-function(p_t,p_c,direction,s,bry_type,alpha_2,theta_cut,samsz,N_max,beta,sim_num){
  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better

  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  # mu=mu_t-mu_c
  comp=length(p_t)
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))

  simu_out=Two_stage_adapt_rej_OC_simu_binary(p_t=p_t,p_c=p_c,direction=direction,
                                              s=s,c_bry_two_stg=c_bry_two_stg,alpha=alpha,beta=beta,theta_cut=theta_cut,samsz=samsz,N_max=N_max,sim_num=sim_num,dist=1)

  p_t_names=paste(paste("p<sub>t</sub>","[",c(1:length(p_t)),"]",sep=""))
  out_1= round(c(p_t,p_c,as.numeric(simu_out)),4)
  names(out_1)=c(p_t_names,"p<sub>c</sub>",names(simu_out))
  # out=data.frame(out_1)
  out=out_1
  return(out)
}
###############################
#two_stg_adapt_power_OC_simu_binary(p_t,p_c,direction,s,bry_type,alpha_2,theta_cut,samsz,N_max,beta,sim_num)
###########################
#two_stg_adapt_power_OC_simu_binary(p_t=c(0.4,0.5,0.6),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=200,beta=0.1,sim_num=1000)

#two_stg_adapt_power_OC_simu_binary(p_t=c(0.6,0.6,0.6),p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=200,beta=0.1,sim_num=100000)


##############################################################

###################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Binary distribution, Optimizing adaptive sequential design
#'
#' @param groups Matrix of assumed event rates for the test arms under study: one row per candidate test-arm scenario and one column per dose, so that the number of doses is the number of columns; each event rate must be > 0 and < 1
#' @param p_c Event rate of the control arm; must be > 0 and < 1
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s1 Vector of information fractions of the first interim analysis schedule; the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule; the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s1; each element must be > 0 and the elements must sum to alpha_2
#' @param alpha_2_alps_s2 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s2; each element must be > 0 and the elements must sum to alpha_2
#' @param theta_cut Lack-of-efficacy threshold: an effect size below which the test drug is not considered to provide clinically meaningful efficacy
#' @param samsz_1 Planned sample size of the control arm for the first sample size scenario
#' @param samsz_2 Planned sample size of the control arm for the second sample size scenario
#' @param N_max_1 Maximum sample size (control arm) allowed after sample size re-estimation for the first maximum sample size scenario
#' @param N_max_2 Maximum sample size (control arm) allowed after sample size re-estimation for the second maximum sample size scenario
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with O'Brien-Fleming boundaries}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with O'Brien-Fleming boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with alpha-spending boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with alpha-spending boundaries}
#' }
#' Each of these four lists has the components "Group sequential design" (data frame with the columns "Looks", "Information fractions" and the critical boundaries, plus "Alpha-spending" and "Alpha-spending boundaries" for the alpha-spending schedules) and the "Simulation results: Information fractions ..., Sample size ..., N_max ..." operating characteristics of the corresponding adaptive design, one data frame per combination of interim analysis schedule, planned sample size and maximum sample size, with one row per assumed event-rate scenario.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' Two_stage_asd_simu_OC_binary(
#' groups = rbind(c(0.3,0.4,0.5),c(0.1,0.2,0.3),c(0.1,0.2,0.3),c(0.1,0.2,0.3)), p_c = 0.2,
#' direction = 1, s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05,
#' alpha_2_alps_s1 = c(0.001,0.001,0.048), alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01,
#' samsz_1 = 100, samsz_2 = 150, N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
Two_stage_asd_simu_OC_binary<-function(groups,p_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
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

  s1_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=dist)[1]))
  s2_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=dist)[1]))

  s1_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s1,alpha=alpha_alsp_s1,comp=comp,dist=dist)[1]))
  s2_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s2,alpha=alpha_alsp_s2,comp=comp,dist=dist)[1]))

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
    simu_out_OF_s1_samsz_1_N_max_1=rbind(simu_out_OF_s1_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_1)=p_t_group_names

  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_1=rbind(simu_out_OF_s2_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_1)=p_t_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_1=rbind(simu_out_OF_s1_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_1)=p_t_group_names

  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_1=rbind(simu_out_OF_s2_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
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
    simu_out_OF_s1_samsz_1_N_max_2=rbind(simu_out_OF_s1_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_2)=p_t_group_names

  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_2=rbind(simu_out_OF_s2_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_2)=p_t_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_2=rbind(simu_out_OF_s1_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_2)=p_t_group_names

  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_2=rbind(simu_out_OF_s2_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
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
    simu_out_alsp_s1_samsz_1_N_max_1=rbind(simu_out_alsp_s1_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_1)=p_t_group_names

  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_1=rbind(simu_out_alsp_s2_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                           bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_1)=p_t_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_1=rbind(simu_out_alsp_s1_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_1)=p_t_group_names

  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_1=rbind(simu_out_alsp_s2_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
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
    simu_out_alsp_s1_samsz_1_N_max_2=rbind(simu_out_alsp_s1_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_2)=p_t_group_names

  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_2=rbind(simu_out_alsp_s2_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
                                                                                                           bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_2)=p_t_group_names
  ## samsz_2###

  ##### s1 ####

  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_2=rbind(simu_out_alsp_s1_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s1,
                                                                                                           bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                           beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_2)=p_t_group_names

  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(p_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_2=rbind(simu_out_alsp_s2_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_binary(p_t=p_t_vec[i,],p_c=p_c,direction=direction,s=s2,
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
#Two_stage_asd_simu_OC_normal(mu_t_vec_1,mu_t_vec_2,mu_c,sigma_t,sigma_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num,dist)
#Two_stage_asd_simu_OC_binary(groups=rbind(c(0.3,0.4,0.5),c(0.1,0.2,0.3),c(0.1,0.2,0.3),c(0.1,0.2,0.3)),p_c=0,direction=1,
#                             s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                            theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1, sim_num=1000)


#############################################
###############################################################


##########################################################################
### II.4 .Simulations  -Poisson distributions                         #####
##########################################################################


##########################################################################
### II.4.1 .Simulations with planned TSSSD -Poisson distributions     #####
##          Type I error, power,  simulations                        #####
#### This will be included in the software                           #####
##########################################################################

#########################################
two_stg_rej_rate_simu_poisson=function(lambda_t,lambda_c,direction,s,c_bry_two_stg,samsz,sim_num){
  # direction==1 means that larger rate lambda is better
  # direction==2 means that smaller rate lambda  is better

  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsaz is sample size

  N=ceiling(samsz*s)
  comp=length(lambda_t)

  if (direction ==1){ #higher rate is better
  mu=lambda_t-lambda_c
  }else if (direction ==2){ #lower rate is better
    mu=lambda_c-lambda_t
  }

  ##########################
  N1=N[1]
  NF=samsz
  #t_two_stg is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "2ASD", function() {

      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      x0=stats::rpois(NF,lambda_c)

      x0_p2=x0[1:N1]
      #x0 is the random sample at the first interim analysis for the control
      x=matrix(nrow=comp,ncol=NF)
      for (i in 1:comp){
      x[i,]=stats::rpois(NF,lambda_t[i])
      }

      x_p2=x[,1:N1]

      #x is the matrix of random samples for each dose from the phase 2 stage
      #for each simulation.
      #x matirx does not include the control sample.
      #The row number is the same as the active dose
      lambda_c_hat_p2=mean(x0_p2)
      lambda_t_hat_p2=NULL
      t_hat_p2=NULL
      theta_hat_p2=NULL

      if (direction ==1){ #higher rate is better
      for (i in 1:comp ){

        lambda_t_hat_p2[i]=mean(x_p2[i,])
        t_hat_p2[i]= N1/(lambda_t_hat_p2[i]+lambda_c_hat_p2)
        theta_hat_p2[i]=lambda_t_hat_p2[i]-lambda_c_hat_p2
      }
      }else  if (direction ==2){ #lower rate is better
        for (i in 1:comp ){

          lambda_t_hat_p2[i]=mean(x_p2[i,])
          t_hat_p2[i]= N1/(lambda_t_hat_p2[i]+lambda_c_hat_p2)
          theta_hat_p2[i]=lambda_c_hat_p2-lambda_t_hat_p2[i]
        }
      }

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################
      S_p2=NULL
      for (i in 1:comp){
        S_p2[i]=theta_hat_p2[i]*t_hat_p2[i]
      }
      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(max(S_p2),S_p2)
     # t_two_stg=N/(p_t[s_slct]*(1-p_t[s_slct])+p_c*(1-p_c))
      #s_slct identifies which dose was selected for the i-th simulation
      ###########################################################
      ### Step 4: generate  random samples for phase 3 stage ####

      x0_two_stg=x0
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_two_stg=x[s_slct,]
      #x_two_stg contains all the random samples(including samples from the phase 2 stage) for the seected dose

      lambda_c_hat_two_stg=NULL
      lambda_t_hat_two_stg=NULL
      theta_hat_two_stg=NULL
      t_two_stg=NULL

      if (direction ==1){ #higher rate is better

      for (i in 1:length(s)){
      lambda_c_hat_two_stg[i]=mean(x0_two_stg[1:N[i]])
      lambda_t_hat_two_stg[i]=mean(x_two_stg[1:N[i]])
      t_two_stg[i]=N[i]/(lambda_c_hat_two_stg[i]+lambda_t_hat_two_stg[i])
      theta_hat_two_stg[i]=lambda_t_hat_two_stg[i]-lambda_c_hat_two_stg[i]
      }

      } else   if (direction ==2){ #lower rate is better

        for (i in 1:length(s)){
          lambda_c_hat_two_stg[i]=mean(x0_two_stg[1:N[i]])
          lambda_t_hat_two_stg[i]=mean(x_two_stg[1:N[i]])
          t_two_stg[i]=N[i]/(lambda_c_hat_two_stg[i]+lambda_t_hat_two_stg[i])
          theta_hat_two_stg[i]=lambda_c_hat_two_stg[i]-lambda_t_hat_two_stg[i]
        }
      }

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      wald_two_stg=NULL
      for (i in 1:length(s) ){

        wald_two_stg[i]=theta_hat_two_stg[i]*sqrt(t_two_stg[i])
        #wald_two_stg containes the vector Wald statistics  for the  simulation
      }
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if (isTRUE(max(test)>0)){

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        last_vt=stg
        control_samsz_simu=N[last_vt]
        tot_samsz_simu=N[1]*comp+control_samsz_simu

        #stg identifies the interim analysis at which the null hypothesis is rejected

      }else{

        rej=0
        stg=0
        last_vt=length(s)
        control_samsz_simu=N[last_vt]
        tot_samsz_simu=N[1]*(comp+1)+2*(control_samsz_simu-N[1])
      }
      #######################################################
      ### Step 7: output and finish the simjulation loop ####
      #######################################################
      out=c(rej,stg,s_slct,control_samsz_simu,tot_samsz_simu)
  })#(end of q loop)

  ############################################
  ### Step 8: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  s_slct= as.numeric(simr[,3])
  control_samsz_simu= as.numeric(simr[,4])
  tot_samsz_simu= as.numeric(simr[,5])
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=sum(rej[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  select_rt=rep(0,comp)
  for (j in 1:comp){
    select_rt[j]=length( s_slct[ s_slct==j])/sim_num
  }

  rej_rt_indiv=rep(0,comp)
  for (j in 1:comp){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/sim_num
  }

  con_rej_rt_indiv=rep(0,comp)
  for (j in 1:comp){
    if (isTRUE(select_rt[j]==0)){
      con_rej_rt_indiv[j]=NA
    }else {
      con_rej_rt_indiv[j]=rej_rt_indiv[j]/select_rt[j]
    }
  }

  results_by_arm=cbind(lambda_t,lambda_c,select_rt, rej_rt_indiv,con_rej_rt_indiv)
  colnames(results_by_arm)=c(LABEL()$lambda_t1,LABEL()$lambda_c,"selection rate","rejection rate by each lambda_t","conditional rejection rate")

  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look")

  tyI_er=sum(rej[mu[s_slct]==0])/sim_num

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/sim_num
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(lambda_t)-1)+2*samsz
  sim_results=c(rej_rt,tyI_er, non_tyI_er_rej_rt,samsz,mean_control_samsz,tot_samsz,mean_tot_samsz)
  names(sim_results)=c("rejection rate","Type I error","Non-Type I error rejection rate",
                       "Planned control sample size","Mean control sample size","Planned total sample size","Mean total sample size")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(seq_design,results_by_arm,sim_results)
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)
}
#######################################
#two_stg_rej_rate_simu_poisson(lambda_c,lambda_t,direction,s,c_bry_two_stg,samsz,sim_num)
#two_stg_rej_rate_simu_poisson(lambda_t=c(0.6,0.7,0.8),lambda_c=1.1,direction,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),samsz=100,sim_num=10000)
#two_stg_rej_rate_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),samsz=100,sim_num=100000)
#two_stg_rej_rate_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),samsz=100,sim_num=100000)


 ###############Estimates#####################


##############################################

#' Phase 2/3 seamless combination Group sequential design: Simulations, Poisson distribution, Type I error and power
#'
#' @param lambda_t Mean of the test arm; each element must be > 0
#' @param lambda_c Mean of the control arm; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param samsz Planned sample size of the control arm used in the simulation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries" and "rejection by look"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "lambda_t1,...", "lambda_c", "selection rate", "rejection rate by each lambda_t" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the overall "rejection rate", the "Type I error", the "Non-Type I error rejection rate", the "Planned control sample size", the "Mean control sample size", the "Planned total sample size" and the "Mean total sample size"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_power_simu_poisson(lambda_t = c(1.1,1.3,1.4), lambda_c = 1.1, direction = 1,
#' s = c(1/3,1), bry_type = 1, alpha_2 = 0.05,
#' samsz = 100, sim_num = 20)
#'
#' @export
two_stg_power_simu_poisson<-function(lambda_t,lambda_c,direction,s,bry_type,alpha_2,samsz,sim_num){
  # direction==1 means that larger rate lambda is better
  # direction==2 means that smaller rate lambda  is better
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=length(lambda_t)
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s,alpha,comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s,alpha,comp,dist=1)
  }
  tot_samsz=ceiling((samsz*s[1])*(comp+1)+2*(samsz*(s[length(s)]-s[1])))
  c_bry_two_stg=as.numeric(unlist(design_out[1]))
  simu_out=two_stg_rej_rate_simu_poisson(lambda_t,lambda_c,direction,s,c_bry_two_stg,samsz,sim_num)
  out= simu_out
 # names(out)=c("simulation summary","planned control sample size ","planned total sample size")
  return(out)
}

##################################
#two_stg_power_simu_poisson(lambda_t,lambda_c,s,bry_type,alpha_2,samsz,sim_num)
#two_stg_power_simu_poisson(lambda_t=c(1.1,1.3,1.4),lambda_c=1.1,direction=1,s=c(1/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)
#two_stg_power_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,1),bry_type=2,alpha_2=c(0.001,0.049),samsz=100,sim_num=100000)

#two_stg_power_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=100000)
#two_stg_power_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),samsz=100,sim_num=100000)


#######################################


###########################################################################
### II.4 .Simulations with adaptive TSSSD -Poisson  distributions     #####
##          Type I error, power, simulations                          #####
##########################################################################

#############################################
Two_stage_adapt_rej_simu_poisson<-function(lambda_t,lambda_c,direction,s,c_bry_two_stg,beta,theta_cut,samsz,N_max,sim_num){

  # direction==1 means that larger rate lambda is better
  # direction==2 means that smaller rate lambda  is better
  ### parameters for two-stage #########

  N=ceiling(samsz*s)

  comp=length(lambda_t)

  N1=N[1]
  #N1 is the sample at the end of phase 2 stage
  NF=N[length(s)]


  if (direction ==1){ #higher rate is better
    mu=lambda_t-lambda_c
  }else if (direction ==2){ #lower rate is better
    mu=lambda_c-lambda_t
  }

  ######################## exit boundaries##########
 # e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined

  ############ the simulation loops ##########
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=N[length(s)]
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    control_sampsz_simu=sampsz_simu; tot_sampsz_simu=(comp+1)*sampsz_simu
    futi=0
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################
      x0=stats::rpois(N_max,lambda_c)
      #x0 is the random sample at the first interim analysis for the control
      x=matrix(nrow=comp,ncol=N_max)
      for (i in 1:comp){
        x[i,]=stats::rpois(N_max,lambda_t[i])
      }

      lambda_c_hat_p2=mean(x0[1:N1])
      lambda_t_hat_p2=NULL
      t_hat_p2=NULL
      theta_hat_p2=NULL

      if (direction ==1){ #higher rate is better
      for (i in 1:comp ){

        lambda_t_hat_p2[i]=mean(x[i,1:N1])
        t_hat_p2[i]= N1/(lambda_t_hat_p2[i]+lambda_c_hat_p2)
        theta_hat_p2[i]=lambda_t_hat_p2[i]-lambda_c_hat_p2
      }
      }else if (direction ==2){ #lower rate is better
        for (i in 1:comp ){

          lambda_t_hat_p2[i]=mean(x[i,1:N1])
          t_hat_p2[i]= N1/(lambda_t_hat_p2[i]+lambda_c_hat_p2)
          theta_hat_p2[i]=lambda_c_hat_p2-lambda_t_hat_p2[i]
        }
      }

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################
      S_p2=NULL
      wald_p2=NULL
      for (i in 1:comp){
        S_p2[i]=theta_hat_p2[i]*t_hat_p2[i]
        wald_p2[i]=theta_hat_p2[i]*sqrt(t_hat_p2[i])
      }
      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(max(S_p2),S_p2)

      x0_two_stg=x0[1:NF]
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_two_stg=x[s_slct,1:NF]
      #x_two_stg contains all the random samples(including samples from the phase 2 stage) for the seected dose

      lambda_c_hat_two_stg=NULL
      lambda_t_hat_two_stg=NULL
      theta_hat_two_stg=NULL
      theta_hat_sd_two_stg=NULL

      if (direction ==1){ #higher rate is better
        for (i in 1:length(s)){
          lambda_c_hat_two_stg[i]=mean(x0_two_stg[1:N[i]])
          lambda_t_hat_two_stg[i]=mean(x_two_stg[1:N[i]])
          theta_hat_two_stg[i]=lambda_t_hat_two_stg[i]-lambda_c_hat_two_stg[i]
        }
      } else  if (direction ==2){ #lower rate is better
        for (i in 1:length(s)){
          lambda_c_hat_two_stg[i]=mean(x0_two_stg[1:N[i]])
          lambda_t_hat_two_stg[i]=mean(x_two_stg[1:N[i]])
          theta_hat_two_stg[i]=lambda_c_hat_two_stg[i]-lambda_t_hat_two_stg[i]
        }
      }
     T_two_stg=N[length(s)]/(lambda_t_hat_two_stg[length(s)]+lambda_c_hat_two_stg[length(s)])

      t_two_stg=N/(lambda_t_hat_two_stg[length(s)]+lambda_c_hat_two_stg[length(s)])
      inter_vt=length(s)-1

      t_inter=(N[inter_vt]/N[length(s)])*T_two_stg
      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      wald_two_stg=NULL
      for (i in 1:length(s) ){

        wald_two_stg[i]=theta_hat_two_stg[i]*sqrt(t_two_stg[i])
        #wald_two_stg containes the vector Wald statistics  for the  simulation
      }
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg

      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test_1=test[-length(s)]  ## test statistics excluding the final analysis

      ##########################################

      ############################################
      if (isTRUE(theta_hat_two_stg[inter_vt]<theta_cut)) {futi=1} else {futi=0}
        if (isTRUE(max(test)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
      }
      if  (isTRUE(max(test_1)>=0)){ #null hypothesis rejected, no more random sampling ot theta_hat at length(s_p3)-1 is almost negative
        #  if ( max(tes_1)>0){ #null hypothesis rejected, no more random sampling
        rej=1
        #If max(tes_1)>0, te nulll hypothesis is rejected
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        sampsz_simu=N[last_vt]
        if (last_vt==1){
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*sampsz_simu
        }else {
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
        }

      }else if (isTRUE(max(test_1)<0)){

        if (futi==1){

          last_vt=length(s)
          theta_hat_last=theta_hat_two_stg[last_vt]
          theta_hat_last_sd=theta_hat_sd_two_stg[last_vt]

          sampsz_simu=N[length(s)]
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
          if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }

         #############  #############

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat_two_stg[inter_vt]
          theta_hat_inter_sd=theta_hat_sd_two_stg[inter_vt]
          score_inter=theta_hat_inter*t_two_stg[inter_vt]
          c_bry=c_bry_two_stg

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t_two_stg[length(s)]-t_two_stg[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t_two_stg[length(s)])-score_inter)/sqrt(t_two_stg[length(s)]-t_two_stg[inter_vt]))

          if (isTRUE(cond_power>=1-beta)){
            last_vt=length(s)
            theta_hat_last=theta_hat_two_stg[last_vt]
            theta_hat_last_sd=theta_hat_sd_two_stg[last_vt]

            sampsz_simu=N[length(s)]
            control_sampsz_simu=sampsz_simu
            tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (isTRUE(cond_power<1-beta)){
            ## proceed to sample size re-estimation #####
            inter_vt=length(s)-1
            theta_hat_inter=theta_hat_two_stg[inter_vt]
            theta_hat_inter_sd=theta_hat_sd_two_stg[inter_vt]

            score_inter=theta_hat_inter*t_two_stg[inter_vt]

            T_new_0=theta_hat_inter^(-2)*((c_bry_two_stg[length(s)]*sqrt(T_two_stg)-score_inter)/sqrt(T_two_stg-t_two_stg[inter_vt])+stats::qnorm(1-beta))^2+t_two_stg[inter_vt]
            new_sampsz=ceiling(T_new_0*(lambda_c_hat_two_stg[length(s)]+lambda_t_hat_two_stg[length(s)]) )

            ###############################################
            if (isTRUE(new_sampsz>N_max)){ ## new sample size is too large, the sample size will not change. Trial ends at the planned sample size
              if (isTRUE(test[length(test)]>0)){ ### the NULL hypothesis is rejected at the final anbalysis, with the planned sample size
                rej=1
                stg=length(s)
              }else{### the NULL hypothesis is NOT rejected at the final anbalysis, with the planned sample size
                rej=0
                stg=0
              }

              ### combined phase 2/3 data ####
              last_vt=length(s)
              sampsz_simu=N[last_vt]
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu


            } else if (isTRUE(new_sampsz<= N_max)) {## new sample size does not exceed the pre-selected maximum. Will change the sample  size to this new sample size

              sampsz_simu=new_sampsz
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              x0_two_stg_new=x0[1:new_sampsz] # combine the previous [combined phase 2/3] data of the control arm (up to inter_vt), and the newly generated additional random sample

              x_two_stg_new=x[s_slct,1:new_sampsz]  # combine the previous [combined phase 2/3] data of the selectd dose arm (up to inter_vt), and the newly generated additional random sample

            ####################

              lambda_c_hat_last_ad=mean(x0_two_stg_new)
              lambda_t_hat_last_ad=mean(x_two_stg_new)
              t_hat_last_ad=new_sampsz/(lambda_c_hat_last_ad+lambda_t_hat_last_ad)

              theta_hat_sd_last_ad=1/sqrt(t_hat_last_ad)

              if (direction ==1){ #higher rate is better
                theta_hat_last_ad=lambda_t_hat_last_ad-lambda_c_hat_last_ad
              }else if (direction ==2){ #lower rate is better
                theta_hat_last_ad=lambda_c_hat_last_ad-lambda_t_hat_last_ad
              }

              wald_last_ad=theta_hat_last_ad*sqrt(t_hat_last_ad)  # Wald statistics at the adjusted final analysis

              T_new=new_sampsz/(lambda_c_hat_last_ad+lambda_t_hat_last_ad)
              T_final_1=NF/(lambda_c_hat_last_ad+lambda_t_hat_last_ad)
              t_inter=N[inter_vt]/(lambda_c_hat_last_ad+lambda_t_hat_last_ad)
              e_bry_new_two_stg=
                (sqrt(T_new-t_inter)/sqrt(T_final_1-t_inter))*(c_bry_two_stg[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter
              c_bry_new_two_stg=e_bry_new_two_stg/sqrt(T_new)

              ##########################################

              test_add=wald_last_ad-c_bry_new_two_stg  # compare the final Wald statistics with the new adjusted critical boundary

              if (isTRUE(test_add>0)){ # the new critical boundary is crossed, the NULL is rejected

                rej=1
                stg=length(s)

              }else{  # the new critical boundary is NOT crossed, the NULL is rejected

                rej=0
                stg=0

              }
              ###################### Parameter estimate with sample size change ##################
            }
          }
        }
      }
      out=c(rej_nossr,rej,stg,s_slct,control_sampsz_simu,tot_sampsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  control_samsz_simu = as.numeric(simr[,5])
  tot_samsz_simu = as.numeric(simr[,6])
  futi=as.numeric(simr[,7])

  #############termination timepoint####
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  rej_time=rep(0,(length(s)))
  for (j in 1:(length(s))){

    count_j=c(1:sim_num)[stg==j]
    rej_time[j]=sum(rej[count_j])/(dim(simr)[1])
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  select_rt=rep(0,comp)
  for (j in 1:comp){
    select_rt[j]=length( s_slct[ s_slct==j])/(dim(simr)[1])
  }

  rej_rt_indiv=rep(0,comp)
  for (j in 1:comp){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/(dim(simr)[1])
  }

  con_rej_rt_indiv=rep(0,comp)
  for (j in 1:comp){
    if (isTRUE(select_rt[j]==0)){
      con_rej_rt_indiv[j]=NA
    }else {
      con_rej_rt_indiv[j]=rej_rt_indiv[j]/select_rt[j]
    }
  }

  samsz_icr=NULL
  for (i in 1:length(control_samsz_simu)){
    if (control_samsz_simu[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_rt_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  futi_rt=mean(futi)

  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= control_samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)
  #############################
  ### Step 9: Final output ####
  #############################

  results_by_arm=cbind(lambda_t,lambda_c,select_rt, rej_rt_indiv,con_rej_rt_indiv)
  colnames(results_by_arm)=c(LABEL()$lambda_t1,LABEL()$lambda_c,"selection rate","rejection rate by each lambda_t","conditional rejection rate")

  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time,samsz_icr_rt_byLook)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look","Sample size increase")

  tyI_er=sum(rej[mu[s_slct]==0])/(dim(simr)[1])
  #  tyI_er=sum(rej_0)/sim_num

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/(dim(simr)[1])
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(lambda_t)-1)+2*samsz
  sim_results=c(rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er, non_tyI_er_rej_rt,futi_rt,samsz,
                N_max,mean_control_samsz,tot_samsz,mean_tot_samsz, mean_samsz_simu_control_samsz_icr)
  names(sim_results)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned control sample size","Maximum control sample size","Mean control sample size",
                       "Planned total sample size","Mean total sample size","mean control sample size conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(seq_design,results_by_arm,sim_results)
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)

}


#################################################
#Two_stage_adapt_rej_simu_poisson(lambda_t,lambda_c,s,c_bry_two_stg,beta,theta_cut,samsz,N_max,sim_num)
#Two_stage_adapt_rej_simu_poisson(lambda_t=c(0.9,0.8,0.7),lambda_c=1.1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000)
#Two_stage_adapt_rej_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000)
#Two_stage_adapt_rej_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000)


###############################################################

#################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Poisson distribution, Operating characteristics
#'
#' @param lambda_t Mean of the test arm; each element must be > 0
#' @param lambda_c Mean of the control arm; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param theta_cut Lack-of-efficacy threshold: an effect size below which the test drug is not considered to provide clinically meaningful efficacy
#' @param N_max Maximum sample size (control arm) allowed after sample size re-estimation
#' @param samsz Planned sample size of the control arm used in the simulation
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned adaptive group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries", "rejection by look" and "Sample size increase"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "lambda_t1,...", "lambda_c", "selection rate", "rejection rate by each lambda_t" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the "Rejection rate without SSR", the "Rejection rate with SSR", the "rate of sample size increase", the "Rejection rate without SSR/(meeting sample size incerase criteria)", the "Rejection rate with SSR/(sample size incerase)", the "Type I error", the "Non-Type I error rejection rate", the "Futility rate", the "Planned control sample size", the "Maximum control sample size", the "Mean control sample size", the "Planned total sample size", the "Mean total sample size" and the "mean control sample size conditional on increase"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_adapt_power_simu_poisson(lambda_t = c(1.1,1.1,1.1), lambda_c = 1.1, direction = 1,
#' s = c(1/3,1), bry_type = 1, alpha_2 = 0.05,
#' theta_cut = 0.05, N_max = 1000, samsz = 100, beta = 0.1, sim_num = 20)
#'
#' @export
two_stg_adapt_power_simu_poisson<-function(lambda_t,lambda_c,direction,s,bry_type,alpha_2,theta_cut,N_max,samsz,beta,sim_num){

  # direction==1 means that larger rate lambda is better
  # direction==2 means that smaller rate lambda  is better
  alpha=alpha_2/2
  comp=as.numeric(length(lambda_t))

  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))

  simu_out=Two_stage_adapt_rej_simu_poisson(lambda_t=lambda_t,lambda_c=lambda_c,direction=direction,s=s,c_bry_two_stg=c_bry_two_stg,beta=beta,theta_cut=theta_cut,samsz=samsz,N_max=N_max,sim_num=sim_num)

  out=simu_out
  return(out)
}

##################################
#two_stg_adapt_power_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,N_max=1000,samsz=100,beta=0.1,sim_num=1000)
#two_stg_adapt_power_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,N_max=1000,samsz=100,beta=0.1,sim_num=100000)
#two_stg_adapt_power_simu_poisson(lambda_t=c(1.1,1.3,1.4),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,N_max=1000,samsz=100,beta=0.1,sim_num=10000)

##############################################################
Two_stage_adapt_rej_OC_simu_poisson<-function(lambda_t,lambda_c,direction,s,c_bry_two_stg,beta,theta_cut,samsz,N_max,sim_num){

  # direction==1 means that larger rate lambda is better
  # direction==2 means that smaller rate lambda  is better
  ### parameters for two-stage #########

  N=ceiling(samsz*s)

  comp=length(lambda_t)

  N1=N[1]
  #N1 is the sample at the end of phase 2 stage
  NF=N[length(s)]


  if (direction ==1){ #higher rate is better
    mu=lambda_t-lambda_c
  }else if (direction ==2){ #lower rate is better
    mu=lambda_c-lambda_t
  }

  ######################## exit boundaries##########
  # e_bry_two_stg=c_bry_two_stg*sqrt(s)# exit boundary for phase 2/3 combined

  ############ the simulation loops ##########
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=N[length(s)]
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    control_sampsz_simu=sampsz_simu; tot_sampsz_simu=(comp+1)*sampsz_simu
    futi=0
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################
      x0=stats::rpois(N_max,lambda_c)
      #x0 is the random sample at the first interim analysis for the control
      x=matrix(nrow=comp,ncol=N_max)
      for (i in 1:comp){
        x[i,]=stats::rpois(N_max,lambda_t[i])
      }

      lambda_c_hat_p2=mean(x0[1:N1])
      lambda_t_hat_p2=NULL
      t_hat_p2=NULL
      theta_hat_p2=NULL

      if (direction ==1){ #higher rate is better
        for (i in 1:comp ){

          lambda_t_hat_p2[i]=mean(x[i,1:N1])
          t_hat_p2[i]= N1/(lambda_t_hat_p2[i]+lambda_c_hat_p2)
          theta_hat_p2[i]=lambda_t_hat_p2[i]-lambda_c_hat_p2
        }
      }else if (direction ==2){ #lower rate is better
        for (i in 1:comp ){

          lambda_t_hat_p2[i]=mean(x[i,1:N1])
          t_hat_p2[i]= N1/(lambda_t_hat_p2[i]+lambda_c_hat_p2)
          theta_hat_p2[i]=lambda_c_hat_p2-lambda_t_hat_p2[i]
        }
      }

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################
      S_p2=NULL
      wald_p2=NULL
      for (i in 1:comp){
        S_p2[i]=theta_hat_p2[i]*t_hat_p2[i]
        wald_p2[i]=theta_hat_p2[i]*sqrt(t_hat_p2[i])
      }
      #S_p2 contains the score function values for each dose vs. control
      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################
      s_slct=match(max(S_p2),S_p2)

      x0_two_stg=x0[1:NF]
      #x0_two_stg contains all the random samples(including samples from the phase 2 stage) for the control
      x_two_stg=x[s_slct,1:NF]
      #x_two_stg contains all the random samples(including samples from the phase 2 stage) for the seected dose

      lambda_c_hat_two_stg=NULL
      lambda_t_hat_two_stg=NULL
      theta_hat_two_stg=NULL
      theta_hat_sd_two_stg=NULL

      if (direction ==1){ #higher rate is better
        for (i in 1:length(s)){
          lambda_c_hat_two_stg[i]=mean(x0_two_stg[1:N[i]])
          lambda_t_hat_two_stg[i]=mean(x_two_stg[1:N[i]])
          theta_hat_two_stg[i]=lambda_t_hat_two_stg[i]-lambda_c_hat_two_stg[i]
        }
      } else  if (direction ==2){ #lower rate is better
        for (i in 1:length(s)){
          lambda_c_hat_two_stg[i]=mean(x0_two_stg[1:N[i]])
          lambda_t_hat_two_stg[i]=mean(x_two_stg[1:N[i]])
          theta_hat_two_stg[i]=lambda_c_hat_two_stg[i]-lambda_t_hat_two_stg[i]
        }
      }
      T_two_stg=N[length(s)]/(lambda_t_hat_two_stg[length(s)]+lambda_c_hat_two_stg[length(s)])

      t_two_stg=N/(lambda_t_hat_two_stg[length(s)]+lambda_c_hat_two_stg[length(s)])
      inter_vt=length(s)-1

      t_inter=(N[inter_vt]/N[length(s)])*T_two_stg
      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      wald_two_stg=NULL
      for (i in 1:length(s) ){

        wald_two_stg[i]=theta_hat_two_stg[i]*sqrt(t_two_stg[i])
        #wald_two_stg containes the vector Wald statistics  for the  simulation
      }
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg

      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test_1=test[-length(s)]  ## test statistics excluding the final analysis

      ##########################################

      ############################################
      if (isTRUE(theta_hat_two_stg[inter_vt]<theta_cut)) {futi=1} else {futi=0}
      ########################
      if (isTRUE(max(test)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
      }
      ######################
      if  (isTRUE(max(test_1)>=0)){ #null hypothesis rejected, no more random sampling ot theta_hat at length(s_p3)-1 is almost negative
        #  if ( max(tes_1)>0){ #null hypothesis rejected, no more random sampling
        rej=1
        #If max(tes_1)>0, te nulll hypothesis is rejected
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        ###########################
        sampsz_simu=N[last_vt]
        if (last_vt==1){
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*sampsz_simu
        }else {
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
        }

      }else if (isTRUE(max(test_1)<0)){

        if (futi==1){

          last_vt=length(s)
          theta_hat_last=theta_hat_two_stg[last_vt]
          theta_hat_last_sd=theta_hat_sd_two_stg[last_vt]

          sampsz_simu=N[length(s)]
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
          if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }

          #############  #############

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat_two_stg[inter_vt]
          theta_hat_inter_sd=theta_hat_sd_two_stg[inter_vt]
          score_inter=theta_hat_inter*t_two_stg[inter_vt]
          c_bry=c_bry_two_stg

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t_two_stg[length(s)]-t_two_stg[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t_two_stg[length(s)])-score_inter)/sqrt(t_two_stg[length(s)]-t_two_stg[inter_vt]))

          if (isTRUE(cond_power>=1-beta)){
            last_vt=length(s)
            theta_hat_last=theta_hat_two_stg[last_vt]
            theta_hat_last_sd=theta_hat_sd_two_stg[last_vt]

            sampsz_simu=N[length(s)]
            control_sampsz_simu=sampsz_simu
            tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (isTRUE(cond_power<1-beta)){
            ## proceed to sample size re-estimation #####
            inter_vt=length(s)-1
            theta_hat_inter=theta_hat_two_stg[inter_vt]
            theta_hat_inter_sd=theta_hat_sd_two_stg[inter_vt]

            score_inter=theta_hat_inter*t_two_stg[inter_vt]

            ###############
            T_new_0=theta_hat_inter^(-2)*((c_bry_two_stg[length(s)]*sqrt(T_two_stg)-score_inter)/sqrt(T_two_stg-t_two_stg[inter_vt])+stats::qnorm(1-beta))^2+t_two_stg[inter_vt]
            new_sampsz=ceiling(T_new_0*(lambda_c_hat_two_stg[length(s)]+lambda_t_hat_two_stg[length(s)]) )

            ###############################################
            if (isTRUE(new_sampsz>N_max)){ ## new sample size is too large, the sample size will not change. Trial ends at the planned sample size
              if (isTRUE(test[length(test)]>0)){ ### the NULL hypothesis is rejected at the final anbalysis, with the planned sample size
                rej=1
                stg=length(s)
              }else{### the NULL hypothesis is NOT rejected at the final anbalysis, with the planned sample size
                rej=0
                stg=0
              }

              ### combined phase 2/3 data ####
              last_vt=length(s)
              sampsz_simu=N[last_vt]
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu


            } else if (isTRUE(new_sampsz<= N_max)) {## new sample size does not exceed the pre-selected maximum. Will change the sample  size to this new sample size

              sampsz_simu=new_sampsz
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              x0_two_stg_new=x0[1:new_sampsz] # combine the previous [combined phase 2/3] data of the control arm (up to inter_vt), and the newly generated additional random sample

              x_two_stg_new=x[s_slct,1:new_sampsz]  # combine the previous [combined phase 2/3] data of the selectd dose arm (up to inter_vt), and the newly generated additional random sample

              ####################

              lambda_c_hat_last_ad=mean(x0_two_stg_new)
              lambda_t_hat_last_ad=mean(x_two_stg_new)
              t_hat_last_ad=new_sampsz/(lambda_c_hat_last_ad+lambda_t_hat_last_ad)

              theta_hat_sd_last_ad=1/sqrt(t_hat_last_ad)

              if (direction ==1){ #higher rate is better
                theta_hat_last_ad=lambda_t_hat_last_ad-lambda_c_hat_last_ad
              }else if (direction ==2){ #lower rate is better
                theta_hat_last_ad=lambda_c_hat_last_ad-lambda_t_hat_last_ad
              }

              wald_last_ad=theta_hat_last_ad*sqrt(t_hat_last_ad)  # Wald statistics at the adjusted final analysis

              T_new=new_sampsz/(lambda_c_hat_last_ad+lambda_t_hat_last_ad)
              T_final_1=NF/(lambda_c_hat_last_ad+lambda_t_hat_last_ad)
              t_inter=N[inter_vt]/(lambda_c_hat_last_ad+lambda_t_hat_last_ad)
              e_bry_new_two_stg=
                (sqrt(T_new-t_inter)/sqrt(T_final_1-t_inter))*(c_bry_two_stg[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter
              c_bry_new_two_stg=e_bry_new_two_stg/sqrt(T_new)

              ##########################################

              test_add=wald_last_ad-c_bry_new_two_stg  # compare the final Wald statistics with the new adjusted critical boundary

              if (isTRUE(test_add>0)){ # the new critical boundary is crossed, the NULL is rejected

                rej=1
                stg=length(s)

              }else{  # the new critical boundary is NOT crossed, the NULL is rejected

                rej=0
                stg=0

              }
              ###################### Parameter estimate with sample size change ##################
            }
          }
        }
      }
      out=c(rej_nossr,rej,stg,s_slct,control_sampsz_simu,tot_sampsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  control_samsz_simu = as.numeric(simr[,5])
  tot_samsz_simu = as.numeric(simr[,6])
  futi=as.numeric(simr[,7])

  #############termination timepoint####
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  samsz_icr=NULL
  for (i in 1:length(control_samsz_simu)){
    if (control_samsz_simu[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  futi_rt=mean(futi)

  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= control_samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)
  #############################
  ### Step 9: Final output ####
  #############################
  tyI_er=sum(rej[mu[s_slct]==0])/(dim(simr)[1])
  #  tyI_er=sum(rej_0)/sim_num

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/(dim(simr)[1])
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(lambda_t)-1)+2*samsz
  sim_results=c(rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er, non_tyI_er_rej_rt,futi_rt,samsz,
                N_max,mean_control_samsz,tot_samsz,mean_tot_samsz, mean_samsz_simu_control_samsz_icr)
  names(sim_results)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned control sample size","Maximum control sample size","Mean control sample size",
                       "Planned total sample size","Mean total sample size","mean control sample size conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=sim_results
  return(out)

}

#################################################
#Two_stage_adapt_rej_OC_simu_poisson(lambda_t,lambda_c,direction,s,c_bry_two_stg,beta,theta_cut,samsz,N_max,sim_num)
#Two_stage_adapt_rej_simu_poisson(lambda_t=c(0.9,0.8,0.7),lambda_c=1.1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000)
#Two_stage_adapt_rej_OC_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000)
#Two_stage_adapt_rej_OC_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000)
#Two_stage_adapt_rej_OC_simu_poisson(lambda_t=c(0.1,0.2,0.3),lambda_c=0.1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938065,2.784632,2.273643),
#                                   theta_cut=0.01,samsz=100,N_max=200,beta=0.1,sim_num=10000)


##################

#############################################
two_stg_adapt_power_OC_simu_poisson<-function(lambda_t,lambda_c,direction,s,bry_type,alpha_2,theta_cut,samsz,N_max,beta,sim_num){
  # direction==1 means that larger mu is better
  # direction==2 means that smaller mu is better

  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=length(lambda_t)
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))

  simu_out=Two_stage_adapt_rej_OC_simu_poisson(lambda_t,lambda_c,direction=direction,s=s,c_bry_two_stg=c_bry_two_stg,
                                               beta=beta,theta_cut=theta_cut,samsz=samsz,N_max=N_max,sim_num=sim_num)
  lambda_t_names=paste(paste(LABEL()$lambda_t,"[",c(1:length(lambda_t)),"]",sep=""))
  out_1= round(c(lambda_t,lambda_c,as.numeric(simu_out)),4)
  names(out_1)=c(lambda_t_names,LABEL()$lambda_c,names(simu_out))
  out=out_1
  return(out)
}
###############################


###########################
#two_stg_adapt_power_OC_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=200,beta=0.1,sim_num=100000)

#two_stg_adapt_power_OC_simu_poisson(lambda_t=c(1.1,1.1,1.1),lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=200,beta=0.1,sim_num=100000)

#two_stg_adapt_power_OC_simu_poisson(lambda_t=c(0.1,0.2,0.3),lambda_c=0.1,direction=1,s=c(1/3,2/3,1),
#                                    bry_type=1,alpha_2=0.05,theta_cut=0.01,samsz=100,N_max=200,
#                                    beta=0.1,sim_num=10000)
###############################################################


###################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Poisson distribution, Optimizing adaptive sequential design
#'
#' @param groups Matrix of assumed means for the test arms under study: one row per candidate test-arm scenario and one column per dose, so that the number of doses is the number of columns; each mean must be > 0
#' @param lambda_c Mean of the control arm; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s1 Vector of information fractions of the first interim analysis schedule; the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule; the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s1; each element must be > 0 and the elements must sum to alpha_2
#' @param alpha_2_alps_s2 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s2; each element must be > 0 and the elements must sum to alpha_2
#' @param theta_cut Lack-of-efficacy threshold: an effect size below which the test drug is not considered to provide clinically meaningful efficacy
#' @param samsz_1 Planned sample size of the control arm for the first sample size scenario
#' @param samsz_2 Planned sample size of the control arm for the second sample size scenario
#' @param N_max_1 Maximum sample size (control arm) allowed after sample size re-estimation for the first maximum sample size scenario
#' @param N_max_2 Maximum sample size (control arm) allowed after sample size re-estimation for the second maximum sample size scenario
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with O'Brien-Fleming boundaries}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with O'Brien-Fleming boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with alpha-spending boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with alpha-spending boundaries}
#' }
#' Each of these four lists has the components "Group sequential design" (data frame with the columns "Looks", "Information fractions" and the critical boundaries, plus "Alpha-spending" and "Alpha-spending boundaries" for the alpha-spending schedules) and the "Simulation results: Information fractions ..., Sample size ..., N_max ..." operating characteristics of the corresponding adaptive design, one data frame per combination of interim analysis schedule, planned sample size and maximum sample size, with one row per assumed mean scenario.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' Two_stage_asd_simu_OC_poisson(
#' groups = rbind(c(0.3,0.4,0.5),c(0.2,0.2,0.5),c(0.2,0.2,0.3),c(0.2,0.3,0.4)), lambda_c = 0.2,
#' direction = 1, s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05,
#' alpha_2_alps_s1 = c(0.001,0.001,0.048), alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01,
#' samsz_1 = 100, samsz_2 = 150, N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
Two_stage_asd_simu_OC_poisson<-function(groups,lambda_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
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

  s1_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=dist)[1]))
  s2_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=dist)[1]))

  s1_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s1,alpha=alpha_alsp_s1,comp=comp,dist=dist)[1]))
  s2_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s2,alpha=alpha_alsp_s2,comp=comp,dist=dist)[1]))

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
    simu_out_OF_s1_samsz_1_N_max_1=rbind(simu_out_OF_s1_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_1)=lambda_t_group_names



  ###############  s2 #######


  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_1=rbind(simu_out_OF_s2_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_1)=lambda_t_group_names


  ## samsz_2###

  ##### s1 ####


  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_1=rbind(simu_out_OF_s1_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_1)=lambda_t_group_names



  ###############  s2 #######


  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_1=rbind(simu_out_OF_s2_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
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
    simu_out_OF_s1_samsz_1_N_max_2=rbind(simu_out_OF_s1_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_2)=lambda_t_group_names



  ###############  s2 #######


  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_2=rbind(simu_out_OF_s2_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_2)=lambda_t_group_names


  ## samsz_2###

  ##### s1 ####


  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_2=rbind(simu_out_OF_s1_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_2)=lambda_t_group_names



  ###############  s2 #######


  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_2=rbind(simu_out_OF_s2_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
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
    simu_out_alsp_s1_samsz_1_N_max_1=rbind(simu_out_alsp_s1_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_1)=lambda_t_group_names



  ###############  s2 #######


  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_1=rbind(simu_out_alsp_s2_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                            bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_1)=lambda_t_group_names


  ## samsz_2###

  ##### s1 ####


  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_1=rbind(simu_out_alsp_s1_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_1)=lambda_t_group_names



  ###############  s2 #######


  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_1=rbind(simu_out_alsp_s2_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
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
    simu_out_alsp_s1_samsz_1_N_max_2=rbind(simu_out_alsp_s1_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_2)=lambda_t_group_names



  ###############  s2 #######


  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_2=rbind(simu_out_alsp_s2_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
                                                                                                            bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_2)=lambda_t_group_names


  ## samsz_2###

  ##### s1 ####


  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_2=rbind(simu_out_alsp_s1_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s1,
                                                                                                            bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                            beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_2)=lambda_t_group_names



  ###############  s2 #######


  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(lambda_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_2=rbind(simu_out_alsp_s2_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_poisson(lambda_t=lambda_t_vec[i,],lambda_c=lambda_c,direction=direction,s=s2,
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
#Two_stage_asd_simu_OC_normal(mu_t_vec_1,mu_t_vec_2,mu_c,sigma_t,sigma_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num,dist)
#Two_stage_asd_simu_OC_poisson(groups=rbind(c(0.3,0.4,0.5),c(0.2,0.2,0.5),c(0.2,0.2,0.3),c(0.2,0.3,0.4)),lambda_c=0.2,direction=1,
#                            s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                             theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1, sim_num=100)


#############################################

##########################################################################
### II.5 .Simulations  -Negative binomial distributions              #####
##########################################################################


#####################################################################################
### II.5.1 .Simulations with planned TSSSD -Negative binomial distributions     #####
##          Type I error, power,  simulations                                   #####
#### This will be included in the software                                      #####
#####################################################################################

#########################################
two_stg_rej_rate_simu_neg_binomial<-function(r_t,r_c,kappa,direction,s,c_bry_two_stg,samsz,sim_num){
  # direction==1 means that larger rate r is better
  # direction==2 means that smaller rate r is better

  #########################
  N=ceiling(samsz*s)
  comp=length(r_t)

  if (direction ==1){ #higher rate r is better
    mu=log(r_t/r_c)
  } else  if (direction ==2){ #lower rate r is better
    mu=-log(r_t/r_c)
  }

  # N is the vetoc of sample sizes at which interim analysis will be conducted

  ##########################
  N1=N[1]
  NF=samsz

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    rej=0; stg=0
    s_slct=0; control_samsz_simu=N[length(s)]
    tot_samsz_simu=N[1]*(comp-1)+2*control_samsz_simu
      #########################################################
      ### Step 1: Generate random samples fro phae 2 stage ####
      #########################################################
      ### geting random sample ################

      x_c=cbind(rep(1,samsz),stats::rnbinom(n=samsz, size=kappa, mu=r_c))

      x=matrix(nrow=samsz,ncol=comp)
      for (i in 1:comp){
        x[,i]=stats::rnbinom(n=samsz, size=kappa, mu=r_t[i])
      }

      x_t=cbind(rep(2,samsz),x)

      t_hat=matrix(nrow=comp,ncol=length(s))
      theta_hat=matrix(nrow=comp,ncol=length(s))
      wald=matrix(nrow=comp,ncol=length(s))
      for (i in 1:length(s)){

        for (j in 1:comp){

          x_c_temp=x_c[1:N[i],]
          x_t_temp=x_t[1:N[i],c(1,j+1)]

          x_temp=rbind(x_c_temp,x_t_temp)

          negglm_inter=tryCatch(glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"), error=function(e) NULL)
          sum_temp=if (is.null(negglm_inter)) NULL else summary(negglm_inter)

          ## a failed negative-binomial fit is reported as NA instead of
          ## aborting the batch inside theta.ml()
          se_temp=if (is.null(sum_temp)) NA_real_ else sum_temp$coefficients[2,2]

          r_c_temp=mean(x_c_temp[,2])
          r_t_temp=mean(x_t_temp[,2])

          if (direction ==1){ #higher rate r is better
          theta_hat[j,i]=log(r_t_temp/r_c_temp)
          } else if (direction ==2){ #lower rate r is better
            theta_hat[j,i]=-log(r_t_temp/r_c_temp)
          }

          t_hat[j,i]=1/se_temp^2
          wald[j,i]=theta_hat[j,i]*sqrt(t_hat[j,i])
        }
      }

      wald_max=NULL

      s_slct=match(max(wald[,1]),wald[,1])
      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      wald_two_stg=wald[s_slct,]

      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg

      # test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if (isTRUE(max(test)>0)){

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        last_vt=stg
        control_samsz_simu=N[last_vt]
        tot_samsz_simu=N[1]*(comp-1)+2*control_samsz_simu

        #stg identifies the interim analysis at which the null hypothesis is rejected

      }else{

        rej=0
        stg=0
        last_vt=length(s)
        control_samsz_simu=N[last_vt]
        tot_samsz_simu=N[1]*(comp-1)+2*control_samsz_simu
      }
      #######################################################
      ### Step 7: output and finish the simjulation loop ####
      #######################################################
      out=c(rej,stg,s_slct,control_samsz_simu,tot_samsz_simu)
  })#(end of q loop)

  ############################################
  ### Step 8: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  s_slct= as.numeric(simr[,3])
  control_samsz_simu= as.numeric(simr[,4])
  tot_samsz_simu= as.numeric(simr[,5])
  rej_rt=mean(rej)

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=sum(rej[stg==j])/(dim(simr)[1])
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  select_rt=rep(0,comp)
  for (j in 1:comp){
    select_rt[j]=length( s_slct[ s_slct==j])/(dim(simr)[1])
  }


  rej_rt_indiv=rep(0,comp)
  for (j in 1:comp){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/(dim(simr)[1])
  }

  con_rej_rt_indiv=rep(0,comp)
  for (j in 1:comp){
    if (isTRUE(select_rt[j]==0)){
      con_rej_rt_indiv[j]=NA
    }else {
      con_rej_rt_indiv[j]=rej_rt_indiv[j]/select_rt[j]
    }
  }



  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  #############################
  ### Step 9: Final output ####
  #############################
  results_by_arm=cbind(r_t, r_c,select_rt, rej_rt_indiv,con_rej_rt_indiv)
  colnames(results_by_arm)=c("r<sub>t1</sub>,...","r<sub>c</sub>","selection rate","rejection rate by each r_t","conditional rejection rate")

  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look")

  tyI_er=sum(rej[mu[s_slct]==0])/(dim(simr)[1])

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/(dim(simr)[1])
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(r_t)-1)+2*samsz
  sim_results=c(rej_rt,tyI_er, non_tyI_er_rej_rt,samsz,mean_control_samsz,tot_samsz,mean_tot_samsz)
  names(sim_results)=c("rejection rate","Type I error","Non-Type I error rejection rate",
                       "Planned control sample size","Mean control sample size","Planned total sample size","Mean total sample size")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(seq_design,results_by_arm,sim_results)
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)
}
#######################################
#two_stg_rej_rate_simu_neg_binomial(r_t,r_c,kappa,nu_t,rand_ratio,s,c_bry_two_stg,samsz,sim_num)
#two_stg_rej_rate_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),samsz=100,sim_num=1000)
#two_stg_rej_rate_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),samsz=100,sim_num=100000)


###############Estimates#####################


##############################################

#' Phase 2/3 seamless combination Group sequential design: Simulations, Negative binomial distribution, Type I error and power
#'
#' @param r_t Mean event rates of the test arms; each element must be > 0
#' @param r_c Mean event rate of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param samsz Planned sample size of the control arm used in the simulation
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries" and "rejection by look"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "r<sub>t1</sub>,...", "r<sub>c</sub>", "selection rate", "rejection rate by each r_t" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the overall "rejection rate", the "Type I error", the "Non-Type I error rejection rate", the "Planned control sample size", the "Mean control sample size", the "Planned total sample size" and the "Mean total sample size"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_power_simu_neg_binomial(r_t = c(1.1,1.1,1.1), r_c = 1.1, kappa = 1, direction = 1,
#' s = c(1/3,1), bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
two_stg_power_simu_neg_binomial<-function(r_t,r_c,kappa,direction,s,bry_type,alpha_2,samsz,sim_num){

  # direction==1 means that larger rate r is better
  # direction==2 means that smaller rate r is better

  alpha=alpha_2/2
  comp=length(r_t)
  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s,alpha,comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s,alpha,comp,dist=1)
  }
  c_bry_two_stg=as.numeric(unlist(design_out[1]))
  simu_out=two_stg_rej_rate_simu_neg_binomial(r_t,r_c,kappa,direction,s,c_bry_two_stg,samsz,sim_num)
  out=simu_out
  return(out)
}

##################################
#two_stg_power_simu_neg_binomial(r_t,r_c,kappa,nu_t,rand_ratio,s,bry_type,alpha_2,samsz,sim_num)
#two_stg_power_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=100000)
#two_stg_power_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,1),bry_type=2,alpha_2=c(0.001,0.049),samsz=100,sim_num=100000)

#two_stg_power_simu_neg_binomial(r_t=c(1.1,1.2,1.3),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)
#two_stg_power_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),samsz=100,sim_num=100000)
#two_stg_power_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),samsz=100,sim_num=100000)
#two_stg_power_simu_neg_binomial(r_t=c(1.3,1.4,1.5),r_c=0.8,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),samsz=129,sim_num=100000)


#######################################


####################################################################################
### II.4 .Simulations with adaptive TSSSD -Negative binomial distributions     #####
##          Type I error, power, simulations                                   #####
####################################################################################

#############################################
Two_stage_adapt_rej_simu_neg_binomial<-function(r_t,r_c,kappa,direction,s,c_bry_two_stg,beta,theta_cut,samsz,N_max,sim_num){

  # direction==1 means that larger rate r is better
  # direction==2 means that smaller rate r is better

  #########################
  N=ceiling(samsz*s)
  comp=length(r_t)


  if (direction ==1){ #higher rate r is better
    mu=log(r_t/r_c)
  } else  if (direction ==2){ #lower rate r is better
    mu=-log(r_t/r_c)
  }

  ##########################
  N1=N[1]
  NF=samsz
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############ the simulation loops ##########
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=N[length(s)]
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    control_sampsz_simu=sampsz_simu; tot_sampsz_simu=(comp+1)*sampsz_simu
    futi=0
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################

      x_c=cbind(rep(1,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_c))

      x=matrix(nrow=N_max,ncol=comp)
      for (i in 1:comp){
        x[,i]=stats::rnbinom(n=N_max, size=kappa, mu=r_t[i])
      }

      x_t=cbind(rep(2,N_max),x)

      t_hat=matrix(nrow=comp,ncol=length(s))
      theta_hat=matrix(nrow=comp,ncol=length(s))
      wald=matrix(nrow=comp,ncol=length(s))
      v_hat=matrix(nrow=comp,ncol=length(s))
      theta_hat_sd=matrix(nrow=comp,ncol=length(s))
      for (i in 1:length(s)){

        for (j in 1:comp){

          x_c_temp=x_c[1:N[i],]
          x_t_temp=x_t[1:N[i],c(1,j+1)]

          x_temp=rbind(x_c_temp,x_t_temp)

          negglm_temp=tryCatch(glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"), error=function(e) NULL)
          sum_temp=if (is.null(negglm_temp)) NULL else summary(negglm_temp)

          ## a failed negative-binomial fit is reported as NA instead of
          ## aborting the batch inside theta.ml()
          se_temp=if (is.null(sum_temp)) NA_real_ else sum_temp$coefficients[2,2]
          r_c_temp=mean(x_c_temp[,2])
          r_t_temp=mean(x_t_temp[,2])

          if (direction ==1){ #higher rate r is better
            theta_hat[j,i]=log(r_t_temp/r_c_temp)
          } else if (direction ==2){ #lower rate r is better
            theta_hat[j,i]=-log(r_t_temp/r_c_temp)
          }

          t_hat[j,i]=1/se_temp^2
          wald[j,i]=theta_hat[j,i]*sqrt(t_hat[j,i])
        }
      }

      s_slct=match(max(wald[,1]),wald[,1])
      theta_hat_two_stg=theta_hat[s_slct,]

      t_two_stg=t_hat[s_slct,]

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      wald_two_stg=wald[s_slct,]

      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg

      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################


      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test_1=test[-length(s)]  ## test statistics excluding the final analysis

      ##########################################
      theta_hat_two_stg=theta_hat[s_slct,]
      inter_vt=length(s)-1
      ############################################
      if (isTRUE(theta_hat_two_stg[inter_vt]<theta_cut)) {futi=1} else {futi=0}
      ##################################

      if (isTRUE(max(test)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
      }
      #################

      if  (isTRUE(max(test_1)>=0)){ #null hypothesis rejected, no more random sampling ot theta_hat at length(s_p3)-1 is almost negative

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        ###########################
        sampsz_simu=N[last_vt]
        if (last_vt==1){
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*sampsz_simu
        }else {
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
        }

      }else if (isTRUE(max(test_1)<0)){

        if (futi==1){

          sampsz_simu=N[length(s)]
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
          if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat_two_stg[inter_vt]
         # theta_hat_inter_sd=theta_hat_sd_two_stg[inter_vt]
          score_inter=theta_hat_inter*t_two_stg[inter_vt]
          T_two_stg=t_two_stg[length(s)]
          c_bry=c_bry_two_stg

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t_two_stg[length(s)]-t_two_stg[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t_two_stg[length(s)])-score_inter)/sqrt(t_two_stg[length(s)]-t_two_stg[inter_vt]))

          if (isTRUE(cond_power>=1-beta)){

            sampsz_simu=N[length(s)]
            control_sampsz_simu=sampsz_simu
            tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (isTRUE(cond_power<1-beta)){

            ## proceed to sample size re-estimation #####
              ###############
            T_new_0=theta_hat_inter^(-2)*((c_bry_two_stg[length(s)]*sqrt(T_two_stg)-score_inter)/sqrt(T_two_stg-t_two_stg[inter_vt])+stats::qnorm(1-beta))^2+t_two_stg[inter_vt]
            new_sampsz=ceiling((T_new_0/t_hat[s_slct,length(s)])*NF )

            ###############################################
            ###############################################
            if (isTRUE(new_sampsz>N_max)){ ## new sample size is too large, the sample size will not change. Trial ends at the planned sample size
              if (isTRUE(test[length(test)]>0)){ ### the NULL hypothesis is rejected at the final anbalysis, with the planned sample size
                rej=1
                stg=length(s)
              }else{### the NULL hypothesis is NOT rejected at the final anbalysis, with the planned sample size
                rej=0
                stg=0
              }

              ### combined phase 2/3 data ####
              last_vt=length(s)
              sampsz_simu=N[last_vt]
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

            } else if (isTRUE(new_sampsz<= N_max)) {## new sample size does not exceed the pre-selected maximum. Will change the sample  size to this new sample size

              sampsz_simu=new_sampsz
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              x0_two_stg_new=x_c[1:new_sampsz,]
              x_t_two_stg_new=x_t[1:new_sampsz,c(1,s_slct+1)]

              x_two_stg_new=rbind(x0_two_stg_new,x_t_two_stg_new)

              negglm_two_stg_new=tryCatch(glm.nb(as.numeric(x_two_stg_new[,2])~as.factor(x_two_stg_new[,1]),link="log"), error=function(e) NULL)
              sum_two_stg_new=if (is.null(negglm_two_stg_new)) NULL else summary(negglm_two_stg_new)

              se_sum_two_stg_new=if (is.null(sum_two_stg_new)) NA_real_ else sum_two_stg_new$coefficients[2,2]
              r_c_two_stg_new=mean(x0_two_stg_new[,2])
              r_t_two_stg_new=mean(x_t_two_stg_new[,2])

              if (direction ==1){ #higher rate r is better
                theta_hat_two_stg_new=log(r_t_two_stg_new/r_c_two_stg_new)
              } else if (direction ==2){ #lower rate r is better
                theta_hat_two_stg_new=-log(r_t_two_stg_new/r_c_two_stg_new)
              }

              t_hat_two_stg_new=1/se_sum_two_stg_new^2
              wald_two_stg_new=theta_hat_two_stg_new*sqrt(t_hat_two_stg_new)
              ####################
              T_new=t_hat_two_stg_new
              T_final_1=(NF/new_sampsz)*t_hat_two_stg_new
              t_inter_1=(N[inter_vt]/new_sampsz)*t_hat_two_stg_new
              e_bry_new_two_stg=
                (sqrt(T_new-t_inter_1)/sqrt(T_final_1-t_inter_1))*(c_bry_two_stg[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter
              c_bry_new_two_stg=e_bry_new_two_stg/sqrt(T_new)

              ##########################################

              test_add= wald_two_stg_new-c_bry_new_two_stg  # compare the final Wald statistics with the new adjusted critical boundary

              if (isTRUE(test_add>0)){ # the new critical boundary is crossed, the NULL is rejected

                rej=1
                stg=length(s)

              }else{  # the new critical boundary is NOT crossed, the NULL is rejected

                rej=0
                stg=0

              }

              ###################### Parameter estimate with sample size change ##################

            }
          }
        }
      }
      out=c(rej_nossr,rej,stg,s_slct,control_sampsz_simu,tot_sampsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  control_samsz_simu = as.numeric(simr[,5])
  tot_samsz_simu = as.numeric(simr[,6])
  futi=as.numeric(simr[,7])
  #############termination timepoint####
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)
  rej_time=rep(0,(length(s)))
  for (j in 1:(length(s))){

    count_j=c(1:sim_num)[stg==j]
    rej_time[j]=sum(rej[count_j])/(dim(simr)[1])

  }
  select_rt=rep(0,comp)
  for (j in 1:comp){
    select_rt[j]=length( s_slct[ s_slct==j])/(dim(simr)[1])
  }


  rej_rt_indiv=rep(0,comp)
  for (j in 1:comp){
    rej_rt_indiv[j]=sum( rej[ s_slct==j])/(dim(simr)[1])
  }

  con_rej_rt_indiv=rep(0,comp)
  for (j in 1:comp){
    if (isTRUE(select_rt[j]==0)){
      con_rej_rt_indiv[j]=NA
    }else {
      con_rej_rt_indiv[j]=rej_rt_indiv[j]/select_rt[j]
    }
  }

  samsz_icr=NULL
  for (i in 1:length(control_samsz_simu)){
    if (control_samsz_simu[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  samsz_icr_rt_byLook=c(rep(0,(length(s)-1)),mean(samsz_icr))
  futi_rt=mean(futi)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= control_samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  #############################
  ### Step 9: Final output ####
  #############################
  results_by_arm=cbind(r_t, r_c,select_rt, rej_rt_indiv,con_rej_rt_indiv)
  colnames(results_by_arm)=c("r<sub>t1</sub>,...","r<sub>c</sub>","selection rate","rejection rate by each r_t","conditional rejection rate")


  seq_design=cbind(c(1:length(s)),s,c_bry_two_stg,rej_time,samsz_icr_rt_byLook)
  colnames(seq_design)=c("Look","Information fraction","Critical boundaries","rejection by look","Sample size increase")

  tyI_er=sum(rej[mu[s_slct]==0])/(dim(simr)[1])

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/(dim(simr)[1])
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(r_t)-1)+2*samsz
  sim_results=c(rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er, non_tyI_er_rej_rt,futi_rt,samsz,
                N_max,mean_control_samsz,tot_samsz,mean_tot_samsz, mean_samsz_simu_control_samsz_icr)
  names(sim_results)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned control sample size","Maximum control sample size","Mean control sample size",
                       "Planned total sample size","Mean total sample size","mean control sample size conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=list(seq_design,results_by_arm,sim_results)
  names(out)=c("Sequential design","Results by arm","Simulation results")
  return(out)

}


#################################################
#Two_stage_adapt_rej_simu_neg_binomial(r_t,r_c,kappa,nu_t,rand_ratio,direction,s,c_bry_two_stg,beta,theta_cut,samsz,N_max,sim_num)
#Two_stage_adapt_rej_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,nu_t=1,rand_ratio=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000)
#Two_stage_adapt_rej_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=10000)
#Two_stage_adapt_rej_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000)
#Two_stage_adapt_rej_simu_neg_binomial(r_t=c(1.2,1.2,1.2),r_c=1.1,kappa=1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=100000)


###############################################################


#################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Negative binomial distribution, Operating characteristics
#'
#' @param r_t Mean event rates of the test arms; each element must be > 0
#' @param r_c Mean event rate of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock; any other value leaves the boundaries undefined
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param theta_cut Lack-of-efficacy threshold: an effect size below which the test drug is not considered to provide clinically meaningful efficacy
#' @param samsz Planned sample size of the control arm used in the simulation
#' @param N_max Maximum sample size (control arm) allowed after sample size re-estimation
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design}{Matrix of the planned adaptive group sequential design, with one row per analysis and the columns "Look", "Information fraction", "Critical boundaries", "rejection by look" and "Sample size increase"}
#'   \item{Results by arm}{Matrix with one row per dose and the columns "r<sub>t1</sub>,...", "r<sub>c</sub>", "selection rate", "rejection rate by each r_t" and "conditional rejection rate"; the selection rate is the probability that the dose is selected to proceed to phase 3, the rejection rate is the probability that the dose is selected and the null hypothesis is rejected, and the conditional rejection rate is the rejection rate divided by the selection rate}
#'   \item{Simulation results}{Named vector of the "Rejection rate without SSR", the "Rejection rate with SSR", the "rate of sample size increase", the "Rejection rate without SSR/(meeting sample size incerase criteria)", the "Rejection rate with SSR/(sample size incerase)", the "Type I error", the "Non-Type I error rejection rate", the "Futility rate", the "Planned control sample size", the "Maximum control sample size", the "Mean control sample size", the "Planned total sample size", the "Mean total sample size" and the "mean control sample size conditional on increase"}
#' }
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stg_adapt_power_simu_neg_binomial(
#' r_t = c(1.1,1.1,1.1), r_c = 1.1, kappa = 1, direction = 1,
#' s = c(1/3,1), bry_type = 1, alpha_2 = 0.05, theta_cut = 0.05, samsz = 100, N_max = 1000,
#' beta = 0.1, sim_num = 20)
#'
#' @export
two_stg_adapt_power_simu_neg_binomial<-function(r_t,r_c,kappa,direction,s,bry_type,alpha_2,theta_cut,samsz,N_max,beta,sim_num){

  # direction==1 means that larger rate r is better
  # direction==2 means that smaller rate r is better

  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2
  comp=as.numeric(length(r_t))

  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))

  simu_out=Two_stage_adapt_rej_simu_neg_binomial(r_t=r_t,r_c=r_c,kappa=kappa,direction=direction,s=s,c_bry_two_stg=c_bry_two_stg,beta=beta,theta_cut=theta_cut,samsz=samsz,N_max=N_max,sim_num=sim_num)

  out=simu_out
  return(out)
}

##################################
#two_stg_adapt_power_simu_neg_binomial(r_t,r_c,kappa,nu_t,rand_ratio,direction,s,bry_type,alpha_2,theta_cut,N_max,samsz,beta,sim_num)

#two_stg_adapt_power_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=1000,beta=0.1,sim_num=1000)
#two_stg_adapt_power_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=1000,beta=0.1,sim_num=100000)
#two_stg_adapt_power_simu_neg_binomial(r_t=c(1.3,1.4,1.5),r_c=0.8,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=129,N_max=1000,beta=0.1,sim_num=100000)
#two_stg_adapt_power_simu_neg_binomial(r_t=c(1.3,1.4,1.5),r_c=0.8,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=300,beta=0.1,sim_num=100000)
#two_stg_adapt_power_simu_neg_binomial(r_t=c(1.2,1.4,1.5),r_c=1.2,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=500,beta=0.1,sim_num=10000)

##############################################################
Two_stage_adapt_rej_OC_simu_neg_binomial<-function(r_t,r_c,kappa,direction,s,c_bry_two_stg,beta,theta_cut,samsz,N_max,sim_num){

  # direction==1 means that larger rate r is better
  # direction==2 means that smaller rate r is better

  #########################
  N=ceiling(samsz*s)
  comp=length(r_t)


  if (direction ==1){ #higher rate r is better
    mu=log(r_t/r_c)
  } else  if (direction ==2){ #lower rate r is better
    mu=-log(r_t/r_c)
  }

  ##########################
  N1=N[1]
  NF=samsz
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############ the simulation loops ##########
  simr <- simuloop(sim_num, "2ASD", function() {
  ## ---- non-evaluable default (added 2026-10-31) ----------------------
  ## A replicate whose statistics cannot be evaluated is recorded as
  ## 'no rejection, no sample size increase, finish at the planned size'.
  ## Without these defaults an NA decision left the output variables
  ## undefined and aborted the whole batch.
    sampsz_simu=N[length(s)]
    rej_nossr=0; rej=0
    stg=0; s_slct=0
    control_sampsz_simu=sampsz_simu; tot_sampsz_simu=(comp+1)*sampsz_simu
    futi=0
      #########################################################
      ### Step 1: Generate random samples fro phase 2 stage ####
      #########################################################

      x_c=cbind(rep(1,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_c))

      x=matrix(nrow=N_max,ncol=comp)
      for (i in 1:comp){
        x[,i]=stats::rnbinom(n=N_max, size=kappa, mu=r_t[i])
      }

      x_t=cbind(rep(2,N_max),x)

      t_hat=matrix(nrow=comp,ncol=length(s))
      theta_hat=matrix(nrow=comp,ncol=length(s))
      wald=matrix(nrow=comp,ncol=length(s))
      v_hat=matrix(nrow=comp,ncol=length(s))
      theta_hat_sd=matrix(nrow=comp,ncol=length(s))
      for (i in 1:length(s)){

        for (j in 1:comp){

          x_c_temp=x_c[1:N[i],]
          x_t_temp=x_t[1:N[i],c(1,j+1)]

          x_temp=rbind(x_c_temp,x_t_temp)

          negglm_temp=tryCatch(glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"), error=function(e) NULL)
          sum_temp=if (is.null(negglm_temp)) NULL else summary(negglm_temp)

          ## a failed negative-binomial fit is reported as NA instead of
          ## aborting the batch inside theta.ml()
          se_temp=if (is.null(sum_temp)) NA_real_ else sum_temp$coefficients[2,2]
          r_c_temp=mean(x_c_temp[,2])
          r_t_temp=mean(x_t_temp[,2])

          if (direction ==1){ #higher rate r is better
            theta_hat[j,i]=log(r_t_temp/r_c_temp)
          } else if (direction ==2){ #lower rate r is better
            theta_hat[j,i]=-log(r_t_temp/r_c_temp)
          }

          t_hat[j,i]=1/se_temp^2
          wald[j,i]=theta_hat[j,i]*sqrt(t_hat[j,i])
        }
      }

      s_slct=match(max(wald[,1]),wald[,1])
      theta_hat_two_stg=theta_hat[s_slct,]

      t_two_stg=t_hat[s_slct,]

      ################################################################################
      ### Step 2: calculate the score statistics for selection  at end of phase 2 ####
      ################################################################################

      #####################################################################
      ### Step 5: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      wald_two_stg=wald[s_slct,]

      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg

      #############################################################################################
      ### Step 3: the dose with the largest S_p2 score is selected to porceed to phase 3 stage ####
      #############################################################################################


      ################# test srtatistics ################
      ###################################################################################################
      ### Step 6: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################
      test=wald_two_stg-c_bry_two_stg  ## comparing the Wald statistics with the critical bouundary  for combined phase 2/3  test

      ################################################

      test_1=test[-length(s)]  ## test statistics excluding the final analysis

      ##########################################
      theta_hat_two_stg=theta_hat[s_slct,]
      inter_vt=length(s)-1
      ############################################
      if (isTRUE(theta_hat_two_stg[inter_vt]<theta_cut)) {futi=1} else {futi=0}
      ##################################

      if (isTRUE(max(test)>=0)){
        rej_nossr=1
      }else if (isTRUE(max(test)<0)){
        rej_nossr=0
      }
      #################

      if  (isTRUE(max(test_1)>=0)){ #null hypothesis rejected, no more random sampling ot theta_hat at length(s_p3)-1 is almost negative

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        ###########################
        sampsz_simu=N[last_vt]
        if (last_vt==1){
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp+1)*sampsz_simu
        }else {
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
        }

      }else if (isTRUE(max(test_1)<0)){

        if (futi==1){

          sampsz_simu=N[length(s)]
          control_sampsz_simu=sampsz_simu
          tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
          if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat_two_stg[inter_vt]
          # theta_hat_inter_sd=theta_hat_sd_two_stg[inter_vt]
          score_inter=theta_hat_inter*t_two_stg[inter_vt]
          T_two_stg=t_two_stg[length(s)]
          c_bry=c_bry_two_stg

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t_two_stg[length(s)]-t_two_stg[inter_vt])-(c_bry_two_stg[length(s)]*sqrt(t_two_stg[length(s)])-score_inter)/sqrt(t_two_stg[length(s)]-t_two_stg[inter_vt]))

          if (isTRUE(cond_power>=1-beta)){

            sampsz_simu=N[length(s)]
            control_sampsz_simu=sampsz_simu
            tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu
            if (isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (isTRUE(cond_power<1-beta)){

            ## proceed to sample size re-estimation #####
            ###############
            T_new_0=theta_hat_inter^(-2)*((c_bry_two_stg[length(s)]*sqrt(T_two_stg)-score_inter)/sqrt(T_two_stg-t_two_stg[inter_vt])+stats::qnorm(1-beta))^2+t_two_stg[inter_vt]
            new_sampsz=ceiling((T_new_0/t_hat[s_slct,length(s)])*NF )

            ###############################################
            ###############################################
            if (isTRUE(new_sampsz>N_max)){ ## new sample size is too large, the sample size will not change. Trial ends at the planned sample size
              if (isTRUE(test[length(test)]>0)){ ### the NULL hypothesis is rejected at the final anbalysis, with the planned sample size
                rej=1
                stg=length(s)
              }else{### the NULL hypothesis is NOT rejected at the final anbalysis, with the planned sample size
                rej=0
                stg=0
              }

              ### combined phase 2/3 data ####
              last_vt=length(s)
              sampsz_simu=N[last_vt]
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

            } else if (isTRUE(new_sampsz<= N_max)) {## new sample size does not exceed the pre-selected maximum. Will change the sample  size to this new sample size

              sampsz_simu=new_sampsz
              control_sampsz_simu=sampsz_simu
              tot_sampsz_simu=(comp-1)*N[1]+2*sampsz_simu

              x0_two_stg_new=x_c[1:new_sampsz,]
              x_t_two_stg_new=x_t[1:new_sampsz,c(1,s_slct+1)]

              x_two_stg_new=rbind(x0_two_stg_new,x_t_two_stg_new)

              negglm_two_stg_new=tryCatch(glm.nb(as.numeric(x_two_stg_new[,2])~as.factor(x_two_stg_new[,1]),link="log"), error=function(e) NULL)
              sum_two_stg_new=if (is.null(negglm_two_stg_new)) NULL else summary(negglm_two_stg_new)

              se_sum_two_stg_new=if (is.null(sum_two_stg_new)) NA_real_ else sum_two_stg_new$coefficients[2,2]
              r_c_two_stg_new=mean(x0_two_stg_new[,2])
              r_t_two_stg_new=mean(x_t_two_stg_new[,2])

              if (direction ==1){ #higher rate r is better
                theta_hat_two_stg_new=log(r_t_two_stg_new/r_c_two_stg_new)
              } else if (direction ==2){ #lower rate r is better
                theta_hat_two_stg_new=-log(r_t_two_stg_new/r_c_two_stg_new)
              }

              t_hat_two_stg_new=1/se_sum_two_stg_new^2
              wald_two_stg_new=theta_hat_two_stg_new*sqrt(t_hat_two_stg_new)
              ####################
              T_new=t_hat_two_stg_new
              T_final_1=(NF/new_sampsz)*t_hat_two_stg_new
              t_inter_1=(N[inter_vt]/new_sampsz)*t_hat_two_stg_new
              e_bry_new_two_stg=
                (sqrt(T_new-t_inter_1)/sqrt(T_final_1-t_inter_1))*(c_bry_two_stg[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter
              c_bry_new_two_stg=e_bry_new_two_stg/sqrt(T_new)

              ##########################################

              test_add= wald_two_stg_new-c_bry_new_two_stg  # compare the final Wald statistics with the new adjusted critical boundary

              if (isTRUE(test_add>0)){ # the new critical boundary is crossed, the NULL is rejected

                rej=1
                stg=length(s)

              }else{  # the new critical boundary is NOT crossed, the NULL is rejected

                rej=0
                stg=0

              }

              ###################### Parameter estimate with sample size change ##################

            }
          }
        }
      }
      out=c(rej_nossr,rej,stg,s_slct,control_sampsz_simu,tot_sampsz_simu,futi)
  })#(end of q loop)
  #parallel computing ends. All resuts are combined and send to simr
  #################
  rej_nossr= as.numeric(simr[,1])
  rej = as.numeric(simr[,2])
  stg = as.numeric(simr[,3])
  s_slct = as.numeric(simr[,4])
  control_samsz_simu = as.numeric(simr[,5])
  tot_samsz_simu = as.numeric(simr[,6])
  futi=as.numeric(simr[,7])
  #############termination timepoint####
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  samsz_icr=NULL
  for (i in 1:length(control_samsz_simu)){
    if (control_samsz_simu[i]>samsz){
      samsz_icr[i]=1
    }else{
      samsz_icr[i]=0
    }
  }
  samsz_icr_rt=mean(samsz_icr)
  futi_rt=mean(futi)
  #### no ssr vs ssr when samsz_icr==1##
  rej_nossr_samsz_icr= rej_nossr[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  samsz_simu_control_samsz_icr= control_samsz_simu[samsz_icr==1]
  rej_samsz_icr= rej[samsz_icr==1]
  rej_nossr_rt_samsz_icr=mean( rej_nossr_samsz_icr)
  rej_rt_samsz_icr=mean( rej_samsz_icr)
  mean_samsz_simu_control_samsz_icr=mean(samsz_simu_control_samsz_icr)

  #############################
  ### Step 9: Final output ####
  #############################

  tyI_er=sum(rej[mu[s_slct]==0])/(dim(simr)[1])

  non_tyI_er_rej_rt=sum(rej[mu[s_slct]>0])/(dim(simr)[1])
  mean_control_samsz=mean(control_samsz_simu)
  mean_tot_samsz=mean(tot_samsz_simu)
  tot_samsz=N1*(length(r_t)-1)+2*samsz
  sim_results=c(rej_nossr_rt,rej_rt,samsz_icr_rt,rej_nossr_rt_samsz_icr,rej_rt_samsz_icr,tyI_er, non_tyI_er_rej_rt,futi_rt,samsz,
                N_max,mean_control_samsz,tot_samsz,mean_tot_samsz, mean_samsz_simu_control_samsz_icr)
  names(sim_results)=c("Rejection rate without SSR","Rejection rate with SSR","rate of sample size increase","Rejection rate without SSR/(meeting sample size incerase criteria)",
                       "Rejection rate with SSR/(sample size incerase)",
                       "Type I error","Non-Type I error rejection rate","Futility rate",
                       "Planned control sample size","Maximum control sample size","Mean control sample size",
                       "Planned total sample size","Mean total sample size","mean control sample size conditional on increase")
  #############################
  ### Step 9: Final output ####
  #############################

  out=sim_results
  return(out)

}

#################################################
#Two_stage_adapt_rej_simu_neg_binomial(r_t,r_c,kappa,nu_t,rand_ratio,direction,s,c_bry_two_stg,beta,theta_cut,samsz,N_max,sim_num)
#Two_stage_adapt_rej_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,nu_t=1,rand_ratio=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000)
#Two_stage_adapt_rej_OC_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry_two_stg=c(3.938244,2.784759,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=10000)
#Two_stage_adapt_rej_OC_simu_neg_binomial(r_t=c(1.1,1.1,1.1),r_c=1.1,kappa=1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000)
#Two_stage_adapt_rej_OC_simu_neg_binomial(r_t=c(1.2,1.2,1.2),r_c=1.1,kappa=1,direction=1,s=c(1/3,1),c_bry_two_stg=c(3.938244,2.273746),beta=0.1,theta_cut=0.05,samsz=100,N_max=200,sim_num=1000)

#################
two_stg_adapt_power_OC_simu_neg_binomial<-function(r_t,r_c,kappa,direction,s,bry_type,alpha_2,theta_cut,samsz,N_max,beta,sim_num){

  # direction==1 means that larger rate r is better
  # direction==2 means that smaller rate r is better


  alpha=alpha_2/2

  comp=length(r_t)
 # tot_samsz=ceiling((samsz*s[1])*(comp+1)+2*(samsz*(s[length(s)]-s[1])))

  if (bry_type==1) {
    design_out=two_stage_OF_bundry(s=s,alpha=alpha,comp=comp,dist=1)

                  #two_stage_power_OF_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2)
  } else if (bry_type==2) {
    design_out=two_stage_alsp_bundry(s=s,alpha=alpha,comp=comp,dist=1)
  }

  c_bry_two_stg=as.numeric(unlist(design_out[1]))


  simu_out=Two_stage_adapt_rej_OC_simu_neg_binomial(r_t=r_t,r_c=r_c,kappa=kappa,direction=direction,s=s,c_bry_two_stg=c_bry_two_stg,beta=beta,theta_cut=theta_cut,samsz=samsz,N_max=N_max,sim_num=sim_num)

  r_t_names=paste(paste("r<sub>t</sub>","[",c(1:length(r_t)),"]",sep=""))
  out_1= round(c(r_t,r_c,as.numeric(simu_out)),4)
  names(out_1)=c(r_t_names,"r<sub>c</sub>", names(simu_out))
  # out=data.frame(out_1)
  out=out_1
  return(out)

}

##################################
#two_stg_adapt_power_OC_simu_neg_binomial(r_t,r_c,kappa,direction,s,bry_type,alpha_2,theta_cut,N_max,samsz,beta,sim_num)

#two_stg_adapt_power_OC_simu_neg_binomial(r_t=c(1.1,1.3,1.4),r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=1000,beta=0.1,sim_num=10000)
#two_stg_adapt_power_OC_simu_neg_binomial(r_t=c(0.6,0.6,0.6),r_c=0.6,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=1000,beta=0.1,sim_num=100000)
#two_stg_adapt_power_OC_simu_neg_binomial(r_t=c(1.3,1.4,1.5),r_c=0.8,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=129,N_max=1000,beta=0.1,sim_num=1000)
#two_stg_adapt_power_OC_simu_neg_binomial(r_t=c(1.3,1.4,1.5),r_c=0.8,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,theta_cut=0.05,samsz=100,N_max=300,beta=0.1,sim_num=1000)
#two_stg_adapt_power_OC_simu_neg_binomial(r_t=c(1.3,1.4,1.5),r_c=0.8,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=c(0.001,0.001,0.048),theta_cut=0.05,samsz=100,N_max=300,beta=0.1,sim_num=1000)

###############################################################


###################
#' Phase 2/3 seamless combination Group sequential design: Simulations, Negative binomial distribution, Optimizing adaptive sequential design
#'
#' @param groups Matrix of assumed mean event rates for the test arms under study: one row per candidate test-arm scenario and one column per dose, so that the number of doses is the number of columns; each mean event rate must be > 0
#' @param r_c Mean event rate of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s1 Vector of information fractions of the first interim analysis schedule; the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule; the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; 0 < alpha_2 < 1
#' @param alpha_2_alps_s1 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s1; each element must be > 0 and the elements must sum to alpha_2
#' @param alpha_2_alps_s2 Vector of two-sided alpha to be spent at each analysis of the alpha-spending schedule s2; each element must be > 0 and the elements must sum to alpha_2
#' @param theta_cut Lack-of-efficacy threshold: an effect size below which the test drug is not considered to provide clinically meaningful efficacy
#' @param samsz_1 Planned sample size of the control arm for the first sample size scenario
#' @param samsz_2 Planned sample size of the control arm for the second sample size scenario
#' @param N_max_1 Maximum sample size (control arm) allowed after sample size re-estimation for the first maximum sample size scenario
#' @param N_max_2 Maximum sample size (control arm) allowed after sample size re-estimation for the second maximum sample size scenario
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample size re-estimation; 0 < beta < 1
#' @param sim_num Number of simulation replicates
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with O'Brien-Fleming boundaries}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with O'Brien-Fleming boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{List holding the group sequential design and the operating characteristics simulated under the schedule s1 with alpha-spending boundaries}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{List holding the group sequential design and the operating characteristics simulated under the schedule s2 with alpha-spending boundaries}
#' }
#' Each of these four lists has the components "Group sequential design" (data frame with the columns "Looks", "Information fractions" and the critical boundaries, plus "Alpha-spending" and "Alpha-spending boundaries" for the alpha-spending schedules) and the "Simulation results: Information fractions ..., Sample size ..., N_max ..." operating characteristics of the corresponding adaptive design, one data frame per combination of interim analysis schedule, planned sample size and maximum sample size, with one row per assumed mean event-rate scenario.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' Two_stage_asd_simu_OC_neg_binomial(
#' groups = rbind(c(0.3,0.4,0.5),c(0.2,0.2,0.5),c(0.2,0.2,0.3),c(0.2,0.3,0.4)), r_c = 0.2,
#' kappa = 1, direction = 1, s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05,
#' alpha_2_alps_s1 = c(0.001,0.001,0.048), alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01,
#' samsz_1 = 100, samsz_2 = 150, N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 2)
#' }
#'
#' @export
Two_stage_asd_simu_OC_neg_binomial<-function(groups,r_c,kappa,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
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

  s1_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s1,alpha=alpha,comp=comp,dist=dist)[1]))
  s2_OF_cbry=as.numeric(unlist(two_stage_OF_bundry(s=s2,alpha=alpha,comp=comp,dist=dist)[1]))

  s1_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s1,alpha=alpha_alsp_s1,comp=comp,dist=dist)[1]))
  s2_alsp_cbry=as.numeric(unlist(two_stage_alsp_bundry(s=s2,alpha=alpha_alsp_s2,comp=comp,dist=dist)[1]))

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
    simu_out_OF_s1_samsz_1_N_max_1=rbind(simu_out_OF_s1_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                  bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                  beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_1)=r_t_group_names



  ###############  s2 #######

  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_1=rbind(simu_out_OF_s2_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                 bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_1)=r_t_group_names


  ## samsz_2###


  ##### s1 ####

  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_1=rbind(simu_out_OF_s1_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                 bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_1)=r_t_group_names



  ###############  s2 #######

  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_1=rbind(simu_out_OF_s2_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
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
    simu_out_OF_s1_samsz_1_N_max_2=rbind(simu_out_OF_s1_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                 bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_1_N_max_2)=r_t_group_names



  ###############  s2 #######

  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_OF_s2_samsz_1_N_max_2=rbind(simu_out_OF_s2_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                 bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s2_samsz_1_N_max_2)=r_t_group_names


  ## samsz_2###


  ##### s1 ####

  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_OF_s1_samsz_2_N_max_2=rbind(simu_out_OF_s1_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                 bry_type=1,alpha_2=alpha_2,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_OF_s1_samsz_2_N_max_2)=r_t_group_names



  ###############  s2 #######

  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_OF_s2_samsz_2_N_max_2=rbind(simu_out_OF_s2_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
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
    simu_out_alsp_s1_samsz_1_N_max_1=rbind(simu_out_alsp_s1_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                 bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_1)=r_t_group_names



  ###############  s2 #######

  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_1=rbind(simu_out_alsp_s2_samsz_1_N_max_1,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                 bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_1,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_1)=r_t_group_names


  ## samsz_2###


  ##### s1 ####

  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_1=rbind(simu_out_alsp_s1_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                 bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_1,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_1)=r_t_group_names



  ###############  s2 #######

  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_1=rbind(simu_out_alsp_s2_samsz_2_N_max_1,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
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
    simu_out_alsp_s1_samsz_1_N_max_2=rbind(simu_out_alsp_s1_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                 bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_1_N_max_2)=r_t_group_names



  ###############  s2 #######

  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_alsp_s2_samsz_1_N_max_2=rbind(simu_out_alsp_s2_samsz_1_N_max_2,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                                                 bry_type=2,alpha_2=alpha_2_alps_s2,theta_cut=theta_cut,samsz=samsz_1,N_max=N_max_2,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s2_samsz_1_N_max_2)=r_t_group_names


  ## samsz_2###


  ##### s1 ####

  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_alsp_s1_samsz_2_N_max_2=rbind(simu_out_alsp_s1_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                                                 bry_type=2,alpha_2=alpha_2_alps_s1,theta_cut=theta_cut,samsz=samsz_2,N_max=N_max_2,
                                                                                                                 beta=beta,sim_num=sim_num))
  }
  rownames(simu_out_alsp_s1_samsz_2_N_max_2)=r_t_group_names



  ###############  s2 #######

  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:dim(r_t_vec)[1]){
    simu_out_alsp_s2_samsz_2_N_max_2=rbind(simu_out_alsp_s2_samsz_2_N_max_2,two_stg_adapt_power_OC_simu_neg_binomial(r_t=r_t_vec[i,],r_c=r_c,kappa=kappa,direction=direction,s=s2,
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
#Two_stage_asd_simu_OC_neg_binomial(r_t_vec_1,r_t_vec_2,r_c,kappa,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num)
#Two_stage_asd_simu_OC_neg_binomial(groups=rbind(c(0.3,0.4,0.5),c(0.2,0.2,0.5),c(0.2,0.2,0.3),c(0.2,0.3,0.4)),r_c=0.2,kappa=1,direction=1,
#                                   s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                                   theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1,sim_num=1000)




