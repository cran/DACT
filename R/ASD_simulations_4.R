# install.packages("glm")
# install.packages("glm2")
# install.packages("MASS")
#library(tidyverse)
#library(ggplot2)



###########################################################
#### 1.   simulations for normal distributions     ########
###########################################################

###########################################################
#### 1.1  GSD simulations for normal distributions  ########
###########################################################

gsd_rej_rate_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,samsz,sim_num){
  #This function simulates the rejection rate, given mu, s,c_bry,sim_num,sigma
  ## for a group sequential design
  ## the random samples are from normal distributions
  #mu is the effect size
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsaz is the sample size
  #sigma is the common stdard deviation
  ## When mu=0, the rejection rate is the type I error
  ## When mu>0, the rejection rate is the power
  if(direction==1){ #larger mu is better
    mu=mu_t-mu_c
  }else if (direction==2){ #smaller mu is better
    mu=mu_c-mu_t
  }
  ######################################
  ###  Preparing the parameters ########
  ######################################
  N=ceiling(samsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted

  t<- N/(sigma_t^2+sigma_c^2)    # t is the vector of information times at each interim anbalysis
  NF=N[length(s)]  ## NF is the sample size at the end of the trial
  #t is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      ########################################
      ### Step 1: Generate random samples ####
      ########################################
      x0=stats::rnorm(NF,0,sigma_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rnorm(NF,mu,sigma_t)

      #####################################################################
      ### Step 2: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      Wald=NULL
      for (j in 1:length(s) ){
        Wald=c(Wald,(mean(x[1:N[j]])-mean(x0[1:N[j]]))*sqrt(t[j]))
      }

      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry
      #at each simulation, test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if ( max(test)>0){

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        #stg identifies the interim analysis at which the null hypothesis is rejected
      }else{
        rej=0
        stg=0
      }
      if (stg>0){
        last_vt=stg
      } else {last_vt=length(s)
      }
      #######################################################
      ### Step 4: output and finish the simjulation loop ####
      #######################################################
      samsz_simu=N[last_vt]
      out=c(rej,stg,samsz_simu)
  })#(end of q loop)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  samsz_simu = as.numeric(simr[,3])
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  rej_num=sum(rej>0)

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  mean_samsz=mean(samsz_simu)

  #############################
  ### Step 6: Final output ####
  #############################

  out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(rej_time)))
  names(out_1)=c("Looks","Information fraction","Critical boundaries","Rejection by look")
  out_2=c(rej_rt,samsz,mean_samsz)
  names(out_2)=c("Rejection rate","Planned control sample size","Mean control sample size")
  out=list(out_1,out_2)
  names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)
}

#gsd_rej_rate_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,samsz,sim_num)

#gsd_rej_rate_normal(mu_t=0.5,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.470381, 2.453930, 2.003625),samsz=86,sim_num=1000)
#gsd_rej_rate_normal(mu_t=0,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.470381, 2.453930, 2.003625),samsz=100,sim_num=100000)
##############################################

#' Group sequential designs / Adaptive sequential designs: Simulations, Normal distribution, Type I error and power
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param sigma_t Common standard deviation of the test arm and the control arm; must be > 0
#' @param sigma_c Common standard deviation of the test arm and the control arm; must be > 0
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param samsz Planned sample size of the control arm used in the simulation; must be > 0
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries and Rejection by look (the
#'     rejection rate at each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Rejection rate (the simulated
#'     power, i.e. the type I error rate when mu_t = mu_c), Planned control sample size (the
#'     sample size used in the simulation) and Mean control sample size (the mean simulated
#'     sample size)}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' gsd_power_simu_normal(mu_t = 0.4, mu_c = 0, sigma_t = 1, sigma_c = 1, direction = 1,
#' s = c(1/4,1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
gsd_power_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,samsz,sim_num){
 # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
 # mu=abs(mu_t-mu_c)
  alpha=alpha_2/2

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  simu_out=gsd_rej_rate_normal(mu_t=mu_t,mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s,c_bry=c_bry,samsz=samsz,sim_num=sim_num)
  #out_1=simu_out[[1]]
  #out_2=c(samsz,as.numeric(unlist(simu_out[[2]])))
  #names(out_2)=c("Planned control sample size","rejection rate","Mean control sample size")
  #out=list(out_1,out_2)
  out=simu_out

  #names(out)=c("Sequential design and rejection by look","Sample size and rejection rate")
  return(out)
}

##################################
#gsd_power_simu_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,samsz,sim_num)
#gsd_power_simu_normal(mu_t=0.4,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/4,1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)
#gsd_power_simu_normal(mu_t=0,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/4,1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=100)
#gsd_power_simu_normal(mu_t=0.5,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=86,sim_num=10000)

#gsd_power_simu_normal(mu_t=0.5,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=86,sim_num=10000)
#gsd_power_simu_normal(mu_t=0.,mu_c=0.5,sigma_t=1,sigma_c=1,direction=2,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=86,sim_num=100000)

#gsd_power_simu_normal(mu_t=0,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=0.05,samsz=100,sim_num=100000)

#gsd_power_simu_normal(mu_t=0,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=3,alpha_2=c(0.001,0.001,0.048),samsz=100,sim_num=100000)


###########################################################################
#### 1.2  Adaptive sequential simulations for normal distributions  ########
###########################################################################

########################################################################
asd_rej_simu_adapt_com_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num){
  #mu is the effect size ##
  #s is the vector of information frations
  #c_bry is the vector of critical values
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsaz is sample size
  #N_max is themaximum sample size
  #sigma is the common stdard deviation
  if(direction==1){ #larger mu is better
    mu=mu_t-mu_c
  }else if (direction==2){ #smaller mu is better
    mu=mu_c-mu_t
  }
  ###############################################################
  ###############################################
  ############## generate the parameters ########
  ###############################################
  T=samsz/(sigma_t^2+sigma_c^2)
  N=ceiling(samsz*s)
  # N is the vetoc of sample sizes at which interim analysis will be conducted
  t<- N/(sigma_t^2+sigma_c^2)
  #t is the vector of information times at the interim analysis
  NF=N[length(s)]


  #t is the vector of information times at the interim analysis
  e_bry=c_bry*sqrt(t)

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {
      #the loop specifies what is generated for each simulation

     #rm(ampsz_simu,case,cond_power,out,futi)
      ########################################
      ### Step 1: Generate random samples ####
      ########################################

      x0=stats::rnorm(N_max,0,sigma_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rnorm(N_max,mu,sigma_t)

      #############################################################################################
      ### Step 2: Obtain the Wald statistics and theta_hat for hypothesis testing and estimate ####
      #############################################################################################

      Wald=NULL
      theta_hat=NULL
      for (j in 1:length(s)){
        Wald= c(Wald,(mean(x[1:N[j]])-mean(x0[1:N[j]]))*sqrt(t[j]))
        theta_hat= c(theta_hat,(mean(x[1:N[j]])-mean(x0[1:N[j]])))
      }

      theta_hat_sd=1/sqrt(t)
      #Wald containes the Wald statistics
      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd=1/sqrt(t[length(s)-1])
      if (theta_hat_inter<=theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry

      test_1=test[1:length(s)-1]
      # test compares the Wald statistics with the critical boundaries
     #############################

      if ( max(test)>=0){
        rej_nossr=1
      }else if ( max(test)<0){
        rej_nossr=0
      }
      ##############
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling
        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        wald_last=Wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        sampsz_simu=N[stg]


        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else  if ( max(test_1)<=0){#


        #################################################################################
        ###                   Step 4.1: Compare theta_hat with theta_cut             ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################


        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation

          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=1/sqrt(t[last_vt])

          if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }
          sampsz_simu=N[length(s)]

          #############  #############

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          score_inter=theta_hat_inter*t[inter_vt]

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))

          if (cond_power>=1-beta){

            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=1/sqrt(t[last_vt])

            sampsz_simu=N[length(s)]
            if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (cond_power<1-beta){

            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##

            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(T-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_sampsz=ceiling(T_new_0*(sigma_t^2+sigma_c^2))
            T_new=new_sampsz/(sigma_t^2+sigma_c^2)
            e_bry_new=
              (sqrt(T_new-t[inter_vt])/sqrt(T-t[inter_vt]))*(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)+ score_inter

            c_bry_new=e_bry_new/sqrt(T_new)

            if  (new_sampsz>N_max){
              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (last(test)>0){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              sampsz_simu=last(N)
            } else if  (new_sampsz<=N_max){
              case=5
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              sampsz_simu=new_sampsz
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              x0_new=x0[1:new_sampsz]

              x_new=x[1:new_sampsz]

              score_last=(mean(x_new)-mean(x0_new))*T_new
              wald_last=(mean(x_new)-mean(x0_new))*sqrt(T_new)

              #################################################################################
              ###                   Step 4.6: Hypothesis testing                           ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              test_add=wald_last-c_bry_new

              if (test_add>0){
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
      out=c(rej_nossr,rej,stg,sampsz_simu,futi)
  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################
  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
    #################
  }
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>samsz){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  sampsz_icr_rt_byLook=c(rep(0,(length(s)-1)),mean(sampsz_icr))
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################


  out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(sampsz_icr_rt_byLook),as.numeric(rej_time)))
  names(out_1)=c("Looks","Information fraction","Critical boundaries","Sample size increase","Rejection by look")
  out_2=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,samsz,N_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out_2)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  out=list(out_1,out_2)
  names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)

}
#################
#asd_rej_simu_adapt_com_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num)
#asd_rej_simu_adapt_com_normal(mu_t=0,mu_c=0,sigma_t=1,sigma_c=1, direction=1,s=c(1/3,2/3,1),c_bry=c(3.470381, 2.453930, 2.003625),samsz=100,N_max=500,theta_cut=0.05,beta=0.1,sim_num=1000)
#asd_rej_simu_adapt_com_normal(mu_t=0.2,mu_c=0,sigma_t=1,sigma_c=1, direction=1,s=c(1/3,2/3,1),c_bry=c(3.470381, 2.453930, 2.003625),samsz=100,N_max=500,theta_cut=0.05,beta=0.1,sim_num=1000)
#asd_rej_simu_adapt_com_normal(mu_t=0.5,mu_c=0,sigma_t=1,sigma_c=1, direction=1,s=c(1/3,2/3,1),c_bry=c(3.470381, 2.453930, 2.003625),samsz=86,N_max=500,theta_cut=0.05,beta=0.1,sim_num=100000)

########################################

#' Group sequential designs / Adaptive sequential designs: Simulations, Normal distribution, Operating characteristics
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param sigma_t Common standard deviation of the test arm and the control arm; must be > 0
#' @param sigma_c Common standard deviation of the test arm and the control arm; must be > 0
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param samsz Planned sample size of the control arm; the trial starts from this sample size and
#'   it is increased only if the conditional power at the interim analysis is below the target
#' @param N_max Maximum sample size of the control arm allowed for the sample size change; the
#'   re-estimated sample size is truncated at N_max
#' @param theta_cut Lack of efficacy threshold on the effect-size scale: if the observed interim
#'   effect size is below theta_cut, the trial continues to the planned sample size with no
#'   sample size increase and no early termination for futility
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries, Sample size increase (the rate of
#'     sample size increase attributed to each look) and Rejection by look (the rejection rate at
#'     each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Rejection rate without SSR,
#'     Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size,
#'     Maximum control sample size, Mean control sample size and
#'     mean control sample size conditional on increase}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_power_simu_normal(mu_t = 0, mu_c = 0, sigma_t = 1, sigma_c = 1, direction = 1,
#' s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, samsz = 100, N_max = 500, theta_cut = 0.05,
#' beta = 0.1, sim_num = 20)
#'
#' @export
asd_power_simu_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,samsz,N_max,theta_cut,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate
  alpha=alpha_2/2

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  # samsz=control_sample_size
  simu_out=asd_rej_simu_adapt_com_normal(mu_t=mu_t,mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s,c_bry=c_bry,samsz=samsz,N_max=N_max,theta_cut=theta_cut,beta=beta,sim_num=sim_num)
 # simu_out_1=simu_out[1]
 # names(simu_out_1)=""


 # simu_out_2_0=c(tot_planned_power,as.numeric(unlist(simu_out[2])))
#  simu_out_2_1=data.frame(simu_out_2_0)
 # colnames(simu_out_2_1)=NULL
#  simu_out_2=data.frame(t(simu_out_2_1))

#  colnames(simu_out_2)=c("Power at planned sample size",names(unlist(simu_out[[2]])))

  out=simu_out
  #names(out)=c("Sequential design","Simulation results")
 # names(out)=c("Sequential design and rejection by look","Power and sample size")

  return(out)
}

##################################
#asd_power_simu_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,bry_type,alpha_2,samsz,N_max,theta_cut,beta,sim_num)

#asd_power_simu_normal(mu_t=0,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,N_max=500,theta_cut=0.05,beta=0.1,sim_num=1000000)
#asd_power_simu_normal(mu_t=0.5,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=86,N_max=150,theta_cut=0.05,beta=0.1,sim_num=1000)
#asd_power_simu_normal(mu_t=0,mu_c=0.5,sigma_t=1,sigma_c=1,direction=2,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=86,N_max=150,theta_cut=0.05,beta=0.1,sim_num=1000000)

#####################################################
########################################################################

########################################################################
asd_rej_simu_adapt_OC_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num){
  #mu is the effect size ##
  #s is the vector of information frations
  #c_bry is the vector of critical values
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsaz is sample size
  #N_max is themaximum sample size
  #sigma is the common stdard deviation
  if(direction==1){ #larger mu is better
    mu=mu_t-mu_c
  }else if (direction==2){ #smaller mu is better
    mu=mu_c-mu_t
  }
  ###############################################################
  ###############################################
  ############## generate the parameters ########
  ###############################################
  T=samsz/(sigma_t^2+sigma_c^2)
  N=ceiling(samsz*s)
  # N is the vetoc of sample sizes at which interim analysis will be conducted
  t<- N/(sigma_t^2+sigma_c^2)
  #t is the vector of information times at the interim analysis
  NF=N[length(s)]


  #t is the vector of information times at the interim analysis
  e_bry=c_bry*sqrt(t)

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {
      #the loop specifies what is generated for each simulation

      #rm(ampsz_simu,case,cond_power,out,futi)
      ########################################
      ### Step 1: Generate random samples ####
      ########################################

      x0=stats::rnorm(N_max,0,sigma_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rnorm(N_max,mu,sigma_t)

      #############################################################################################
      ### Step 2: Obtain the Wald statistics and theta_hat for hypothesis testing and estimate ####
      #############################################################################################

      Wald=NULL
      theta_hat=NULL
      for (j in 1:length(s)){
        Wald= c(Wald,(mean(x[1:N[j]])-mean(x0[1:N[j]]))*sqrt(t[j]))
        theta_hat= c(theta_hat,(mean(x[1:N[j]])-mean(x0[1:N[j]])))
      }

      theta_hat_sd=1/sqrt(t)
      #Wald containes the Wald statistics
      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd=1/sqrt(t[length(s)-1])
      if (theta_hat_inter<=theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry

      test_1=test[1:length(s)-1]
      # test compares the Wald statistics with the critical boundaries
      #############################

      if ( max(test)>=0){
        rej_nossr=1
      }else if ( max(test)<0){
        rej_nossr=0
      }
      ##############
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling
        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        wald_last=Wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        sampsz_simu=N[stg]


        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else  if ( max(test_1)<=0){#


        #################################################################################
        ###                   Step 4.1: Compare theta_hat with theta_cut             ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################


        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation

          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=1/sqrt(t[last_vt])

          if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }
          sampsz_simu=N[length(s)]

          #############  #############

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          score_inter=theta_hat_inter*t[inter_vt]

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))

          if (cond_power>=1-beta){

            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=1/sqrt(t[last_vt])

            sampsz_simu=N[length(s)]
            if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (cond_power<1-beta){

            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##

            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(T-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_sampsz=ceiling(T_new_0*(sigma_t^2+sigma_c^2))
            T_new=new_sampsz/(sigma_t^2+sigma_c^2)
            e_bry_new=
              (sqrt(T_new-t[inter_vt])/sqrt(T-t[inter_vt]))*(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)+ score_inter

            c_bry_new=e_bry_new/sqrt(T_new)

            if  (new_sampsz>N_max){
              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (last(test)>0){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              sampsz_simu=last(N)
            } else if  (new_sampsz<=N_max){
              case=5
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              sampsz_simu=new_sampsz
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              x0_new=x0[1:new_sampsz]

              x_new=x[1:new_sampsz]

              score_last=(mean(x_new)-mean(x0_new))*T_new
              wald_last=(mean(x_new)-mean(x0_new))*sqrt(T_new)

              #################################################################################
              ###                   Step 4.6: Hypothesis testing                           ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              test_add=wald_last-c_bry_new

              if (test_add>0){
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
      out=c(rej_nossr,rej,stg,sampsz_simu,futi)
  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>samsz){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################

  out=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,samsz,N_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  return(out)

}

######################
#asd_rej_simu_adapt_OC_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num)
####################
#asd_rej_simu_adapt_OC_normal(mu_t=0.5,mu_c=0,sigma_t=1,sigma_c=1,direction=1,s=c(1/3,2/3,1),
#                            c_bry=c(3.470381, 2.453930, 2.003625),
 #                        samsz=100, N_max=200, theta_cut=0.01,beta=0.1,sim_num=100)
########################################

##########################

#' Group sequential designs / Adaptive sequential designs: Simulations, Normal distribution, Optimizing adaptive sequential design
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param sigma_t Common standard deviation of the test arm and the control arm; must be > 0
#' @param sigma_c Common standard deviation of the test arm and the control arm; must be > 0
#' @param direction Test direction: 1 = larger mu is better, otherwise smaller mu is better
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g.
#'   c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g.
#'   c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param alpha_2_alps_s1 Two-sided type I error spent at the interim analyses of schedule s1 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error spent at the interim analyses of schedule s2 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param theta_cut Lack of efficacy threshold on the effect-size scale: if the observed interim
#'   effect size is below theta_cut, the trial continues to the planned sample size with no
#'   sample size increase and no early termination for futility
#' @param samsz_1 Planned sample size of the control arm for the first stage of the simulation grid
#' @param samsz_2 Planned sample size of the control arm for the second stage of the simulation grid
#' @param N_max_1 Maximum sample size of the control arm allowed for the sample size change in the
#'   first scenario of the simulation grid
#' @param N_max_2 Maximum sample size of the control arm allowed for the sample size change in the
#'   second scenario of the simulation grid
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with six components, one for each combination of critical boundary type
#' (O'Brien-Fleming, Pocock, alpha-spending) and interim analysis schedule (s1, s2):
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{Named list with the group
#'     sequential design (Looks, Information fractions, O'Brien-Fleming boundaries) and the
#'     simulation results for the four combinations of planned sample size (samsz_1, samsz_2) and
#'     maximum sample size (N_max_1, N_max_2), each a data frame with one row per test-arm mean and
#'     the columns mu_t1,..., mu_c, sample size for fixed design, Rejection rate without SSR,
#'     Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size,
#'     Maximum control sample size, Mean control sample size,
#'     mean control sample size conditional on increase and
#'     Ratio of adaptive sample size over fixed sample size}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with the information fractions s2}
#'   \item{Pocock boundaries, interim analysis schedule=s1}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s1}
#'   \item{Pocock boundaries, interim analysis schedule=s2}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s2}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s1}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s2}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' asd_simu_OC_normal(mu_t = c(0.3,0.4,0.5), mu_c = 0.2, sigma_t = 1, sigma_c = 1, direction = 1,
#' s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05, alpha_2_alps_s1 = c(0.001,0.001,0.048),
#' alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01, samsz_1 = 100, samsz_2 = 150,
#' N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 20)
#' }
#'
#' @export
asd_simu_OC_normal<-function(mu_t,mu_c,sigma_t,sigma_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)


  ###### sample size calculations for fixed design ####

  if(direction==1){ #larger mu is better
    mu=mu_t-mu_c
  }else if (direction==2){ #smaller mu is better
    mu=mu_c-mu_t
  }
  sampsz_mu=NULL

  for (i in 1:length(mu)){
    sampsz_mu[i]=ceiling((((a+b)/mu[i])**2)*(sigma_t**2+sigma_c**2))
  }

  ############ critical boundary calculations #####

  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(OF_bundry(s1,alpha)[1]))
  s2_OF_cbry=as.numeric(unlist(OF_bundry(s2,alpha)[1]))

  s1_pocock_cbry=as.numeric(unlist(pocock_bundry(s1,alpha)[1]))
  s2_pocock_cbry=as.numeric(unlist(pocock_bundry(s2,alpha)[1]))

  s1_alsp_cbry=as.numeric(unlist(alsp_bundry(s1,alpha_alsp_s1)[1]))
  s2_alsp_cbry=as.numeric(unlist(alsp_bundry(s2,alpha_alsp_s2)[1]))

  ##### GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_pocock_design=data.frame(cbind(c(1:length(s1)),s1,s1_pocock_cbry))
  names(s1_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s2_pocock_design=data.frame(cbind(c(1:length(s2)),s2,s2_pocock_cbry))
  names(s2_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")

  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_result_title = c(LABEL()$lambda_t1,LABEL()$lambda_c,"sample size  for fixed design", "Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                        "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                        "Futility rate","Planned control sample size","Maximum control sample size",
                        "Mean control sample size","mean control sample size conditional on increase",
                        "Ratio of adaptive sample size over fixed sample size")
  simu_out_OF_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_OF_s1_samsz_1_N_max_1= rbind(simu_out_OF_s1_samsz_1_N_max_1,
                                          asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                       c_bry=s1_OF_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }


  samsz_ratio_OF_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_OF_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_OF_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_1_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_OF_s1_samsz_1_N_max_1,samsz_ratio_OF_s1_samsz_1_N_max_1)))
  names(simu_out_OF_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_OF_s2_samsz_1_N_max_1= rbind(simu_out_OF_s2_samsz_1_N_max_1,
                                          asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                       c_bry=s2_OF_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_OF_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_OF_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s2_samsz_1_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_OF_s2_samsz_1_N_max_1,samsz_ratio_OF_s2_samsz_1_N_max_1)))
  names(simu_out_OF_s2_samsz_1_N_max_1)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_OF_s1_samsz_2_N_max_1= rbind(simu_out_OF_s1_samsz_2_N_max_1,
                                          asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                       c_bry=s1_OF_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_OF_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_OF_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_2_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_OF_s1_samsz_2_N_max_1,samsz_ratio_OF_s1_samsz_2_N_max_1)))
  names(simu_out_OF_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_OF_s2_samsz_2_N_max_1= rbind(simu_out_OF_s2_samsz_2_N_max_1,
                                          asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                       c_bry=s2_OF_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_OF_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_OF_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s2_samsz_2_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_OF_s2_samsz_2_N_max_1,samsz_ratio_OF_s2_samsz_2_N_max_1)))
  names(simu_out_OF_s2_samsz_2_N_max_1)=simu_result_title

  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####
  simu_out_OF_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_OF_s1_samsz_1_N_max_2= rbind(simu_out_OF_s1_samsz_1_N_max_2,
                                          asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                       c_bry=s1_OF_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_OF_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_OF_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_1_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_OF_s1_samsz_1_N_max_2,samsz_ratio_OF_s1_samsz_1_N_max_2)))
  names(simu_out_OF_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_OF_s2_samsz_1_N_max_2= rbind(simu_out_OF_s2_samsz_1_N_max_2,
                                          asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                       c_bry=s2_OF_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_OF_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_OF_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s2_samsz_1_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_OF_s2_samsz_1_N_max_2,samsz_ratio_OF_s2_samsz_1_N_max_2)))
  names(simu_out_OF_s2_samsz_1_N_max_2)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_OF_s1_samsz_2_N_max_2= rbind(simu_out_OF_s1_samsz_2_N_max_2,
                                          asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                       c_bry=s1_OF_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_OF_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_OF_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_2_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_OF_s1_samsz_2_N_max_2,samsz_ratio_OF_s1_samsz_2_N_max_2)))
  names(simu_out_OF_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_OF_s2_samsz_2_N_max_2= rbind(simu_out_OF_s2_samsz_2_N_max_2,
                                          asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                       c_bry=s2_OF_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_OF_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_OF_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s2_samsz_2_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_OF_s2_samsz_2_N_max_2,samsz_ratio_OF_s2_samsz_2_N_max_2)))
  names(simu_out_OF_s2_samsz_2_N_max_2)=simu_result_title

  #############################################
  #### OC simulations for Pocock boundaries ###
  #############################################


  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_pocock_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_pocock_s1_samsz_1_N_max_1= rbind(simu_out_pocock_s1_samsz_1_N_max_1,
                                              asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                           c_bry=s1_pocock_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_pocock_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_pocock_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_1_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_pocock_s1_samsz_1_N_max_1,samsz_ratio_pocock_s1_samsz_1_N_max_1)))
  names(simu_out_pocock_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_pocock_s2_samsz_1_N_max_1= rbind(simu_out_pocock_s2_samsz_1_N_max_1,
                                              asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                           c_bry=s2_pocock_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_pocock_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_pocock_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_1_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_pocock_s2_samsz_1_N_max_1,samsz_ratio_pocock_s2_samsz_1_N_max_1)))
  names(simu_out_pocock_s2_samsz_1_N_max_1)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_pocock_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_pocock_s1_samsz_2_N_max_1= rbind(simu_out_pocock_s1_samsz_2_N_max_1,
                                              asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                           c_bry=s1_pocock_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_pocock_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_pocock_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_2_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_pocock_s1_samsz_2_N_max_1,samsz_ratio_pocock_s1_samsz_2_N_max_1)))
  names(simu_out_pocock_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_pocock_s2_samsz_2_N_max_1= rbind(simu_out_pocock_s2_samsz_2_N_max_1,
                                              asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                           c_bry=s2_pocock_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_pocock_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_pocock_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_2_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_pocock_s2_samsz_2_N_max_1,samsz_ratio_pocock_s2_samsz_2_N_max_1)))
  names(simu_out_pocock_s2_samsz_2_N_max_1)=simu_result_title

  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####
  simu_out_pocock_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_pocock_s1_samsz_1_N_max_2= rbind(simu_out_pocock_s1_samsz_1_N_max_2,
                                              asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                           c_bry=s1_pocock_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_pocock_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_pocock_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_1_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_pocock_s1_samsz_1_N_max_2,samsz_ratio_pocock_s1_samsz_1_N_max_2)))
  names(simu_out_pocock_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_pocock_s2_samsz_1_N_max_2= rbind(simu_out_pocock_s2_samsz_1_N_max_2,
                                              asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                           c_bry=s2_pocock_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_pocock_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_pocock_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_1_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_pocock_s2_samsz_1_N_max_2,samsz_ratio_pocock_s2_samsz_1_N_max_2)))
  names(simu_out_pocock_s2_samsz_1_N_max_2)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_pocock_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_pocock_s1_samsz_2_N_max_2= rbind(simu_out_pocock_s1_samsz_2_N_max_2,
                                              asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                           c_bry=s1_pocock_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_pocock_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_pocock_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_2_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_pocock_s1_samsz_2_N_max_2,samsz_ratio_pocock_s1_samsz_2_N_max_2)))
  names(simu_out_pocock_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_pocock_s2_samsz_2_N_max_2= rbind(simu_out_pocock_s2_samsz_2_N_max_2,
                                              asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                           c_bry=s2_pocock_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_pocock_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_pocock_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_2_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_pocock_s2_samsz_2_N_max_2,samsz_ratio_pocock_s2_samsz_2_N_max_2)))
  names(simu_out_pocock_s2_samsz_2_N_max_2)=simu_result_title

  ##############

  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_alsp_s1_samsz_1_N_max_1= rbind(simu_out_alsp_s1_samsz_1_N_max_1,
                                            asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                         c_bry=s1_alsp_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_alsp_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_alsp_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_1_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_alsp_s1_samsz_1_N_max_1,samsz_ratio_alsp_s1_samsz_1_N_max_1)))
  names(simu_out_alsp_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_alsp_s2_samsz_1_N_max_1= rbind(simu_out_alsp_s2_samsz_1_N_max_1,
                                            asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                         c_bry=s2_alsp_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_alsp_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_alsp_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_1_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_alsp_s2_samsz_1_N_max_1,samsz_ratio_alsp_s2_samsz_1_N_max_1)))
  names(simu_out_alsp_s2_samsz_1_N_max_1)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_alsp_s1_samsz_2_N_max_1= rbind(simu_out_alsp_s1_samsz_2_N_max_1,
                                            asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                         c_bry=s1_alsp_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_alsp_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_alsp_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_2_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_alsp_s1_samsz_2_N_max_1,samsz_ratio_alsp_s1_samsz_2_N_max_1)))
  names(simu_out_alsp_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(mu_t)){
    simu_out_alsp_s2_samsz_2_N_max_1= rbind(simu_out_alsp_s2_samsz_2_N_max_1,
                                            asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                         c_bry=s2_alsp_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_alsp_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_alsp_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_2_N_max_1=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_alsp_s2_samsz_2_N_max_1,samsz_ratio_alsp_s2_samsz_2_N_max_1)))
  names(simu_out_alsp_s2_samsz_2_N_max_1)=simu_result_title

  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_alsp_s1_samsz_1_N_max_2= rbind(simu_out_alsp_s1_samsz_1_N_max_2,
                                            asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                         c_bry=s1_alsp_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_alsp_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_alsp_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_1_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_alsp_s1_samsz_1_N_max_2,samsz_ratio_alsp_s1_samsz_1_N_max_2)))
  names(simu_out_alsp_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_alsp_s2_samsz_1_N_max_2= rbind(simu_out_alsp_s2_samsz_1_N_max_2,
                                            asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                         c_bry=s2_alsp_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_alsp_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_alsp_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_1_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_alsp_s2_samsz_1_N_max_2,samsz_ratio_alsp_s2_samsz_1_N_max_2)))
  names(simu_out_alsp_s2_samsz_1_N_max_2)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_alsp_s1_samsz_2_N_max_2= rbind(simu_out_alsp_s1_samsz_2_N_max_2,
                                            asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s1,
                                                                         c_bry=s1_alsp_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_alsp_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_alsp_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_2_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_alsp_s1_samsz_2_N_max_2,samsz_ratio_alsp_s1_samsz_2_N_max_2)))
  names(simu_out_alsp_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(mu_t)){
    simu_out_alsp_s2_samsz_2_N_max_2= rbind(simu_out_alsp_s2_samsz_2_N_max_2,
                                            asd_rej_simu_adapt_OC_normal(mu_t=mu_t[i],mu_c=mu_c,sigma_t=sigma_t,sigma_c=sigma_c,direction=direction,s=s2,
                                                                         c_bry=s2_alsp_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(mu_t))){
    samsz_ratio_alsp_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_alsp_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_2_N_max_2=data.frame((cbind(mu_t,mu_c,sampsz_mu,simu_out_alsp_s2_samsz_2_N_max_2,samsz_ratio_alsp_s2_samsz_2_N_max_2)))
  names(simu_out_alsp_s2_samsz_2_N_max_2)=simu_result_title

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


  list_pocock_1=list(s1_pocock_design,
                     simu_out_pocock_s1_samsz_1_N_max_1,
                     simu_out_pocock_s1_samsz_2_N_max_1,
                     simu_out_pocock_s1_samsz_1_N_max_2,
                     simu_out_pocock_s1_samsz_2_N_max_2
                     )
  names(list_pocock_1)=c("Group sequential design",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 2")


  list_pocock_2=list(s2_pocock_design,
                     simu_out_pocock_s2_samsz_1_N_max_1,
                     simu_out_pocock_s2_samsz_2_N_max_1,
                     simu_out_pocock_s2_samsz_1_N_max_2,
                     simu_out_pocock_s2_samsz_2_N_max_2
  )
  names(list_pocock_2)=c("Group sequential design",
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
           list_pocock_1,list_pocock_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Pocock boundaries, interim analysis schedule=s1",
               "Pocock boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}

##################################
#asd_simu_OC_normal(mu_t,mu_c,sigma_t,sigma_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num)
################################################

#####################
#asd_simu_OC_normal(mu_t=c(0.3,0.4,0.5), mu_c=0,sigma_t=1,sigma_c=1,direction=1,s1=c(1/3,2/3,1),s2=c(0.5,1),
#                  alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#               theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1,sim_num=100)



#########################################################
#### 2.   simulations for survival analysis      ########
#########################################################

###########################################################
#### 2.1  GSD simulations for survival analysis    ########
###########################################################

#####################################################
gsd_rej_rate_survival<-function(HR,lambda,s,c_bry,surv_rate,evt_num,samp_multip,sim_num){


    ###############################################################
    ###############################################
    ############## generate the parameters ########
    ###############################################

    evt_num_s=ceiling(evt_num*s)
    sampsz=ceiling(samp_multip*evt_num/surv_rate)

    ############################################
    ####  begining of the simulation loops  ####
    ############################################
    simr <- simuloop(sim_num, "ASD", function() {

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
        # for (i in 1:sampsz){
        #   evt_count[i]=sum(x0_status[1:i])+sum(x_status[1:i]) #number of events after 2*i observations
        # }
        total_evt_count = 0
        for (i in 1:sampsz){
          total_evt_count = total_evt_count + x0_status[i] + x_status[i]
          evt_count[i]=total_evt_count
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

        ###################################################################################################
        ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
        ###################################################################################################

        test=wald[1:length(s)]-c_bry

        # test compares the Wald statistics with the critical boundaries
        #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
        if ( max(test)>0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

          rej=1
          test_t=test[test>0]
          stg=match(test_t[1],test)
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=stg
          #  last_vt=2
          wald_last=wald[last_vt]
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=theta_hat_sd[last_vt]
          samsz_simu=evt_num_s[last_vt]

          #################################################################################
          ###                   Step 4: Compare theta_hat with theta_cut               ####
          ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
        }else {
              rej=0
              stg=0
              last_vt=length(s)
              samsz_simu=evt_num_s[last_vt]
            }

        #######################################################
        ### Step 5: output and finish the simjulation loop ####
        #######################################################
        out=c(rej,stg,samsz_simu)
    })#(end of simulation loops)
    #############termination timepoint####
    ############################################
    ### Step 6: Generate summary statistics ####
    ############################################
    rej= as.numeric(simr[,1])
    stg= as.numeric(simr[,2])
    samsz_simu=as.numeric(simr[,3])
    ###########################
    rej_time=rep(0,length(s))
    for (j in 1:length(s)){
      rej_time[j]=length(stg[stg==j])/(dim(simr)[1])
      #rej_time[j] is the proportion
      #that the null hypothesis were rejected at the j-th interim analysis
      #################
    }
    power=mean(rej)

    av_samsz=mean(samsz_simu)
    rej_rt=mean(rej)
    ###################
    #############################
    ### Step 7: Final output ####
    #############################
    out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(rej_time)))
    names(out_1)=c("Looks","Information fraction","Critical boundaries","Rejection by look")
    out_2=c(rej_rt,evt_num,av_samsz)
    names(out_2)=c("Rejection rate","Planned sample size","Mean sample size")
    out=list(out_1,out_2)
    #names(out)=c("Sequential design","power and sample size")
    names(out)=c("Sequential design and rejection by look","Power and sample size")
    return(out)

}

#######################################
#gsd_rej_rate_survival(HR,lambda,s,c_bry,surv_rate,evt_num,samp_multip,sim_num)

#gsd_rej_rate_survival(HR=0.8,lambda=4,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),surv_rate=0.5,evt_num=845,samp_multip=5,sim_num=100)
#gsd_rej_rate_survival(HR=1,lambda=4,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),surv_rate=0.5,evt_num=845,samp_multip=5,sim_num=10000)

########################################################################

##############################################


#' Group sequential designs / Adaptive sequential designs: Simulations, Survival analysis, Type I error and power
#'
#' @param HR Hazards ratio of the test arm vs. the control arm; HR < 1 favours the test arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param surv_rate Survival rate at the end of the trial, used together with the number of events
#'   to derive the sample size; must be between 0 and 1
#' @param evt_num Planned total number of events at the final analysis; must be > 0
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries and Rejection by look (the
#'     rejection rate at each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Rejection rate (the simulated
#'     power, i.e. the type I error rate when HR = 1), Planned sample size (the planned number of
#'     events) and Mean sample size (the mean simulated number of events)}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' gsd_power_simu_survival(HR = 0.8, s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05,
#' surv_rate = 0.5, evt_num = 845, sim_num = 20)
#'
#' @export
gsd_power_simu_survival<-function(HR,s,bry_type,alpha_2,surv_rate,evt_num,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  lambda=1

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  samp_multip=5
  simu_out=gsd_rej_rate_survival(HR=HR,lambda=1,s=s,c_bry=c_bry,surv_rate=surv_rate,evt_num=evt_num,samp_multip=samp_multip,sim_num=sim_num)

 # out_1=simu_out[[1]]
 # out_2=c(evt_num,as.numeric(unlist(simu_out[[2]])))
 # names(out_2)=c("Planned sample size","Rejection rate","Mean sample size")
  out=simu_out
  #names(out)=c("Sequential design","Sample size and rejection rate")
  #names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)
}

##################################
#gsd_power_simu_survival(HR,s,bry_type,alpha_2,surv_rate,evt_num,sim_num)
#gsd_power_simu_survival(HR=0.8,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=845,sim_num=100)
#gsd_power_simu_survival(HR=0.8,s=c(1/3,2/3,1),bry_type=2,alpha_2=0.05,surv_rate=0.5,evt_num=845,sim_num=1000)

#gsd_power_simu_survival(HR=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,surv_rate=0.5,evt_num=845,sim_num=10000)

##############################################


###########################################################################
#### 2.1  Adaptive sequential simulations for survival analysis    ########
###########################################################################

asd_rej_simu_adapt_com_survival<-function(HR,lambda,s,c_bry,surv_rate,evt_num,evt_max,HR_cut,samp_multip,sim_num,beta){
   #s is the vector of information frations
  #c_bry is the vector of critical values
  ###############################################################
  ###############################################
  ############## generate the parameters ########
  ###############################################

  evt_num_s=c(ceiling(evt_num*s),evt_max)
  sampsz=ceiling(samp_multip*evt_max/surv_rate)
  theta_cut=-log(HR_cut)

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

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
      # for (i in 1:sampsz){
      #   evt_count[i]=sum(x0_status[1:i])+sum(x_status[1:i]) #number of events after 2*i observations
      # }
      total_evt_count = 0
      for (i in 1:sampsz){
        total_evt_count = total_evt_count + x0_status[i] + x_status[i]
        evt_count[i]=total_evt_count
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

      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd=theta_hat_sd[length(s)-1]
      if (theta_hat_inter<theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test_all=wald[1:length(s)]-c_bry

      test_1=test_all[1:length(s)-1]

      #############################

      if ( max(test_all>=0)){
        rej_nossr=1
      }else if ( max(test_all)<0){
        rej_nossr=0
      }
      ##############
      # test_all compares the Wald statistics with the critical boundaries
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        wald_last=wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        samsz_simu=evt_num_s[last_vt]

        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else {


        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation

          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=theta_hat_sd[last_vt]


          samsz_simu=evt_num_s[length(s)]
          if (test_all[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
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
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          event_inter=evt_num_s[inter_vt]
          score_inter=theta_hat_inter*t[inter_vt]

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))

          if (cond_power>=1-beta){
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=theta_hat_sd[last_vt]


            samsz_simu=evt_num_s[length(s)]
            if (test_all[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else{

            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##
            inter_vt=length(s)-1
            event_inter=evt_num_s[inter_vt]
            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_sampsz=ceiling(event_inter*T_new_0/t[inter_vt])

            if  (new_sampsz>evt_max){
              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (last(test_all)>0){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              samsz_simu=evt_num_s[length(s)]
            } else  if  (new_sampsz<=evt_max){
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              samsz_simu=new_sampsz
              case=5
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              evt_pos=evt_count[evt_count>new_sampsz]#substring of  evt_count that are larger than evt_num
              evt_cut=match(evt_pos[1],evt_count) #the first incidence that evt_count is larger than evt_num

              #  vector of zeros, length of the vector is the same as the original smaple size
              time=c(time_x0[1:evt_cut],time_x[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
              event=c(x0_status[1:evt_cut],x_status[1:evt_cut]) #censoring status before total number of events exceed evt_num
              grp=c(grp_x0[1:evt_cut],grp_x[1:evt_cut]) #vector of zero, with length=2*evt_cut


              test= list(time=time,status=event,grp)#survival data
              coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
              # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
              # test is the same as data=test
              theta_hat_last= coxcf$coefficients[1] #estimated log(HR)
              theta_hat_last_sd= coxcf$coefficients[3] #estimated standard error of log(HR)
              t_last=(theta_hat_last_sd)^(-2)
              T_new=t_last

              e_bry_new=
                (sqrt(T_new-t[inter_vt])/sqrt(t[length(s)]-t[inter_vt]))*(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)+ score_inter

              c_bry_new=e_bry_new/sqrt(T_new)

              wald_last=theta_hat_last/ theta_hat_last_sd
              score_last=theta_hat_last*t_last

              #################################################################################
              ###                   Step 4.6: Hypothesis testing                           ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              test_add=wald_last-c_bry_new
              samsz_simu=new_sampsz
              if (test_add>0){
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
        out=c(rej_nossr,rej,stg,samsz_simu,futi)

  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################
  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
    #################
  }
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>evt_num){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  sampsz_icr_rt_byLook=c(rep(0,(length(s)-1)),mean(sampsz_icr))
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################


  out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(sampsz_icr_rt_byLook),as.numeric(rej_time)))
  names(out_1)=c("Looks","Information fraction","Critical boundaries","Sample size increase","Rejection by look")
  out_2=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,evt_num,evt_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out_2)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  out=list(out_1,out_2)
  names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)

}
#######################################
#asd_rej_simu_adapt_com_survival(HR,lambda,s,c_bry,surv_rate,evt_num,evt_max,HR_cut,samp_multip,sim_num,beta)

#asd_rej_simu_adapt_com_survival(HR=1,lambda=4,s=c(1/3,2/3,1),c_bry=c(3.471095, 2.454434,2.004037),surv_rate=0.5,evt_num=845,evt_max=2000,HR_cut=0.9,samp_multip=5,sim_num=1000,beta=0.1)
##################################################################################################
#asd_rej_simu_adapt_com_survival(HR=0.8,lambda=4,s=c(0.5,1),c_bry=c(2^32,1.96),surv_rate=0.5,evt_num=845,evt_max=2000,HR_cut=0.9,samp_multip=5,sim_num=100,beta=0.1)


#######################################################################
##############################################

#' Group sequential designs / Adaptive sequential designs: Simulations, Survival analysis, Operating characteristics
#'
#' @param HR Hazards ratio of the test arm vs. the control arm; HR < 1 favours the test arm
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param evt_num Planned total number of events at the final analysis; the trial starts from this
#'   number of events and it is increased only if the conditional power at the interim analysis is
#'   below the target
#' @param evt_max Maximum number of events allowed for the sample size change; the re-estimated
#'   number of events is truncated at evt_max
#' @param HR_cut Lack of efficacy threshold on the hazards-ratio scale; the corresponding
#'   threshold on the log-hazards-ratio scale is theta_cut = -log(HR_cut). If the observed interim
#'   log-hazards-ratio is below theta_cut (i.e. the observed HR is greater than HR_cut), the trial
#'   continues to the planned number of events with no sample size increase and no early
#'   termination for futility
#' @param surv_rate Survival rate at the end of the trial, used together with the number of events
#'   to derive the sample size; must be between 0 and 1
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries, Sample size increase (the rate of
#'     sample size increase attributed to each look) and Rejection by look (the rejection rate at
#'     each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Rejection rate without SSR,
#'     Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size
#'     (the planned number of events), Maximum control sample size (the maximum number of events
#'     allowed), Mean control sample size (the mean simulated number of events) and
#'     mean control sample size conditional on increase}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_power_simu_survival(HR = 0.8, s = c(1/3,2/3,1), bry_type = 1, alpha_2 = 0.05, evt_num = 845,
#' evt_max = 1500, HR_cut = 0.95, surv_rate = 0.5, beta = 0.1, sim_num = 20)
#'
#' @export
asd_power_simu_survival<-function(HR,s,bry_type,alpha_2,evt_num,evt_max,HR_cut,surv_rate,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  lambda=1

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  simu_out=asd_rej_simu_adapt_com_survival(HR=HR,lambda=1,s=s,c_bry=c_bry,surv_rate=surv_rate,evt_num=evt_num,evt_max=evt_max,HR_cut=HR_cut,samp_multip=5,sim_num=sim_num,beta=beta)

  out=  simu_out

  return(out)
}

##################################
#asd_power_simu_survival(HR,s,bry_type,alpha_2,evt_num,evt_max,HR_cut,surv_rate,beta,sim_num)
#asd_power_simu_survival(HR=0.8,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,evt_num=845,evt_max=1500,HR_cut=0.95,surv_rate=0.5,beta=0.1,sim_num=100)
#asd_power_simu_survival(HR=0.78,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,evt_num=845,evt_max=1500,HR_cut=0.95,surv_rate=0.5,beta=0.1,sim_num=100)


#######################################
##################
asd_rej_simu_adapt_OC_survival<-function(HR,lambda,s,c_bry,surv_rate,evt_num,evt_max,HR_cut,samp_multip,sim_num,beta){
  #s is the vector of information frations
  #c_bry is the vector of critical values
  ###############################################################
  ###############################################
  ############## generate the parameters ########
  ###############################################

  evt_num_s=c(ceiling(evt_num*s),evt_max)
  sampsz=ceiling(samp_multip*evt_max/surv_rate)
  theta_cut=-log(HR_cut)

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

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
      # for (i in 1:sampsz){
      #   evt_count[i]=sum(x0_status[1:i])+sum(x_status[1:i]) #number of events after 2*i observations
      # }
      total_evt_count = 0
      for (i in 1:sampsz){
        total_evt_count = total_evt_count + x0_status[i] + x_status[i]
        evt_count[i]=total_evt_count
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

      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd=theta_hat_sd[length(s)-1]
      if (theta_hat_inter<theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test_all=wald[1:length(s)]-c_bry

      test_1=test_all[1:length(s)-1]

      #############################

      if ( max(test_all>=0)){
        rej_nossr=1
      }else if ( max(test_all)<0){
        rej_nossr=0
      }
      ##############
      # test_all compares the Wald statistics with the critical boundaries
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        wald_last=wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        samsz_simu=evt_num_s[last_vt]

        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else {


        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation

          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=theta_hat_sd[last_vt]


          samsz_simu=evt_num_s[length(s)]
          if (test_all[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
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
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          event_inter=evt_num_s[inter_vt]
          score_inter=theta_hat_inter*t[inter_vt]

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))

          if (cond_power>=1-beta){
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=theta_hat_sd[last_vt]


            samsz_simu=evt_num_s[length(s)]
            if (test_all[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else{

            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##
            inter_vt=length(s)-1
            event_inter=evt_num_s[inter_vt]
            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_sampsz=ceiling(event_inter*T_new_0/t[inter_vt])

            if  (new_sampsz>evt_max){
              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (last(test_all)>0){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              samsz_simu=evt_num_s[length(s)]
            } else  if  (new_sampsz<=evt_max){
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              samsz_simu=new_sampsz
              case=5
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              evt_pos=evt_count[evt_count>new_sampsz]#substring of  evt_count that are larger than evt_num
              evt_cut=match(evt_pos[1],evt_count) #the first incidence that evt_count is larger than evt_num

              #  vector of zeros, length of the vector is the same as the original smaple size
              time=c(time_x0[1:evt_cut],time_x[1:evt_cut]) #combined samples befroe total number of events exceed evt_num
              event=c(x0_status[1:evt_cut],x_status[1:evt_cut]) #censoring status before total number of events exceed evt_num
              grp=c(grp_x0[1:evt_cut],grp_x[1:evt_cut]) #vector of zero, with length=2*evt_cut


              test= list(time=time,status=event,grp)#survival data
              coxcf=summary(coxph(Surv(time,status)~grp,test)) #run the coxph function and get the estimates
              # ~grp means that coxph compares the HR of levels in grp (x vs. x0)
              # test is the same as data=test
              theta_hat_last= coxcf$coefficients[1] #estimated log(HR)
              theta_hat_last_sd= coxcf$coefficients[3] #estimated standard error of log(HR)
              t_last=(theta_hat_last_sd)^(-2)
              T_new=t_last

              e_bry_new=
                (sqrt(T_new-t[inter_vt])/sqrt(t[length(s)]-t[inter_vt]))*(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)+ score_inter

              c_bry_new=e_bry_new/sqrt(T_new)

              wald_last=theta_hat_last/ theta_hat_last_sd
              score_last=theta_hat_last*t_last

              #################################################################################
              ###                   Step 4.6: Hypothesis testing                           ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              test_add=wald_last-c_bry_new
              samsz_simu=new_sampsz
              if (test_add>0){
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
      out=c(rej_nossr,rej,stg,samsz_simu,futi)

  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################

  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>evt_num){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################

  out=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,evt_num,evt_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  return(out)

}

####################
#asd_rej_simu_adapt_OC_survival(HR,lambda,s,c_bry,surv_rate,evt_num,evt_max,HR_cut,samp_multip,sim_num,beta)

#asd_rej_simu_adapt_OC_survival(HR=0.8,lambda=4,s=c(1/3,2/3,1),c_bry=c(3.471095, 2.454434,2.004037),surv_rate=0.5,evt_num=845,evt_max=2000,HR_cut=0.9,samp_multip=5,sim_num=100,beta=0.1)
#asd_rej_simu_adapt_OC_survival(HR=0.7,lambda=1,s=c(1/3,2/3,1),c_bry=c(3.471095, 2.454434,2.004037),surv_rate=0.5,evt_num=600,evt_max=1000,HR_cut=0.95,samp_multip=5,sim_num=1000,beta=0.1)

#asd_rej_simu_adapt_com_survival(HR=1,lambda=4,s=c(1/3,2/3,1),c_bry=c(3.471095, 2.454434,2.004037),surv_rate=0.5,evt_num=845,evt_max=2000,HR_cut=0.9,samp_multip=5,sim_num=10000,beta=0.1)
##################################################################################################

##########################

#' Group sequential designs / Adaptive sequential designs: Simulations, Survival analysis, Optimizing adaptive sequential design
#'
#' @param HR Hazards ratio of the test arm vs. the control arm; HR < 1 favours the test arm. The
#'   simulation grid varies HR over the values supplied
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g.
#'   c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g.
#'   c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param alpha_2_alps_s1 Two-sided type I error spent at the interim analyses of schedule s1 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error spent at the interim analyses of schedule s2 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param surv_rate Survival rate at the end of the trial, used together with the number of events
#'   to derive the sample size; must be between 0 and 1
#' @param evt_num_1 Planned number of events for the first scenario of the simulation grid
#' @param evt_num_2 Planned number of events for the second scenario of the simulation grid
#' @param evt_max_1 Maximum number of events allowed for the sample size change in the first
#'   scenario of the simulation grid
#' @param evt_max_2 Maximum number of events allowed for the sample size change in the second
#'   scenario of the simulation grid
#' @param HR_cut Lack of efficacy threshold on the hazards-ratio scale; the corresponding
#'   threshold on the log-hazards-ratio scale is theta_cut = -log(HR_cut). If the observed interim
#'   log-hazards-ratio is below theta_cut (i.e. the observed HR is greater than HR_cut), the trial
#'   continues to the planned number of events with no sample size increase and no early
#'   termination for futility
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with six components, one for each combination of critical boundary type
#' (O'Brien-Fleming, Pocock, alpha-spending) and interim analysis schedule (s1, s2):
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{Named list with the group
#'     sequential design (Looks, Information fractions, O'Brien-Fleming boundaries) and the
#'     simulation results for the four combinations of planned number of events (evt_num_1,
#'     evt_num_2) and maximum number of events (evt_max_1, evt_max_2), each a data frame with one
#'     row per hazards ratio and the columns HR<sub>1</sub>,..., sample size for fixed design,
#'     Rejection rate without SSR, Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned sample size
#'     (the planned number of events), Maximum sample size (the maximum number of events allowed),
#'     Mean sample size (the mean simulated number of events),
#'     mean sample size conditional on increase and
#'     Ratio of adaptive sample size over fixed sample size}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with the information fractions s2}
#'   \item{Pocock boundaries, interim analysis schedule=s1}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s1}
#'   \item{Pocock boundaries, interim analysis schedule=s2}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s2}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s1}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s2}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' asd_simu_OC_survival(HR = c(0.75,0.8,0.85), s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05,
#' alpha_2_alps_s1 = c(0.001,0.001,0.048), alpha_2_alps_s2 = c(0.001,0.049), surv_rate = 0.5,
#' evt_num_1 = 200, evt_num_2 = 300, evt_max_1 = 2000, evt_max_2 = 3000, HR_cut = 0.95, beta = 0.1,
#' sim_num = 20)
#' }
#'
#' @export
asd_simu_OC_survival<-function(HR,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,surv_rate,evt_num_1,evt_num_2,evt_max_1,evt_max_2,HR_cut,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)

  mu=abs(log(HR))


  ###### sample size calculations for fixed design ####

  sampsz_HR=NULL

  for (i in 1:length(HR)){
    sampsz_HR[i]=ceiling(((a+b)/(0.5*mu[i]))**2)
  }

  ############ critical boundary calculations #####

  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(OF_bundry(s1,alpha)[1]))
  s2_OF_cbry=as.numeric(unlist(OF_bundry(s2,alpha)[1]))

  s1_pocock_cbry=as.numeric(unlist(pocock_bundry(s1,alpha)[1]))
  s2_pocock_cbry=as.numeric(unlist(pocock_bundry(s2,alpha)[1]))

  s1_alsp_cbry=as.numeric(unlist(alsp_bundry(s1,alpha_alsp_s1)[1]))
  s2_alsp_cbry=as.numeric(unlist(alsp_bundry(s2,alpha_alsp_s2)[1]))


  ##### GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_pocock_design=data.frame(cbind(c(1:length(s1)),s1,s1_pocock_cbry))
  names(s1_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s2_pocock_design=data.frame(cbind(c(1:length(s2)),s2,s2_pocock_cbry))
  names(s2_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")


  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## evt_max_1 ###
  ##############

  ## evt_num_1###

  ##### s1 ####
  simu_result_title = c("HR<sub>1</sub>,...","sample size for fixed design","Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                        "Rejection rate without SSR/(meeting sample size incerase criteria)",
                        "Rejection rate with SSR/(sample size incerase)",
                        "Futility rate","Planned sample size","Maximum sample size",
                        "Mean sample size","mean sample size conditional on increase","Ratio of adaptive sample size over fixed sample size")
  simu_out_OF_s1_EVTnum_1_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_OF_s1_EVTnum_1_EVTmax_1= rbind(simu_out_OF_s1_EVTnum_1_EVTmax_1,
                                            asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_OF_cbry,
                                                                           surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }



  samsz_ratio_OF_s1_EVTnum_1_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_OF_s1_EVTnum_1_EVTmax_1[i]=as.numeric(simu_out_OF_s1_EVTnum_1_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_OF_s1_EVTnum_1_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_OF_s1_EVTnum_1_EVTmax_1,samsz_ratio_OF_s1_EVTnum_1_EVTmax_1)))
  names(simu_out_OF_s1_EVTnum_1_EVTmax_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_EVTnum_1_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_OF_s2_EVTnum_1_EVTmax_1= rbind(simu_out_OF_s2_EVTnum_1_EVTmax_1,
                                            asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_OF_cbry,
                                                                           surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_EVTnum_1_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_OF_s2_EVTnum_1_EVTmax_1[i]=as.numeric(simu_out_OF_s2_EVTnum_1_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_OF_s2_EVTnum_1_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_OF_s2_EVTnum_1_EVTmax_1,samsz_ratio_OF_s2_EVTnum_1_EVTmax_1)))
  names(simu_out_OF_s2_EVTnum_1_EVTmax_1)=simu_result_title

  ## evt_num_2###


  ##### s1 ####
  simu_out_OF_s1_EVTnum_2_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_OF_s1_EVTnum_2_EVTmax_1= rbind(simu_out_OF_s1_EVTnum_2_EVTmax_1,
                                            asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_OF_cbry,
                                                                           surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_EVTnum_2_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_OF_s1_EVTnum_2_EVTmax_1[i]=as.numeric(simu_out_OF_s1_EVTnum_2_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_OF_s1_EVTnum_2_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_OF_s1_EVTnum_2_EVTmax_1,samsz_ratio_OF_s1_EVTnum_2_EVTmax_1)))
  names(simu_out_OF_s1_EVTnum_2_EVTmax_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_EVTnum_2_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_OF_s2_EVTnum_2_EVTmax_1= rbind(simu_out_OF_s2_EVTnum_2_EVTmax_1,
                                            asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_OF_cbry,
                                                                           surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_EVTnum_2_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_OF_s2_EVTnum_2_EVTmax_1[i]=as.numeric(simu_out_OF_s2_EVTnum_2_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_OF_s2_EVTnum_2_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_OF_s2_EVTnum_2_EVTmax_1,samsz_ratio_OF_s2_EVTnum_2_EVTmax_1)))
  names(simu_out_OF_s2_EVTnum_2_EVTmax_1)=simu_result_title

  ##############
  ## evt_max_2 ###
  ##############


  ## evt_num_1###

  ##### s1 ####
  simu_out_OF_s1_EVTnum_1_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_OF_s1_EVTnum_1_EVTmax_2= rbind(simu_out_OF_s1_EVTnum_1_EVTmax_2,
                                            asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_OF_cbry,
                                                                           surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_EVTnum_1_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_OF_s1_EVTnum_1_EVTmax_2[i]=as.numeric(simu_out_OF_s1_EVTnum_1_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_OF_s1_EVTnum_1_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_OF_s1_EVTnum_1_EVTmax_2,samsz_ratio_OF_s1_EVTnum_1_EVTmax_2)))
  names(simu_out_OF_s1_EVTnum_1_EVTmax_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_EVTnum_1_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_OF_s2_EVTnum_1_EVTmax_2= rbind(simu_out_OF_s2_EVTnum_1_EVTmax_2,
                                            asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_OF_cbry,
                                                                           surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_EVTnum_1_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_OF_s2_EVTnum_1_EVTmax_2[i]=as.numeric(simu_out_OF_s2_EVTnum_1_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_OF_s2_EVTnum_1_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_OF_s2_EVTnum_1_EVTmax_2,samsz_ratio_OF_s2_EVTnum_1_EVTmax_2)))
  names(simu_out_OF_s2_EVTnum_1_EVTmax_2)=simu_result_title

  ## evt_num_2###


  ##### s1 ####
  simu_out_OF_s1_EVTnum_2_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_OF_s1_EVTnum_2_EVTmax_2= rbind(simu_out_OF_s1_EVTnum_2_EVTmax_2,
                                            asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_OF_cbry,
                                                                           surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_EVTnum_2_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_OF_s1_EVTnum_2_EVTmax_2[i]=as.numeric(simu_out_OF_s1_EVTnum_2_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_OF_s1_EVTnum_2_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_OF_s1_EVTnum_2_EVTmax_2,samsz_ratio_OF_s1_EVTnum_2_EVTmax_2)))
  names(simu_out_OF_s1_EVTnum_2_EVTmax_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_EVTnum_2_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_OF_s2_EVTnum_2_EVTmax_2= rbind(simu_out_OF_s2_EVTnum_2_EVTmax_2,
                                            asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_OF_cbry,
                                                                           surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_EVTnum_2_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_OF_s2_EVTnum_2_EVTmax_2[i]=as.numeric(simu_out_OF_s2_EVTnum_2_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_OF_s2_EVTnum_2_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_OF_s2_EVTnum_2_EVTmax_2,samsz_ratio_OF_s2_EVTnum_2_EVTmax_2)))
  names(simu_out_OF_s2_EVTnum_2_EVTmax_2)=simu_result_title

  ##############
  #############################################
  #### OC simulations for Pocock boundaries ###
  #############################################


  ##############
  ## evt_max_1 ###
  ##############

  ## evt_num_1###

  ##### s1 ####
  simu_out_pocock_s1_EVTnum_1_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_pocock_s1_EVTnum_1_EVTmax_1= rbind(simu_out_pocock_s1_EVTnum_1_EVTmax_1,
                                                asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_pocock_cbry,
                                                                               surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_EVTnum_1_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_pocock_s1_EVTnum_1_EVTmax_1[i]=as.numeric(simu_out_pocock_s1_EVTnum_1_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_pocock_s1_EVTnum_1_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_pocock_s1_EVTnum_1_EVTmax_1,samsz_ratio_pocock_s1_EVTnum_1_EVTmax_1)))
  names(simu_out_pocock_s1_EVTnum_1_EVTmax_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_EVTnum_1_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_pocock_s2_EVTnum_1_EVTmax_1= rbind(simu_out_pocock_s2_EVTnum_1_EVTmax_1,
                                                asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_pocock_cbry,
                                                                               surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_EVTnum_1_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_pocock_s2_EVTnum_1_EVTmax_1[i]=as.numeric(simu_out_pocock_s2_EVTnum_1_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_pocock_s2_EVTnum_1_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_pocock_s2_EVTnum_1_EVTmax_1,samsz_ratio_pocock_s2_EVTnum_1_EVTmax_1)))
  names(simu_out_pocock_s2_EVTnum_1_EVTmax_1)=simu_result_title

  ## evt_num_2###


  ##### s1 ####
  simu_out_pocock_s1_EVTnum_2_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_pocock_s1_EVTnum_2_EVTmax_1= rbind(simu_out_pocock_s1_EVTnum_2_EVTmax_1,
                                                asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_pocock_cbry,
                                                                               surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_EVTnum_2_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_pocock_s1_EVTnum_2_EVTmax_1[i]=as.numeric(simu_out_pocock_s1_EVTnum_2_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_pocock_s1_EVTnum_2_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_pocock_s1_EVTnum_2_EVTmax_1,samsz_ratio_pocock_s1_EVTnum_2_EVTmax_1)))
  names(simu_out_pocock_s1_EVTnum_2_EVTmax_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_EVTnum_2_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_pocock_s2_EVTnum_2_EVTmax_1= rbind(simu_out_pocock_s2_EVTnum_2_EVTmax_1,
                                                asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_pocock_cbry,
                                                                               surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_EVTnum_2_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_pocock_s2_EVTnum_2_EVTmax_1[i]=as.numeric(simu_out_pocock_s2_EVTnum_2_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_pocock_s2_EVTnum_2_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_pocock_s2_EVTnum_2_EVTmax_1,samsz_ratio_pocock_s2_EVTnum_2_EVTmax_1)))
  names(simu_out_pocock_s2_EVTnum_2_EVTmax_1)=simu_result_title

  ##############
  ## evt_max_2 ###
  ##############


  ## evt_num_1###

  ##### s1 ####
  simu_out_pocock_s1_EVTnum_1_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_pocock_s1_EVTnum_1_EVTmax_2= rbind(simu_out_pocock_s1_EVTnum_1_EVTmax_2,
                                                asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_pocock_cbry,
                                                                               surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_EVTnum_1_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_pocock_s1_EVTnum_1_EVTmax_2[i]=as.numeric(simu_out_pocock_s1_EVTnum_1_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_pocock_s1_EVTnum_1_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_pocock_s1_EVTnum_1_EVTmax_2,samsz_ratio_pocock_s1_EVTnum_1_EVTmax_2)))
  names(simu_out_pocock_s1_EVTnum_1_EVTmax_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_EVTnum_1_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_pocock_s2_EVTnum_1_EVTmax_2= rbind(simu_out_pocock_s2_EVTnum_1_EVTmax_2,
                                                asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_pocock_cbry,
                                                                               surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num))
  }

  samsz_ratio_pocock_s2_EVTnum_1_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_pocock_s2_EVTnum_1_EVTmax_2[i]=as.numeric(simu_out_pocock_s2_EVTnum_1_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_pocock_s2_EVTnum_1_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_pocock_s2_EVTnum_1_EVTmax_2,samsz_ratio_pocock_s2_EVTnum_1_EVTmax_2)))
  names(simu_out_pocock_s2_EVTnum_1_EVTmax_2)=simu_result_title

  ## evt_num_2###


  ##### s1 ####
  simu_out_pocock_s1_EVTnum_2_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_pocock_s1_EVTnum_2_EVTmax_2= rbind(simu_out_pocock_s1_EVTnum_2_EVTmax_2,
                                                asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_pocock_cbry,
                                                                               surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_EVTnum_2_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_pocock_s1_EVTnum_2_EVTmax_2[i]=as.numeric(simu_out_pocock_s1_EVTnum_2_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_pocock_s1_EVTnum_2_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_pocock_s1_EVTnum_2_EVTmax_2,samsz_ratio_pocock_s1_EVTnum_2_EVTmax_2)))
  names(simu_out_pocock_s1_EVTnum_2_EVTmax_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_EVTnum_2_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_pocock_s2_EVTnum_2_EVTmax_2= rbind(simu_out_pocock_s2_EVTnum_2_EVTmax_2,
                                                asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_pocock_cbry,
                                                                               surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_EVTnum_2_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_pocock_s2_EVTnum_2_EVTmax_2[i]=as.numeric(simu_out_pocock_s2_EVTnum_2_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_pocock_s2_EVTnum_2_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_pocock_s2_EVTnum_2_EVTmax_2,samsz_ratio_pocock_s2_EVTnum_2_EVTmax_2)))
  names(simu_out_pocock_s2_EVTnum_2_EVTmax_2)=simu_result_title

  ##############

  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################

  ##############
  ## evt_max_1 ###
  ##############

  ## evt_num_1###

  ##### s1 ####
  simu_out_alsp_s1_EVTnum_1_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_alsp_s1_EVTnum_1_EVTmax_1= rbind(simu_out_alsp_s1_EVTnum_1_EVTmax_1,
                                              asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_alsp_cbry,
                                                                             surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_EVTnum_1_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_alsp_s1_EVTnum_1_EVTmax_1[i]=as.numeric(simu_out_alsp_s1_EVTnum_1_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_alsp_s1_EVTnum_1_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_alsp_s1_EVTnum_1_EVTmax_1,samsz_ratio_alsp_s1_EVTnum_1_EVTmax_1)))
  names(simu_out_alsp_s1_EVTnum_1_EVTmax_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_1_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_alsp_s2_EVTnum_1_EVTmax_1= rbind(simu_out_alsp_s2_EVTnum_1_EVTmax_1,
                                              asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_alsp_cbry,
                                                                             surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_EVTnum_1_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_alsp_s2_EVTnum_1_EVTmax_1[i]=as.numeric(simu_out_alsp_s2_EVTnum_1_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_alsp_s2_EVTnum_1_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_alsp_s2_EVTnum_1_EVTmax_1,samsz_ratio_alsp_s2_EVTnum_1_EVTmax_1)))
  names(simu_out_alsp_s2_EVTnum_1_EVTmax_1)=simu_result_title

  ## evt_num_2###


  ##### s1 ####
  simu_out_alsp_s1_EVTnum_2_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_alsp_s1_EVTnum_2_EVTmax_1= rbind(simu_out_alsp_s1_EVTnum_2_EVTmax_1,
                                              asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_alsp_cbry,
                                                                             surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_EVTnum_2_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_alsp_s1_EVTnum_2_EVTmax_1[i]=as.numeric(simu_out_alsp_s1_EVTnum_2_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_alsp_s1_EVTnum_2_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_alsp_s1_EVTnum_2_EVTmax_1,samsz_ratio_alsp_s1_EVTnum_2_EVTmax_1)))
  names(simu_out_alsp_s1_EVTnum_2_EVTmax_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_2_EVTmax_1=NULL
  for (i in 1:length(HR)){
    simu_out_alsp_s2_EVTnum_2_EVTmax_1= rbind(simu_out_alsp_s2_EVTnum_2_EVTmax_1,
                                              asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_alsp_cbry,
                                                                             surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_1,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_EVTnum_2_EVTmax_1=NULL
  for (i in (1:length(HR))){
    samsz_ratio_alsp_s2_EVTnum_2_EVTmax_1[i]=as.numeric(simu_out_alsp_s2_EVTnum_2_EVTmax_1[i,9])/sampsz_HR[i]
  }


  simu_out_alsp_s2_EVTnum_2_EVTmax_1=data.frame((cbind(HR,sampsz_HR,simu_out_alsp_s2_EVTnum_2_EVTmax_1,samsz_ratio_alsp_s2_EVTnum_2_EVTmax_1)))
  names(simu_out_alsp_s2_EVTnum_2_EVTmax_1)=simu_result_title
  ##############
  ## evt_max_2 ###
  ##############


  ## evt_num_1###

  ##### s1 ####
  simu_out_alsp_s1_EVTnum_1_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_alsp_s1_EVTnum_1_EVTmax_2= rbind(simu_out_alsp_s1_EVTnum_1_EVTmax_2,
                                              asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_alsp_cbry,
                                                                             surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_EVTnum_1_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_alsp_s1_EVTnum_1_EVTmax_2[i]=as.numeric(simu_out_alsp_s1_EVTnum_1_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_alsp_s1_EVTnum_1_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_alsp_s1_EVTnum_1_EVTmax_2,samsz_ratio_alsp_s1_EVTnum_1_EVTmax_2)))
  names(simu_out_alsp_s1_EVTnum_1_EVTmax_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_1_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_alsp_s2_EVTnum_1_EVTmax_2= rbind(simu_out_alsp_s2_EVTnum_1_EVTmax_2,
                                              asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_alsp_cbry,
                                                                             surv_rate=surv_rate,evt_num=evt_num_1,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_EVTnum_1_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_alsp_s2_EVTnum_1_EVTmax_2[i]=as.numeric(simu_out_alsp_s2_EVTnum_1_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_alsp_s2_EVTnum_1_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_alsp_s2_EVTnum_1_EVTmax_2,samsz_ratio_alsp_s2_EVTnum_1_EVTmax_2)))
  names(simu_out_alsp_s2_EVTnum_1_EVTmax_2)=simu_result_title

  ## evt_num_2###


  ##### s1 ####
  simu_out_alsp_s1_EVTnum_2_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_alsp_s1_EVTnum_2_EVTmax_2= rbind(simu_out_alsp_s1_EVTnum_2_EVTmax_2,
                                              asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s1,c_bry=s1_alsp_cbry,
                                                                             surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_EVTnum_2_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_alsp_s1_EVTnum_2_EVTmax_2[i]=as.numeric(simu_out_alsp_s1_EVTnum_2_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_alsp_s1_EVTnum_2_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_alsp_s1_EVTnum_2_EVTmax_2,samsz_ratio_alsp_s1_EVTnum_2_EVTmax_2)))
  names(simu_out_alsp_s1_EVTnum_2_EVTmax_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_EVTnum_2_EVTmax_2=NULL
  for (i in 1:length(HR)){
    simu_out_alsp_s2_EVTnum_2_EVTmax_2= rbind(simu_out_alsp_s2_EVTnum_2_EVTmax_2,
                                              asd_rej_simu_adapt_OC_survival(HR=HR[i],lambda=1,s=s2,c_bry=s2_alsp_cbry,
                                                                             surv_rate=surv_rate,evt_num=evt_num_2,evt_max=evt_max_2,HR_cut=HR_cut,samp_multip=5,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_EVTnum_2_EVTmax_2=NULL
  for (i in (1:length(HR))){
    samsz_ratio_alsp_s2_EVTnum_2_EVTmax_2[i]=as.numeric(simu_out_alsp_s2_EVTnum_2_EVTmax_2[i,9])/sampsz_HR[i]
  }


  simu_out_alsp_s2_EVTnum_2_EVTmax_2=data.frame((cbind(HR,sampsz_HR,simu_out_alsp_s2_EVTnum_2_EVTmax_2,samsz_ratio_alsp_s2_EVTnum_2_EVTmax_2)))
  names(simu_out_alsp_s2_EVTnum_2_EVTmax_2)=simu_result_title

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

  list_pocock_1=list(s1_pocock_design,
                     simu_out_pocock_s1_EVTnum_1_EVTmax_1,
                     simu_out_pocock_s1_EVTnum_2_EVTmax_1,
                     simu_out_pocock_s1_EVTnum_1_EVTmax_2,
                     simu_out_pocock_s1_EVTnum_2_EVTmax_2
                     )

  names(list_pocock_1)=c("Group sequential design",
                     "Simulation results: Information fractions 1, Sample size 1, evt_max 1",
                     "Simulation results: Information fractions 1, Sample size 2, evt_max 1",
                     "Simulation results: Information fractions 1, Sample size 1, evt_max 2",
                     "Simulation results: Information fractions 1, Sample size 2, evt_max 2")

  list_pocock_2=list(s2_pocock_design,
                     simu_out_pocock_s2_EVTnum_1_EVTmax_1,
                     simu_out_pocock_s2_EVTnum_2_EVTmax_1,
                     simu_out_pocock_s2_EVTnum_1_EVTmax_2,
                     simu_out_pocock_s2_EVTnum_2_EVTmax_2
  )

  names(list_pocock_2)=c("Group sequential design",
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
           list_pocock_1,list_pocock_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Pocock boundaries, interim analysis schedule=s1",
               "Pocock boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}

##################################
#asd_simu_OC_survival(HR=c(0.7,0.78,0.8),s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                   surv_rate=0.5,evt_num_1=600,evt_num_2=800,evt_max_1=1000,evt_max_2=1500,HR_cut=0.95,beta=0.1,sim_num=100)
################################################


#########################################################
#### 3.   simulations for binary distributions   ########
#########################################################

##############################################################
#### 3.1  GSD simulations for binary distributions    ########
##############################################################

################################################
gsd_rej_rate_binary<-function(p_t,p_c,direction,s,c_bry,samsz,sim_num){
  #This function simulates the rejection rate, given p_t, p_c,s,c_bry,sim_num
  ## for a group sequential design
  ## the random samples are from normal distributions
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is the sample size

  ######################################
  ###  Preparing the parameters ########
  ######################################
  N=ceiling(samsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted

  t<- N/(p_t*(1-p_t)+p_c*(1-p_c))  # t is the vector of information times at each interim anbalysis
  NF=N[length(s)]  ## NF is the sample size at the end of the trial
  #t is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {


      ########################################
      ### Step 1: Generate random samples ####
      ########################################
      x0=stats::rbinom(NF,1,p_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rbinom(NF,1,p_t)

      #####################################################################
      ### Step 2: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      Wald=NULL
      theta=NULL

      if (direction==1){ #higher rate p is better
      for (j in 1:length(s) ){

        theta[j]=mean(x[1:N[j]])-mean(x0[1:N[j]])
        Wald[j]=theta[j]*sqrt(t[j])
      }
      }else if (direction==2){ #lower rate p is better
        for (j in 1:length(s) ){

          theta[j]=mean(x0[1:N[j]])-mean(x[1:N[j]])
          Wald[j]=theta[j]*sqrt(t[j])

        }
      }
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry
      #at each simulation, test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if ( max(test)>0){

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        #stg identifies the interim analysis at which the null hypothesis is rejected
      }else{
        rej=0
        stg=0
      }
      if (stg>0){
        last_vt=stg
      } else {last_vt=length(s)
      }
      #######################################################
      ### Step 4: output and finish the simjulation loop ####
      #######################################################
      samsz_simu=N[last_vt]
      out=c(rej,stg,samsz_simu)
  })#(end of q loop)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  samsz_simu = as.numeric(simr[,3])
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  rej_num=sum(rej>0)

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  mean_samsz=mean(samsz_simu)

  #############################
  ### Step 6: Final output ####
  #############################
  out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(rej_time)))
  names(out_1)=c("Looks","Information fraction","Critical boundaries","Rejection by look")
  out_2=c(samsz,rej_rt,mean_samsz)
  names(out_2)=c("Planned control sample size","Rejection rate","Mean control sample size")
  out=list(out_1,out_2)
  names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)
}
#######################################

#gsd_rej_rate_binary(p_t,p_c,direction,s,c_bry,samsz,sim_num)


#gsd_rej_rate_binary(p_t=0.7,p_c=0.2,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,sim_num=100)
##############################################

#' Group sequential designs / Adaptive sequential designs: Simulations, Binary distribution, Type I error and power
#'
#' @param p_t Event rate of the test arm; must be between 0 and 1
#' @param p_c Event rate of the control arm; must be between 0 and 1
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param samsz Planned sample size of the control arm used in the simulation; must be > 0
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries and Rejection by look (the
#'     rejection rate at each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Planned control sample size (the
#'     sample size used in the simulation), Rejection rate (the simulated power, i.e. the type I
#'     error rate when p_t = p_c) and Mean control sample size (the mean simulated sample size)}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' gsd_power_simu_binary(p_t = 0.4, p_c = 0.4, direction = 1, s = c(1/3,2/3,1), bry_type = 1,
#' alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
gsd_power_simu_binary<-function(p_t,p_c,direction,s,bry_type,alpha_2,samsz,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  simu_out=gsd_rej_rate_binary(p_t=p_t,p_c=p_c,direction=direction, s=s,c_bry=c_bry,samsz=samsz,sim_num=sim_num)

  # out_1=simu_out[[1]]
 #  out_2=c(samsz,as.numeric(unlist(simu_out[[2]])))
#   names(out_2)=c("Planned sample size","rejection rate","Mean sample size")
   out=simu_out
  # names(out)=c("Sequential design","Sample size and rejection rate")
   return(out)
}

##################################
#gsd_power_simu_binary(p_t,p_c,direction,s,bry_type,alpha_2,sample_size,sim_num)
#gsd_power_simu_binary(p_t=0.4,p_c=0.4,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000000)
#gsd_power_simu_binary(p_t=0.4,p_c=0.2,direction=1,s=c(1/4,1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=107,sim_num=1000)


###########################################################################
#### 3.2  Adaptive sequential simulations for binary distributions  #######
###########################################################################
##############
asd_rej_simu_adapt_com_binary<-function(p_t,p_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num){
  # p_t is the event rate for the test arm
  # p_c is the evnet rate for the control arm
  #s is the vector of information frations
  #c_bry is the vector of critical values
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size
  #N_max is themaximum sample size

  ######################################
  ###  Preparing the parameters ########
  ######################################
  N=ceiling(samsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted

 t<- N/(p_t*(1-p_t)+p_c*(1-p_c))  # t is the vector of information times at each interim anbalysis
  NF=N[length(s)]  ## NF is the sample size at the end of the trial
  #t is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      ########################################
      ### Step 1: Generate random samples ####
      ########################################
      x0=stats::rbinom(N_max,1,p_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rbinom(N_max,1,p_t)
      #############################################################################################
      ### Step 2: Obtain the Wald statistics and theta_hat for hypothesis testing and estimate ####
      #############################################################################################

      Wald=NULL
      theta_hat=NULL

      if (direction==1){ #higher rate p is better
        for (j in 1:length(s) ){

          theta_hat[j]=mean(x[1:N[j]])-mean(x0[1:N[j]])
          Wald[j]=theta_hat[j]*sqrt(t[j])
        }
      }else if (direction==2){ #lower rate p is better
        for (j in 1:length(s) ){

          theta_hat[j]=mean(x0[1:N[j]])-mean(x[1:N[j]])
          Wald[j]=theta_hat[j]*sqrt(t[j])

        }
      }

      theta_hat_sd=1/sqrt(t)
      #Wald containes the Wald statistics

      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd=theta_hat_sd[length(s)-1]
      if (theta_hat_inter<theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry

      test_1=test[1:length(s)-1]

      ##################

      if ( max(test)>=0){
        rej_nossr=1
      }else if ( max(test)<0){
        rej_nossr=0
      }
      ##############
      # test compares the Wald statistics with the critical boundaries
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        #  last_vt=2
        wald_last=Wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        sampsz_simu=N[stg]

        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else {

        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation
          case=2
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=1/sqrt(t[last_vt])


         if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }
          sampsz_simu=N[length(s)]

          #############  #############

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          score_inter=theta_hat_inter*t[inter_vt]

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))

        #  cond_power=asd_c_power_theta(theta=theta_hat_inter,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)

          if (cond_power>=1-beta){
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=1/sqrt(t[last_vt])


            sampsz_simu=N[length(s)]
            if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else{

          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_sampsz=ceiling(T_new_0*(p_t*(1-p_t)+p_c*(1-p_c)))
            #        T_new=new_sampsz/2/sigma^2
            T_new=new_sampsz/(p_t*(1-p_t)+p_c*(1-p_c))
            e_bry_new=
              (sqrt(T_new-t[inter_vt])/sqrt(t[length(s)]-t[inter_vt]))*(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)+ score_inter

            c_bry_new=e_bry_new/sqrt(T_new)


         # T_new_0=new_info_time_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,beta)
   #       new_sampsz=ceiling(T_new_0*2*sigma^2)

         # new_c_bry=new_c_bry_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,T_new)
          if  (new_sampsz>N_max){
            #################################################################################
            ###                   Step 4.3: Compare theta_hat with theta_cut             ####
            ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
            ########                   and calculate confidenc einterval                 ####
            #################################################################################

            if (last(test)>0){
              rej=1
              stg=length(s)
            }else{
              rej=0
              stg=0
            }
            sampsz_simu=last(N)
          } else {
            #################################################################################
            ###                   Step 4.4: Compare theta_hat with theta_cut             ####
            ###                     if new_sampsz<=N_max, proceed to new sample size     ####
            ########                   and calculate confidenc einterval                 ####
            #################################################################################
            sampsz_simu=new_sampsz
            ##### Generate addtional random samples #######
            inter_vt=length(s)-1
            #################################################################################
            ###                   Step 4.5: generate new random samples                  ####
            ###                         and calculate statistics                         ####
            #################################################################################

            x0_new=x0[1:new_sampsz]
            x_new= x[1:new_sampsz]

            if (direction==1){ #higher rate p is better
              theta_hat_last=mean(x_new)-mean(x0_new)
            }else if (direction==2){ #hlower rate p is better
              theta_hat_last=mean(x0_new)-mean(x_new)
            }

            score_last=theta_hat_last*T_new
            wald_last=theta_hat_last*sqrt(T_new)

            #################################################################################
            ###                   Step 4.6: Hypothesis testing                           ####
            ########                   and calculate confidenc einterval                 ####
            #################################################################################
            test_add=wald_last-c_bry_new

            if (test_add>0){
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
      out=c(rej_nossr,rej,stg,sampsz_simu,futi)
  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################
  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
    #################
  }
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>samsz){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  sampsz_icr_rt_byLook=c(rep(0,(length(s)-1)),mean(sampsz_icr))
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################


  out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(sampsz_icr_rt_byLook),as.numeric(rej_time)))
  names(out_1)=c("Looks","Information fraction","Critical boundaries","Sample size increase","Rejection by look")
  out_2=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,samsz,N_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out_2)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  out=list(out_1,out_2)
  names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)

}
#######################################
#asd_rej_simu_adapt_com_binary(p_t,p_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num)

#asd_rej_simu_adapt_com_binary(p_t=0.6,p_c=0.6,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_rej_simu_adapt_com_binary(p_t=0.8,p_c=0.6,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)

########################################################################

#' Group sequential designs / Adaptive sequential designs: Simulations, Binary distribution, Operating characteristics
#'
#' @param p_t Event rate of the test arm; must be between 0 and 1
#' @param p_c Event rate of the control arm; must be between 0 and 1
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param samsz Planned sample size of the control arm; the trial starts from this sample size and
#'   it is increased only if the conditional power at the interim analysis is below the target
#' @param N_max Maximum sample size of the control arm allowed for the sample size change; the
#'   re-estimated sample size is truncated at N_max
#' @param theta_cut Lack of efficacy threshold on the effect-size scale: if the observed interim
#'   effect size is below theta_cut, the trial continues to the planned sample size with no sample
#'   size increase and no early termination for futility
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries, Sample size increase (the rate of
#'     sample size increase attributed to each look) and Rejection by look (the rejection rate at
#'     each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Rejection rate without SSR,
#'     Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size,
#'     Maximum control sample size, Mean control sample size and
#'     mean control sample size conditional on increase}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_power_simu_binary(p_t = 0.4, p_c = 0.2, direction = 1, s = c(1/3,2/3,1), bry_type = 1,
#' alpha_2 = 0.05, samsz = 107, N_max = 500, theta_cut = 0.05, beta = 0.1, sim_num = 20)
#'
#' @export
asd_power_simu_binary<-function(p_t,p_c,direction,s,bry_type,alpha_2,samsz,N_max,theta_cut,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }


  c_bry=as.numeric(unlist(boundary_out[1]))

  # samsz=control_sample_size
  simu_out=asd_rej_simu_adapt_com_binary(p_t=p_t,p_c=p_c,direction=direction,s=s,c_bry=c_bry,samsz=samsz,N_max=N_max,theta_cut=theta_cut,beta=beta,sim_num=sim_num)

  out= simu_out


  return(out)
}

##################################
#asd_power_simu_binary(p_t,p_c,direction,s,bry_type,alpha_2,samsz,N_max,theta_cut,beta,sim_num)
#asd_power_simu_binary(p_t=0.4,p_c=0.2,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=107,N_max=500,theta_cut=0.05,beta=0.1,sim_num=1000)
#asd_power_simu_binary(p_t=0.35,p_c=0.2,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=107,N_max=500,theta_cut=0.05,beta=0.1,sim_num=10000)

#
#####################################################

##################
asd_rej_simu_adapt_OC_binary<-function(p_t,p_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num){
  # p_t is the event rate for the test arm
  # p_c is the evnet rate for the control arm
  #s is the vector of information frations
  #c_bry is the vector of critical values
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is sample size
  #N_max is themaximum sample size

  ######################################
  ###  Preparing the parameters ########
  ######################################
  N=ceiling(samsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted

  t<- N/(p_t*(1-p_t)+p_c*(1-p_c))  # t is the vector of information times at each interim anbalysis
  NF=N[length(s)]  ## NF is the sample size at the end of the trial
  #t is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      ########################################
      ### Step 1: Generate random samples ####
      ########################################
      x0=stats::rbinom(N_max,1,p_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rbinom(N_max,1,p_t)
      #############################################################################################
      ### Step 2: Obtain the Wald statistics and theta_hat for hypothesis testing and estimate ####
      #############################################################################################

      Wald=NULL
      theta_hat=NULL

      if (direction==1){ #higher rate p is better
        for (j in 1:length(s) ){

          theta_hat[j]=mean(x[1:N[j]])-mean(x0[1:N[j]])
          Wald[j]=theta_hat[j]*sqrt(t[j])
        }
      }else if (direction==2){ #lower rate p is better
        for (j in 1:length(s) ){

          theta_hat[j]=mean(x0[1:N[j]])-mean(x[1:N[j]])
          Wald[j]=theta_hat[j]*sqrt(t[j])

        }
      }

      theta_hat_sd=1/sqrt(t)
      #Wald containes the Wald statistics

      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd=theta_hat_sd[length(s)-1]
      if (theta_hat_inter<theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry

      test_1=test[1:length(s)-1]

      ##################

      if ( max(test)>=0){
        rej_nossr=1
      }else if ( max(test)<0){
        rej_nossr=0
      }
      ##############
      # test compares the Wald statistics with the critical boundaries
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        #  last_vt=2
        wald_last=Wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        sampsz_simu=N[stg]

        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else {

        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation
          case=2
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=1/sqrt(t[last_vt])


          if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }
          sampsz_simu=N[length(s)]

          #############  #############

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          score_inter=theta_hat_inter*t[inter_vt]

          cond_power=stats::pnorm(theta_hat_inter*sqrt(t[length(s)]-t[inter_vt])-(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt]))

          #  cond_power=asd_c_power_theta(theta=theta_hat_inter,s,c_bry,inter_vt,theta_hat_inter,theta_hat_inter_sd)

          if (cond_power>=1-beta){
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=1/sqrt(t[last_vt])


            sampsz_simu=N[length(s)]
            if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else{

            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##
            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(t[length(s)])-score_inter)/sqrt(t[length(s)]-t[inter_vt])+stats::qnorm(1-beta))^2+t[inter_vt]
            new_sampsz=ceiling(T_new_0*(p_t*(1-p_t)+p_c*(1-p_c)))
            #        T_new=new_sampsz/2/sigma^2
            T_new=new_sampsz/(p_t*(1-p_t)+p_c*(1-p_c))
            e_bry_new=
              (sqrt(T_new-t[inter_vt])/sqrt(t[length(s)]-t[inter_vt]))*(c_bry[length(s)]*sqrt(t[length(s)])-score_inter)+ score_inter

            c_bry_new=e_bry_new/sqrt(T_new)


            # T_new_0=new_info_time_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,beta)
            #       new_sampsz=ceiling(T_new_0*2*sigma^2)

            # new_c_bry=new_c_bry_com(s,c_bry,theta_hat_inter,theta_hat_inter_sd,T_new)
            if  (new_sampsz>N_max){
              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (last(test)>0){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              sampsz_simu=last(N)
            } else {
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              sampsz_simu=new_sampsz
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              x0_new=x0[1:new_sampsz]
              x_new= x[1:new_sampsz]

              if (direction==1){ #higher rate p is better
                theta_hat_last=mean(x_new)-mean(x0_new)
              }else if (direction==2){ #hlower rate p is better
                theta_hat_last=mean(x0_new)-mean(x_new)
              }

              score_last=theta_hat_last*T_new
              wald_last=theta_hat_last*sqrt(T_new)

              #################################################################################
              ###                   Step 4.6: Hypothesis testing                           ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              test_add=wald_last-c_bry_new

              if (test_add>0){
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
      out=c(rej_nossr,rej,stg,sampsz_simu,futi)
  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################

  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>samsz){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################
  out=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,samsz,N_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  return(out)

}

#######################################
#asd_rej_simu_adapt_OC_binary(p_t,p_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num)
#asd_rej_simu_adapt_OC_binary(p_t=0.6,p_c=0.6,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
########################################################################

##############################
#' Group sequential designs / Adaptive sequential designs: Simulations, Binary distribution, Optimizing adaptive sequential design
#'
#' @param p_t Event rate of the test arm; must be between 0 and 1. The simulation grid varies p_t
#'   over the values supplied
#' @param p_c Event rate of the control arm; must be between 0 and 1
#' @param direction Test direction: 1 = larger p is better, otherwise smaller p is better
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g.
#'   c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g.
#'   c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param alpha_2_alps_s1 Two-sided type I error spent at the interim analyses of schedule s1 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error spent at the interim analyses of schedule s2 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param theta_cut Lack of efficacy threshold on the effect-size scale: if the observed interim
#'   effect size is below theta_cut, the trial continues to the planned sample size with no sample
#'   size increase and no early termination for futility
#' @param samsz_1 Planned sample size of the control arm for the first scenario of the simulation
#'   grid
#' @param samsz_2 Planned sample size of the control arm for the second scenario of the simulation
#'   grid
#' @param N_max_1 Maximum sample size of the control arm allowed for the sample size change in the
#'   first scenario of the simulation grid
#' @param N_max_2 Maximum sample size of the control arm allowed for the sample size change in the
#'   second scenario of the simulation grid
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with six components, one for each combination of critical boundary type
#' (O'Brien-Fleming, Pocock, alpha-spending) and interim analysis schedule (s1, s2):
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{Named list with the group
#'     sequential design (Looks, Information fractions, O'Brien-Fleming boundaries) and the
#'     simulation results for the four combinations of planned sample size (samsz_1, samsz_2) and
#'     maximum sample size (N_max_1, N_max_2), each a data frame with one row per test-arm event
#'     rate and the columns p<sub>t1</sub>,..., p<sub>c</sub>, sample size for fixed design,
#'     Rejection rate without SSR, Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size,
#'     Maximum control sample size, Mean control sample size and
#'     mean control sample size conditional on increase}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with the information fractions s2}
#'   \item{Pocock boundaries, interim analysis schedule=s1}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s1}
#'   \item{Pocock boundaries, interim analysis schedule=s2}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s2}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s1}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s2}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' asd_simu_OC_binary(p_t = c(0.3,0.35,0.4), p_c = 0.2, direction = 1, s1 = c(1/3,2/3,1),
#' s2 = c(0.5,1), alpha_2 = 0.05, alpha_2_alps_s1 = c(0.001,0.001,0.048),
#' alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01, samsz_1 = 100, samsz_2 = 150,
#' N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 20)
#' }
#'
#' @export
asd_simu_OC_binary<-function(p_t,p_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)


  ###### sample size calculations for fixed design ####

  if(direction==1){ #larger mu is better
    mu=p_t-p_c
  }else if (direction==2){ #smaller mu is better
    mu=p_c-p_t
  }
  var_pc=p_c*(1-p_c)
  var_pt=p_t*(1-p_t)
  sampsz_mu=NULL

  for (i in 1:length(p_t)){
    sampsz_mu[i]=ceiling((((a+b)/mu[i])**2)*(var_pc+var_pt[i]))
  }


  ############ critical boundary calculations #####

  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(OF_bundry(s1,alpha)[1]))
  s2_OF_cbry=as.numeric(unlist(OF_bundry(s2,alpha)[1]))

  s1_pocock_cbry=as.numeric(unlist(pocock_bundry(s1,alpha)[1]))
  s2_pocock_cbry=as.numeric(unlist(pocock_bundry(s2,alpha)[1]))

  s1_alsp_cbry=as.numeric(unlist(alsp_bundry(s1,alpha_alsp_s1)[1]))
  s2_alsp_cbry=as.numeric(unlist(alsp_bundry(s2,alpha_alsp_s2)[1]))


  ##### GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_pocock_design=data.frame(cbind(c(1:length(s1)),s1,s1_pocock_cbry))
  names(s1_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s2_pocock_design=data.frame(cbind(c(1:length(s2)),s2,s2_pocock_cbry))
  names(s2_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")


  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_result_title = c("p<sub>t1</sub>,...","p<sub>c</sub>","sample size  for fixed design", "Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                        "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                        "Futility rate","Planned control sample size","Maximum control sample size",
                        "Mean control sample size","mean control sample size conditional on increase","Ratio of adaptive sample size over fixed sample size")
  simu_out_OF_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_OF_s1_samsz_1_N_max_1= rbind(simu_out_OF_s1_samsz_1_N_max_1,
                                          asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                       c_bry=s1_OF_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_OF_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_OF_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_1_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_OF_s1_samsz_1_N_max_1,samsz_ratio_OF_s1_samsz_1_N_max_1)))
  names(simu_out_OF_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_OF_s2_samsz_1_N_max_1= rbind(simu_out_OF_s2_samsz_1_N_max_1,
                                          asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                       c_bry=s2_OF_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_OF_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_OF_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s2_samsz_1_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_OF_s2_samsz_1_N_max_1,samsz_ratio_OF_s2_samsz_1_N_max_1)))
  names(simu_out_OF_s2_samsz_1_N_max_1)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_OF_s1_samsz_2_N_max_1= rbind(simu_out_OF_s1_samsz_2_N_max_1,
                                          asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                       c_bry=s1_OF_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_OF_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_OF_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_2_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_OF_s1_samsz_2_N_max_1,samsz_ratio_OF_s1_samsz_2_N_max_1)))
  names(simu_out_OF_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_OF_s2_samsz_2_N_max_1= rbind(simu_out_OF_s2_samsz_2_N_max_1,
                                          asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                       c_bry=s2_OF_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_OF_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_OF_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s2_samsz_2_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_OF_s2_samsz_2_N_max_1,samsz_ratio_OF_s2_samsz_2_N_max_1)))
  names(simu_out_OF_s2_samsz_2_N_max_1)=simu_result_title



  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####
  simu_out_OF_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_OF_s1_samsz_1_N_max_2= rbind(simu_out_OF_s1_samsz_1_N_max_2,
                                          asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                       c_bry=s1_OF_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_OF_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_OF_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_1_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_OF_s1_samsz_1_N_max_2,samsz_ratio_OF_s1_samsz_1_N_max_2)))
  names(simu_out_OF_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_OF_s2_samsz_1_N_max_2= rbind(simu_out_OF_s2_samsz_1_N_max_2,
                                          asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                       c_bry=s2_OF_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_OF_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_OF_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s2_samsz_1_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_OF_s2_samsz_1_N_max_2,samsz_ratio_OF_s2_samsz_1_N_max_2)))
  names(simu_out_OF_s2_samsz_1_N_max_2)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_OF_s1_samsz_2_N_max_2= rbind(simu_out_OF_s1_samsz_2_N_max_2,
                                          asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                       c_bry=s1_OF_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_OF_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_OF_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_2_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_OF_s1_samsz_2_N_max_2,samsz_ratio_OF_s1_samsz_2_N_max_2)))
  names(simu_out_OF_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_OF_s2_samsz_2_N_max_2= rbind(simu_out_OF_s2_samsz_2_N_max_2,
                                          asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                       c_bry=s2_OF_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_OF_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_OF_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s2_samsz_2_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_OF_s2_samsz_2_N_max_2,samsz_ratio_OF_s2_samsz_2_N_max_2)))
  names(simu_out_OF_s2_samsz_2_N_max_2)=simu_result_title



  #############################################
  #### OC simulations for Pocock boundaries ###
  #############################################


  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_pocock_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_pocock_s1_samsz_1_N_max_1= rbind(simu_out_pocock_s1_samsz_1_N_max_1,
                                              asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                           c_bry=s1_pocock_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_pocock_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_pocock_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_1_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_pocock_s1_samsz_1_N_max_1,samsz_ratio_pocock_s1_samsz_1_N_max_1)))
  names(simu_out_pocock_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_pocock_s2_samsz_1_N_max_1= rbind(simu_out_pocock_s2_samsz_1_N_max_1,
                                              asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                           c_bry=s2_pocock_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_pocock_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_pocock_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_1_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_pocock_s2_samsz_1_N_max_1,samsz_ratio_pocock_s2_samsz_1_N_max_1)))
  names(simu_out_pocock_s2_samsz_1_N_max_1)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_pocock_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_pocock_s1_samsz_2_N_max_1= rbind(simu_out_pocock_s1_samsz_2_N_max_1,
                                              asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                           c_bry=s1_pocock_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_pocock_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_pocock_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_2_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_pocock_s1_samsz_2_N_max_1,samsz_ratio_pocock_s1_samsz_2_N_max_1)))
  names(simu_out_pocock_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_pocock_s2_samsz_2_N_max_1= rbind(simu_out_pocock_s2_samsz_2_N_max_1,
                                              asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                           c_bry=s2_pocock_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_pocock_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_pocock_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_2_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_pocock_s2_samsz_2_N_max_1,samsz_ratio_pocock_s2_samsz_2_N_max_1)))
  names(simu_out_pocock_s2_samsz_2_N_max_1)=simu_result_title



  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####
  simu_out_pocock_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_pocock_s1_samsz_1_N_max_2= rbind(simu_out_pocock_s1_samsz_1_N_max_2,
                                              asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                           c_bry=s1_pocock_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_pocock_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_pocock_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_1_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_pocock_s1_samsz_1_N_max_2,samsz_ratio_pocock_s1_samsz_1_N_max_2)))
  names(simu_out_pocock_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_pocock_s2_samsz_1_N_max_2= rbind(simu_out_pocock_s2_samsz_1_N_max_2,
                                              asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                           c_bry=s2_pocock_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_pocock_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_pocock_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_1_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_pocock_s2_samsz_1_N_max_2,samsz_ratio_pocock_s2_samsz_1_N_max_2)))
  names(simu_out_pocock_s2_samsz_1_N_max_2)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_pocock_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_pocock_s1_samsz_2_N_max_2= rbind(simu_out_pocock_s1_samsz_2_N_max_2,
                                              asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                           c_bry=s1_pocock_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_pocock_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_pocock_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_2_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_pocock_s1_samsz_2_N_max_2,samsz_ratio_pocock_s1_samsz_2_N_max_2)))
  names(simu_out_pocock_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_pocock_s2_samsz_2_N_max_2= rbind(simu_out_pocock_s2_samsz_2_N_max_2,
                                              asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                           c_bry=s2_pocock_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_pocock_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_pocock_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_2_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_pocock_s2_samsz_2_N_max_2,samsz_ratio_pocock_s2_samsz_2_N_max_2)))
  names(simu_out_pocock_s2_samsz_2_N_max_2)=simu_result_title



  ##############

  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_alsp_s1_samsz_1_N_max_1= rbind(simu_out_alsp_s1_samsz_1_N_max_1,
                                            asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                         c_bry=s1_alsp_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_alsp_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_alsp_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_1_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_alsp_s1_samsz_1_N_max_1,samsz_ratio_alsp_s1_samsz_1_N_max_1)))
  names(simu_out_alsp_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_alsp_s2_samsz_1_N_max_1= rbind(simu_out_alsp_s2_samsz_1_N_max_1,
                                            asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                         c_bry=s2_alsp_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_alsp_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_alsp_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_1_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_alsp_s2_samsz_1_N_max_1,samsz_ratio_alsp_s2_samsz_1_N_max_1)))
  names(simu_out_alsp_s2_samsz_1_N_max_1)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_alsp_s1_samsz_2_N_max_1= rbind(simu_out_alsp_s1_samsz_2_N_max_1,
                                            asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                         c_bry=s1_alsp_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_alsp_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_alsp_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_2_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_alsp_s1_samsz_2_N_max_1,samsz_ratio_alsp_s1_samsz_2_N_max_1)))
  names(simu_out_alsp_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(p_t)){
    simu_out_alsp_s2_samsz_2_N_max_1= rbind(simu_out_alsp_s2_samsz_2_N_max_1,
                                            asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                         c_bry=s2_alsp_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_alsp_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_alsp_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_2_N_max_1=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_alsp_s2_samsz_2_N_max_1,samsz_ratio_alsp_s2_samsz_2_N_max_1)))
  names(simu_out_alsp_s2_samsz_2_N_max_1)=simu_result_title



  ##############
  ## N_max_2 ###
  ##############


  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_alsp_s1_samsz_1_N_max_2= rbind(simu_out_alsp_s1_samsz_1_N_max_2,
                                            asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                         c_bry=s1_alsp_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_alsp_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_alsp_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_1_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_alsp_s1_samsz_1_N_max_2,samsz_ratio_alsp_s1_samsz_1_N_max_2)))
  names(simu_out_alsp_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_alsp_s2_samsz_1_N_max_2= rbind(simu_out_alsp_s2_samsz_1_N_max_2,
                                            asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                         c_bry=s2_alsp_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_alsp_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_alsp_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_1_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_alsp_s2_samsz_1_N_max_2,samsz_ratio_alsp_s2_samsz_1_N_max_2)))
  names(simu_out_alsp_s2_samsz_1_N_max_2)=simu_result_title
  ## samsz_2###

  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_alsp_s1_samsz_2_N_max_2= rbind(simu_out_alsp_s1_samsz_2_N_max_2,
                                            asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s1,
                                                                         c_bry=s1_alsp_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_alsp_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_alsp_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_2_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_alsp_s1_samsz_2_N_max_2,samsz_ratio_alsp_s1_samsz_2_N_max_2)))
  names(simu_out_alsp_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(p_t)){
    simu_out_alsp_s2_samsz_2_N_max_2= rbind(simu_out_alsp_s2_samsz_2_N_max_2,
                                            asd_rej_simu_adapt_OC_binary(p_t=p_t[i],p_c=p_c,direction=direction,s=s2,
                                                                         c_bry=s2_alsp_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(p_t))){
    samsz_ratio_alsp_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_alsp_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_2_N_max_2=data.frame((cbind(p_t,p_c,sampsz_mu,simu_out_alsp_s2_samsz_2_N_max_2,samsz_ratio_alsp_s2_samsz_2_N_max_2)))
  names(simu_out_alsp_s2_samsz_2_N_max_2)=simu_result_title


  ##############
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


  list_pocock_1=list(s1_pocock_design,
                     simu_out_pocock_s1_samsz_1_N_max_1,
                     simu_out_pocock_s1_samsz_2_N_max_1,
                     simu_out_pocock_s1_samsz_1_N_max_2,
                     simu_out_pocock_s1_samsz_2_N_max_2
  )
  names(list_pocock_1)=c("Group sequential design",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 2")


  list_pocock_2=list(s2_pocock_design,
                     simu_out_pocock_s2_samsz_1_N_max_1,
                     simu_out_pocock_s2_samsz_2_N_max_1,
                     simu_out_pocock_s2_samsz_1_N_max_2,
                     simu_out_pocock_s2_samsz_2_N_max_2
  )
  names(list_pocock_2)=c("Group sequential design",
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
           list_pocock_1,list_pocock_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Pocock boundaries, interim analysis schedule=s1",
               "Pocock boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}
#################################
#asd_simu_OC_binary(p_t,p_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num)
##########################
#asd_simu_OC_binary(p_t=c(0.3,0.4,0.5),p_c=0.2,direction=1,
#                  s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                 theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1,sim_num=100)
##########################################################
#### 4.   simulations for Poisson distributions   ########
##########################################################

###########################################################
#### 4.1  GSD simulations for Poisson distributions  ######
###########################################################

##############

gsd_rej_rate_poisson<-function(lambda_c,lambda_t,direction,s,c_bry,samsz,sim_num){
  #This function simulates the rejection rate, given mu, s,c_bry,sim_num,sigma
  ## for a group sequential design
  ## the random samples are from normal distributions
  #mu is the effect size
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsaz is the sample size
  #sigma is the common stdard deviation
  ## When mu=0, the rejection rate is the type I error
  ## When mu>0, the rejection rate is the power

  ######################################
  ###  Preparing the parameters ########
  ######################################
  N=ceiling(samsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted

  #t<- N/2/sigma^2   # t is the vector of information times at each interim anbalysis
  NF=samsz  ## NF is the sample size at the end of the trial
  #t is the vector of information times at the interim analysis
  ############################################
  ####  begining of the simulation loops  ####
  ############################################

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      ########################################
      ### Step 1: Generate random samples ####
      ########################################
      #the loop specifies what is generated for each simulation

      ########################################
      x0=stats::rpois(samsz,lambda_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rpois(samsz,lambda_t)

      #####################################################################
      ### Step 2: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      lambda_c_hat=NULL
      lambda_t_hat=NULL
      t_hat=NULL
      theta_hat=NULL

      if (direction==1){ #higher rate p is better
      for (i in 1:length(s) ){
        lambda_c_hat[i]=mean(x0[1:N[i]])
        lambda_t_hat[i]=mean(x[1:N[i]])
        t_hat[i]= N[i]/(lambda_t_hat[i]+lambda_c_hat[i])
        theta_hat[i]=mean(x[1:N[i]])-mean(x0[1:N[i]])
      }
      } else  if (direction==2){ #lower rate p is better
      for (i in 1:length(s) ){
        lambda_c_hat[i]=mean(x0[1:N[i]])
        lambda_t_hat[i]=mean(x[1:N[i]])
        t_hat[i]= N[i]/(lambda_t_hat[i]+lambda_c_hat[i])
        theta_hat[i]=mean(x0[1:N[i]])-mean(x[1:N[i]])
      }
    }

      ###############################################
      ############## generate the parameters ########
      ###############################################



      Wald=NULL
      for (i in 1:length(s) ){
        Wald[i]=theta_hat[i]*sqrt(t_hat[i])
      }


      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry
      #at each simulation, test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if ( max(test)>0){

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        #stg identifies the interim analysis at which the null hypothesis is rejected
      }else{
        rej=0
        stg=0
      }
      if (stg>0){
        last_vt=stg
      } else {last_vt=length(s)
      }
      #######################################################
      ### Step 4: output and finish the simjulation loop ####
      #######################################################
      samsz_simu=N[last_vt]
      out=c(rej,stg,samsz_simu)
  })#(end of q loop)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  samsz_simu = as.numeric(simr[,3])
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
  rej_num=sum(rej>0)

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
  }
  mean_samsz=mean(samsz_simu)

  #############################
  ### Step 6: Final output ####
  #############################
  out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(rej_time)))
  names(out_1)=c("Looks","Information fraction","Critical boundaries","Rejection by look")
  out_2=c(samsz,rej_rt,mean_samsz)
  names(out_2)=c("Planned control sample size","Rejection rate","Mean control sample size")
  out=list(out_1,out_2)
  names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)
}
#######################################
#gsd_rej_rate_poisson(lambda_c,lambda_t,direction,s,c_bry,samsz,sim_num)
#gsd_rej_rate_poisson(lambda_c=0.8,lambda_t=0.8,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,sim_num=1000)

#####################
#############
#' Group sequential designs / Adaptive sequential designs: Simulations, Poisson distribution, Type I error and power
#'
#' @param lambda_c Mean event rate of the control arm; must be > 0
#' @param lambda_t Mean event rate of the test arm; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param samsz Planned sample size of the control arm used in the simulation; must be > 0
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries and Rejection by look (the
#'     rejection rate at each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Planned control sample size (the
#'     sample size used in the simulation), Rejection rate (the simulated power, i.e. the type I
#'     error rate when lambda_t = lambda_c) and Mean control sample size (the mean simulated
#'     sample size)}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' gsd_power_simu_poisson(lambda_c = 0.5, lambda_t = 0.8, direction = 1, s = c(1/3,2/3,1),
#' bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
gsd_power_simu_poisson<-function(lambda_c,lambda_t,direction,s,bry_type,alpha_2,samsz,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  simu_out=gsd_rej_rate_poisson(lambda_c,lambda_t,direction,s=s,c_bry=c_bry,samsz=samsz,sim_num=sim_num)
 # out_1=simu_out[[1]]
 # out_2=c(samsz,as.numeric(unlist(simu_out[[2]])))
 # names(out_2)=c("Planned sample size","rejection rate","Mean sample size")
  out=simu_out
# names(out)=c("Sequential design","Sample size and rejection rate")
  return(out)
}

##################################
#gsd_power_simu_poisson(lambda_c,lambda_t,direction,s,bry_type,alpha_2,samsz,sim_num)
#gsd_power_simu_poisson(lambda_c=0.5,lambda_t=0.8,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)
#gsd_power_simu_poisson(lambda_c=0.8,lambda_t=0.5,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000)

#gsd_power_simu_poisson(lambda_c=0.3,lambda_t=0.7,direction=1,s=c(1/3,2/3,1),bry_type=3,alpha_2=c(0.001,0.001,0.048),samsz=67,sim_num=10000)
#gsd_power_simu_poisson(lambda_c=0.7,lambda_t=0.3,direction=2,s=c(1/3,2/3,1),bry_type=3,alpha_2=c(0.001,0.001,0.048),samsz=67,sim_num=10000)

#gsd_power_simu_poisson(lambda_c=0.3,lambda_t=0.3,direction=1,s=c(1/3,2/3,1),bry_type=3,alpha_2=c(0.001,0.001,0.048),samsz=67,sim_num=100000)

#######################################

###########################################################################
#### 4.2  Adaptive sequential simulations for Poissonb distributions  #####
###########################################################################
######################################################################
asd_rej_simu_adapt_com_poisson<-function(lambda_t,lambda_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num){

  ###############################################################
  ###############################################
  ############## generate the parameters ########
  ###############################################

  N=ceiling(samsz*s)

  NF=N[length(s)]


  #t is the vector of information times at the interim analysis
  # e_bry=c_bry*sqrt(t)

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      #the loop specifies what is generated for each simulation

      ########################################
      ### Step 1: Generate random samples ####
      ########################################

      ########################################
      x0=stats::rpois(N_max,lambda_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rpois(N_max,lambda_t)

      #####################################################################
      ### Step 2: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      lambda_c_hat=NULL
      lambda_t_hat=NULL
      t_hat=NULL
      theta_hat=NULL

      if (direction==1){ #higher rate p is better
        for (i in 1:length(s) ){
          lambda_c_hat[i]=mean(x0[1:N[i]])
          lambda_t_hat[i]=mean(x[1:N[i]])
          t_hat[i]= N[i]/(lambda_t_hat[i]+lambda_c_hat[i])
          theta_hat[i]=mean(x[1:N[i]])-mean(x0[1:N[i]])
        }
      } else  if (direction==2){ #lower rate p is better
        for (i in 1:length(s) ){
          lambda_c_hat[i]=mean(x0[1:N[i]])
          lambda_t_hat[i]=mean(x[1:N[i]])
          t_hat[i]= N[i]/(lambda_t_hat[i]+lambda_c_hat[i])
          theta_hat[i]=mean(x0[1:N[i]])-mean(x[1:N[i]])
        }
      }

      ###############################################
      ############## generate the parameters ########
      ###############################################
      Wald=NULL
      for (i in 1:length(s) ){
        Wald[i]=theta_hat[i]*sqrt(t_hat[i])
      }

      theta_hat_sd=1/sqrt(t_hat)
      #############################################################################################
      ### Step 2: Obtain the Wald statistics and theta_hat for hypothesis testing and estimate ####
      #############################################################################################

      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd=theta_hat_sd[length(s)-1]
      if (theta_hat_inter<=theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry

      test_1=test[1:length(s)-1]

      ###################

      if ( max(test)>=0){
        rej_nossr=1
      }else if ( max(test)<0){
        rej_nossr=0
      }
      #####################
      # test compares the Wald statistics with the critical boundaries
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>=0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        wald_last=Wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        sampsz_simu=N[stg]
        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else  if ( max(test_1)<0){

        #################################################################################
        ###                   Step 4.1: Compare theta_hat with theta_cut             ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation
          #stg identifies the interim analysis at which the null hypothesis is rejected

          if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }
          sampsz_simu=N[length(s)]
          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          t_hat_inter_0= N[inter_vt]/(lambda_t_hat[length(s)-1]+lambda_c_hat[length(s)-1])
          T_final_0=NF/(lambda_t_hat[length(s)-1]+lambda_c_hat[length(s)-1])
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          score_inter=theta_hat_inter*t_hat_inter_0

          cond_power=stats::pnorm(theta_hat_inter*sqrt(T_final_0-t_hat_inter_0)-(c_bry[length(s)]*sqrt(T_final_0)-score_inter)/sqrt(T_final_0-t_hat_inter_0))

          if (cond_power>=1-beta){
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=theta_hat_sd[last_vt]


            sampsz_simu=N[length(s)]
            if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else{


            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##
            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(T_final_0)-score_inter)/sqrt(T_final_0-t_hat_inter_0)+stats::qnorm(1-beta))^2+t_hat_inter_0
            new_sampsz_0=ceiling(N[inter_vt]*T_new_0/t_hat_inter_0)
            new_sampsz_1=max(new_sampsz_0,NF)

            new_sampsz=min(N_max,new_sampsz_1)

            #################################################################################
            ###                   Step 4.4: Compare theta_hat with theta_cut             ####
            ###                     if new_sampsz<=N_max, proceed to new sample size     ####
            ########                   and calculate confidenc einterval                 ####
            #################################################################################
            sampsz_simu=new_sampsz
            ##### Generate addtional random samples #######
            inter_vt=length(s)-1
            #################################################################################
            ###                   Step 4.5: generate new random samples                  ####
            ###                         and calculate statistics                         ####
            #################################################################################

            x0_new=x0[1:new_sampsz]

            x_new=x[1:new_sampsz]

            lambda_c_hat_new=mean(x0_new)
            lambda_t_hat_new=mean(x_new)
            t_hat_new= new_sampsz/(lambda_c_hat_new+lambda_t_hat_new)

            if (direction==1){ #higher rate p is better
            theta_hat_new=mean(x_new)-mean(x0_new)
            }else  if (direction==2){ #lower rate p is better
              theta_hat_new=mean(x0_new)-mean(x_new)
            }

            theta_hat_sd_new=1/sqrt(t_hat_new)


            T_new=t_hat_new
            t_inter_1=N[inter_vt]/(lambda_c_hat_new+lambda_t_hat_new)
            T_final_1=NF/(lambda_c_hat_new+lambda_t_hat_new)

            if (T_new<=t_inter_1){
              e_bry_new=c_bry[length(s)]*sqrt(T_new)
            }else if (T_new>t_inter_1){
            e_bry_new=
              (sqrt(T_new-t_inter_1)/sqrt(T_final_1-t_inter_1))*(c_bry[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter
            }
            c_bry_new=e_bry_new/sqrt(T_new)

            ###############################################
            ############## generate the parameters ########
            ###############################################

            wald_last=theta_hat_new*sqrt(t_hat_new)

            #################################################################################
            ###                   Step 4.6: Hypothesis testing                           ####
            ########                   and calculate confidenc einterval                 ####
            #################################################################################

            test_add=wald_last-c_bry_new

            if (test_add>0){
              rej=1
              stg=length(s)
            }else{
              rej=0
              stg=0
            }

          }
        }
      }

      #######################################################
      ### Step 5: output and finish the simjulation loop ####
      #######################################################
      out=c(rej_nossr,rej,stg,sampsz_simu,futi)
  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################
  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(stg[stg==j])/sim_num
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
    #################
  }
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>samsz){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  sampsz_icr_rt_byLook=c(rep(0,(length(s)-1)),mean(sampsz_icr))
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################


  out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(sampsz_icr_rt_byLook),as.numeric(rej_time)))
  names(out_1)=c("Looks","Information fraction","Critical boundaries","Sample size increase","Rejection by look")
  out_2=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,samsz,N_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out_2)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  out=list(out_1,out_2)
  names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)

}
#######################################
#asd_rej_simu_adapt_com_poisson(lambda_c,lambda_t,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num)
#asd_rej_simu_adapt_com_poisson(lambda_c=0.3,lambda_t=0.7,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=67,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_rej_simu_adapt_com_poisson(lambda_t=1.1,lambda_c=1.1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=100000)
########################################################################

#' Group sequential designs / Adaptive sequential designs: Simulations, Poisson distribution, Operating characteristics
#'
#' @param lambda_t Mean event rate of the test arm; must be > 0
#' @param lambda_c Mean event rate of the control arm; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param samsz Planned sample size of the control arm; the trial starts from this sample size and
#'   it is increased only if the conditional power at the interim analysis is below the target
#' @param N_max Maximum sample size of the control arm allowed for the sample size change; the
#'   re-estimated sample size is truncated at N_max
#' @param theta_cut Lack of efficacy threshold on the effect-size scale: if the observed interim
#'   effect size is below theta_cut, the trial continues to the planned sample size with no
#'   sample size increase and no early termination for futility
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries, Sample size increase (the rate of
#'     sample size increase attributed to each look) and Rejection by look (the rejection rate at
#'     each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Rejection rate without SSR,
#'     Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size,
#'     Maximum control sample size, Mean control sample size and
#'     mean control sample size conditional on increase}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_power_simu_poisson(lambda_t = 1.1, lambda_c = 0.8, direction = 1, s = c(1/3,2/3,1),
#' bry_type = 1, alpha_2 = 0.05, samsz = 100, N_max = 200, theta_cut = 0.01,
#' beta = 0.1, sim_num = 20)
#'
#' @export
asd_power_simu_poisson<-function(lambda_t,lambda_c,direction,s,bry_type,alpha_2,samsz,N_max,theta_cut,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  # samsz=control_sample_size
  simu_out=asd_rej_simu_adapt_com_poisson(lambda_t=lambda_t,lambda_c=lambda_c,direction=direction,s=s,c_bry=c_bry,samsz=samsz,N_max=N_max,theta_cut=theta_cut,beta=beta,sim_num=sim_num)

  out=simu_out
  #names(out)=c("Sequential design","Simulation results")
 # names(out)=c("Sequential design and rejection by look","Power and sample size")

  return(out)
}

##################################
#asd_power_simu_poisson(lambda_t,lambda_c,direction=1,s,bry_type,alpha_2,samsz,N_max,theta_cut,beta,sim_num)
#asd_power_simu_poisson(lambda_t=1.1,lambda_c=0.8,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_power_simu_poisson(lambda_t=1.1,lambda_c=1.1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=1000)
#asd_power_simu_poisson(lambda_t=0.3,lambda_c=0.7,direction=2,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=67,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)

########################################################################
#########################
asd_rej_simu_adapt_OC_poisson<-function(lambda_t,lambda_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num){

  ###############################################################
  ###############################################
  ############## generate the parameters ########
  ###############################################
  N=ceiling(samsz*s)

  NF=N[length(s)]
  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      #the loop specifies what is generated for each simulation

      ########################################
      ### Step 1: Generate random samples ####
      ########################################

      ########################################
      x0=stats::rpois(N_max,lambda_c)
      #x0 is the random sample at the first interim analysis for the control
      x=stats::rpois(N_max,lambda_t)

      #####################################################################
      ### Step 2: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      lambda_c_hat=NULL
      lambda_t_hat=NULL
      t_hat=NULL
      theta_hat=NULL

      if (direction==1){ #higher rate p is better
        for (i in 1:length(s) ){
          lambda_c_hat[i]=mean(x0[1:N[i]])
          lambda_t_hat[i]=mean(x[1:N[i]])
          t_hat[i]= N[i]/(0.01+lambda_t_hat[i]+lambda_c_hat[i])
          theta_hat[i]=mean(x[1:N[i]])-mean(x0[1:N[i]])
        }
      } else  if (direction==2){ #lower rate p is better
        for (i in 1:length(s) ){
          lambda_c_hat[i]=mean(x0[1:N[i]])
          lambda_t_hat[i]=mean(x[1:N[i]])
          t_hat[i]= N[i]/(0.01+lambda_t_hat[i]+lambda_c_hat[i])
          theta_hat[i]=mean(x0[1:N[i]])-mean(x[1:N[i]])
        }
      }

      ###############################################
      ############## generate the parameters ########
      ###############################################
      Wald=NULL
      for (i in 1:length(s) ){
        Wald[i]=theta_hat[i]*sqrt(t_hat[i])
      }

      theta_hat_sd=1/sqrt(t_hat)
      #############################################################################################
      ### Step 2: Obtain the Wald statistics and theta_hat for hypothesis testing and estimate ####
      #############################################################################################

      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd=theta_hat_sd[length(s)-1]
      if (theta_hat_inter<=theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry

      test_1=test[1:length(s)-1]

      ###################

      if ( max(test)>=0){
        rej_nossr=1
      }else if ( max(test)<0){
        rej_nossr=0
      }
      #####################
      # test compares the Wald statistics with the critical boundaries
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>=0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg
        wald_last=Wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        sampsz_simu=N[stg]
        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else  if ( max(test_1)<0){

        #################################################################################
        ###                   Step 4.1: Compare theta_hat with theta_cut             ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################
        if (futi==1){## (theta_hat_two_stg[inter_vt]>theta_cut) indicated no efficacy, the trial stops at the planned final analysis, no sample size re-estimation
          #stg identifies the interim analysis at which the null hypothesis is rejected

          if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }
          sampsz_simu=N[length(s)]
          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          t_hat_inter_0= N[inter_vt]/(lambda_t_hat[length(s)-1]+lambda_c_hat[length(s)-1])
          T_final_0=NF/(lambda_t_hat[length(s)-1]+lambda_c_hat[length(s)-1])
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          score_inter=theta_hat_inter*t_hat_inter_0

          cond_power=stats::pnorm(theta_hat_inter*sqrt(T_final_0-t_hat_inter_0)-(c_bry[length(s)]*sqrt(T_final_0)-score_inter)/sqrt(T_final_0-t_hat_inter_0))

          if (cond_power>=1-beta){
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=theta_hat_sd[last_vt]


            sampsz_simu=N[length(s)]
            if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else{


            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##
            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(T_final_0)-score_inter)/sqrt(T_final_0-t_hat_inter_0)+stats::qnorm(1-beta))^2+t_hat_inter_0
            new_sampsz_0=ceiling(N[inter_vt]*T_new_0/t_hat_inter_0)
            new_sampsz_1=max(new_sampsz_0,NF)

            new_sampsz=min(N_max,new_sampsz_1)

            #################################################################################
            ###                   Step 4.4: Compare theta_hat with theta_cut             ####
            ###                     if new_sampsz<=N_max, proceed to new sample size     ####
            ########                   and calculate confidenc einterval                 ####
            #################################################################################
            sampsz_simu=new_sampsz
            ##### Generate addtional random samples #######
            inter_vt=length(s)-1
            #################################################################################
            ###                   Step 4.5: generate new random samples                  ####
            ###                         and calculate statistics                         ####
            #################################################################################

            x0_new=x0[1:new_sampsz]

            x_new=x[1:new_sampsz]

            lambda_c_hat_new=mean(x0_new)
            lambda_t_hat_new=mean(x_new)
            t_hat_new= new_sampsz/(lambda_c_hat_new+lambda_t_hat_new)

            if (direction==1){ #higher rate p is better
              theta_hat_new=mean(x_new)-mean(x0_new)
            }else  if (direction==2){ #lower rate p is better
              theta_hat_new=mean(x0_new)-mean(x_new)
            }

            theta_hat_sd_new=1/sqrt(t_hat_new)


            T_new=t_hat_new
            t_inter_1=N[inter_vt]/(lambda_c_hat_new+lambda_t_hat_new)
            T_final_1=NF/(lambda_c_hat_new+lambda_t_hat_new)

            if (T_new<=t_inter_1){
              e_bry_new=c_bry[length(s)]*sqrt(T_new)
            }else if (T_new>t_inter_1){
              e_bry_new=
                (sqrt(T_new-t_inter_1)/sqrt(T_final_1-t_inter_1))*(c_bry[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter
            }
            c_bry_new=e_bry_new/sqrt(T_new)

            ###############################################
            ############## generate the parameters ########
            ###############################################

            wald_last=theta_hat_new*sqrt(t_hat_new)

            #################################################################################
            ###                   Step 4.6: Hypothesis testing                           ####
            ########                   and calculate confidenc einterval                 ####
            #################################################################################

            test_add=wald_last-c_bry_new

            if (test_add>0){
              rej=1
              stg=length(s)
            }else{
              rej=0
              stg=0
            }

          }
        }
      }

      #######################################################
      ### Step 5: output and finish the simjulation loop ####
      #######################################################
      out=c(rej_nossr,rej,stg,sampsz_simu,futi)
  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>samsz){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################

  out=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,samsz,N_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  return(out)

}

############
#asd_rej_simu_adapt_OC_poisson(lambda_t,lambda_c,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num)
#asd_rej_simu_adapt_OC_poisson(lambda_t=1.1,lambda_c=0.8,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_rej_simu_adapt_OC_poisson(lambda_t=c(0.7,0.78,0.8),lambda_c=0.2,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)

###############################

##########################
#' Group sequential designs / Adaptive sequential designs: Simulations, Poisson distribution, Optimizing adaptive sequential design
#'
#' @param lambda_t Mean event rate of the test arm; must be > 0. The simulation grid varies
#'   lambda_t over the values supplied
#' @param lambda_c Mean event rate of the control arm; must be > 0
#' @param direction Test direction: 1 = larger lambda is better, otherwise smaller lambda is better
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g.
#'   c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g.
#'   c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param alpha_2_alps_s1 Two-sided type I error spent at the interim analyses of schedule s1 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error spent at the interim analyses of schedule s2 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param theta_cut Lack of efficacy threshold on the effect-size scale: if the observed interim
#'   effect size is below theta_cut, the trial continues to the planned sample size with no
#'   sample size increase and no early termination for futility
#' @param samsz_1 Planned sample size of the control arm for the first scenario of the simulation
#'   grid
#' @param samsz_2 Planned sample size of the control arm for the second scenario of the simulation
#'   grid
#' @param N_max_1 Maximum sample size of the control arm allowed for the sample size change in the
#'   first scenario of the simulation grid
#' @param N_max_2 Maximum sample size of the control arm allowed for the sample size change in the
#'   second scenario of the simulation grid
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with six components, one for each combination of critical boundary type
#' (O'Brien-Fleming, Pocock, alpha-spending) and interim analysis schedule (s1, s2):
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{Named list with the group
#'     sequential design (Looks, Information fractions, O'Brien-Fleming boundaries) and the
#'     simulation results for the four combinations of planned sample size (samsz_1, samsz_2) and
#'     maximum sample size (N_max_1, N_max_2), each a data frame with one row per test-arm event
#'     rate and the columns lambda_t1,..., lambda_c, sample size for fixed design,
#'     Rejection rate without SSR, Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size,
#'     Maximum control sample size, Mean control sample size and
#'     mean control sample size conditional on increase}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with the information fractions s2}
#'   \item{Pocock boundaries, interim analysis schedule=s1}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s1}
#'   \item{Pocock boundaries, interim analysis schedule=s2}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s2}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s1}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s2}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_simu_OC_poisson(lambda_t = c(0.3,0.4,0.5), lambda_c = 0.2, direction = 1, s1 = c(1/3,2/3,1),
#' s2 = c(0.5,1), alpha_2 = 0.05, alpha_2_alps_s1 = c(0.001,0.001,0.048),
#' alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01, samsz_1 = 50, samsz_2 = 80,
#' N_max_1 = 120, N_max_2 = 150, beta = 0.1, sim_num = 20)
#'
#' @export
asd_simu_OC_poisson<-function(lambda_t,lambda_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)


  ###### sample size calculations for fixed design ####

  if(direction==1){ #larger mu is better
    mu=lambda_t-lambda_c
  }else if (direction==2){ #smaller mu is better
    mu=lambda_c-lambda_t
  }

  sampsz_mu=NULL

  for (i in 1:length(lambda_t)){
    sampsz_mu[i]=ceiling((((a+b)/mu[i])**2)*(lambda_t[i]+lambda_c))
  }


  ############ critical boundary calculations #####

  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(OF_bundry(s1,alpha)[1]))
  s2_OF_cbry=as.numeric(unlist(OF_bundry(s2,alpha)[1]))

  s1_pocock_cbry=as.numeric(unlist(pocock_bundry(s1,alpha)[1]))
  s2_pocock_cbry=as.numeric(unlist(pocock_bundry(s2,alpha)[1]))

  s1_alsp_cbry=as.numeric(unlist(alsp_bundry(s1,alpha_alsp_s1)[1]))
  s2_alsp_cbry=as.numeric(unlist(alsp_bundry(s2,alpha_alsp_s2)[1]))

  ##### GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_pocock_design=data.frame(cbind(c(1:length(s1)),s1,s1_pocock_cbry))
  names(s1_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s2_pocock_design=data.frame(cbind(c(1:length(s2)),s2,s2_pocock_cbry))
  names(s2_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")


  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_result_title = c(LABEL()$lambda_t1,LABEL()$lambda_c,"sample size  for fixed design", "Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                        "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                        "Futility rate","Planned control sample size","Maximum control sample size",
                        "Mean control sample size","mean control sample size conditional on increase","Ratio of adaptive sample size over fixed sample size")
  simu_out_OF_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_OF_s1_samsz_1_N_max_1= rbind(simu_out_OF_s1_samsz_1_N_max_1,
                                          asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                        c_bry=s1_OF_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_OF_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_OF_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s1_samsz_1_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_OF_s1_samsz_1_N_max_1,samsz_ratio_OF_s1_samsz_1_N_max_1)))
  names(simu_out_OF_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_OF_s2_samsz_1_N_max_1= rbind(simu_out_OF_s2_samsz_1_N_max_1,
                                          asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                        c_bry=s2_OF_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_OF_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_OF_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s2_samsz_1_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_OF_s2_samsz_1_N_max_1,samsz_ratio_OF_s2_samsz_1_N_max_1)))
  names(simu_out_OF_s2_samsz_1_N_max_1)=simu_result_title

  ## samsz_2###

  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_OF_s1_samsz_2_N_max_1= rbind(simu_out_OF_s1_samsz_2_N_max_1,
                                          asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                        c_bry=s1_OF_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_OF_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_OF_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s1_samsz_2_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_OF_s1_samsz_2_N_max_1,samsz_ratio_OF_s1_samsz_2_N_max_1)))
  names(simu_out_OF_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_OF_s2_samsz_2_N_max_1= rbind(simu_out_OF_s2_samsz_2_N_max_1,
                                          asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                        c_bry=s2_OF_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_OF_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_OF_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s2_samsz_2_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_OF_s2_samsz_2_N_max_1,samsz_ratio_OF_s2_samsz_2_N_max_1)))
  names(simu_out_OF_s2_samsz_2_N_max_1)=simu_result_title

  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_OF_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_OF_s1_samsz_1_N_max_2= rbind(simu_out_OF_s1_samsz_1_N_max_2,
                                          asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                        c_bry=s1_OF_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_OF_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_OF_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s1_samsz_1_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_OF_s1_samsz_1_N_max_2,samsz_ratio_OF_s1_samsz_1_N_max_2)))
  names(simu_out_OF_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_OF_s2_samsz_1_N_max_2= rbind(simu_out_OF_s2_samsz_1_N_max_2,
                                          asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                        c_bry=s2_OF_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_OF_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_OF_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s2_samsz_1_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_OF_s2_samsz_1_N_max_2,samsz_ratio_OF_s2_samsz_1_N_max_2)))
  names(simu_out_OF_s2_samsz_1_N_max_2)=simu_result_title

  ## samsz_2###

  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_OF_s1_samsz_2_N_max_2= rbind(simu_out_OF_s1_samsz_2_N_max_2,
                                          asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                        c_bry=s1_OF_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_OF_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_OF_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_OF_s1_samsz_2_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_OF_s1_samsz_2_N_max_2,samsz_ratio_OF_s1_samsz_2_N_max_2)))
  names(simu_out_OF_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_OF_s2_samsz_2_N_max_2= rbind(simu_out_OF_s2_samsz_2_N_max_2,
                                          asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                        c_bry=s2_OF_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_OF_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_OF_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s2_samsz_2_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_OF_s2_samsz_2_N_max_2,samsz_ratio_OF_s2_samsz_2_N_max_2)))
  names(simu_out_OF_s2_samsz_2_N_max_2)=simu_result_title


  #############################################
  #### OC simulations for Pocock boundaries ###
  #############################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_pocock_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_pocock_s1_samsz_1_N_max_1= rbind(simu_out_pocock_s1_samsz_1_N_max_1,
                                              asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                            c_bry=s1_pocock_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_pocock_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_pocock_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s1_samsz_1_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_pocock_s1_samsz_1_N_max_1,samsz_ratio_pocock_s1_samsz_1_N_max_1)))
  names(simu_out_pocock_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_pocock_s2_samsz_1_N_max_1= rbind(simu_out_pocock_s2_samsz_1_N_max_1,
                                              asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                            c_bry=s2_pocock_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_pocock_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_pocock_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s2_samsz_1_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_pocock_s2_samsz_1_N_max_1,samsz_ratio_pocock_s2_samsz_1_N_max_1)))
  names(simu_out_pocock_s2_samsz_1_N_max_1)=simu_result_title

  ## samsz_2###

  ##### s1 ####
  simu_out_pocock_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_pocock_s1_samsz_2_N_max_1= rbind(simu_out_pocock_s1_samsz_2_N_max_1,
                                              asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                            c_bry=s1_pocock_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_pocock_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_pocock_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s1_samsz_2_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_pocock_s1_samsz_2_N_max_1,samsz_ratio_pocock_s1_samsz_2_N_max_1)))
  names(simu_out_pocock_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_pocock_s2_samsz_2_N_max_1= rbind(simu_out_pocock_s2_samsz_2_N_max_1,
                                              asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                            c_bry=s2_pocock_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_pocock_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_pocock_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s2_samsz_2_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_pocock_s2_samsz_2_N_max_1,samsz_ratio_pocock_s2_samsz_2_N_max_1)))
  names(simu_out_pocock_s2_samsz_2_N_max_1)=simu_result_title

  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_pocock_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_pocock_s1_samsz_1_N_max_2= rbind(simu_out_pocock_s1_samsz_1_N_max_2,
                                              asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                            c_bry=s1_pocock_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_pocock_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_pocock_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_1_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_pocock_s1_samsz_1_N_max_2,samsz_ratio_pocock_s1_samsz_1_N_max_2)))
  names(simu_out_pocock_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_pocock_s2_samsz_1_N_max_2= rbind(simu_out_pocock_s2_samsz_1_N_max_2,
                                              asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                            c_bry=s2_pocock_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_pocock_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_pocock_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_1_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_pocock_s2_samsz_1_N_max_2,samsz_ratio_pocock_s2_samsz_1_N_max_2)))
  names(simu_out_pocock_s2_samsz_1_N_max_2)=simu_result_title


  ## samsz_2###

  ##### s1 ####
  simu_out_pocock_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_pocock_s1_samsz_2_N_max_2= rbind(simu_out_pocock_s1_samsz_2_N_max_2,
                                              asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                            c_bry=s1_pocock_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_pocock_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_pocock_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s1_samsz_2_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_pocock_s1_samsz_2_N_max_2,samsz_ratio_pocock_s1_samsz_2_N_max_2)))
  names(simu_out_pocock_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_pocock_s2_samsz_2_N_max_2= rbind(simu_out_pocock_s2_samsz_2_N_max_2,
                                              asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                            c_bry=s2_pocock_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_pocock_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_pocock_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_pocock_s2_samsz_2_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_pocock_s2_samsz_2_N_max_2,samsz_ratio_pocock_s2_samsz_2_N_max_2)))
  names(simu_out_pocock_s2_samsz_2_N_max_2)=simu_result_title

  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_alsp_s1_samsz_1_N_max_1= rbind(simu_out_alsp_s1_samsz_1_N_max_1,
                                            asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                          c_bry=s1_alsp_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_alsp_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_alsp_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_1_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_alsp_s1_samsz_1_N_max_1,samsz_ratio_alsp_s1_samsz_1_N_max_1)))
  names(simu_out_alsp_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_alsp_s2_samsz_1_N_max_1= rbind(simu_out_alsp_s2_samsz_1_N_max_1,
                                            asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                          c_bry=s2_alsp_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_alsp_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_alsp_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_1_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_alsp_s2_samsz_1_N_max_1,samsz_ratio_alsp_s2_samsz_1_N_max_1)))
  names(simu_out_alsp_s2_samsz_1_N_max_1)=simu_result_title


  ## samsz_2###

  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_alsp_s1_samsz_2_N_max_1= rbind(simu_out_alsp_s1_samsz_2_N_max_1,
                                            asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                          c_bry=s1_alsp_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_alsp_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_alsp_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_2_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_alsp_s1_samsz_2_N_max_1,samsz_ratio_alsp_s1_samsz_2_N_max_1)))
  names(simu_out_alsp_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(lambda_t)){
    simu_out_alsp_s2_samsz_2_N_max_1= rbind(simu_out_alsp_s2_samsz_2_N_max_1,
                                            asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                          c_bry=s2_alsp_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_alsp_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_alsp_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_2_N_max_1=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_alsp_s2_samsz_2_N_max_1,samsz_ratio_alsp_s2_samsz_2_N_max_1)))
  names(simu_out_alsp_s2_samsz_2_N_max_1)=simu_result_title


  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_alsp_s1_samsz_1_N_max_2= rbind(simu_out_alsp_s1_samsz_1_N_max_2,
                                            asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                          c_bry=s1_alsp_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_alsp_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_alsp_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_1_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_alsp_s1_samsz_1_N_max_2,samsz_ratio_alsp_s1_samsz_1_N_max_2)))
  names(simu_out_alsp_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_alsp_s2_samsz_1_N_max_2= rbind(simu_out_alsp_s2_samsz_1_N_max_2,
                                            asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                          c_bry=s2_alsp_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_alsp_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_alsp_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_1_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_alsp_s2_samsz_1_N_max_2,samsz_ratio_alsp_s2_samsz_1_N_max_2)))
  names(simu_out_alsp_s2_samsz_1_N_max_2)=simu_result_title


  ## samsz_2###

  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_alsp_s1_samsz_2_N_max_2= rbind(simu_out_alsp_s1_samsz_2_N_max_2,
                                            asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s1,
                                                                          c_bry=s1_alsp_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_alsp_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_alsp_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s1_samsz_2_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_alsp_s1_samsz_2_N_max_2,samsz_ratio_alsp_s1_samsz_2_N_max_2)))
  names(simu_out_alsp_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(lambda_t)){
    simu_out_alsp_s2_samsz_2_N_max_2= rbind(simu_out_alsp_s2_samsz_2_N_max_2,
                                            asd_rej_simu_adapt_OC_poisson(lambda_t=lambda_t[i],lambda_c,direction=direction,s=s2,
                                                                          c_bry=s2_alsp_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(lambda_t))){
    samsz_ratio_alsp_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_alsp_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }


  simu_out_alsp_s2_samsz_2_N_max_2=data.frame((cbind(lambda_t,lambda_c,sampsz_mu,simu_out_alsp_s2_samsz_2_N_max_2,samsz_ratio_alsp_s2_samsz_2_N_max_2)))
  names(simu_out_alsp_s2_samsz_2_N_max_2)=simu_result_title

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


  list_pocock_1=list(s1_pocock_design,
                     simu_out_pocock_s1_samsz_1_N_max_1,
                     simu_out_pocock_s1_samsz_2_N_max_1,
                     simu_out_pocock_s1_samsz_1_N_max_2,
                     simu_out_pocock_s1_samsz_2_N_max_2
  )
  names(list_pocock_1)=c("Group sequential design",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 2")


  list_pocock_2=list(s2_pocock_design,
                     simu_out_pocock_s2_samsz_1_N_max_1,
                     simu_out_pocock_s2_samsz_2_N_max_1,
                     simu_out_pocock_s2_samsz_1_N_max_2,
                     simu_out_pocock_s2_samsz_2_N_max_2
  )
  names(list_pocock_2)=c("Group sequential design",
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
           list_pocock_1,list_pocock_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Pocock boundaries, interim analysis schedule=s1",
               "Pocock boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}
#####################
##########################
#asd_simu_OC_poisson(lambda_t,lambda_c,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num)
###########
# asd_simu_OC_poisson(lambda_t=c(0.3,0.4,0.5),lambda_c=0.2,direction=1,
#                     s1=c(1/3,2/3,1),s2=c(0.5,1),alpha_2=0.05, alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                     theta_cut=0.01,samsz_1=50,samsz_2=80,N_max_1=120,N_max_2=150,beta=0.1,sim_num=1000)
# asd_simu_OC_poisson(lambda_t=c(0.3,0.4,0.5),lambda_c=0.2,direction=1,
#                     s1=c(0.3,0.7,1),s2=c(0.5,1),alpha_2=0.05, alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                    theta_cut=0.01,samsz_1=50,samsz_2=80,N_max_1=120,N_max_2=150,beta=0.1,sim_num=1000)

####################################################################
#### 5.   simulations for negative binomial distributions   ########
####################################################################

#########################################################################
#### 5.1  GSD simulations for negative binomial distributions    ########
#########################################################################

################################################
gsd_rej_rate_negbinom<-function(r_t,r_c,kappa,direction,s,c_bry,samsz,sim_num){
  #This function simulates the rejection rate, given p_t, p_c,s,c_bry,sim_num
  ## for a group sequential design
  ## the random samples are from normal distributions
  #c_bry is the vector of critical values
  #s is the vector of information frations
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsz is the sample size

  #####################################################


  ######################################
  ###  Preparing the parameters ########
  ######################################
  N=ceiling(samsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  NF=N[length(s)] # is the sample size at the final analysis
  mu=-log(r_t/r_c)

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      ########################################
      ### Step 1: Generate random samples ####
      ########################################
      ###########################################
      ### geting random sample ################
      x_c=cbind(rep(1,samsz),stats::rnbinom(n=samsz, size=kappa, mu=r_c))
      x_t=cbind(rep(2,samsz),stats::rnbinom(n=samsz, size=kappa, mu=r_t))
      ####  analysis parameters ###
      r_c_hat=NULL
      r_t_hat=NULL
      t_hat=NULL
      theta_hat=NULL

      for (i in 1:length(s)){

        x_c_inter=x_c[1:N[i],]
        x_t_inter=x_t[1:N[i],]

        x_inter=rbind(x_c_inter,x_t_inter)

        negglm_inter=glm.nb(as.numeric(x_inter[,2])~as.factor(x_inter[,1]),link="log")
        sum_inter=summary(negglm_inter)

        se_inter=sum_inter$coefficients[2,2]
        r_c_inter=mean(x_c_inter[,2])
        r_t_inter=mean(x_t_inter[,2])

        if (direction ==1){ #higher rate r is better
          theta_hat_inter=log(r_t_inter/r_c_inter)
        } else if (direction ==2){ #lower rate r is better
          theta_hat_inter=-log(r_t_inter/r_c_inter)
        }

        theta_hat[i]=theta_hat_inter

        t_hat[i]=1/se_inter^2

      }

      #####################################################################
      ### Step 2: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      Wald=NULL

      for (i in 1:length(s) ){
        Wald=theta_hat[i]*sqrt(t_hat[i])

      }

      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry
      #at each simulation, test compares the Wald statistics with the critical boundaries
      #any value of test>0 indicates that the wald statistics crossed the critical boundaries
      if ( max(test)>0){

        rej=1
        #If max(test)>0, te nulll hypothesis is rejected
        test_t=test[test>0]
        stg=match(test_t[1],test)
        #stg identifies the interim analysis at which the null hypothesis is rejected
      }else{
        rej=0
        stg=0
      }
      if (stg>0){
        last_vt=stg
      } else {last_vt=length(s)
      }
      #######################################################
      ### Step 4: output and finish the simjulation loop ####
      #######################################################
      samsz_simu=N[last_vt]
      out=c(rej,stg,samsz_simu)
  })#(end of q loop)

  ############################################
  ### Step 5: Generate summary statistics ####
  ############################################
  rej= as.numeric(simr[,1])
  stg = as.numeric(simr[,2])
  samsz_simu = as.numeric(simr[,3])
  rej_rt=mean(rej)
  #rej_rt is the rejection rate
 # rej_num=sum(rej>0)

  rej_time=rep(0,length(s))
  for (j in 1:length(s)){
    rej_time[j]=length(rej[stg==j])/(dim(simr)[1])
  }
    #rej_time[j] is the proportion
    #that the null hypothesis were rejected at the j-th interim analysis
    #################
  mean_samsz=mean(samsz_simu)

  #############################
  ### Step 6: Final output ####
  #############################
  out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(rej_time)))
  names(out_1)=c("Looks","Information fraction","Critical boundaries","Rejection by look")
  out_2=c(samsz,rej_rt,mean_samsz)
  names(out_2)=c("Planned  control sample size","Rejection rate","Mean control sample size")
  out=list(out_1,out_2)
  names(out)=c("Sequential design and rejection by look","Power and sample size")
  return(out)
}
#######################################
#gsd_rej_rate_negbinom(r_t,r_c,kappa,direction,s,c_bry,samsz,sim_num)
#gsd_rej_rate_negbinom(r_t=1,r_c=1.4,kappa=1,direction=2,s=c(1/3,2/3,1),c_bry=c(6.441526,3.540084,1.960064),samsz=100,sim_num=100)
#gsd_rej_rate_negbinom(r_t=1,r_c=1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(6.441526,3.540084,1.960064),samsz=100,sim_num=100000)


##############################################

#' Group sequential designs / Adaptive sequential designs: Simulations, Negative binomial distribution, Type I error and power
#'
#' @param r_t Mean of the test arm; must be > 0
#' @param r_c Mean of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param samsz Planned sample size of the control arm used in the simulation; must be > 0
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries and Rejection by look (the
#'     rejection rate at each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Planned  control sample size
#'     (the sample size used in the simulation), Rejection rate (the simulated power, i.e. the
#'     type I error rate when r_t = r_c) and Mean control sample size (the mean simulated sample
#'     size)}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' gsd_power_simu_negbinom(r_t = 1.1, r_c = 1.1, kappa = 1, direction = 1, s = c(1/3,2/3,1),
#' bry_type = 1, alpha_2 = 0.05, samsz = 100, sim_num = 20)
#'
#' @export
gsd_power_simu_negbinom<-function(r_t,r_c,kappa,direction,s,bry_type,alpha_2,samsz,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  simu_out=gsd_rej_rate_negbinom(r_t=r_t,r_c=r_c,kappa=kappa,direction=direction,s=s,c_bry=c_bry,samsz=samsz,sim_num=sim_num)
#  out_1=simu_out[[1]]
 # out_2=c(samsz,as.numeric(unlist(simu_out[[2]])))
 # names(out_2)=c("Planned sample size","rejection rate","Mean sample size")
  out=simu_out
  #names(out)=c("Sequential design","Sample size and rejection rate")
  return(out)
}

##################################
#gsd_power_simu_negbinom(r_t,r_c,kappa,direction,s,bry_type,alpha_2,samsz,sim_num)
#gsd_power_simu_negbinom(r_t=1.1,r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,sim_num=1000000)
#gsd_power_simu_negbinom(r_t=1.1,r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=2,alpha_2=0.05,samsz=100,sim_num=1000000)

#gsd_power_simu_negbinom(r_t=1.1,r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=3,alpha_2=c(0.001,0.001,0.048),samsz=100,sim_num=1000000)

#gsd_power_simu_negbinom(r_t=1,r_c=0.5,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=112,sim_num=10000)
#gsd_power_simu_negbinom(r_t=0.3,r_c=0.2,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=112,sim_num=100)



######################################################################################
#### 5.2  Adaptive sequential simulations for negative binomial distributions  #######
######################################################################################
asd_rej_simu_adapt_com_neg_binomial<-function(r_t,r_c,kappa,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num){
  #s is the vector of information frations
  #c_bry is the vector of critical values
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsaz is sample size
  #N_max is themaximum sample size

  ######################################
  ###  Preparing the parameters ########
  ######################################
  N=ceiling(samsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  NF=samsz # is the sample size at the final analysis
  mu=abs(log(r_t/r_c))

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      ########################################
      ### Step 1: Generate random samples ####
      ########################################
      ### geting random sample ################
      x_c=cbind(rep(1,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_c))
      x_t=cbind(rep(2,N_max),stats::rnbinom(n=N_max, size=kappa, mu=r_t))
      ####  analysis parameters ###
      t_hat=NULL
      theta_hat=NULL
      theta_hat_sd=NULL

      for (i in 1:length(s)){

        x_c_temp=x_c[1:N[i],]
        x_t_temp=x_t[1:N[i],]

        x_temp=rbind(x_c_temp,x_t_temp)

        negglm_temp=glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log")
        sum_temp=summary(negglm_temp)

        se_temp=sum_temp$coefficients[2,2]
        r_c_temp=mean(x_c_temp[,2])
        r_t_temp=mean(x_t_temp[,2])

        if (direction ==1){ #higher rate r is better
          theta_hat_temp=log(r_t_temp/r_c_temp)
        } else if (direction ==2){ #lower rate r is better
          theta_hat_temp=-log(r_t_temp/r_c_temp)
        }

        theta_hat[i]=theta_hat_temp
        t_hat[i]=1/se_temp^2
        theta_hat_sd[i]=1/sqrt(t_hat[i])
      }

      #####################################################################
      ### Step 2: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      Wald=NULL

      for (i in 1:length(s) ){
        Wald[i]=theta_hat[i]*sqrt(t_hat[i])

      }

      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd= theta_hat_sd[length(s)-1]
      if (theta_hat_inter<=theta_cut) {futi=1}else{futi=0}
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry

      test_1=test[1:length(s)-1]
      ##################################

      if ( max(test)>=0){
        rej_nossr=1
      }else if ( max(test)<0){
        rej_nossr=0
      }
      #################

      # test compares the Wald statistics with the critical boundaries
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if ( max(test_1)>=0){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        wald_last=Wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        sampsz_simu=N[stg]
        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else if ( max(test_1)<0){

        #################################################################################
        ###                   Step 4.1: Compare theta_hat with theta_cut             ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (futi==1){##  the trial stops at the planned final analysis, no sample size re-estimation
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=theta_hat_sd[last_vt]


          if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }
          sampsz_simu=N[length(s)]

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          score_inter=theta_hat_inter*t_hat[inter_vt]
          T_f_1=max(t_hat[length(s)],(NF/N[inter_vt])*t_hat[inter_vt])
          cond_power=stats::pnorm(theta_hat_inter*sqrt(T_f_1-t_hat[inter_vt])-(c_bry[length(s)]*sqrt(t_hat[length(s)])-score_inter)/sqrt(T_f_1-t_hat[inter_vt]))

          if (cond_power>=1-beta){

            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=1/sqrt(t_hat[last_vt])

            sampsz_simu=N[length(s)]
            if (test[length(s)]>0){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else if (cond_power<1-beta){

            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##
           # T_f_1=max(t_hat[length(s)],(NF/N[inter_vt])*t_hat[inter_vt])
            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(t_hat[length(s)])-score_inter)/sqrt(T_f_1-t_hat[inter_vt])+stats::qnorm(1-beta))^2+t_hat[inter_vt]
            new_sampsz=max(ceiling(N[inter_vt]*T_new_0/t_hat[inter_vt]),NF)

            if  (new_sampsz>N_max){
              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (last(test)>0){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              sampsz_simu=last(N)
            } else if  (new_sampsz<=N_max){
              case=5
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              sampsz_simu=new_sampsz
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

            x_c_new=x_c[1:new_sampsz,]
            x_t_new=x_t[1:new_sampsz,]

            x_new=rbind(x_c_new,x_t_new)

            negglm_new=glm.nb(as.numeric(x_new[,2])~as.factor(x_new[,1]),link="log")
            sum_new=summary(negglm_new)

            se_new=sum_new$coefficients[2,2]
            r_c_new=mean(x_c_new[,2])
            r_t_new=mean(x_t_new[,2])

            if (direction ==1){ #higher rate r is better
              theta_hat_new=log(r_t_new/r_c_new)
            } else if (direction ==2){ #lower rate r is better
              theta_hat_new=-log(r_t_new/r_c_new)
            }

            t_hat_new=1/se_new^2
            T_new=t_hat_new

            t_inter_1=(N[inter_vt]/new_sampsz)*t_hat_new
            T_final_1= (NF/new_sampsz)*t_hat_new

           # e_bry_new=
           #   (sqrt(T_new-t_inter_1)/sqrt(T_final_1-t_inter_1))*(c_bry[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter
            e_bry_new=
              (sqrt(new_sampsz-N[inter_vt])/sqrt(NF-N[inter_vt]))*(c_bry[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter

            c_bry_new=e_bry_new/sqrt(T_new)
            ###############################################
            ############## generate the parameters ########
            ###############################################
            score_new=theta_hat_new*t_hat_new
            wald_new=theta_hat_new*sqrt(t_hat_new)

            #################################################################################
            ###                   Step 4.6: Hypothesis testing                           ####
            ########                   and calculate confidenc einterval                 ####
            #################################################################################

            test_add=wald_new-c_bry_new

            if (test_add>0){
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
      out=c(rej_nossr,rej,stg,sampsz_simu,futi)

  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
      rej_nossr= as.numeric(simr[,1])
      rej= as.numeric(simr[,2])
      stg= as.numeric(simr[,3])
      sampsz_simu=as.numeric(simr[,4])
      futi=as.numeric(simr[,5])
      ###########################
      rej_time=rep(0,length(s))
      for (j in 1:length(s)){
        rej_time[j]=length(stg[stg==j])/dim(simr)[1]
        #rej_time[j] is the proportion
        #that the null hypothesis were rejected at the j-th interim analysis
        #################
      }
      rej_nossr_rt=mean(rej_nossr)
      rej_rt=mean(rej)

      av_samsz=mean(sampsz_simu)
      futi_rt=mean(futi)
      sampsz_icr=NULL
      for (i in 1:length(sampsz_simu)){
        if (sampsz_simu[i]>samsz){
          sampsz_icr[i]=1
        }else{
          sampsz_icr[i]=0
        }
      }
      sampsz_icr_rt=mean(sampsz_icr)
      sampsz_icr_rt_byLook=c(rep(0,(length(s)-1)),mean(sampsz_icr))
      #### no ssr vs ssr when sampsz_icr==1##
      rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
      rej_sampsz_icr= rej[sampsz_icr==1]
      samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
      rej_sampsz_icr= rej[sampsz_icr==1]
      rej_nossr_rt_sampsz_icr=mean( rej_nossr_sampsz_icr)
      rej_rt_sampsz_icr=mean(rej_sampsz_icr)
      mean_samsz_simu_control_sampsz_icr=mean(samsz_simu_sampsz_icr)
      #############################
      ### Step 7: Final output ####
      #############################


      out_1=data.frame(cbind(c(1:length(s)),as.numeric(s),as.numeric(c_bry),as.numeric(sampsz_icr_rt_byLook),as.numeric(rej_time)))
      names(out_1)=c("Looks","Information fraction","Critical boundaries","Sample size increase","Rejection by look")
      out_2=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
              futi_rt,samsz,N_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
      names(out_2)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                     "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                     "Futility rate","Planned control sample size","Maximum control sample size",
                     "Mean control sample size","mean control sample size conditional on increase")
      out=list(out_1,out_2)
      names(out)=c("Sequential design and rejection by look","Power and sample size")
      return(out)

}
#######################################
#asd_rej_simu_adapt_com_neg_binomial(r_t,r_c,kappa,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num)
#asd_rej_simu_adapt_com_neg_binomial(r_c=0.5,r_t=1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=112,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_rej_simu_adapt_com_neg_binomial(r_c=1.1,r_t=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=100000)
#asd_rej_simu_adapt_com_neg_binomial(r_t=0.3,r_c=0.2,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_rej_simu_adapt_com_neg_binomial(r_t=0.4,r_c=0.2,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=1000)

########################################################################
#' Group sequential designs / Adaptive sequential designs: Simulations, Negative binomial distribution, Operating characteristics
#'
#' @param r_t Mean of the test arm; must be > 0
#' @param r_c Mean of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param bry_type Critical boundary type: 1 = O'Brien-Fleming, 2 = Pocock, 3 = alpha-spending
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param samsz Planned sample size of the control arm; the trial starts from this sample size and
#'   it is increased only if the conditional power at the interim analysis is below the target
#' @param N_max Maximum sample size of the control arm allowed for the sample size change; the
#'   re-estimated sample size is truncated at N_max
#' @param theta_cut Lack of efficacy threshold on the effect-size scale: if the observed interim
#'   effect size is below theta_cut, the trial continues to the planned sample size with no
#'   sample size increase and no early termination for futility
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Sequential design and rejection by look}{Data frame with one row per analysis and the
#'     columns Looks, Information fraction, Critical boundaries, Sample size increase (the rate of
#'     sample size increase attributed to each look) and Rejection by look (the rejection rate at
#'     each interim analysis)}
#'   \item{Power and sample size}{Named vector with the components Rejection rate without SSR,
#'     Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size,
#'     Maximum control sample size, Mean control sample size and
#'     mean control sample size conditional on increase}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' asd_power_simu_neg_binomial(r_t = 0.3, r_c = 0.2, kappa = 1, direction = 1, s = c(1/3,2/3,1),
#' bry_type = 1, alpha_2 = 0.05, samsz = 112, N_max = 200, theta_cut = 0.01,
#' beta = 0.1, sim_num = 20)
#'
#' @export
asd_power_simu_neg_binomial<-function(r_t,r_c,kappa,direction,s,bry_type,alpha_2,samsz,N_max,theta_cut,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  alpha=alpha_2/2

  if (bry_type==1) {
    boundary_out=OF_bundry(s,alpha)
  } else if (bry_type==2) {
    boundary_out=pocock_bundry(s,alpha)
  }else if (bry_type==3) {
    boundary_out=alsp_bundry(s,alpha)
  }

  c_bry=as.numeric(unlist(boundary_out[1]))

  # samsz=control_sample_size
  simu_out=asd_rej_simu_adapt_com_neg_binomial(r_t=r_t,r_c=r_c,kappa=kappa,direction=direction,s=s,c_bry=c_bry,samsz=samsz,N_max=N_max,theta_cut=theta_cut,beta=beta,sim_num=sim_num)


  out=simu_out


  return(out)
}

##################################
#asd_power_simu_neg_binomial(r_t,r_c,kappa,direction,s,bry_type,alpha_2,samsz,N_max,theta_cut,beta,sim_num)
#asd_power_simu_neg_binomial(r_t=1.1,r_c=1.1,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=100000)
#asd_power_simu_neg_binomial(r_t=1,r_c=0.5,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=112,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_power_simu_neg_binomial(r_t=0.3,r_c=0.2,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=112,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_power_simu_neg_binomial(r_t=0.4,r_c=0.2,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=112,N_max=500,theta_cut=0.01,beta=0.1,sim_num=10000)
#asd_power_simu_neg_binomial(r_t=0.5,r_c=0.2,kappa=1,direction=1,s=c(1/3,2/3,1),bry_type=1,alpha_2=0.05,samsz=112,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)

########################################################################
asd_rej_simu_adapt_OC_neg_binomial<-function(r_t,r_c,kappa,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num){
  #s is the vector of information frations
  #c_bry is the vector of critical values
  #c_bry and s have the same length
  #sim_num=simulation number
  #samsaz is sample size
  #N_max is themaximum sample size
  #
  # Numerical safeguards (2026-10-31):
  # When a simulated replicate gives zero counts in one or both arms, the
  # effect estimate is +/-Inf or NaN and the information may be NA. The
  # original code then evaluated if (max(test) >= 0) etc. on an NA condition,
  # which raises "missing value where TRUE/FALSE needed" and aborts the whole
  # batch of simulations. All conditions below are therefore evaluated with
  # base::isTRUE(), and a replicate that cannot be evaluated is recorded as
  # "no rejection, no sample size increase, finish at the planned sample size".
  # The arithmetic on well behaved replicates is unchanged.

  ######################################
  ###  Preparing the parameters ########
  ######################################
  N=ceiling(samsz*s)  # N is the vetoc of sample sizes at which interim analysis will be conducted
  NF=samsz # is the sample size at the final analysis
  mu=abs(log(r_t/r_c))
  # number of observations that must be generated: N_max only caps the
  # re-estimated sample size; if N_max < planned sample size we still need the
  # planned sample size, otherwise x_c[1:N[i],] would return NA rows
  N_gen=max(N_max,max(N))

  ############################################
  ####  begining of the simulation loops  ####
  ############################################
  simr <- simuloop(sim_num, "ASD", function() {

      ########################################
      ### Step 1: Generate random samples ####
      ########################################
      ### geting random sample ################
      x_c=cbind(rep(1,N_gen),stats::rnbinom(n=N_gen, size=kappa, mu=r_c))
      x_t=cbind(rep(2,N_gen),stats::rnbinom(n=N_gen, size=kappa, mu=r_t))
      ####  analysis parameters ###
      t_hat=NULL
      theta_hat=NULL
      theta_hat_sd=NULL

      # defaults: a replicate that cannot be evaluated finishes at the planned
      # sample size without rejection (keeps every output field well defined)
      rej=0
      stg=0
      sampsz_simu=N[length(s)]
      futi=0

      for (i in 1:length(s)){

        x_c_temp=x_c[1:N[i],]
        x_t_temp=x_t[1:N[i],]

        x_temp=rbind(x_c_temp,x_t_temp)

        negglm_temp=tryCatch(MASS::glm.nb(as.numeric(x_temp[,2])~as.factor(x_temp[,1]),link="log"),
                             error=function(e) NULL)

        if (is.null(negglm_temp)) {
          theta_hat[i]=NA
          t_hat[i]=NA
          theta_hat_sd[i]=NA
        } else {
          sum_temp=summary(negglm_temp)

          se_temp=sum_temp$coefficients[2,2]
          r_c_temp=mean(x_c_temp[,2])
          r_t_temp=mean(x_t_temp[,2])

          if (direction ==1){ #higher rate r is better
            theta_hat_temp=log(r_t_temp/r_c_temp)
          } else if (direction ==2){ #lower rate r is better
            theta_hat_temp=-log(r_t_temp/r_c_temp)
          }

          theta_hat[i]=theta_hat_temp
          t_hat[i]=1/se_temp^2
          theta_hat_sd[i]=1/sqrt(t_hat[i])
        }
      }

      #####################################################################
      ### Step 2: calculate the Wald statistics for hypothesis testing ####
      #####################################################################
      Wald=NULL

      for (i in 1:length(s) ){
        Wald[i]=theta_hat[i]*sqrt(t_hat[i])

      }

      theta_hat_inter=theta_hat[length(s)-1]
      theta_hat_inter_sd= theta_hat_sd[length(s)-1]
      futi=if (base::isTRUE(theta_hat_inter<=theta_cut)) 1 else 0
      ###################################################################################################
      ### Step 3: Hypothesis testing: comparing the Wald statistics with the critical boundary c_bry ####
      ###################################################################################################

      test=Wald-c_bry

      test_1=test[1:length(s)-1]
      ##################################

      if (base::isTRUE(max(test)>=0)){
        rej_nossr=1
      }else{
        rej_nossr=0
      }
      #################

      # test compares the Wald statistics with the critical boundaries
      #any value of test_1>0 indicates that the Wald statistics crossed the critical boundaries at up to length(s)-1
      if (base::isTRUE(max(test_1)>=0)){#  if ( max(test_1)>0){ #null hypothesis rejected, no more random sampling

        rej=1
        test_t=test_1[test_1>0]
        stg=match(test_t[1],test_1)
        #stg identifies the interim analysis at which the null hypothesis is rejected
        last_vt=stg

        wald_last=Wald[last_vt]
        theta_hat_last=theta_hat[last_vt]
        theta_hat_last_sd=theta_hat_sd[last_vt]
        sampsz_simu=N[stg]
        #################################################################################
        ###                   Step 4: Compare theta_hat with theta_cut               ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ###              if theta_hat>theta_cut, re-estimate sampe size              ####
        #################################################################################
      }else{

        #################################################################################
        ###                   Step 4.1: Compare theta_hat with theta_cut             ####
        ### if theta_hat<=theta_cut, conplete trial to planned sample size and stop  ####
        ########                   and calculate confidenc einterval                 ####
        #################################################################################

        if (futi==1){##  the trial stops at the planned final analysis, no sample size re-estimation
          #stg identifies the interim analysis at which the null hypothesis is rejected
          last_vt=length(s)
          theta_hat_last=theta_hat[last_vt]
          theta_hat_last_sd=theta_hat_sd[last_vt]


          if (base::isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
            rej=1
            stg=length(s)
          }else { # the NULL hypothesis is NOT rejected at the final analysis
            rej=0
            stg=0
          }
          sampsz_simu=N[length(s)]

          ####################### this completes the estimation when there is no saple size change ########

        } else if (futi==0){ # No futility
          #################################################################################
          ###                   Step 4.2: Compare theta_hat with theta_cut             ####
          ###              if theta_hat>theta_cut, re-estimate sampe size              ####
          #################################################################################
          #sample size re-estimation procedure##
          inter_vt=length(s)-1
          theta_hat_inter=theta_hat[inter_vt]
          theta_hat_inter_sd=theta_hat_sd[inter_vt]
          score_inter=theta_hat_inter*t_hat[inter_vt]
          T_f_1=max(t_hat[length(s)],(NF/N[inter_vt])*t_hat[inter_vt])
          cond_power=stats::pnorm(theta_hat_inter*sqrt(T_f_1-t_hat[inter_vt])-(c_bry[length(s)]*sqrt(t_hat[length(s)])-score_inter)/sqrt(T_f_1-t_hat[inter_vt]))

          if (!is.finite(cond_power)) {
            # conditional power cannot be evaluated (degenerate counts or
            # unavailable information): do not increase, finish as planned
            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=theta_hat_sd[last_vt]

            sampsz_simu=N[length(s)]
            if (base::isTRUE(test[length(s)]>0)){
              rej=1
              stg=length(s)
            }else{
              rej=0
              stg=0
            }

          } else if (cond_power>=1-beta){

            last_vt=length(s)
            theta_hat_last=theta_hat[last_vt]
            theta_hat_last_sd=1/sqrt(t_hat[last_vt])

            sampsz_simu=N[length(s)]
            if (base::isTRUE(test[length(s)]>0)){ # the NULL hypothesis is rejected at the final analysis
              rej=1
              stg=length(s)
            }else { # the NULL hypothesis is NOT rejected at the final analysis
              rej=0
              stg=0
            }

          }else{

            #################################################################################
            ###                   Step 4.2: Compare theta_hat with theta_cut             ####
            ###              if theta_hat>theta_cut, re-estimate sampe size              ####
            #################################################################################
            #sample size re-estimation procedure##
            # T_f_1=max(t_hat[length(s)],(NF/N[inter_vt])*t_hat[inter_vt])
            T_new_0=theta_hat_inter^(-2)*((c_bry[length(s)]*sqrt(t_hat[length(s)])-score_inter)/sqrt(T_f_1-t_hat[inter_vt])+stats::qnorm(1-beta))^2+t_hat[inter_vt]
            new_sampsz=max(ceiling(N[inter_vt]*T_new_0/t_hat[inter_vt]),NF)

            if  (!is.finite(new_sampsz) || new_sampsz>N_max){
              #################################################################################
              ###                   Step 4.3: Compare theta_hat with theta_cut             ####
              ###     if new_sampsz>N_max, conplete trial to planned sample size and stop  ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################

              if (base::isTRUE(test[length(s)]>0)){
                rej=1
                stg=length(s)
              }else{
                rej=0
                stg=0
              }
              sampsz_simu=N[length(s)]
            } else {
              case=5
              #################################################################################
              ###                   Step 4.4: Compare theta_hat with theta_cut             ####
              ###                     if new_sampsz<=N_max, proceed to new sample size     ####
              ########                   and calculate confidenc einterval                 ####
              #################################################################################
              sampsz_simu=new_sampsz
              ##### Generate addtional random samples #######
              inter_vt=length(s)-1
              #################################################################################
              ###                   Step 4.5: generate new random samples                  ####
              ###                         and calculate statistics                         ####
              #################################################################################

              x_c_new=x_c[1:new_sampsz,]
              x_t_new=x_t[1:new_sampsz,]

              x_new=rbind(x_c_new,x_t_new)

              negglm_new=tryCatch(MASS::glm.nb(as.numeric(x_new[,2])~as.factor(x_new[,1]),link="log"),
                                  error=function(e) NULL)

              if (is.null(negglm_new)) {
                # the refitted model is not available: report no rejection
                rej=0
                stg=0
              } else {
                sum_new=summary(negglm_new)

                se_new=sum_new$coefficients[2,2]
                r_c_new=mean(x_c_new[,2])
                r_t_new=mean(x_t_new[,2])

                if (direction ==1){ #higher rate r is better
                  theta_hat_new=log(r_t_new/r_c_new)
                } else if (direction ==2){ #lower rate r is better
                  theta_hat_new=-log(r_t_new/r_c_new)
                }

                t_hat_new=1/se_new^2
                T_new=t_hat_new

                t_inter_1=(N[inter_vt]/new_sampsz)*t_hat_new
                T_final_1= (NF/new_sampsz)*t_hat_new

                # e_bry_new=
                #   (sqrt(T_new-t_inter_1)/sqrt(T_final_1-t_inter_1))*(c_bry[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter
                e_bry_new=
                  (sqrt(new_sampsz-N[inter_vt])/sqrt(NF-N[inter_vt]))*(c_bry[length(s)]*sqrt(T_final_1)-score_inter)+ score_inter

                c_bry_new=e_bry_new/sqrt(T_new)
                ###############################################
                ############## generate the parameters ########
                ###############################################
                score_new=theta_hat_new*t_hat_new
                wald_new=theta_hat_new*sqrt(t_hat_new)

                #################################################################################
                ###                   Step 4.6: Hypothesis testing                           ####
                ########                   and calculate confidenc einterval                 ####
                #################################################################################

                test_add=wald_new-c_bry_new

                if (base::isTRUE(test_add>0)){
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
      }

      #######################################################
      ### Step 5: output and finish the simjulation loop ####
      #######################################################
      out=c(rej_nossr,rej,stg,sampsz_simu,futi)

  })#(end of simulation loops)
  #############termination timepoint####
  ############################################
  ### Step 6: Generate summary statistics ####
  ############################################
  rej_nossr= as.numeric(simr[,1])
  rej= as.numeric(simr[,2])
  stg= as.numeric(simr[,3])
  sampsz_simu=as.numeric(simr[,4])
  futi=as.numeric(simr[,5])
  ###########################
  rej_nossr_rt=mean(rej_nossr)
  rej_rt=mean(rej)

  av_samsz=mean(sampsz_simu)
  futi_rt=mean(futi)
  sampsz_icr=NULL
  for (i in 1:length(sampsz_simu)){
    if (sampsz_simu[i]>samsz){
      sampsz_icr[i]=1
    }else{
      sampsz_icr[i]=0
    }
  }
  sampsz_icr_rt=mean(sampsz_icr)
  #### no ssr vs ssr when sampsz_icr==1##
  rej_nossr_sampsz_icr= rej_nossr[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  samsz_simu_sampsz_icr= sampsz_simu[sampsz_icr==1]
  rej_sampsz_icr= rej[sampsz_icr==1]
  rej_nossr_rt_sampsz_icr=if (length(rej_nossr_sampsz_icr)==0) NA_real_ else mean( rej_nossr_sampsz_icr)
  rej_rt_sampsz_icr=if (length(rej_sampsz_icr)==0) NA_real_ else mean(rej_sampsz_icr)
  mean_samsz_simu_control_sampsz_icr=if (length(samsz_simu_sampsz_icr)==0) NA_real_ else mean(samsz_simu_sampsz_icr)
  #############################
  ### Step 7: Final output ####
  #############################
  out=c(rej_nossr_rt,rej_rt,sampsz_icr_rt,rej_nossr_rt_sampsz_icr,rej_rt_sampsz_icr,
          futi_rt,samsz,N_max,av_samsz,mean_samsz_simu_control_sampsz_icr)
  names(out)=c("Rejection rate without SSR","Rejection rate with SSR","Rate of sample size increase",
                 "Rejection rate without SSR/(meeting sample size incerase criteria)","Rejection rate with SSR/(sample size incerase)",
                 "Futility rate","Planned control sample size","Maximum control sample size",
                 "Mean control sample size","mean control sample size conditional on increase")
  return(out)

}

##################
#asd_rej_simu_adapt_OC_neg_binomial(r_t,r_c,kappa,direction,s,c_bry,samsz,N_max,theta_cut,beta,sim_num)
#asd_rej_simu_adapt_OC_neg_binomial(r_t=0.4,r_c=0.2,kappa=1,direction=1,s=c(1/3,2/3,1),c_bry=c(3.471095,2.454434,2.004037),samsz=100,N_max=200,theta_cut=0.01,beta=0.1,sim_num=10000)
################
#######################
#' Group sequential designs / Adaptive sequential designs: Simulations, Negative binomial distribution, Optimizing adaptive sequential design
#'
#' @param r_t Mean of the test arm; must be > 0. The simulation grid varies r_t over the values
#'   supplied
#' @param r_c Mean of the control arm; must be > 0
#' @param kappa Dispersion parameter of the negative binomial model; must be > 0
#' @param direction Test direction: 1 = larger r is better, otherwise smaller r is better
#' @param s1 Vector of information fractions of the first interim analysis schedule, e.g.
#'   c(1/3, 2/3, 1); the last element must be 1
#' @param s2 Vector of information fractions of the second interim analysis schedule, e.g.
#'   c(0.5, 1); the last element must be 1
#' @param alpha_2 Type I error for the two-sided test; must be between 0 and 1, and alpha_2/2 is the
#'   one-sided type I error used for the critical boundaries
#' @param alpha_2_alps_s1 Two-sided type I error spent at the interim analyses of schedule s1 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param alpha_2_alps_s2 Two-sided type I error spent at the interim analyses of schedule s2 with
#'   the alpha-spending boundary; each element must be between 0 and 1 and the elements must sum to
#'   alpha_2
#' @param theta_cut Lack of efficacy threshold on the effect-size scale: if the observed interim
#'   effect size is below theta_cut, the trial continues to the planned sample size with no
#'   sample size increase and no early termination for futility
#' @param samsz_1 Planned sample size of the control arm for the first scenario of the simulation
#'   grid
#' @param samsz_2 Planned sample size of the control arm for the second scenario of the simulation
#'   grid
#' @param N_max_1 Maximum sample size of the control arm allowed for the sample size change in the
#'   first scenario of the simulation grid
#' @param N_max_2 Maximum sample size of the control arm allowed for the sample size change in the
#'   second scenario of the simulation grid
#' @param beta Conditional type II error; 1 - beta is the target conditional power for the sample
#'   size re-estimation; must be between 0 and 1
#' @param sim_num Number of simulation replicates; must be > 0
#'
#' @return A named list with six components, one for each combination of critical boundary type
#' (O'Brien-Fleming, Pocock, alpha-spending) and interim analysis schedule (s1, s2):
#' \describe{
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s1}{Named list with the group
#'     sequential design (Looks, Information fractions, O'Brien-Fleming boundaries) and the
#'     simulation results for the four combinations of planned sample size (samsz_1, samsz_2) and
#'     maximum sample size (N_max_1, N_max_2), each a data frame with one row per test-arm mean and
#'     the columns r<sub>t1</sub>,..., r<sub>c</sub>, sample size  for fixed design,
#'     Rejection rate without SSR, Rejection rate with SSR, Rate of sample size increase,
#'     Rejection rate without SSR/(meeting sample size incerase criteria),
#'     Rejection rate with SSR/(sample size incerase), Futility rate, Planned control sample size,
#'     Maximum control sample size, Mean control sample size and
#'     mean control sample size conditional on increase}
#'   \item{O'Brien-Fleming boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with the information fractions s2}
#'   \item{Pocock boundaries, interim analysis schedule=s1}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s1}
#'   \item{Pocock boundaries, interim analysis schedule=s2}{The same group sequential design and
#'     simulation results computed with Pocock boundaries and the information fractions s2}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s1}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s1}
#'   \item{Alpha-spending boundaries, interim analysis schedule=s2}{The same group sequential
#'     design and simulation results computed with alpha-spending boundaries and the information
#'     fractions s2}
#' }
#'
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}.
#' @references H.H. Müller, H. Schäfer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#' @references Gao, P., L. Liu, and C. Mehta. 2013 Oct 15. Exact inference for adaptive group sequential designs. Statistics in Medicine 32(23):3991–4005. doi:10.1002/sim.5847.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' asd_simu_OC_neg_binomial(r_t = c(0.3,0.4,0.5), r_c = 0.2, kappa = 1, direction = 1,
#' s1 = c(1/3,2/3,1), s2 = c(0.5,1), alpha_2 = 0.05, alpha_2_alps_s1 = c(0.001,0.001,0.048),
#' alpha_2_alps_s2 = c(0.001,0.049), theta_cut = 0.01, samsz_1 = 100, samsz_2 = 150,
#' N_max_1 = 1000, N_max_2 = 1500, beta = 0.1, sim_num = 20)
#' }
#'
#' @export
asd_simu_OC_neg_binomial<-function(r_t,r_c,kappa,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num){
  # This simulation function checks the actua power and sample size for different scenerios##
  # alpha_s is two sided type I error rate

  ## Input checks: fail early with an explicit message instead of producing
  ## undefined boundaries or NA rows deep inside the simulation loop.
  if (length(r_t) < 1 || any(!is.finite(r_t)) || any(r_t <= 0) ||
      !is.finite(r_c) || r_c <= 0) {
    stop("r_t and r_c must be positive finite values")
  }
  if (!is.finite(kappa) || kappa <= 0) {
    stop("kappa must be a positive finite value")
  }
  if (!(length(direction) == 1 && direction %in% c(1, 2))) {
    stop("direction must be 1 (larger r is better) or 2 (smaller r is better)")
  }
  if (length(s1) != length(alpha_2_alps_s1)) {
    stop("length(s1) must equal length(alpha_2_alps_s1)")
  }
  if (length(s2) != length(alpha_2_alps_s2)) {
    stop("length(s2) must equal length(alpha_2_alps_s2)")
  }
  if (any(diff(c(0, s1)) <= 0) || abs(s1[length(s1)] - 1) > 1e-8) {
    stop("s1 must be strictly increasing and end at 1")
  }
  if (any(diff(c(0, s2)) <= 0) || abs(s2[length(s2)] - 1) > 1e-8) {
    stop("s2 must be strictly increasing and end at 1")
  }
  if (any(alpha_2_alps_s1 < 0) || any(alpha_2_alps_s2 < 0)) {
    stop("alpha spending values must be non-negative")
  }
  if (any(samsz_1 > N_max_1) || any(samsz_2 > N_max_2)) {
    warning("planned sample size exceeds N_max in at least one scenario; ",
            "no sample size increase is possible in that scenario")
  }

  alpha=alpha_2/2
  alpha_alsp_s1=alpha_2_alps_s1/2
  alpha_alsp_s2=alpha_2_alps_s2/2

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)

  #####################################################
  ## varaince function per Zhu, harris, 2012 ####
  v_fun=function(r_0,r_1,theta,nu_t,kappa){
    V_01=(1+theta)/(nu_t*theta*r_0)+(1+theta)*kappa/theta
    V_02=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
    V_03=(1+theta)**2/(nu_t*theta*(r_0+theta*r_1))+(1+theta)*kappa/theta
    V_max=max(V_01,V_02,V_02)
    return(V_max)
  }
  ###########
  rand_ratio=1
  nu_t=1

  var_0=NULL

  for (i in 1:length(r_t)){

    var_0[i]=v_fun(r_0=r_c,r_1=r_t[i],theta=rand_ratio,nu_t=nu_t,kappa=kappa)
  }
  ###### sample size calculations for fixed design ####

  if(direction==1){ #larger r is better
    mu=log(r_t/r_c)
  }else if (direction==2){ #smaller r is better
    mu=-log(r_t/r_c)
  }

  sampsz_mu=NULL

  for (i in 1:length(r_t)){
    sampsz_mu[i]=ceiling((((a+b)/mu[i])**2)*var_0[i])
  }
  ############ critical boundary calculations #####

  ### critical boundaries for s1 ############

  s1_OF_cbry=as.numeric(unlist(OF_bundry(s1,alpha)[1]))
  s2_OF_cbry=as.numeric(unlist(OF_bundry(s2,alpha)[1]))

  s1_pocock_cbry=as.numeric(unlist(pocock_bundry(s1,alpha)[1]))
  s2_pocock_cbry=as.numeric(unlist(pocock_bundry(s2,alpha)[1]))

  s1_alsp_cbry=as.numeric(unlist(alsp_bundry(s1,alpha_alsp_s1)[1]))
  s2_alsp_cbry=as.numeric(unlist(alsp_bundry(s2,alpha_alsp_s2)[1]))


  ##### GSD designs #####
  s1_OF_design=data.frame(cbind(c(1:length(s1)),s1,s1_OF_cbry))
  names(s1_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s2_OF_design=data.frame(cbind(c(1:length(s2)),s2,s2_OF_cbry))
  names(s2_OF_design)=c("Looks","Information fractions","O'Brien-Fleming boundaries")

  s1_pocock_design=data.frame(cbind(c(1:length(s1)),s1,s1_pocock_cbry))
  names(s1_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s2_pocock_design=data.frame(cbind(c(1:length(s2)),s2,s2_pocock_cbry))
  names(s2_pocock_design)=c("Looks","Information fractions","Pocock boundaries")

  s1_alsp_design=data.frame(cbind(c(1:length(s1)),s1,alpha_alsp_s1,s1_alsp_cbry))
  names(s1_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending boundaries")

  s2_alsp_design=data.frame(cbind(c(1:length(s2)),s2,alpha_alsp_s2,s2_alsp_cbry))
  names(s2_alsp_design)=c("Looks","Information fractions","Alpha-spending","Alpha-spending  boundaries")


  #########################################
  #### OC simulations for OF boundaries ###
  #########################################

  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_result_title = c("r<sub>t1</sub>,...","r<sub>c</sub>","sample size  for fixed design","Rejection rate without SSR",
                        "Rejection rate with SSR","Rate of sample size increase",
                        "Rejection rate without SSR/(meeting sample size incerase criteria)",
                        "Rejection rate with SSR/(sample size incerase)",
                        "Futility rate","Planned control sample size","Maximum control sample size",
                        "Mean control sample size","mean control sample size conditional on increase",
                        "Ratio of adaptive sample size over fixed sample size")
  simu_out_OF_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_OF_s1_samsz_1_N_max_1= rbind(simu_out_OF_s1_samsz_1_N_max_1,
                                          asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                             c_bry=s1_OF_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_OF_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_OF_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s1_samsz_1_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_OF_s1_samsz_1_N_max_1,samsz_ratio_OF_s1_samsz_1_N_max_1)))
  names(simu_out_OF_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_OF_s2_samsz_1_N_max_1= rbind(simu_out_OF_s2_samsz_1_N_max_1,
                                          asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                             c_bry=s2_OF_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_OF_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_OF_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s2_samsz_1_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_OF_s2_samsz_1_N_max_1,samsz_ratio_OF_s2_samsz_1_N_max_1)))
  names(simu_out_OF_s2_samsz_1_N_max_1)=simu_result_title

  ## samsz_2###
  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_OF_s1_samsz_2_N_max_1= rbind(simu_out_OF_s1_samsz_2_N_max_1,
                                          asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                             c_bry=s1_OF_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_OF_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_OF_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s1_samsz_2_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_OF_s1_samsz_2_N_max_1,samsz_ratio_OF_s1_samsz_2_N_max_1)))
  names(simu_out_OF_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_OF_s2_samsz_2_N_max_1= rbind(simu_out_OF_s2_samsz_2_N_max_1,
                                          asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                             c_bry=s2_OF_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_OF_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_OF_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s2_samsz_2_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_OF_s2_samsz_2_N_max_1,samsz_ratio_OF_s2_samsz_2_N_max_1)))
  names(simu_out_OF_s2_samsz_2_N_max_1)=simu_result_title
  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_OF_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_OF_s1_samsz_1_N_max_2= rbind(simu_out_OF_s1_samsz_1_N_max_2,
                                          asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                             c_bry=s1_OF_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_OF_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_OF_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s1_samsz_1_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_OF_s1_samsz_1_N_max_2,samsz_ratio_OF_s1_samsz_1_N_max_2)))
  names(simu_out_OF_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_OF_s2_samsz_1_N_max_2= rbind(simu_out_OF_s2_samsz_1_N_max_2,
                                          asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                             c_bry=s2_OF_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_OF_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_OF_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s2_samsz_1_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_OF_s2_samsz_1_N_max_2,samsz_ratio_OF_s2_samsz_1_N_max_2)))
  names(simu_out_OF_s2_samsz_1_N_max_2)=simu_result_title

  ## samsz_2###
  ##### s1 ####
  simu_out_OF_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_OF_s1_samsz_2_N_max_2= rbind(simu_out_OF_s1_samsz_2_N_max_2,
                                          asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                             c_bry=s1_OF_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_OF_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_OF_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s1_samsz_2_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_OF_s1_samsz_2_N_max_2,samsz_ratio_OF_s1_samsz_2_N_max_2)))
  names(simu_out_OF_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_OF_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_OF_s2_samsz_2_N_max_2= rbind(simu_out_OF_s2_samsz_2_N_max_2,
                                          asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                             c_bry=s2_OF_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_OF_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_OF_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_OF_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_OF_s2_samsz_2_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_OF_s2_samsz_2_N_max_2,samsz_ratio_OF_s2_samsz_2_N_max_2)))
  names(simu_out_OF_s2_samsz_2_N_max_2)=simu_result_title

  ##############

  #############################################
  #### OC simulations for Pocock boundaries ###
  #############################################


  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_pocock_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_pocock_s1_samsz_1_N_max_1= rbind(simu_out_pocock_s1_samsz_1_N_max_1,
                                              asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                 c_bry=s1_pocock_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_pocock_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_pocock_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s1_samsz_1_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_pocock_s1_samsz_1_N_max_1,samsz_ratio_pocock_s1_samsz_1_N_max_1)))
  names(simu_out_pocock_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_pocock_s2_samsz_1_N_max_1= rbind(simu_out_pocock_s2_samsz_1_N_max_1,
                                              asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                 c_bry=s2_pocock_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_pocock_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_pocock_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s2_samsz_1_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_pocock_s2_samsz_1_N_max_1,samsz_ratio_pocock_s2_samsz_1_N_max_1)))
  names(simu_out_pocock_s2_samsz_1_N_max_1)=simu_result_title

  ## samsz_2###
  ##### s1 ####
  simu_out_pocock_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_pocock_s1_samsz_2_N_max_1= rbind(simu_out_pocock_s1_samsz_2_N_max_1,
                                              asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                 c_bry=s1_pocock_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_pocock_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_pocock_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s1_samsz_2_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_pocock_s1_samsz_2_N_max_1,samsz_ratio_pocock_s1_samsz_2_N_max_1)))
  names(simu_out_pocock_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_pocock_s2_samsz_2_N_max_1= rbind(simu_out_pocock_s2_samsz_2_N_max_1,
                                              asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                 c_bry=s2_pocock_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_pocock_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_pocock_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s2_samsz_2_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_pocock_s2_samsz_2_N_max_1,samsz_ratio_pocock_s2_samsz_2_N_max_1)))
  names(simu_out_pocock_s2_samsz_2_N_max_1)=simu_result_title
  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_pocock_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_pocock_s1_samsz_1_N_max_2= rbind(simu_out_pocock_s1_samsz_1_N_max_2,
                                              asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                 c_bry=s1_pocock_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_pocock_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_pocock_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s1_samsz_1_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_pocock_s1_samsz_1_N_max_2,samsz_ratio_pocock_s1_samsz_1_N_max_2)))
  names(simu_out_pocock_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_pocock_s2_samsz_1_N_max_2= rbind(simu_out_pocock_s2_samsz_1_N_max_2,
                                              asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                 c_bry=s2_pocock_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_pocock_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_pocock_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s2_samsz_1_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_pocock_s2_samsz_1_N_max_2,samsz_ratio_pocock_s2_samsz_1_N_max_2)))
  names(simu_out_pocock_s2_samsz_1_N_max_2)=simu_result_title

  ## samsz_2###
  ##### s1 ####
  simu_out_pocock_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_pocock_s1_samsz_2_N_max_2= rbind(simu_out_pocock_s1_samsz_2_N_max_2,
                                              asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                                 c_bry=s1_pocock_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_pocock_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_pocock_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s1_samsz_2_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_pocock_s1_samsz_2_N_max_2,samsz_ratio_pocock_s1_samsz_2_N_max_2)))
  names(simu_out_pocock_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_pocock_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_pocock_s2_samsz_2_N_max_2= rbind(simu_out_pocock_s2_samsz_2_N_max_2,
                                              asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                                 c_bry=s2_pocock_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_pocock_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_pocock_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_pocock_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_pocock_s2_samsz_2_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_pocock_s2_samsz_2_N_max_2,samsz_ratio_pocock_s2_samsz_2_N_max_2)))
  names(simu_out_pocock_s2_samsz_2_N_max_2)=simu_result_title
  ##############

  ###########################################
  #### OC simulations for alsp boundaries ###
  ###########################################


  ##############
  ## N_max_1 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_alsp_s1_samsz_1_N_max_1= rbind(simu_out_alsp_s1_samsz_1_N_max_1,
                                            asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                               c_bry=s1_alsp_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_1_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_alsp_s1_samsz_1_N_max_1[i]=as.numeric(simu_out_alsp_s1_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_alsp_s1_samsz_1_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_alsp_s1_samsz_1_N_max_1,samsz_ratio_alsp_s1_samsz_1_N_max_1)))
  names(simu_out_alsp_s1_samsz_1_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_alsp_s2_samsz_1_N_max_1= rbind(simu_out_alsp_s2_samsz_1_N_max_1,
                                            asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                               c_bry=s2_alsp_cbry,samsz=samsz_1,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_1_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_alsp_s2_samsz_1_N_max_1[i]=as.numeric(simu_out_alsp_s2_samsz_1_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_alsp_s2_samsz_1_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_alsp_s2_samsz_1_N_max_1,samsz_ratio_alsp_s2_samsz_1_N_max_1)))
  names(simu_out_alsp_s2_samsz_1_N_max_1)=simu_result_title

  ## samsz_2###
  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_alsp_s1_samsz_2_N_max_1= rbind(simu_out_alsp_s1_samsz_2_N_max_1,
                                            asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                               c_bry=s1_alsp_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_2_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_alsp_s1_samsz_2_N_max_1[i]=as.numeric(simu_out_alsp_s1_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_alsp_s1_samsz_2_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_alsp_s1_samsz_2_N_max_1,samsz_ratio_alsp_s1_samsz_2_N_max_1)))
  names(simu_out_alsp_s1_samsz_2_N_max_1)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_1=NULL
  for (i in 1:length(r_t)){
    simu_out_alsp_s2_samsz_2_N_max_1= rbind(simu_out_alsp_s2_samsz_2_N_max_1,
                                            asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                               c_bry=s2_alsp_cbry,samsz=samsz_2,N_max=N_max_1,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_2_N_max_1=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_alsp_s2_samsz_2_N_max_1[i]=as.numeric(simu_out_alsp_s2_samsz_2_N_max_1[i,9])/sampsz_mu[i]
  }

  simu_out_alsp_s2_samsz_2_N_max_1=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_alsp_s2_samsz_2_N_max_1,samsz_ratio_alsp_s2_samsz_2_N_max_1)))
  names(simu_out_alsp_s2_samsz_2_N_max_1)=simu_result_title
  ##############
  ## N_max_2 ###
  ##############

  ## samsz_1###

  ##### s1 ####
  simu_out_alsp_s1_samsz_1_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_alsp_s1_samsz_1_N_max_2= rbind(simu_out_alsp_s1_samsz_1_N_max_2,
                                            asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                               c_bry=s1_alsp_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_1_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_alsp_s1_samsz_1_N_max_2[i]=as.numeric(simu_out_alsp_s1_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_alsp_s1_samsz_1_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_alsp_s1_samsz_1_N_max_2,samsz_ratio_alsp_s1_samsz_1_N_max_2)))
  names(simu_out_alsp_s1_samsz_1_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_1_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_alsp_s2_samsz_1_N_max_2= rbind(simu_out_alsp_s2_samsz_1_N_max_2,
                                            asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                               c_bry=s2_alsp_cbry,samsz=samsz_1,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_1_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_alsp_s2_samsz_1_N_max_2[i]=as.numeric(simu_out_alsp_s2_samsz_1_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_alsp_s2_samsz_1_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_alsp_s2_samsz_1_N_max_2,samsz_ratio_alsp_s2_samsz_1_N_max_2)))
  names(simu_out_alsp_s2_samsz_1_N_max_2)=simu_result_title

  ## samsz_2###
  ##### s1 ####
  simu_out_alsp_s1_samsz_2_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_alsp_s1_samsz_2_N_max_2= rbind(simu_out_alsp_s1_samsz_2_N_max_2,
                                            asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s1,
                                                                               c_bry=s1_alsp_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s1_samsz_2_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_alsp_s1_samsz_2_N_max_2[i]=as.numeric(simu_out_alsp_s1_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_alsp_s1_samsz_2_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_alsp_s1_samsz_2_N_max_2,samsz_ratio_alsp_s1_samsz_2_N_max_2)))
  names(simu_out_alsp_s1_samsz_2_N_max_2)=simu_result_title
  ###############  s2 #######
  simu_out_alsp_s2_samsz_2_N_max_2=NULL
  for (i in 1:length(r_t)){
    simu_out_alsp_s2_samsz_2_N_max_2= rbind(simu_out_alsp_s2_samsz_2_N_max_2,
                                            asd_rej_simu_adapt_OC_neg_binomial(r_t=r_t[i],r_c=r_c,kappa=kappa,direction=direction,s=s2,
                                                                               c_bry=s2_alsp_cbry,samsz=samsz_2,N_max=N_max_2,theta_cut=theta_cut,beta=beta,sim_num=sim_num) )
  }

  samsz_ratio_alsp_s2_samsz_2_N_max_2=NULL
  for (i in (1:length(r_t))){
    samsz_ratio_alsp_s2_samsz_2_N_max_2[i]=as.numeric(simu_out_alsp_s2_samsz_2_N_max_2[i,9])/sampsz_mu[i]
  }

  simu_out_alsp_s2_samsz_2_N_max_2=data.frame((cbind(r_t,r_c,sampsz_mu,simu_out_alsp_s2_samsz_2_N_max_2,samsz_ratio_alsp_s2_samsz_2_N_max_2)))
  names(simu_out_alsp_s2_samsz_2_N_max_2)=simu_result_title

  ##############
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


  list_pocock_1=list(s1_pocock_design,
                     simu_out_pocock_s1_samsz_1_N_max_1,
                     simu_out_pocock_s1_samsz_2_N_max_1,
                     simu_out_pocock_s1_samsz_1_N_max_2,
                     simu_out_pocock_s1_samsz_2_N_max_2
  )
  names(list_pocock_1)=c("Group sequential design",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 1",
                         "Simulation results: Information fractions 1, Sample size 1, N_max 2",
                         "Simulation results: Information fractions 1, Sample size 2, N_max 2")


  list_pocock_2=list(s2_pocock_design,
                     simu_out_pocock_s2_samsz_1_N_max_1,
                     simu_out_pocock_s2_samsz_2_N_max_1,
                     simu_out_pocock_s2_samsz_1_N_max_2,
                     simu_out_pocock_s2_samsz_2_N_max_2
  )
  names(list_pocock_2)=c("Group sequential design",
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
           list_pocock_1,list_pocock_2,
           list_alsp_1,list_alsp_2)



  names(out)=c("O'Brien-Fleming boundaries, interim analysis schedule=s1",
               "O'Brien-Fleming boundaries, interim analysis schedule=s2",
               "Pocock boundaries, interim analysis schedule=s1",
               "Pocock boundaries, interim analysis schedule=s2",
               "Alpha-spending boundaries, interim analysis schedule=s1",
               "Alpha-spending boundaries, interim analysis schedule=s2"
  )
  return(out)
}
################################


#asd_simu_OC_neg_binomial_1(r_t,r_c,kappa,direction,s1,s2,alpha_2,alpha_2_alps_s1,alpha_2_alps_s2,theta_cut,samsz_1,samsz_2,N_max_1,N_max_2,beta,sim_num)
#asd_simu_OC_neg_binomial(r_t=c(0.3,0.4,0.5),r_c=0.2,kappa=1,direction=1,s1=c(1/3,2/3,1),s2=c(0.5,1),
#                         alpha_2=0.05,alpha_2_alps_s1=c(0.001,0.001,0.048),alpha_2_alps_s2=c(0.001,0.049),
#                         theta_cut=0.01,samsz_1=100,samsz_2=150,N_max_1=200,N_max_2=300,beta=0.1,sim_num=100)

