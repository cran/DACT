#install.packages("nph")
#install.packages("mdir.logrank")
#rm(list=ls())
#library(mdir.logrank)
#library(survival)
#library(nph)

##############################################################################
#### Part I.1. Sample size and power clculation: Normal distributions   ####
##############################################################################


###################################
####    supeiority trials   #######
###################################

#' Fixed sample designs: Normal distribution, Superiority trial, Sample size calculation
#'
#' @param mu_t The mean of the test arm.
#' @param mu_c The mean of the control arm.
#' @param sigma_t The standard deviation for the test arm.
#' @param sigma_c The standard deviation for the control arm.
#' @param alpha_2 The type I error for the two-sided test.
#' @param beta Type II error.
#' @param rand_ratio The randomization ratio of test arm:control arm.
#'
#' @return The output will show the sample size for the control arm and the test arm.
#'
#' @examples
#' sample_size_two_sample_normal(mu_t = 2, mu_c = 0, sigma_t = 9.79, sigma_c = 9.79,
#' alpha_2 = 0.05, beta = 0.2, rand_ratio = 1)
#'
#' @export
sample_size_two_sample_normal<-function(mu_t,mu_c,sigma_t, sigma_c,alpha_2, beta, rand_ratio){
  ### This function calculates the sample size for given mu_t,mu_c,sigma_t, sigma_c,alpha_2, beta, rand_ratio
  #mu_t is the mean of the tesrt arm
  #mu_c is the mean of the control arm
  #mu=abs(mu_t-mu-c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #sigma_c>0 is the standard deviation for control arm
  #sigma_t>0 is the standard deviation for test arm
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  r = rand_ratio
  mu = abs(mu_t - mu_c)
  if (mu<1e-4){return(NA)}

  a = stats::qnorm(1 - alpha_2 / 2)
  b = stats::qnorm(1 - beta)
  N_c=(((a+b)/mu)**2)*(sigma_c**2+sigma_t**2/r);
  N_c=ceiling(N_c)
  out=list(N_c,N_c*r)
  names(out)=c("Sample size for control arm ","Sample size for test arm ")
  return(out)
}
################################
#sample_size_two_sample_normal(mu_t,mu_c,sigma_t, sigma_c,alpha_2, beta, rand_ratio)
#sample_size_two_sample_normal(2,0,9.79,9.79,0.05, 0.2, 1)
#sample_size_two_sample_normal(2,0,9.79,9.79,0.05, 0.1, 1)
#sample_size_two_sample_normal(2,0,9.79, 9.79,0.05, 0.2, 1)
####################################
#' Fixed sample designs: Normal distribution, Superiority trial, Power calculation
#'
#' @param mu_t The mean of the test arm.
#' @param mu_c The mean of the control arm.
#' @param sigma_t The standard deviation for the test arm.
#' @param sigma_c The standard deviation for the control arm.
#' @param alpha_2 The type I error for the two-sided test.
#' @param n_control The sample size for the control arm.
#' @param rand_ratio The randomization ratio of test arm:control arm.
#'
#' @return The output will show the power.
#'
#' @examples
#' power_two_sample_normal(mu_t = 2.6, mu_c = 0.6, sigma_t = 9.43, sigma_c = 9.86,
#' alpha_2 = 0.05, n_control = 53, rand_ratio = 1)
#'
#' @export
power_two_sample_normal<-function(mu_t,mu_c,sigma_t, sigma_c,alpha_2, n_control, rand_ratio){
  ##This function calculate the power for given  mu_t,mu_c,sigma_t, sigma_c,alpha_2, n_control, rand_ratio
  #mu_t is the mean of the tesrt arm
  #mu_c is the mean of the control arm
  #mu=abs(mu_t-mu-c) is the absolute value of the mean difference
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #sigma_c>0 is the standard deviation for control arm
  #sigma_t>0 is the standard deviation for test arm
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #n_control is the sample size for the control arm
  #n_test=r*n_control
  mu=abs(mu_t-mu_c)
  r=rand_ratio
  a=stats::qnorm(1-alpha_2/2);
  d=mu*sqrt(n_control)/sqrt(sigma_c**2+sigma_t**2/r);
  p=stats::pnorm(d-a);
  return(p)
}
###################################
#power_two_sample_normal(mu_t,mu_c,sigma_t, sigma_c,alpha_2, n_control, rand_ratio)
#power_two_sample_normal(0.3,0,1, 1,0.05, 234, 1)
#power_two_sample_normal(2,0,9.79, 9.79,0.05, 377, 1)
#power_two_sample_normal(2,0,9.79, 9.79,0.05, 504, 1)

##########################################
####      non-inferiority trials   #######
##########################################
#' Fixed sample designs: Normal distribution, Non-inferiority trial, Sample size calculation
#'
#' @param mu_t The mean of the test arm.
#' @param mu_c The mean of the control arm.
#' @param sigma_t The standard deviation for the test arm.
#' @param sigma_c The standard deviation for the control arm.
#' @param margin The non-inferiority margin
#' @param direction If larger mu is better
#' @param alpha The type I error for the two-sided test.
#' @param beta Type II error.
#' @param rand_ratio The randomization ratio of test arm:control arm.
#'
#' @return The output will show the sample size for the control arm and the test arm.
#'
#' @examples
#' sample_size_two_sample_normal_NI(mu_t = 2, mu_c = 2.5, sigma_t = 9.79, sigma_c = 9.79,
#' margin = 0.2, direction = 2, alpha = 0.05, beta = 0.2, rand_ratio = 1)
#'
#' @export
sample_size_two_sample_normal_NI<-function(mu_t,mu_c,sigma_t, sigma_c,margin,direction,alpha, beta, rand_ratio){
  ### This function calculates the sample size for given mu_t,mu_c,sigma_t, sigma_c,alpha, beta, rand_ratio for non-inferiority trial
  #mu_t is the mean of the tesrt arm
  #mu_c is the mean of the control arm
  #mu=mu_t-mu-c is the differenc of the mean differences between the test arm and the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #sigma_c>0 is the standard deviation for control arm
  #sigma_t>0 is the standard deviation for test arm
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
  a=stats::qnorm(1-alpha);
  b=stats::qnorm(1-beta);

  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  if (mu<1e-4){return(NA)}

  N_c=(((a+b)/mu)**2)*(sigma_c**2+sigma_t**2/r);
  N_c=ceiling(N_c)
  out=list(N_c,N_c*r)
  names(out)=c("Sample size for control arm ","Sample size for test arm ")
  return(out)
}
################################
#sample_size_two_sample_normal_NI(mu_t,mu_c,sigma_c, sigma_t,margin,direction,alpha, beta, rand_ratio)
#sample_size_two_sample_normal_NI(0,1.6,10.25, 10.25,4.5,1,0.025, 0.1, 1)
#sample_size_two_sample_normal_NI(0,1.8,10.25, 10.25,4.5,1,0.025, 0.1, 1)
#sample_size_two_sample_normal_NI(2,2.5,9.79, 9.79,0.2,2,0.025, 0.2, 1)

####################################

#######################################
#' Fixed sample designs: Normal distribution, Non-inferiority trial, Power calculation
#'
#' @param mu_t The mean of the test arm.
#' @param mu_c The mean of the control arm.
#' @param sigma_t The standard deviation for the test arm.
#' @param sigma_c The standard deviation for the control arm.
#' @param margin The non-inferiority margin
#' @param direction If larger mu is better
#' @param alpha The type I error for the two-sided test.
#' @param n_control The sample size for the control arm.
#' @param rand_ratio The randomization ratio of test arm:control arm.
#'
#' @return The output will show the power.
#'
#' @examples
#' power_two_sample_normal_NI(mu_t = 0.3, mu_c = 2, sigma_t = 1, sigma_c = 1, margin = 0.4,
#' direction = 1, alpha = 0.025, n_control = 100, rand_ratio = 1)
#'
#' @export
power_two_sample_normal_NI<-function(mu_t,mu_c,sigma_t, sigma_c,margin,direction,alpha,n_control, rand_ratio){
  ### This function calculates the power for given mu_t,mu_c,sigma_t, sigma_c,alpha, n_control, rand_ratio for non-inferiority trial
  #mu_t is the mean of the tesrt arm
  #mu_c is the mean of the control arm
  #mu=mu_t-mu-c is the differenc of the mean differences between the test arm and the control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #sigma_c>0 is the standard deviation for control arm
  #sigma_t>0 is the standard deviation for test arm
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  #n_control is the sample size for the control arm
  #n_test=r*n_control
  #direction=1 if larger mu is better
  #direction ne 1 if smaller mu is better
  r=rand_ratio
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  if (mu<=0){return(NA)}

  a=stats::qnorm(1-alpha);
  d=mu*sqrt(n_control)/sqrt(sigma_c**2+sigma_t**2/r);
  p=stats::pnorm(d-a);

  return(p)
}
###################################
#power_two_sample_normal_NI(mu_t,mu_c,sigma_c, sigma_t,margin,direction,alpha,n_control,rand_ratio)
#power_two_sample_normal_NI(0.3,2,1, 1,0.4,1,0.025, 100, 1)
#power_two_sample_normal_NI(0,1.6,10.25, 10.25,4.5,1,0.025, 263, 1)


######################################################################################
####     Part I.2. Sample size and power calculations: binary distributions    #######
######################################################################################

#####################################
####     superiority trials   #######
#####################################

#' Fixed sample designs: Binary distribution, Superiority trial, Sample size calculation
#'
#' @param p_t The event rate of the test arm
#' @param p_c The event rate of the control arm
#' @param alpha_2 One-sided type I error rate used for boundary calibration
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Randomization ratio is the ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' sample_size_binary_diff(p_t = 0.2, p_c = 0.4, alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
sample_size_binary_diff<-function(p_t,p_c,alpha_2,beta,rand_ratio){
  ### This function calculates the sample size for given p_t,p_c,alpha, beta, rand_ratio
  #p_t is the event rate for the experimental arm
  # 0<p_t<1
  #p_c is the event rate for the control arm
  # 0<p_c<1
  #rand_ratio is the randomization ratio of test arm:control arm
  # rand_ratio>0
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #0<beta<1
  #1-beta is power
  r=rand_ratio
  a=stats::qnorm(1-alpha_2/2)
  b=stats::qnorm(1-beta)
  mu=abs(p_t-p_c)
  if (mu<1e-4){return(NA)}

  sigma=sqrt(p_t*(1-p_t)/r+p_c*(1-p_c))
  n_c=((a+b)/mu)^2*sigma^2
  out=list(ceiling(n_c),r*ceiling(n_c))
  names(out)=c("Sample size for control arm","sample size for experimental arm")
  return(out)
}
######################
#sample_size_binaryl_diff<-function(p_t,p_c,alpha_2,beta,rand_ratio)
#sample_size_binary_diff(0.2,0.4,0.05,0.1,1)

############################################
#' Fixed sample designs: Binary distribution, Superiority trial, Power calculation
#' @param p_t The event rate of the test arm
#' @param p_c The event rate of the control arm
#' @param alpha_2 One-sided type I error rate used for boundary calibration
#' @param n_control Sample size of the control arm
#' @param rand_ratio Randomization ratio is the ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the power
#'
#' @examples
#' power_binary_diff(p_t = 0.4, p_c = 0.3, alpha_2 = 0.05, n_control = 53, rand_ratio = 1)
#'
#' @export
power_binary_diff<-function(p_t,p_c,alpha_2,n_control,rand_ratio){
  ##This function calculate the power for given p_t,p_c,alpha, n_control, rand_ratio
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #p_t is the event rate for the experimental arm
  # 0<p_t<1
  #p_c is the event rate for the control arm
  # 0<p_c<1
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #n_control is the samples size per arm
  r=rand_ratio
  mu=abs(p_t-p_c)
  sigma_t=sqrt(p_t*(1-p_t))
  sigma_c=sqrt(p_c*(1-p_c))
  a=stats::qnorm(1-alpha_2/2);
  sigma=sqrt(p_t*(1-p_t)/r+p_c*(1-p_c))
  d=mu*sqrt(n_control)/sigma;
  p=stats::pnorm(d-a);
  return(p)
}
####################################
#power_binary_diff(p_t,p_c,alpha_2, n_control, rand_ratio)
#power_binary_diff(0.2,0.4,0.05,106,1)

##########################################
##########################################
####      non-inferiority trials   #######
##########################################

#' Fixed sample designs: Binary distribution, Non-inferiority trial, Sample size calculation
#'
#' @param p_t The event rate of the test arm
#' @param p_c The event rate of the control arm
#' @param margin the non-inferiority margin
#' @param direction Test direction: 1 = larger p_t,p_c are better, 0 = smaller p_t,p_c are better
#' @param alpha The type I error for the one-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Randomization ratio is the ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' sample_size_binary_diff_NI(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0, alpha = 0.025,
#' beta = 0.1, rand_ratio = 1)
#'
#' @export
sample_size_binary_diff_NI<-function(p_t,p_c,margin,direction,alpha,beta,rand_ratio){
  ### This function calculates the sample size for a non-inferiroity trial,given p_t,p_c,alpha, beta, rand_ratio for non-inferiority trial
  #p_t is the event rate for the experimental arm
  # 0<p_t<1
  #p_c is the event rate for the control arm
  # 0<p_c<1
  #rand_ratio is the randomization ratio of test arm:control arm
  # rand_ratio >0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #0<beta<1
  #1-beta is power
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger p_t,p_c are better
  #direction ne 1 if smaller p_t,p_c are better
  r=rand_ratio
  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  if (mu<1e-4){return(NA)}

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)
  sigma=sqrt(p_t*(1-p_t)/r+p_c*(1-p_c))
  n_c=((a+b)/mu)^2*sigma^2
  out=list(ceiling(n_c),r*ceiling(n_c))
  names(out)=c("Sample size for control arm","sample size for experimental arm")
  return(out)
}
######################
#sample_size_binary_diff_NI(p_t,p_c,margin,direction,alpha,beta,rand_ratio)
#sample_size_binary_diff_NI(0.4,0.2,0.03,1,0.025,0.1,1)
#sample_size_binary_diff_NI(0.2,0.5,0.5,2,0.025,0.1,1)
#sample_size_binary_diff_NI(0.2,0.6,0.5,1,0.025,0.1,1)
#sample_size_binary_diff_NI(0.3,0.6,0.5,2,0.025,0.1,1)
#sample_size_binary_diff_NI(0.3,0.6,0.5,1,0.025,0.1,1)
############################################

#' Fixed sample designs: Binary distribution, Non-inferiority trial, Power calculation
#'
#' @param p_t The event rate of the test arm
#' @param p_c The event rate of the control arm
#' @param margin the non-inferiority margin
#' @param direction Test direction: 1 = larger p_t,p_c are better, 0 = smaller p_t,p_c are better
#' @param alpha The type I error for the one-sided test
#' @param n_control Sample size of the control arm
#' @param rand_ratio Randomization ratio is the ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the power
#'
#' @examples
#' power_binary_diff_NI(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0, alpha = 0.025,
#' n_control = 53, rand_ratio = 1)
#'
#' @export
power_binary_diff_NI<-function(p_t,p_c,margin,direction, alpha, n_control, rand_ratio){
  ##This function calculate the power for given p_t,p_c,alpha, n_control, rand_ratio for non-inferiority trial
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #p_t is the event rate for the experimental arm
  # 0<p_t<1
  #p_c is the event rate for the control arm
  # 0<p_c<1
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #n_control is the samples size for the control arm
  #margin is the non-inferiority margin
  #margin>0
  #direction=1 if larger p_t,p_c are better
  #direction ne 1 if smaller p_t,p_c are better
  r=rand_ratio

  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }

  sigma_t=sqrt(p_t*(1-p_t))
  sigma_c=sqrt(p_c*(1-p_c))
  a=stats::qnorm(1-alpha);
  sigma=sqrt(p_t*(1-p_t)/r+p_c*(1-p_c))
  d=mu*sqrt(n_control)/sigma;
  p=stats::pnorm(d-a);
  return(p)
}
####################################
#power_binary_diff_NI(p_t,p_c,margin,direction, alpha, n_control, rand_ratio){
#power_binary_diff_NI(0.4,0.2,0.03,1,0.025,80,1)
#power_binary_diff_NI(0.2,0.4,0.03,2,0.025,80,1)


###########################################################
####    Part I.3. sample size and power calculations ######
##  survival analysis -- exponential distributions   ######
####          proportional hazards                   ######
###########################################################

#####################################
####     superiority trials   #######
#####################################

#################################################
#' Fixed sample designs: Survival analysis, Superiority trial, Sample size calculation
#'
#' @param HR The hazards ratio of test arm vs. the control arm
#' @param alpha_2 The type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the total number of events for the control arm and the test arm
#'
#' @examples
#' sample_size_two_sample_surv_Schoenfeld(HR = 0.7, alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
sample_size_two_sample_surv_Schoenfeld<-function(HR,alpha_2,beta,rand_ratio){
  ### This function calculates the sample size for given sigma_x, sigma_y,mu,alpha, beta
  # The formula is from John M. Lachin (1981)
  #h_r is the hazrds ratio
  #h_r>0
  #alpha is two-sided type I error
  #beta -s type II error
  #1-beta is power

  a=stats::qnorm(1-alpha_2/2);
  b=stats::qnorm(1-beta);
  mu=abs(log(HR))

  if (mu<1e-4){return(NA)}

  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  sigma=1/sqrt((p_c*p_t))
  N=((a+b)/mu)^2*sigma^2
  N=ceiling(N)
  out=list(N, " Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. ")
  names(out)=c("Total number of events","Reference")
  return(out)
}
##############################################
#sample_size_two_sample_surv_Schoenfeld(HR,alpha_2,beta, rand_ratio)
#sample_size_two_sample_surv_Schoenfeld(0.8,0.05,0.1,1)
#sample_size_two_sample_surv_Schoenfeld(0.6,0.05,0.1,1)
#sample_size_two_sample_surv_Schoenfeld(0.667,0.1,0.2,1)
#sample_size_two_sample_surv_Schoenfeld(0.7,0.1,0.2,1)
#sample_size_two_sample_surv_Schoenfeld(HTR=0.7,alpha_2=0.05,beta=0.1,rand_ratio=1)
################################################################

#' Fixed sample designs: Survival analysis, Superiority trial, Power calculation
#'
#' @param HR The hazards ratio of test arm vs. the control arm
#' @param alpha_2 The type I error for the two-sided test
#' @param event_tot Total number of events
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the power
#'
#' @examples
#' power_two_sample_surv_Schoenfeld(HR = 0.7, alpha_2 = 0.05, event_tot = 100, rand_ratio = 1)
#'
#' @export
power_two_sample_surv_Schoenfeld<-function(HR,alpha_2,event_tot,rand_ratio){
  ### This function calculates the power for given hazards ratio,, alpha_2, n_control and rand_ratio
  #HR is the hazards ratio of test arm:control arm
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #n_control is the samples size in the control arm
  a=stats::qnorm(1-alpha_2/2);
  mu=abs(log(HR))
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  sigma=1/sqrt(p_c*p_t)
  d=mu*sqrt(event_tot)/sigma;
  p=stats::pnorm(d-a);
  ref=" Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(p,ref )
  names(out)=c("Power","Reference")
  return(out)
}
####################################
#power_two_sample_surv_Schoenfeld(HR,alpha_2,event_tot,rand_ratio)
#power_two_sample_surv_Schoenfeld(0.8,0.05,845,1)
#power_two_sample_surv_Schoenfeld(0.6,0.05,162,1)
#power_two_sample_surv_Schoenfeld(0.667,0.2,151,1)
#power_two_sample_surv_Schoenfeld(HR=0.7,alpha_2=0.0101,166,1)
#power_two_sample_surv_Schoenfeld(0.7,0.05,276,3)
#power_two_sample_surv_Schoenfeld(0.7,0.05,138,1)
#################################################

##########################################
####      non-inferiority trials   #######
##########################################

#################################################
#' Fixed sample designs: Survival analysis, Non-inferiority trial, Sample size calculation
#'
#' @param HR The hazards ratio of test arm vs. the control arm
#' @param margin The non-inferiority margin
#' @param direction Test direction: 1 = smaller HR is better, 0 = larger HR is better
#' @param alpha_2 The type I error for the one-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the total number of events for the control arm and the test arm
#'
#' @examples
#' sample_size_two_sample_surv_Schoenfeld_NI(HR = 0.7, margin = 0.1, direction = 0,
#' alpha_2 = 0.05, beta = 0.1, rand_ratio = 1)
#'
#' @export
sample_size_two_sample_surv_Schoenfeld_NI<-function(HR,margin,direction,alpha_2,beta, rand_ratio){
  ### This function calculates the sample size for given hazards ratio, alpha, beta and rand_ratio for a non-inferiority trial
  #HR is the hazards ratio of test arm:control arm
  #margin is the non-inferiority margin
  #margin>0
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #1-beta is power
  #direction=1 if smaller HR is better
  #direction=2 if larger HR is better
  if (direction==1) {
    mu=log(HR/margin)
  } else{
    mu=-log(HR/margin)
  }
  if (mu<1e-4){return(NA)}

  a=stats::qnorm(1-alpha_2/2);
  b=stats::qnorm(1-beta);
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  sigma=1/sqrt((p_c*p_t))
  N=((a+b)/mu)^2*sigma^2
  N=ceiling(N)
  ref=" Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(N,ref )
  names(out)=c("Total number of events","Reference")
  return(out)
}
##############################################
#sample_size_two_sample_surv_Schoenfeld_NI(HR,margin,direction,alpha,beta, rand_ratio)
#sample_size_two_sample_surv_Schoenfeld_NI(0.9,1.2,2,0.025,0.1,1)
#sample_size_two_sample_surv_Schoenfeld_NI(0.9,0.85,1,0.025,0.1,1)
#sample_size_two_sample_surv_Schoenfeld_NI(1.1,0.85,1,0.025,0.1,1)
################################################################

#' Fixed sample designs: Survival analysis, Non-inferiority trial, Power calculation
#'
#' @param HR The hazards ratio of test arm vs. the control arm
#' @param margin The non-inferiority margin
#' @param direction Test direction: 1 = smaller HR is better, 0 = larger HR is better
#' @param alpha The type I error for the one-sided test
#' @param event_tot Total number of events
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the power
#'
#' @examples
#' power_two_sample_surv_Schoenfeld_NI(HR = 0.7, margin = 0.1, direction = 0, alpha = 0.025,
#' event_tot = 100, rand_ratio = 1)
#'
#' @export
power_two_sample_surv_Schoenfeld_NI<-function(HR,margin,direction,alpha,event_tot, rand_ratio){
  ### This function calculates the power for given hazards ratio, alpha, n_control and rand_ratio for a non-inferiority trial
  #HR is the hazards ratio of test arm:control arm
  #margin is the non-inferiority margin
  #margin>0
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #direction=1 if smaller HR is better
  #direction=2 if larger HR is better
  #n_control is the samples size in the control arm
  if (direction==1) {
    mu=log(HR/margin)
  } else{
    mu=-log(HR/margin)
  }

  a=stats::qnorm(1-alpha/2)
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  sigma=1/sqrt(p_c*p_t)
  d=mu*sqrt(event_tot)/sigma;
  p=stats::pnorm(d-a);
  ref=" Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(p,ref )
  names(out)=c("Power","Reference")
  return(out)
}
####################################
#power_two_sample_surv_Schoenfeld_NI(HR,margin,direction,alpha,event_tot, rand_ratio)
#power_two_sample_surv_Schoenfeld_NI(0.9,1.2,1,0.025,600,1)
#power_two_sample_surv_Schoenfeld_NI(0.9,1.2,2,0.025,600,1)
#power_two_sample_surv_Schoenfeld_NI(1.1,1.2,2,0.025,347,1)
#power_two_sample_surv_Schoenfeld_NI(0.9,1.2,2,0.025,600,1)
#power_two_sample_surv_Schoenfeld_NI(1.1,0.85,1,0.025,747,1)

#########################################################################
###############  Part I.4. sample size and power calculation    #########
#########    two sample Poisson distribution -superiority trials     ####
#########################################################################

#####################################
####     superiority trials   #######
#####################################

#################################################
#' Fixed sample designs: Poisson distribution, Superiority trial, Sample size calculation
#'
#' @param lambda_t The event rate of the test arm
#' @param lambda_c The event rate of the control arm
#' @param alpha_2 The type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' sample_size_two_sample_Poisson(lambda_t = 1.1, lambda_c = 1.0, alpha_2 = 0.05,
#' beta = 0.1, rand_ratio = 1)
#'
#' @export
sample_size_two_sample_Poisson<-function(lambda_t,lambda_c,alpha_2, beta, rand_ratio){
  ### This function calculates the sample size for given lambda_t,lambda_c,alpha, beta, rand_ratio
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #1-beta is power
  #lambda_t is the mean event rates for arm t
  #lambda_c is the mean event rates for arm t
  #beta -s type II error
  #1-beta is power

  r=rand_ratio
  a=stats::qnorm(1-alpha_2/2);
  b=stats::qnorm(1-beta);
  mu=abs(lambda_t-lambda_c)
  if (mu<1e-4){return(NA)}

  sigma=sqrt(lambda_t/r+lambda_c)
  n_c=((a+b)/mu)^2*sigma^2
  n_c=ceiling(n_c)
  n_t=r*n_c
  out=list(n_t,n_c)
  names(out)=c("Sample size for experimental arm","Sample size for control arm")
  return(out)
}
################################
#sample_size_two_sample_Poisson(lambda_t, lambda_c, alpha_2,beta, rand_ratio)
#sample_size_two_sample_Poisson(0.8,1.2,0.05,0.1,1)


####################################
#' Fixed sample designs: Poisson distribution, Superiority trial, Power calculation
#'
#' @param lambda_t The event rate of the test arm
#' @param lambda_c The event rate of the control arm
#' @param alpha_2 The type I error for the two-sided test
#' @param n_control Sample size of the control arm
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the power
#'
#' @examples
#' power_two_sample_Poisson(lambda_t = 1.1, lambda_c = 1.0, alpha_2 = 0.05,
#' n_control = 53, rand_ratio = 1)
#'
#' @export
power_two_sample_Poisson<-function(lambda_t, lambda_c,alpha_2, n_control, rand_ratio){
  ### This function calculates the powerfor given lambda_t, lambda_c,alpha_2, n_control, rand_ratio
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #1-beta is power
  #lambda_t is the mean event rates for arm t
  #lambda_c is the mean event rates for arm t
  #n_control is the samples size for control arm

  r=rand_ratio
  a=stats::qnorm(1-alpha_2/2);
  mu=abs(lambda_t-lambda_c)
  sigma=sqrt(lambda_t/r+lambda_c)
  d=mu*sqrt(n_control)/sigma;
  p=stats::pnorm(d-a);
  return(p)
}
####################################
#power_two_sample_Poisson(lambda_t, lambda_c,alpha_2, n_control, rand_ratio)
#power_two_sample_Poisson(9,11,0.05,40,1)
#power_two_sample_Poisson(0.8,1.2,0.05,132,1)

######################################################################################

##########################################
####      non-inferiority trials   #######
##########################################

#################################################
#' Fixed sample designs: Poisson distribution, Non-inferiority trial, Sample size calculation
#'
#' @param lambda_t The event rate of the test arm
#' @param lambda_c The event rate of the control arm
#' @param margin The non-inferiority margin
#' @param direction Test direction: 1 = smaller lambda is better, 0 = larger lambda is better
#' @param alpha The type I error for the one-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' sample_size_two_sample_Poisson_NI(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1, direction = 0,
#' alpha = 0.025, beta = 0.1, rand_ratio = 1)
#'
#' @export
sample_size_two_sample_Poisson_NI<-function(lambda_t,lambda_c,margin,direction,alpha, beta, rand_ratio){
  ### This function calculates the sample size for given lambda_t,lambda_c,margin,direction,alpha, beta, rand_ratio for a non-inferiority trial
  #margin is the non-inferiority margin
  #margin>0
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #1-beta is power
  #direction=1 if smaller lambda is better
  #direction=2 if larger lambda is better
  #lambda_t is the mean event rates for test arm
  #lambda_c is the mean event rates for control arm

  r=rand_ratio

  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  if (mu<1e-4){return(NA)}

  a=stats::qnorm(1-alpha)
  b=stats::qnorm(1-beta)

  n_c=(((a+b)/mu)**2)*(lambda_t/r+lambda_c)
  n_c=ceiling(n_c)
  n_t=r*n_c
  out=list(n_t,n_c)
  names(out)=c("Sample size for experimental arm","Sample size for control arm")
  return(out)
}
################################
#sample_size_two_sample_Poisson_NI(lambda_t, lambda_c, margin,direction,alpha,beta, rand_ratio)
#sample_size_two_sample_Poisson_NI(9,11,2.5,1,0.025,0.1,1)
#sample_size_two_sample_Poisson_NI(9,11,2.5,2,0.025,0.1,1)

####################################
#' Fixed sample designs: Poisson distribution, Non-inferiority trial, Power calculation
#'
#' @param lambda_t The event rate of the test arm
#' @param lambda_c The event rate of the control arm
#' @param margin The non-inferiority margin
#' @param direction Test direction: 1 = smaller lambda is better, 0 = larger lambda is better
#' @param alpha The type I error for the one-sided test
#' @param n_control Sample size of the control arm
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#'
#' @return the output will show the power
#'
#' @examples
#' power_two_sample_Poisson_NI(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1, direction = 0,
#' alpha = 0.025, n_control = 53, rand_ratio = 1)
#'
#' @export
power_two_sample_Poisson_NI<-function(lambda_t, lambda_c,margin,direction, alpha, n_control, rand_ratio){
  ### This function calculates the power for given lambda_t,lambda_c,margin,direction,alpha, n_control, rand_ratio for a non-inferiority trial
  #margin is the non-inferiority margin
  #margin>0
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #beta is type II error
  #1-beta is power
  #direction=1 if smaller lambda is better
  #direction=2 if larger lambda is better
  #lambda_t is the mean event rates for test arm
  #lambda_c is the mean event rates for control arm
  #n_control is the samples size for control arm

  r=rand_ratio
  if (direction==1){
    mu=lambda_t-lambda_c+margin
  }else{
    mu=-(lambda_t-lambda_c-margin)
  }
  a=stats::qnorm(1-alpha);

  sigma=sqrt(lambda_t/r+lambda_c)
  d=mu*sqrt(n_control)/sigma;
  p=stats::pnorm(d-a);
  return(p)
}
####################################
#power_two_sample_Poisson_NI(lambda_t, lambda_c,margin,direction,alpha, n, r)
#power_two_sample_Poisson_NI(9,11,2.5,1,0.025,841,1)
#power_two_sample_Poisson_NI(9,11,2.5,2,0.05,11,1)


############################################################################
####     Part I.5. Two sample Negative binomial distribtions         #######
#####        sample size and power calculations                      #######
############################################################################

###################################
####    supeiority trials   #######
###################################

#' Fixed sample designs: Negative binomial distribution, Superiority trial, Sample size calculation
#'
#' @param r_t The event rate of the test arm
#' @param r_c The event rate of the control arm
#' @param alpha_2 The type I error for the two-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t The exposure time
#' @param kappa The dispersion parameter
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' sample_size_neg_binomial(r_t = 0.85, r_c = 1.25, alpha_2 = 0.05, beta = 0.2, rand_ratio = 1,
#' nu_t = 1, kappa = 1)
#'
#' @export
sample_size_neg_binomial<-function(r_t,r_c,alpha_2,beta,rand_ratio,nu_t,kappa){
  ### This function calculates the sample size for given r_t,r_c,alpha_2,beta,rand_ratio,nu_t,kappa
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #1-beta is power
  #r_c is the event rate in the control arm
  #r_c>0
  #r_t is the event rate in the test/experimental arm arm
  #r_t>0
  #nu_t is the exposure time
  #nu_t>0
  #kappa is the dispersion parameter
  #kappa>0
  r_0=r_c
  r_1=r_t
  theta=rand_ratio
  a=stats::qnorm(1-alpha_2/2);
  b=stats::qnorm(1-beta);
  mu=abs(log(r_1/r_0))
  if (mu<1e-4){return(NA)}

  #  n_c=(((a+b)/mu)**2)*(lambda_t/r+lambda_c);

  V<-(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  n0_c=ceiling(((a+b)/mu)^2*V)
  n0_t=ceiling(rand_ratio*n0_c)
  mu_stzd=mu/sqrt(V)
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(mu_stzd,n0_t,n0_c,ref)
  names(out)=c("Standardized effect size","Sample size for test arm","Sample size for control arm","Reference")
  return( out)
}
####################################
#sample_size_neg_binomial(r_t,r_c,alpha_2,beta,rand_ratio,nu_t,kappa)
#sample_size_neg_binomial(r_t=0.85,r_c=1.25,alpha_2=0.05,beta=0.2,rand_ratio=1,nu_t=1,kappa=1)
#sample_size_neg_binomial(r_t=0.85,r_c=1.25,alpha_2=0.05,beta=0.2,rand_ratio=2,nu_t=1,kappa=1)

####################################################
#' Fixed sample designs: Negative binomial distribution, Superiority trial, Power calculation
#'
#' @param r_t The event rate of the test arm
#' @param r_c The event rate of the control arm
#' @param alpha_2 The type I error for the two-sided test
#' @param n_control Sample size of the control arm
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t The exposure time
#' @param kappa The dispersion parameter
#'
#' @return The output will show the power
#'
#' @examples
#' power_two_sample_neg_binomial(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, alpha_2 = 0.05,
#' n_control = 53, rand_ratio = 1, nu_t = 6, kappa = 1)
#'
#' @export
power_two_sample_neg_binomial<-function(r_t,r_c,alpha_2,n_control,rand_ratio,nu_t,kappa){
  ### This function calculates the sample size for given r_t,r_c,alpha_2,beta,rand_ratio,nu_t,kappa
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha_2 is two-sided type I error
  #0<alpha_2<1
  #typically, alpha_2=0.05
  #beta is type II error
  #1-beta is power
  #r_c is the event rate in the control arm
  #r_c>0
  #r_t is the event rate in the test/experimental arm arm
  #r_t>0
  #nu_t is the exposure time
  #nu_t>0
  #kappa is the dispersion parameter
  #kappa>0
  r_0=r_c
  r_1=r_t
  theta=rand_ratio
  mu=abs(log(r_1/r_0))
  a=stats::qnorm(1-alpha_2/2);
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(V)
  d=mu*sqrt(n_control)/sigma;
  p=stats::pnorm(d-a);
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(p,ref)
  names(out)=c( "power","Reference")
  return(out)
}
####################################
#power_two_sample_neg_binomial(r_t,r_c,alpha_2,n_control,rand_ratio,nu_t,kappa)
#power_two_sample_neg_binomial(r_t=0.85,r_c=1.25,alpha_2=0.05,n_control=153,rand_ratio=2,nu_t=1,kappa=1)

############################################################################

#########################################
####     non-inferiority trials   #######
#########################################

##############################################################################
#' Fixed sample designs: Negative binomial distribution, Non-inferiority trial, Sample size calculation
#'
#' @param r_t The event rate of the test arm
#' @param r_c The event rate of the control arm
#' @param margin The non-inferiority margin
#' @param direction Test direction: 1 = larger r_t,r_c are better, 0 = smaller r_t,r_c are better
#' @param alpha The type I error for the one-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t The exposure time
#' @param kappa The dispersion parameter
#'
#' @return The output will show the sample size for the control arm and the test arm
#'
#' @examples
#' sample_size_neg_binomial_NI(r_t = 0.9, r_c = 1, margin = 0.3, direction = 1, alpha = 0.025,
#' beta = 0.1, rand_ratio = 1, nu_t = 1, kappa = 1)
#'
#' @export
sample_size_neg_binomial_NI<-function(r_t,r_c,margin,direction,alpha,beta,rand_ratio,nu_t,kappa){
  ##This function calculates the sample size  for given r_t,r_c,margin,direction,alpha,beta,rand_ratio,nu_t,kappa fro a non-inferiority trial
  #r_c is the event rate in the control arm
  #r_c>0
  #r_t is the event rate in the test/experimental arm arm
  #r_t>0
  #rand_ratio is the randomization ratio of test arm:control arm
  #theta>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #nu_t is the exposure time.
  #nu_t>0
  #kappa is the dispersion parameter
  #kappa>0
  #margin is the non-inferiroity margin
  #margin>0
  #direction=1 if larger r_t,r_c are better
  #direction ne 1 if smaller r_t,r_c are better
  r=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
    }
  if (mu<1e-4){return(NA)}

  r_0=r_c
  theta=rand_ratio
  a=stats::qnorm(1-alpha);
  b=stats::qnorm(1-beta);

  # n_c=(((a+b)/mu)**2)*(lambda_t/r+lambda_c);

  V<-(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  n0_c=ceiling(((a+b)/mu)^2*V)
  n0_t=ceiling(rand_ratio*n0_c)
  mu_stzd=mu/sqrt(V)
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(mu_stzd,n0_t,n0_c,ref)
  names(out)=c("Standardized effect size","Sample size for test arm","Sample size for control arm","Reference")
  return( out)
}
###################################################################
#sample_size_neg_binomial_NI(r_t,r_c,margin,direction,alpha,beta,r_ratio,nu_t,kappa){
#sample_size_neg_binomial_NI(0.9,1,0.3,1,0.025,0.1,1,1,1)
#sample_size_neg_binomial_NI(1,0.8,0.3,1,0.05,0.1,1,1,1)
#sample_size_neg_binomial_NI(1.1,1,0.3,1,0.05,0.1,1,1,1)
#sample_size_neg_binomial_NI(1.1,1,0.3,2,0.05,0.1,1,1,1)

####################################

#' Fixed sample designs: Negative binomial distribution, Non-inferiority trial, Power calculation
#'
#' @param r_t The event rate of the test arm
#' @param r_c The event rate of the control arm
#' @param margin The non-inferiority margin
#' @param direction Test direction: 1 = larger r_t,r_c are better, 0 = smaller r_t,r_c are better
#' @param alpha The type I error for the one-sided test
#' @param n_control Sample size of the control arm
#' @param rand_ratio The ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t The exposure time
#' @param kappa The dispersion parameter
#'
#' @return The output will show the power
#'
#' @examples
#' power_two_sample_neg_binomial_NI(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, margin = 0.1, direction = 0,
#' alpha = 0.025, n_control = 53, rand_ratio = 1, nu_t = 6, kappa = 1)
#'
#' @export
power_two_sample_neg_binomial_NI<-function(r_t,r_c,margin,direction,alpha,n_control,rand_ratio,nu_t,kappa){
  ##This function calculates the power  for given r_t,r_c,margin,direction,alpha,n_control,rand_ratio,nu_t,kappa fro a non-inferiority trial
  #r_c is the event rate in the control arm
  #r_c>0
  #r_t is the event rate in the test/experimental arm arm
  #r_t>0
  #rand_ratio is the randomization ratio of test arm:control arm
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #nu_t is the exposure time.
  #nu_t>0
  #kappa is the dispersion parameter
  #kappa>0
  #margin is the non-inferiroity margin
  #margin>0
  #direction=1 if larger r_t,r_c are better
  #direction ne 1 if smaller r_t,r_c are better
  #n_control is the sample size for the control arm
  r=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  r_0=r_c
  theta=rand_ratio
  a=stats::qnorm(1-alpha);
   V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(V)
  d=mu*sqrt(n_control)/sigma;
  p=stats::pnorm(d-a);
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(p,ref)
  names(out)=c( "power","Reference")
  return(out)
}
####################################
#power_two_sample_neg_binomial_NI(r_t,r_c,margin,direction,alpha,n_control,rand_ratio,nu_t,kappa){
#power_two_sample_neg_binomial_NI(1,0.8,0.3,1,0.05,147,1,1,1)

