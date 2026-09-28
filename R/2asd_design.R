#rm(list=ls())
#library(cubature)
#library(foreach)
#library(parallel)
#library(doParallel)
#library(bigstatsr)

integTol = 1e-5
unirTol = 1e-6

###################################################
##########   Part I                    ############
###  Critical bundary calculations     ############
###################################################
###########################################
## Calculates the O'Brien-Fleming type  boundaries
##It also calculate the alpha-spending  at each interim analysis###
## length(s)=<4 ###
## In practice,  length(s)=3 may be the most common###

#########################################
two_stage_OF_bundry<-function(s,alpha,comp,dist){
  #this function selects O'Brien-Fleming type critical boundaries for multiple dose or endpoint comparison#
  #s is the vector of information fractions. The length of s can be any integer#
  #alpha is th eone-sided type I errro level
  #dist is distribution. Normal distirbution ; dist-=1; Nonnomral distributions: dist =2
  #M=comp is the number of comparisons, it can be for dose, or endpoints#
  #### the next section constructs the covariance matrix###
  s_0=s[1] #information fraction at the end of phase 2 stage
  r=(s[-1]-s[1])
  #information fraction conditonal on phase 2 sthage
  s_p3=(s[-1]-s[1])/(last(s)-s[1])
  ## s_1 is the scaled information fraction for selection the ##
  ## critical boundary for phase 3 stage, this boundary will be only used for point estiante and upper limit of confiedence intgerval##

  ########## non-normal distribution- phase 2 stage matrix###
  cov_mtx_0_nonn=diag(comp)
  cov_mtx_0_nonn=s_0*cov_mtx_0_nonn

  ##############non-normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_nonn=diag(comp-1)
  cov_mtx_0c_nonn=s_0*cov_mtx_01_nonn

  ########## normal distribution- phase 2 stage matrix###
  cov_mtx_0_n=diag(comp)
  cov_mtx_0_n[lower.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n[upper.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n=s_0*cov_mtx_0_n

  ##############normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_n=diag(comp-1)
  cov_mtx_01_n[lower.tri(cov_mtx_01_n)]=0.5
  cov_mtx_01_n[upper.tri(cov_mtx_01_n)]=0.5

  a_mtx_01=(1/4)*matrix(rep(1,(comp-1)*(comp-1)),(comp-1),(comp-1))

  con_mtx_01_n=cov_mtx_01_n-a_mtx_01
  cov_mtx_0c_n=s_0*con_mtx_01_n

  ###### Matirx for phase 3 stage #######

  cov_mtx_1= matrix(rep(r,length(r)),nrow=length(r), ncol=length(r))
  cov_mtx_1[lower.tri(cov_mtx_1)]=t(cov_mtx_1)[lower.tri(cov_mtx_1)]

  #######################Civariance matrix for phase 3 ##########
  cov_mtx_2= matrix(rep(s_p3,length(s_p3)),nrow=length(s_p3), ncol=length(s_p3))
  cov_mtx_2[lower.tri(cov_mtx_2)]=t(cov_mtx_2)[lower.tri(cov_mtx_2)]

  ##############################
  ### Normal distribution ####
  ##############################
  tp1 <-function(e_0,s_0,r,sigma0,sigma1,sigma2,comp){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(1/2*x,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
         (1-pmvnorm(upper=rep((e_0-x),length(r)),lower=-Inf,mean=rep(0,length(r)),sigma=sigma2,algorithm=Miwa())[1])
      return(rs)
    }

    rs= 1-pmvnorm(lower=-Inf, upper=rep(e_0,comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]+
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper = e_0,subdivisions = 1000L,abs.tol = integTol)$value
    # rs here is the exit probability (P(E_m) in the manuscript)
    return(rs)
  }
  ##############################
  ### non-Normal distribution ####
  ##############################
  tp2 <-function(e_0,s_0,r,sigma0,sigma1,sigma2,comp){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(0,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
         (1-pmvnorm(upper=rep((e_0-x),length(r)),lower=-Inf,mean=rep(0,length(r)),sigma=sigma2,algorithm=Miwa())[1])
      return(rs)
    }

    rs= 1-pmvnorm(lower=-Inf, upper=rep(e_0,comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]+
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper = e_0,subdivisions = 1000L,abs.tol = integTol)$value
    # rs here is the exit probability (P(E_m) in the manuscript)
    return(rs)
  }
  ###################################
  tp3 <-function(e_0,s_p3,sigma3){
    #This function is used to calculate the critical bpundary##
    # for phase 3 stage as if it was an independent trial with information fraction s_p3##

    rs= 1-pmvnorm(lower=-Inf,upper=rep(e_0,length(s_p3)),mean=rep(0,length(s_p3)),sigma=sigma3,algorithm=Miwa())[1]
    return(rs)
  }
  #################################
  #Normal distribution, off-diagnal elements are known to be 1/2
  ####################################
  if (dist==1){
    e_bry_two_stg=as.numeric(stats::uniroot(function(x) tp1(x,s_0,r,cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_1,comp)-alpha,c(-2^32,2^32),extendInt="yes",tol = unirTol)[1])

  }else{
    e_bry_two_stg=as.numeric(stats::uniroot(function(x) tp2(x,s_0,r,cov_mtx_0_nonn, cov_mtx_0c_nonn,cov_mtx_1,comp)-alpha,c(-2^32,2^32),extendInt="yes",tol = unirTol)[1])
    #non-Normal distribution, off-diagnal elements are unknown and taken to be zero
  }

  c_bry_two_stg= e_bry_two_stg/sqrt(s)
  ############################################
  ## calculating alpha spending with the selected boundaries######
  #### alpha spend at first interim analysis####
  tp01 <-function( e_bry_two_stg,sigma0){
    rs= 1-pmvnorm(lower=-Inf, upper=rep( e_bry_two_stg,comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]
    return(rs)
  }
  #################################################
  #############################################################
  ### Normal distributions --phase 3 stage exit probability ####
  #############################################################
  ## Calculating alpha spending with the selected boundaries#####
  #############################################################


  #############################################################
  #### cumulative alpha spending at inter_vt-th interim analysis####
  ###########################################################

  ###############################################
  tp_norm_cumu <-function( e_bry_two_stg,inter_vt,sigma0,sigma1,sigma2,comp){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2
    #   inter_vt>=2
    #e_bry_two_stg is a single number, not a vector
    rank=inter_vt
    intg <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s[1]))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(1/2*x,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
         (1-pmvnorm(upper= rep(e_bry_two_stg,(rank-1))-x,lower=-Inf,mean=rep(0,(rank-1)),sigma=sigma2,algorithm=Miwa())[1]  )
      return(rs)
    }

    rs= 1-pmvnorm(lower=-Inf, upper=rep( e_bry_two_stg,comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]+
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper = e_bry_two_stg,subdivisions = 1000L,abs.tol = integTol)$value
    # rs here is the exit probability (P(E_m) in the manuscript)
    return(rs)
  }

  #### alpha spending at second interim analysis####
  ################################################


  #############################################################
  ### Non-Normal distributions --phase 3 stage exit probability ####
  #############################################################
  ## Calculating alpha spending with the selected boundaries#####
  #############################################################
  #### alpha spending at second interim analysis####
  ################################################

  tp_nonnorm_cumu <-function( e_bry_two_stg,inter_vt,sigma0,sigma1,sigma2,comp){

    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 0
    #   inter_vt>=2
    rank=inter_vt
    intg <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s[1]))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(0,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
        (1-pmvnorm(upper= rep(e_bry_two_stg,(rank-1))-x,lower=-Inf,mean=rep(0,(rank-1)),sigma=sigma2,algorithm=Miwa())[1])

      rs=stats::dnorm(x,mean=0,sd=sqrt(s[1]))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(1/2*x,comp-1),sigma=sigma1,algorithm=Miwa())[1]*

      return(rs)
    }

    rs= 1-pmvnorm(lower=-Inf, upper=rep( e_bry_two_stg,comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]+
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper = e_bry_two_stg,subdivisions = 1000L,abs.tol = integTol)$value
    # rs here is the exit probability (P(E_m) in the manuscript)
    return(rs)
  }

  ##############################################################
  #### cumulative alpha spending at each interim analysis####
  ###############################################################
  cumu_alpha=c(rep(0,length(s)))
  if (dist==1){
    cumu_alpha[1]=tp01( e_bry_two_stg,cov_mtx_0_n)
    for (i in 2:length(s)){
      cumu_alpha[i]=tp_norm_cumu(e_bry_two_stg,i,cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_1[1:(i-1),1:(i-1)],comp)
    }
  } else {
    cumu_alpha[1]=tp01( e_bry_two_stg,cov_mtx_0_nonn)
    for (i in 2:length(s)){
      cumu_alpha[i]=tp_nonnorm_cumu(e_bry_two_stg,i,cov_mtx_0_nonn, cov_mtx_0c_nonn,cov_mtx_1[1:(i-1),1:(i-1)],comp)
    }
  }

  ############### alpha spending######
  alpha_s=c(rep(0,length(s)))
  if(length(s)==1){
    alpha_s=cumu_alpha
  } else {
    alpha_s[1]=cumu_alpha[1]
    for (i in 2:length(s)){
      alpha_s[i]=cumu_alpha[i]-cumu_alpha[i-1]
    }
  }
  #################################################
  # return(c_bry)
  e_bry_p3=as.numeric(stats::uniroot(function(x) tp3(x,s_p3,cov_mtx_2)-alpha,c(-2^32,2^32),extendInt="yes", tol = unirTol)[1])
  c_bry_p3=e_bry_p3/s_p3
  nom_alpha=1-stats::pnorm(c_bry_two_stg)
  nom_alpha_p3=1-stats::pnorm(c_bry_p3)
  #    sum_nom_alpha_p3=sum(nom_alpha_p3)

  ph3=list(s_p3,c_bry_p3,nom_alpha_p3)
  names(ph3)=c("phase 3 information fraction", "pahse 3 critical boundary","nominal p-value for phase 3")

  out=list(c_bry_two_stg,alpha_s,nom_alpha,cumu_alpha,ph3)
  names(out)=c("Critical boundary","1-sided alpha spending at each interim analysis","1-sided nominal alpha  at each interim analysis","1-sided accumulative alpha","phase 3 ")
  return(out)
}

##############################################
#two_stage_OF_bundry(s,alpha,comp,dist)
#two_stage_OF_bundry(c(1/3,2/3,1),0.025,2,2)
#two_stage_OF_bundry(c(1/3,1),0.025,3,1)
#two_stage_OF_bundry(c(1/3,1),0.025,3,2)
#two_stage_OF_bundry(c(1/3,2/3,1),0.025,3,1)
#two_stage_OF_bundry(s=c(1/3,2/3,1),alpha=0.025,comp=2,dist=1)
#two_stage_OF_bundry(c(1/3,2/3,1),0.025,3,2)
#two_stage_OF_bundry(c(1/4,1/3,2/3,1),0.025,3,1)
#two_stage_OF_bundry(c(1/4,1/3,2/3,1),0.025,3,2)
###########################################
###########################################
## Calculates the boundaries according to the alpha-spending plan.####
###alpha is a vector and length(alpha)=length(s)###
##  sum(alpha)=0.025####
## length(s)=<4 ###
## In practice,  length(s)=3 may be the most common###
###########################################

#######################################
#######################################

two_stage_alsp_bundry<-function(s,alpha,comp,dist){
  #this function selects O'Brien-Fleming type critical boundaries for multiple dose or endpoint comparison#
  #the boundary selection is based on alphs spending at eatch interm analysis#
  #s is the vector of information fractions. The length of s =<4#
  #dist is distribution. Normal distirbution ; dist-=1; Nonnomral distributions: dist =2
  #alpha is the vector of how the alpha will be spent at the interim analyses. The length of alpha is the same as that of s#
  #M=comp is the number of comparisons, it can be foe dose, or endpoints#
  #comp>1
  #### the next section constructs the covariance matrix###
  if (comp==1){return(NA)}

  s_0=s[1] #information fraction at the end of phase 2 stage
  r=(s[-1]-s[1])
  #information fraction conditonal on phase 2 sthage
  #alpha=c(0.001,0.001,0.0248)
  ########## non-normal distribution- phase 2 stage matrix###
  cov_mtx_0_nonn=diag(comp)
  cov_mtx_0_nonn=s_0*cov_mtx_0_nonn

  ##############non-normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_nonn=diag(comp-1)
  cov_mtx_0c_nonn=s_0*cov_mtx_01_nonn

  ########## normal distribution- phase 2 stage matrix###
  cov_mtx_0_n=diag(comp)
  cov_mtx_0_n[lower.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n[upper.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n=s_0*cov_mtx_0_n

  ##############normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_n=diag(comp-1)
  cov_mtx_01_n[lower.tri(cov_mtx_01_n)]=0.5
  cov_mtx_01_n[upper.tri(cov_mtx_01_n)]=0.5

  a_mtx_01=(1/4)*matrix(rep(1,(comp-1)*(comp-1)),(comp-1),(comp-1))

  con_mtx_01_n=cov_mtx_01_n-a_mtx_01
  cov_mtx_0c_n=s_0*con_mtx_01_n

  ###### Matirx for phase 3 stage #######

  cov_mtx_1= matrix(rep(r,length(r)),nrow=length(r), ncol=length(r))
  cov_mtx_1[lower.tri(cov_mtx_1)]=t(cov_mtx_1)[lower.tri(cov_mtx_1)]

  #######################Civariance matrix for phase 3 ##########
  s_p3=(s[-1]-s[1])/(last(s)-s[1])
  cov_mtx_p3= matrix(rep(s_p3,length(s_p3)),nrow=length(s_p3), ncol=length(s_p3))
  cov_mtx_p3[lower.tri(cov_mtx_p3)]=t(cov_mtx_p3)[lower.tri(cov_mtx_p3)]
  ##############################
  ## Phase 2 exit probability  ##
  ###########################################################################
  ### calculating exit boundary at first interim (end of phase 2) analysis##
  ###########################################################################

  tp01 <-function(e_0,sigma0){
    rs= 1-pmvnorm(lower=-Inf, upper=rep(e_0,comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]
    return(rs)
  }

  ##################################################################
  ## Calculating exit boundary at inter_vt-th interim  analysis#####
  ##################################################################
  tp_norm_2 <-function(e_1,e_bry,inter_vt,sigma0,sigma1,sigma2,comp){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2
    rank=inter_vt
    u_bry=c(e_bry,e_1)

    intg <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s[1]))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(1/2*x,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
         (1-pmvnorm(upper=u_bry[2:i]-x,lower=-Inf,mean=rep(0,(rank-1)),sigma=sigma2,algorithm=Miwa())[1]  )
      return(rs)
    }

    rs= 1-pmvnorm(lower=-Inf, upper=rep(e_bry[1],comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]+
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper =e_bry[1],subdivisions = 1000L,abs.tol = integTol)$value
    # rs here is the exit probability (P(E_m) in the manuscript)
    return(rs)
  }
  ##################################################################
  ### Non-Normal distributions --phase 3 stage exit probability ####
  ##################################################################

  tp_nonnorm_2 <-function(e_1,e_bry,inter_vt,sigma0,sigma1,sigma2,comp){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 0
    rank=inter_vt
    u_bry=c(e_bry,e_1)

    intg <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s[1]))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(0,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
        (1-pmvnorm(upper=u_bry[2:i]-x,lower=-Inf,mean=rep(0,(rank-1)),sigma=sigma2,algorithm=Miwa())[1]  )
      return(rs)
    }

    rs= 1-pmvnorm(lower=-Inf, upper=rep(e_bry[1],comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]+
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper =e_bry[1],subdivisions = 1000L,abs.tol = integTol)$value
    # rs here is the exit probability (P(E_m) in the manuscript)
    return(rs)
  }
  ####################################################
  ### functions for alpha-spending calculation #######
  ####################################################

  #############################################################
  tp_norm_2_cumu <-function(e_bry,inter_vt,sigma0,sigma1,sigma2,comp){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2
    rank=inter_vt
    u_bry=e_bry

    intg <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s[1]))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(1/2*x,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
         (1-pmvnorm(upper=u_bry[2:i]-x,lower=-Inf,mean=rep(0,(rank-1)),sigma=sigma2,algorithm=Miwa())[1])
      return(rs)
    }

    rs= 1-pmvnorm(lower=-Inf, upper=rep(e_bry[1],comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]+
      	comp* stats::integrate( Vectorize(intg),lower= -Inf,upper = e_bry[1],subdivisions = 1000L,abs.tol = integTol)$value
    # rs here is the exit probability (P(E_m) in the manuscript)
    return(rs)
  }
  #############################################################
  #############################################################
  ### Non-Normal distributions --phase 3 stage exit probability ####
  #############################################################

  tp_nonnorm_2_cumu <-function(e_bry,inter_vt,sigma0,sigma1,sigma2,comp){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 0

    rank=inter_vt
    u_bry=e_bry

    intg <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s[1]))*pmvnorm(lower=-Inf,upper=rep(x,comp-1),mean=rep(0,comp-1),sigma=sigma1,algorithm=Miwa())[1]*
        (1-pmvnorm(upper=u_bry[2:i]-x,lower=-Inf,mean=rep(0,(rank-1)),sigma=sigma2,algorithm=Miwa())[1]  )

      return(rs)
    }

    rs= 1-pmvnorm(lower=-Inf, upper=rep(e_bry[1],comp),mean=rep(0,comp),sigma=sigma0,algorithm=Miwa())[1]+
        comp* stats::integrate( Vectorize(intg),lower= -Inf,upper = e_bry[1],subdivisions = 1000L,abs.tol = integTol)$value
    # rs here is the exit probability (P(E_m) in the manuscript)
    return(rs)
  }

  ##############################
  ####  exit boundary  #########
  ##############################
  if (dist==1){
    if (length(s)==1){
      e_bry=as.numeric(stats::uniroot(function(x) tp01(x,cov_mtx_0_n)-alpha[1],c(-2^8,2^8),extendInt="yes",tol = 1e-4)[1])
    } else {
      #Normal distribution, off-diagnal elements are known to be 1/2
      e_bry=as.numeric(stats::uniroot(function(x) tp01(x,cov_mtx_0_n)-alpha[1],c(-2^8,2^8),extendInt="yes",tol = unirTol)[1])
      for (i in (2:length(s))) {
        e_bry_i=as.numeric(stats::uniroot(function(x) tp_norm_2(x,e_bry,i,cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_1[1:(i-1),1:(i-1)],comp)-sum(alpha[1:i]),c(-2^8,2^8),extendInt="yes",tol = unirTol)[1])
        e_bry=c(e_bry,e_bry_i)
      }
    }
  }else{
    if (length(s)==1){
      e_bry=as.numeric(stats::uniroot(function(x) tp01(x,cov_mtx_0_nonn)-alpha[1],c(-2^8,2^8),extendInt="yes",tol = unirTol)[1])
    } else {
      #Normal distribution, off-diagnal elements are known to be 0
      e_bry=as.numeric(stats::uniroot(function(x) tp01(x,cov_mtx_0_nonn)-alpha[1],c(-2^8,2^8),extendInt="yes",tol = unirTol)[1])
      for (i in (2:length(s))) {
        e_bry_i=as.numeric(stats::uniroot(function(x)  tp_nonnorm_2(x,e_bry,i,cov_mtx_0_nonn, cov_mtx_0c_nonn,cov_mtx_1[1:(i-1),1:(i-1)],comp)-sum(alpha[1:i]),c(-2^8,2^8),extendInt="yes",tol = unirTol)[1])
        e_bry=c(e_bry,e_bry_i)
      }
    }
  }
  ###########################################
  ####  alpha-spending  calcualtion #########
  ###########################################

  cumu_alpha=c(rep(0,length(s)))

  if (dist==1){
    if (length(s)==1){
      cumu_alpha=tp01(e_bry[1],cov_mtx_0_n)
    } else {
      cumu_alpha[1]=tp01(e_bry[1],cov_mtx_0_n)
      for (i in (2:length(s))){
        cumu_alpha[i]=tp_norm_2_cumu(e_bry[1:i],i,cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_1[1:(i-1),1:(i-1)],comp)
      }
    }
  }else {
    if (length(s)==1){
      cumu_alpha=tp01(e_bry[1],cov_mtx_0_nonn)
    } else {
      cumu_alpha[1]=tp01(e_bry[1],cov_mtx_0_nonn)
      for (i in (2:length(s))){
        cumu_alpha[i]= tp_nonnorm_2_cumu(e_bry[1:i],i,cov_mtx_0_nonn,cov_mtx_0c_nonn,cov_mtx_1[1:(i-1),1:(i-1)],comp)
      }
    }
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

  c_bry_two_stg=e_bry/sqrt(s)

  ###################################################################
  #### Getting the critical bundaries for the phase 3 stage only ####
  ###################################################################

  #### calculate the alpha spending for ppahse 3 only#########
  if (length(s)>2) {
    alpha_p3=c(alpha[2:(length(s)-1)],alpha[length(s)]+alpha[1])
  } else if (length(s)==2){
    alpha_p3=sum(alpha)
  }
  ############################################################################

  c_p3_cal1<-function(e0_1){
    #e0_1 is the  exit boundary for the first interim analysis#
    up_bry=e0_1
    tmp = 1-stats::pnorm(up_bry)
    return(tmp)
  }
  ##########################################
  c_p3_cal2<-function(bry,e0_2,sigma1){
    #bry1 is the  exit boundary for the first interim analysis#
    #e0_2 is the  exit boundary for the second interim analysis#
    #sigma2 is the covariance matrix at the second interim analysis#
    up_bry_0=c(bry,e0_2)
    rank=length(up_bry_0)
    up_bry=up_bry_0*sqrt(s_p3[1:rank])
    sigma2=sigma1[1:rank,1:rank]
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rank),sigma=sigma2,algorithm=Miwa())[1]
    return(tmp)
  }
  ###########################################
  c_p3_cal3<-function(up_bry_0,rank,sigma){
    #e0 is the  exit boundary #
    #sigma is the covariance matrix, it has the rank of =rep_time x reptime#
    sigma1=sigma[1:rank,1:rank]
    up_bry=up_bry_0*sqrt(s_p3[1:rank])
    tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,rank),sigma=sigma1,algorithm=Miwa())[1]
    return(tmp)
  }

  if(length(s_p3)==1){
    bry_p3=as.numeric(stats::uniroot(function(x) c_p3_cal1(x)-alpha_p3,c(-2^8,2^8),extendInt="yes",tol = unirTol)[1])
    #bry1,bry2, bry3 are exit boundaries. c_bry are the ctritical bundaries#
    c_bry_p3=bry_p3
    cum_alpha_0_p3=1-stats::pnorm(bry_p3)
  }

  ###################################################################
  if(length(s_p3)>1){
    bry_p3=as.numeric(stats::uniroot(function(x) c_p3_cal1(x)-alpha_p3[1],c(-2^8,2^8),extendInt="yes",tol = unirTol)[1])

    for (i in (2:length(s_p3))){
      bry_p3_i=as.numeric(stats::uniroot(function(x) c_p3_cal2(bry_p3,x,cov_mtx_p3)-sum(alpha_p3[1:i]),c(-2^8,2^8),extendInt="yes",tol = unirTol)[1])
      bry_p3=c(bry_p3,bry_p3_i)
    }
  }
  c_bry_p3=bry_p3
  ## the next loop canculate the exit bpundaries bry1, bry2,bry3 and critical boundaries c_bry, when there are a to of three interim analyses##
  #####################################

  cumu_alpha_p3=c(rep(0,length(s_p3)))
  for (i in 1:length(s_p3)){
    cumu_alpha_p3[i]=c_p3_cal3(c_bry_p3[1:i],i,cov_mtx_p3)
  }
  alpha_s_p3=c(rep(0,length(s_p3)))
  if(length(s_p3)==1){
    alpha_s_p3=cumu_alpha_p3
  } else {
    alpha_s_p3[1]=cumu_alpha_p3[1]
    for (i in 2:length(s_p3)){
      alpha_s_p3[i]=cumu_alpha_p3[i]-cumu_alpha_p3[i-1]
    }
  }

  nom_alpha=1-stats::pnorm(c_bry_two_stg)
  nom_alpha_p3=1-stats::pnorm(c_bry_p3)
  #    sum_nom_alpha_p3=sum(nom_alpha_p3)

  ph3=list(s_p3,c_bry_p3,nom_alpha_p3)
  names(ph3)=c("phase 3 information fraction", "pahse 3 critical boundary","nominal p-value for phase 3")

  out=list(c_bry_two_stg,alpha_s,nom_alpha,cumu_alpha,ph3)
  names(out)=c("Critical boundary","1-sided alpha spending at each interim analysis","1-sided nominal alpha  at each interim analysis","1-sided accumulative alpha","phase 3 ")
  return(out)
}
################################

#two_stage_alsp_bundry(s,alpha,comp,dist)
#two_stage_alsp_bundry(c(0.5,1),c(0.0001,0.0249),2,1)
#two_stage_alsp_bundry(c(0.5,1),c(0.0001,0.0249),2,2)
#two_stage_alsp_bundry(c(1/3,1),c(0.001,0.0249),3,1)
#two_stage_alsp_bundry(c(1/3,1),c(0.001,0.0249),3,2)
#two_stage_alsp_bundry(c(1/3,2/3,1),c(0.0001,0.0001,0.0248),2,2)
#two_stage_alsp_bundry(c(1/3,2/3,1),c(0.00005,0.00005,0.0249),3,1)
#two_stage_alsp_bundry(c(1/3,2/3,1),c(0.0001,0.0001,0.0248),3,2)
#two_stage_alsp_bundry(c(1/4,1/3,2/3,1),c(0.00005,0.00005,0.00005,0.02485),3,1)
#two_stage_alsp_bundry(c(1/4,1/3,2/3,1),c(0.00005,0.00005,0.00005,0.02485),3,2)

###################################################
####                 Part II                  #####
###  Information time and power  calculation  #####
###################################################

###################################################
# This function calculates the required information time for the desired power

two_stage_info_time<-function(mu,s,c_bry,beta,alpha,dist){
  #This function calculates the infromation required for 1-beta power
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #c_bry is the critial booundary
  #1-beta is the desired power
  #mu is the vector of effct sizes
  #alpha  is the target one-sided type I error#
  #mu=c(0.2,0.3,0.4)
  #s=c(1/3,2/3,1)
  # c_bry=c(4.115415, 2.910038, 2.376036)
  #beta=0.1
  # alpha=0.025
  #dist=1

  ########################################
  #### the next section constructs the covariance matrix###
  s_0=s[1] #information fraction at the end of phase 2 stage
  r=(s[-1]-s[1])
  s_1=(s[-1]-s[1])/(last(s)-s[1])

  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  #information fraction conditonal on phase 2 sthage
  #alpha=c(0.001,0.001,0.0248)
  ########## non-normal distribution- phase 2 stage matrix###
  cov_mtx_0_nonn=diag(comp)
  cov_mtx_0_nonn=s_0*cov_mtx_0_nonn

  ##############non-normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_nonn=diag(comp-1)
  cov_mtx_0c_nonn=s_0*cov_mtx_01_nonn

  ########## normal distribution- phase 2 stage matrix###
  cov_mtx_0_n=diag(comp)
  cov_mtx_0_n[lower.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n[upper.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n=s_0*cov_mtx_0_n

  ##############normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_n=diag(comp-1)
  cov_mtx_01_n[lower.tri(cov_mtx_01_n)]=0.5
  cov_mtx_01_n[upper.tri(cov_mtx_01_n)]=0.5

  a_mtx_01=(1/4)*matrix(rep(1,(comp-1)*(comp-1)),(comp-1),(comp-1))

  con_mtx_01_n=cov_mtx_01_n-a_mtx_01
  cov_mtx_0c_n=s_0*con_mtx_01_n

  ###### Matirx for phase 3 stage #######
  cov_mtx_1= matrix(nrow=length(r), ncol=length(r))

  for (i in 1:length(r)){
    cov_mtx_1[i,1:length(r)]=r[i]
  }

  cov_mtx_1[lower.tri(cov_mtx_1)]=t(cov_mtx_1)[lower.tri(cov_mtx_1)]

  #####################################################

  ### the function calculates the  probability of a dose being at the end of pahse 2 ###
  tp_norm_0 <-function(low_bry,up_bry,sigma1){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg_norm_0  <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=up_bry+x,mean=rep(1/2*x,comp-1),sigma=sigma1,abseps=1e-6,algorithm=Miwa())[1]
      return(rs)
    }
    ## rs1 is th eprobability being selected as the largest###
    ### rs2 is the probabiity of being selected and exit the upper boundary##
    rs1= stats::integrate( Vectorize(intg_norm_0),lower= -Inf,upper =Inf,subdivisions = 1000L,abs.tol = integTol)$value
    rs2= stats::integrate( Vectorize(intg_norm_0),lower= low_bry,upper =Inf,subdivisions = 1000L,abs.tol = integTol)$value
    rs3=rs2/rs1
    return(c(rs1,rs2,rs3))
  }

  tp_norm <-function(up_bry_1,up_bry_2,up_bry_3,sigma1,sigma2){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg_norm  <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0)) * pmvnorm(lower=-Inf,upper=up_bry_2+x,mean=rep(1/2*x,comp-1),sigma=sigma1,abseps=1e-6,algorithm=Miwa())[1] *
         (1-pmvnorm(upper=up_bry_3-x,lower=-Inf,mean=rep(0,length(r)),sigma=sigma2,abseps=1e-6,algorithm=Miwa())[1])
      return(rs)
    }

    rs= stats::integrate( Vectorize(intg_norm),lower= -Inf,upper =up_bry_1,subdivisions = 1000L,abs.tol = integTol)$value
    return(rs)
  }

  tp_norm_all<-function(time,mu,sigma0,sigma1,sigma2){

    #########exiting probability at the end of phase 2 stage#####
    up_bry_0=c_bry[1]*sqrt(s[1])-mu*sqrt(time)*s[1]
    exit_prob_0= 1-pmvnorm(lower=-Inf, upper=up_bry_0,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]

    ######### Exiting probability in the phase 3 stage
    up_bry_1=NULL
    up_bry_2=NULL
    up_bry_3=NULL
    exit_prob=NULL
    for (i in 1:length(mu)) {
      up_bry_1[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(time)*s[1]
      up_bry_2=rbind(up_bry_2,c((mu[i]-mu[-i])*sqrt(time)*s[1]))
      up_bry_3=rbind(up_bry_3,c(c_bry[-1]*sqrt(s[-1])-mu[i]*sqrt(time)*s[-1]))
      exit_prob[i]=tp_norm(up_bry_1[i],up_bry_2[i,],up_bry_3[i,],sigma1,sigma2)
    }
    temp= exit_prob_0+sum(exit_prob)
  }

  tp_norm_all_vector<-function(time,mu,sigma0,sigma1,sigma2){

    #########exiting probability at the end of phase 2 stage#####
    up_bry_0=c_bry[1]*sqrt(s[1])-mu*sqrt(time)*s[1]
    exit_prob_0= 1-pmvnorm(lower=-Inf, upper=up_bry_0,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]

    ######### Exiting probability in the phase 3 stage
    up_bry_1=NULL
    up_bry_2=NULL
    up_bry_3=NULL
    exit_prob_vector=NULL
    for (i in 1:length(mu)) {
      up_bry_1[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(time)*s[1]
      up_bry_2=rbind(up_bry_2,c((mu[i]-mu[-i])*sqrt(time)*s[1]))
      up_bry_3=rbind(up_bry_3,c(c_bry[-1]*sqrt(s[-1])-mu[i]*sqrt(time)*s[-1]))
      exit_prob_vector[i]=tp_norm(up_bry_1[i],up_bry_2[i,],up_bry_3[i,],sigma1,sigma2)
    }
     temp=exit_prob_vector
  }

  tp_non_norm <-function(up_bry_1,up_bry_2,up_bry_3,sigma1,sigma2){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg_non_norm  <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=up_bry_2+x,mean=rep(0,comp-1),sigma=sigma1,abseps=1e-6,algorithm=Miwa())[1]*
        (1-pmvnorm(upper=up_bry_3-x,lower=-Inf,mean=rep(0,length(r)),sigma=sigma2,abseps=1e-6,algorithm=Miwa())[1]  )
      return(rs)
    }

    rs= stats::integrate( Vectorize(intg_non_norm),lower= -Inf,upper =up_bry_1,subdivisions = 1000L,abs.tol = integTol)$value
    return(rs)
  }

  tp_non_norm_0 <-function(low_bry,up_bry,sigma1){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg_non_norm_0  <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=up_bry+x,mean=rep(0,comp-1),sigma=sigma1,abseps=1e-6,algorithm=Miwa())[1]
      return(rs)
    }
    ## rs1 is th eprobability being selected as the largest###
    ### rs2 is the probabiity of being selected and exit the upper boundary##
    rs1= stats::integrate( Vectorize(intg_non_norm_0),lower= -Inf,upper =Inf,subdivisions = 1000L,abs.tol = integTol)$value
    rs2= stats::integrate( Vectorize(intg_non_norm_0),lower= low_bry,upper =Inf,subdivisions = 1000L,abs.tol = integTol)$value
    rs3=rs2/rs1
    return(c(rs1,rs2,rs3))
  }

  tp_non_norm_all<-function(time,mu,sigma0,sigma1,sigma2){

    #########exiting probability at the end of phase 2 stage#####
    up_bry_0=c_bry[1]*sqrt(s[1])-mu*sqrt(time)*s[1]
    exit_prob_0= 1-pmvnorm(lower=-Inf, upper=up_bry_0,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]

    ######### Exiting probability in the phase 3 stage
    up_bry_1=NULL
    up_bry_2=NULL
    up_bry_3=NULL
    exit_prob=NULL
    for (i in 1:length(mu)) {
      up_bry_1[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(time)*s[1]
      up_bry_2=rbind(up_bry_2,c((mu[i]-mu[-i])*sqrt(time)*s[1]))
      up_bry_3=rbind(up_bry_3,c(c_bry[-1]*sqrt(s[-1])-mu[i]*sqrt(time)*s[-1]))
      exit_prob[i]=tp_non_norm(up_bry_1[i],up_bry_2[i,],up_bry_3[i,],sigma1,sigma2)
    }
    temp= exit_prob_0+sum(exit_prob)
  }

  tp_non_norm_all_vector<-function(time,mu,sigma0,sigma1,sigma2){

    #########exiting probability at the end of phase 2 stage#####
    up_bry_0=c_bry[1]*sqrt(s[1])-mu*sqrt(time)*s[1]
    exit_prob_0= 1-pmvnorm(lower=-Inf, upper=up_bry_0,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]

    ######### Exiting probability in the phase 3 stage
    up_bry_1=NULL
    up_bry_2=NULL
    up_bry_3=NULL
    exit_prob_vector=NULL
    for (i in 1:length(mu)) {
      up_bry_1[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(time)*s[1]
      up_bry_2=rbind(up_bry_2,c((mu[i]-mu[-i])*sqrt(time)*s[1]))
      up_bry_3=rbind(up_bry_3,c(c_bry[-1]*sqrt(s[-1])-mu[i]*sqrt(time)*s[-1]))
      exit_prob_vector[i]=tp_non_norm(up_bry_1[i],up_bry_2[i,],up_bry_3[i,],sigma1,sigma2)
    }
    temp=exit_prob_vector
  }

  if (dist==1) {
    total_time=as.numeric(stats::uniroot(function(x) tp_norm_all(x,mu, cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_1)-1+beta,c(0,2^32),extendInt="yes",tol = unirTol)[1])
    exit_prob_vector=tp_norm_all_vector(total_time,mu, cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_1)
    #   return(total_time)
  }else{
    total_time=as.numeric(stats::uniroot(function(x) tp_non_norm_all(x,mu, cov_mtx_0_nonn,cov_mtx_0c_nonn,cov_mtx_1)-1+beta,c(0,2^32),extendInt="yes",tol = unirTol)[1])
    exit_prob_vector=tp_non_norm_all_vector(total_time,mu, cov_mtx_0_nonn,cov_mtx_0c_nonn,cov_mtx_1)
    #  return(total_time)
  }

  if (dist==1) {
    low_bry=NULL
    up_bry=NULL
    up_bry_0=NULL
    exit_prob_0=NULL
    for (i in 1:length(mu)) {
      low_bry[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(total_time)*s[1]
      up_bry_0=rbind(up_bry_0,c(mu[i]-mu[-i]))
      up_bry=rbind(up_bry,c((mu[i]-mu[-i])*sqrt(total_time)*s[1]))
      exit_prob_0=rbind(exit_prob_0,tp_norm_0(low_bry[i],up_bry[i,],cov_mtx_0c_n))
      #tp_norm_0 <-function(low_bry,up_bry,sigma1)
    }
    exit_prob_0_tot=sum( exit_prob_0[,1])
  } else {
    low_bry=NULL
    up_bry=NULL
    up_bry_0=NULL
    exit_prob_0=NULL
    for (i in 1:length(mu)) {
      low_bry[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(total_time)*s[1]
      up_bry_0=rbind(up_bry_0,c(mu[i]-mu[-i]))
      up_bry=rbind(up_bry,c((mu[i]-mu[-i])*sqrt(total_time)*s[1]))
      exit_prob_0=rbind(exit_prob_0,tp_non_norm_0(low_bry[i],up_bry[i,],cov_mtx_0c_nonn))
    }
    exit_prob_0_tot=sum( exit_prob_0[,1])
  }

  #################################################################
  two_asd_info_time_s<-function(mu,s,c_bry,beta,alpha){
    #This function calculates the infromation required for 1-beta power
    #s is the vector of information fractions. The length of s can be either 2 or 3#
    #c_bry is the critial booundary
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
      tmp = 1-pmvnorm(lower=-Inf,upper=up_bry,mean=rep(0,length(s)),sigma=sigma1,abseps=1e-8,algorithm=Miwa())[1]
      return(tmp)
    }

    if (length(s)==1){
      total_time=as.numeric(stats::uniroot(function(x) cal1(x)-1+beta,c(0,2^32),extendInt="yes",tol = unirTol)[1])
    }else {
      total_time=as.numeric(stats::uniroot(function(x) cal2(x,cov_mtx)-1+beta,c(0,2^32),extendInt="yes",tol = unirTol)[1])
    }
    return(total_time)
  }

  info_time_ind=NULL
  for (i in 1:length(mu)) {
    if (mu[i]<1e-4) {
      info_time_ind[i]=NA
    }else{
      info_time_ind[i]=two_asd_info_time_s(mu[i],s,c_bry,beta,alpha)
    }
  }
  #################################################################
  exit_prob_vector_tot=exit_prob_vector+exit_prob_0[,2]
  cond_exit_prob_vector=exit_prob_vector_tot/exit_prob_0[,1]
  ##################################################################
  out=list(total_time,info_time_ind,exit_prob_0[,1], exit_prob_vector_tot,sum(exit_prob_0[,1]),sum(exit_prob_vector_tot),cond_exit_prob_vector)
  names(out)=c("total information time","Individual information time","selection probabiltiy for each mu","exit probabiltiy for each mu","total prob-phase 2","tot exit prob","conditional exit prob if mu selected")
  return(out)
}
########################
#two_stage_info_time(mu,s,c_bry,beta,alpha,dist)
#two_stage_info_time(c(-0.2,0,0.01),c(1/3,2/3,1),c(3.938306,2.784803,2.273782),0.1,0.025,1)
#two_stage_info_time(c(0.2,0.3,0.4),c(1/3,2/3,1),c(3.938306,2.784803,2.273782),0.1,0.025,1)
#two_stage_info_time(c(0.0001,0.0001,0.0001),c(1/3,2/3,1),c(3.938306,2.784803,2.273782),0.1,0.025,1)
#two_stage_info_time(c(1e-5,1e-5,1e-5),c(1/3,2/3,1),c(3.938306,2.784803,2.273782),0.1,0.025,1)
#two_stage_info_time(c(0.00001,0.00001,0.00001),c(1/3,2/3,1),c(3.938306,2.784803,2.273782),0.1,0.025,1)
#two_stage_info_time(c(0.3,0.3,0.3),c(1/3,2/3,1),c(3.938306,2.784803,2.273782),0.1,0.025,1)
#two_stage_info_time(c(0.2,0.3,0.4),c(1/3,2/3,1),c( 4.057782, 2.869285,2.342762),0.1,0.025,2)
#two_stage_info_time(c(0.3,0.3,0.3),c(1/3,2/3,1),c( 4.057782, 2.869285,2.342762),0.1,0.025,2)
#two_stage_info_time(c(0.2,0.3,0.3),c(1/3,2/3,1),c( 4.057782, 2.869285,2.342762),0.1,0.025,2)
#two_stage_info_time(c(0,0.3,0.6),c(1/3,2/3,1),c( 4.057782, 2.869285,2.342762),0.1,0.025,2)

two_stage_time_power<-function(mu,s,c_bry,T,alpha,dist){
  #This function calculates the infromation required for 1-beta power
  #s is the vector of information fractions. The length of s can be either 2 or 3#
  #c_bry is the critial boundary
  #1-beta is the desired power
  #mu is the vector of effct sizes
  #alpha  is the target one-sided type I error#
  #mu=c(0.2,0.3,0.4)
  #s=c(1/3,2/3,1)
  # c_bry=c(4.115415, 2.910038, 2.376036)
  #beta=0.1
  # alpha=0.025
  #dist=1

  ########################################
  #### the next section constructs the covariance matrix###

  s_0=s[1] #information fraction at the end of phase 2 stage

  r=(s[-1]-s[1])
  s_1=(s[-1]-s[1])/(last(s)-s[1])
  comp=length(mu)
  #information fraction conditonal on phase 2 sthage
  #alpha=c(0.001,0.001,0.0248)
  ########## non-normal distribution- phase 2 stage matrix###
  cov_mtx_0_nonn=diag(comp)
  cov_mtx_0_nonn=s_0*cov_mtx_0_nonn

  ##############non-normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_nonn=diag(comp-1)
  cov_mtx_0c_nonn=s_0*cov_mtx_01_nonn

  ########## normal distribution- phase 2 stage matrix###
  cov_mtx_0_n=diag(comp)
  cov_mtx_0_n[lower.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n[upper.tri(cov_mtx_0_n)]=0.5
  cov_mtx_0_n=s_0*cov_mtx_0_n

  ##############normal distribution- phase 2 stage, conditional distribution matrix###

  cov_mtx_01_n=diag(comp-1)
  cov_mtx_01_n[lower.tri(cov_mtx_01_n)]=0.5
  cov_mtx_01_n[upper.tri(cov_mtx_01_n)]=0.5

  a_mtx_01=(1/4)*matrix(rep(1,(comp-1)*(comp-1)),(comp-1),(comp-1))

  con_mtx_01_n=cov_mtx_01_n-a_mtx_01
  cov_mtx_0c_n=s_0*con_mtx_01_n

  ###### Matirx for phase 3 stage #######
  cov_mtx_1= matrix(nrow=length(r), ncol=length(r))
  for (i in 1:length(r)){
    cov_mtx_1[i,1:length(r)]=r[i]
  }
  cov_mtx_1[lower.tri(cov_mtx_1)]=t(cov_mtx_1)[lower.tri(cov_mtx_1)]

  #####################################################

  ### the function calculates the  probability of a dose being at the end of pahse 2 ###
  tp_norm_0 <-function(low_bry,up_bry,sigma1){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg_norm_0  <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)
      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=up_bry+x,mean=rep(1/2*x,comp-1),sigma=sigma1,abseps=1e-6,algorithm=Miwa())[1]
      return(rs)
    }
    ## rs1 is th eprobability being selected as the largest###
    ### rs2 is the probabiity of being selected and exit the upper boundary##
    rs1= stats::integrate( Vectorize(intg_norm_0),lower= -Inf,upper =Inf,subdivisions = 1000L,abs.tol = integTol)$value
    rs2= stats::integrate( Vectorize(intg_norm_0),lower= low_bry,upper =Inf,subdivisions = 1000L,abs.tol = integTol)$value
    rs3=rs2/rs1
    return(c(rs1,rs2,rs3))
  }

  tp_norm <-function(up_bry_1,up_bry_2,up_bry_3,sigma1,sigma2){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg_norm  <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=up_bry_2+x,mean=rep(1/2*x,comp-1),sigma=sigma1,abseps=1e-6,algorithm=Miwa())[1]*
         (1-pmvnorm(upper=up_bry_3-x,lower=-Inf,mean=rep(0,length(r)),sigma=sigma2,abseps=1e-6,algorithm=Miwa())[1])
      return(rs)
    }

    rs= stats::integrate( Vectorize(intg_norm),lower= -Inf,upper =up_bry_1,subdivisions = 1000L,abs.tol = integTol)$value
    return(rs)
  }

  tp_norm_all<-function(time,mu,sigma0,sigma1,sigma2){

    #########exiting probability at the end of phase 2 stage#####
    up_bry_0=c_bry[1]*sqrt(s[1])-mu*sqrt(time)*s[1]
    exit_prob_0= 1-pmvnorm(lower=-Inf, upper=up_bry_0,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]

    ######### Exiting probability in the phase 3 stage
    up_bry_1=NULL
    up_bry_2=NULL
    up_bry_3=NULL
    exit_prob=NULL
    for (i in 1:length(mu)) {
      up_bry_1[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(time)*s[1]
      up_bry_2=rbind(up_bry_2,c((mu[i]-mu[-i])*sqrt(time)*s[1]))
      up_bry_3=rbind(up_bry_3,c(c_bry[-1]*sqrt(s[-1])-mu[i]*sqrt(time)*s[-1]))
      exit_prob[i]=tp_norm(up_bry_1[i],up_bry_2[i,],up_bry_3[i,],sigma1,sigma2)
    }
    temp= exit_prob_0+sum(exit_prob)
  }

  tp_norm_all_vector<-function(time,mu,sigma0,sigma1,sigma2){

    #########exiting probability at the end of phase 2 stage#####
    up_bry_0=c_bry[1]*sqrt(s[1])-mu*sqrt(time)*s[1]
    exit_prob_0= 1-pmvnorm(lower=-Inf, upper=up_bry_0,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]

    ######### Exiting probability in the phase 3 stage
    up_bry_1=NULL
    up_bry_2=NULL
    up_bry_3=NULL
    exit_prob_vector=NULL
    for (i in 1:length(mu)) {
      up_bry_1[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(time)*s[1]
      up_bry_2=rbind(up_bry_2,c((mu[i]-mu[-i])*sqrt(time)*s[1]))
      up_bry_3=rbind(up_bry_3,c(c_bry[-1]*sqrt(s[-1])-mu[i]*sqrt(time)*s[-1]))
      exit_prob_vector[i]=tp_norm(up_bry_1[i],up_bry_2[i,],up_bry_3[i,],sigma1,sigma2)
    }
    temp=exit_prob_vector
  }

  tp_non_norm <-function(up_bry_1,up_bry_2,up_bry_3,sigma1,sigma2){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg_non_norm  <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)

      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=up_bry_2+x,mean=rep(0,comp-1),sigma=sigma1,abseps=1e-6,algorithm=Miwa())[1]*
         (1-pmvnorm(upper=up_bry_3-x,lower=-Inf,mean=rep(0,length(r)),sigma=sigma2,abseps=1e-6,algorithm=Miwa())[1]  )
      return(rs)
    }

    rs= stats::integrate( Vectorize(intg_non_norm),lower= -Inf,upper =up_bry_1,subdivisions = 1000L,abs.tol = integTol)$value
    return(rs)
  }

  tp_non_norm_0 <-function(low_bry,up_bry,sigma1){
    #This function is used when the off-diagnal elements
    #of the covariance matrix is  known to be 1/2

    intg_non_norm_0  <- function(x){
      #this function is the product of three parts:
      #stats::dnorm.(..) is the  density of the selected dose;
      # pmvnorm(..) is the conditional probability of the dtat in the phase 2 stage
      #(1-pmvnorm(...) is the conditional probabiltity of events happend in the phase 3 stage)
      rs=stats::dnorm(x,mean=0,sd=sqrt(s_0))*pmvnorm(lower=-Inf,upper=up_bry+x,mean=rep(0,comp-1),sigma=sigma1,abseps=1e-6,algorithm=Miwa())[1]
      return(rs)
    }
    ## rs1 is th eprobability being selected as the largest###
    ### rs2 is the probabiity of being selected and exit the upper boundary##
    rs1= stats::integrate( Vectorize(intg_non_norm_0),lower= -Inf,upper =Inf,subdivisions = 1000L,abs.tol = integTol)$value
    rs2= stats::integrate( Vectorize(intg_non_norm_0),lower= low_bry,upper =Inf,subdivisions = 1000L,abs.tol = integTol)$value
    rs3=rs2/rs1
    return(c(rs1,rs2,rs3))
  }

  tp_non_norm_all<-function(time,mu,sigma0,sigma1,sigma2){
    #########exiting probability at the end of phase 2 stage#####
    up_bry_0=c_bry[1]*sqrt(s[1])-mu*sqrt(time)*s[1]
    exit_prob_0= 1-pmvnorm(lower=-Inf, upper=up_bry_0,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]

    ######### Exiting probability in the phase 3 stage
    up_bry_1=NULL
    up_bry_2=NULL
    up_bry_3=NULL
    exit_prob=NULL
    for (i in 1:length(mu)) {
      up_bry_1[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(time)*s[1]
      up_bry_2=rbind(up_bry_2,c((mu[i]-mu[-i])*sqrt(time)*s[1]))
      up_bry_3=rbind(up_bry_3,c(c_bry[-1]*sqrt(s[-1])-mu[i]*sqrt(time)*s[-1]))
      exit_prob[i]=tp_non_norm(up_bry_1[i],up_bry_2[i,],up_bry_3[i,],sigma1,sigma2)
    }
    temp= exit_prob_0+sum(exit_prob)
  }

  tp_non_norm_all_vector<-function(time,mu,sigma0,sigma1,sigma2){
    #########exiting probability at the end of phase 2 stage#####
    up_bry_0=c_bry[1]*sqrt(s[1])-mu*sqrt(time)*s[1]
    exit_prob_0= 1-pmvnorm(lower=-Inf, upper=up_bry_0,mean=rep(0,comp),sigma=sigma0,abseps=1e-6,algorithm=Miwa())[1]

    ######### Exiting probability in the phase 3 stage
    up_bry_1=NULL
    up_bry_2=NULL
    up_bry_3=NULL
    exit_prob_vector=NULL
    for (i in 1:length(mu)) {
      up_bry_1[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(time)*s[1]
      up_bry_2=rbind(up_bry_2,c((mu[i]-mu[-i])*sqrt(time)*s[1]))
      up_bry_3=rbind(up_bry_3,c(c_bry[-1]*sqrt(s[-1])-mu[i]*sqrt(time)*s[-1]))
      exit_prob_vector[i]=tp_non_norm(up_bry_1[i],up_bry_2[i,],up_bry_3[i,],sigma1,sigma2)
    }
    temp=exit_prob_vector
  }

  if (dist==1) {
    power_tot=tp_norm_all(T,mu, cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_1)
    exit_prob_vector=tp_norm_all_vector(T,mu, cov_mtx_0_n, cov_mtx_0c_n,cov_mtx_1)
    #   return(total_time)
    # temp=c( power_tot,exit_prob_vector)
  }else{
    power_tot=tp_norm_all(T,mu, cov_mtx_0_n, cov_mtx_0c_nonn,cov_mtx_1)
    exit_prob_vector=tp_non_norm_all_vector(T,mu, cov_mtx_0_nonn,cov_mtx_0c_nonn,cov_mtx_1)
    #  return(total_time)
    # temp=c( power_tot,exit_prob_vector)
  }

  if (dist==1) {
    low_bry=NULL
    up_bry=NULL
    up_bry_0=NULL
    exit_prob_0=NULL
    for (i in 1:length(mu)) {
      low_bry[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(T)*s[1]
      up_bry_0=rbind(up_bry_0,c(mu[i]-mu[-i]))
      up_bry=rbind(up_bry,c((mu[i]-mu[-i])*sqrt(T)*s[1]))
      exit_prob_0=rbind(exit_prob_0,tp_norm_0(low_bry[i],up_bry[i,],cov_mtx_0c_n))
      #tp_norm_0 <-function(low_bry,up_bry,sigma1)
    }
    exit_prob_0_tot=sum( exit_prob_0[,1])
  } else {
    low_bry=NULL
    up_bry=NULL
    up_bry_0=NULL
    exit_prob_0=NULL
    for (i in 1:length(mu)) {
      low_bry[i]=c_bry[1]*sqrt(s[1])-mu[i]*sqrt(T)*s[1]
      up_bry_0=rbind(up_bry_0,c(mu[i]-mu[-i]))
      up_bry=rbind(up_bry,c((mu[i]-mu[-i])*sqrt(T)*s[1]))
      exit_prob_0=rbind(exit_prob_0,tp_non_norm_0(low_bry[i],up_bry[i,],cov_mtx_0c_nonn))
    }
    exit_prob_0_tot=sum( exit_prob_0[,1])
  }

  #################################################################

  #################################################################
  exit_prob_vector_tot=exit_prob_vector+exit_prob_0[,2]
  cond_exit_prob_vector=exit_prob_vector_tot/exit_prob_0[,1]
  ##################################################################
  out=list(power_tot,exit_prob_0[,1], exit_prob_vector_tot,sum(exit_prob_0[,1]),sum(exit_prob_vector_tot),cond_exit_prob_vector)
  names(out)=c("total power","selection probabiltiy for each mu","exit probabiltiy for each mu","total prob-phase 2","tot exit prob","conditional exit prob if mu selected")
  return(out)
}
########################

#two_stage_time_power(mu,s,c_bry,T,alpha,dist)
#two_stage_time_power(c(0.2,0.3,0.4),c(1/3,2/3,1),c(3.938306,2.784803,2.273782),89.95164,0.025,1)
#two_stage_time_power(c(0.3,0.3,0.3),c(1/3,2/3,1),c(3.938306,2.784803,2.273782),109.4869,0.025,1)
#two_stage_time_power(c(0.2,0.3,0.4),c(1/3,2/3,1),c( 4.057782, 2.869285,2.342762),91.62886,0.025,2)
#two_stage_time_power(c(0.3,0.3,0.3),c(1/3,2/3,1),c( 4.057782, 2.869285,2.342762),100.578,0.025,2)
#two_stage_time_power(c(0.2,0.3,0.3),c(1/3,2/3,1),c( 4.057782, 2.869285,2.342762),122.0311,0.025,2)
#two_stage_time_power(c(0.1,0.3,0.6),c(1/3,2/3,1),c( 4.057782, 2.869285,2.342762),51.34378,0.025,2)


###################################################
####                 Part III                  #####
###  sample size  and power  calculation       #####
###################################################

################################################################
#### Part III. Two stage Sample size and power clculation   ####
################################################################

#############################################################################
#### Part III.1 power and Sample size clculation: Normal distributions   ####
#############################################################################

###############################################################################
#### Part III.1.1 Sample size and power clculation: Normal distributions   ####
###############################################################################
#######################################################
####    supeiority and non-inferiority trials   #######
#######################################################


################################
### O'Brien_Fleming boundary ###
################################

#' Phase 2/3 seamless combination Group sequential design: Normal distribution, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param mu_t Means of the test arm
#' @param mu_c Means of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'test arm sample size', 'control arm sample size',
#' 'sample size for each dose in a two arm study with same critical boudary',
#' 'sample size for control arm in a two arm study with same critical boudary', 'type I error control', 'exit probabilities'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_OF_boundary_normal(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' margin = 0.1, direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1,
#' dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_OF_boundary_normal<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given mu_t,mu_c,std_t,std_c,s,alpha,beta,rand_ratio,dist1,dist2
  #for superiority and non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #mu_t is the vector of means for the test arms
  #length(mu_t)>=2
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
  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  ##########################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))

  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  ########## two stage sample size###
  n_c=ceiling(total_time*ratio)
  n_t=n_c*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ###################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,info_time_00)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm in a two arm study with same critical boudary","type I error control","exit probabilities")
  return(out)
}

####################################################################
#two_stage_Sample_size_OF_boundary_normal(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
#two_stage_Sample_size_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0,1,c(1/3,2/3,1),0.025,0.1,1,1,1)

#two_stage_Sample_size_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_normal(c(0.2,0.3,0.4),0.2,1,1,0,1,c(1/3,2/3,1),0.025,0.1,1,2,1)
#two_stage_Sample_size_OF_boundary_normal(c(0.2,0.3,0.4),0.4,1,1,0,2,c(1/3,2/3,1),0.025,0.1,1,2,1)
#two_stage_Sample_size_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.15,1,c(1/3,2/3,1),0.025,0.1,1,2,2)
#two_stage_Sample_size_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.15,2,c(1/3,2/3,1),0.025,0.1,1,2,2)
#######################################################################
######################################################################
#' Phase 2/3 seamless combination Group sequential design: Normal distribution, Power calculation, O'Brien-Fleming boundary
#'
#' @param mu_t Means of the test arm
#' @param mu_c Means of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'type I error control'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_OF_boundary_normal(mu_t = c(0.2,0.3,0.4), mu_c = 0, std_t = 1, std_c = 1,
#' margin = 0.01, direction = 1, s = c(1/3,2/3,1), alpha = 0.025, sampsz_control = 180,
#' rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_OF_boundary_normal<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for given mu_t,mu_c,std_t,std_c,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #mu_t is the vector of means for the test arms
  #length(mu_t)>=2
  #mu_c is the mean of the control arm
  #std_t>0 is the common standard deviation of the test arms
  #std_c>0 is the standard deviation of the control arm
  #alpha_2 is the two-sided tye I error
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

  comp=length(mu)
  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################################
  out=list(power_out,bry_out)
  names(out)=c("power","type I error control")
  return(out)
}

####################################################################
#two_stage_power_OF_boundary_normal(mu_t,mu_c,std_t,std_c,margin,direcction,s,alpha,sampsz_control,rand_ratio,dist1,dist2)
#two_stage_power_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.01,1,c(1/3,2/3,1),0.025,180,1,1,1)
#two_stage_power_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.01,2,c(1/3,2/3,1),0.025,180,1,1,1)
#two_stage_power_OF_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),0.025,110,1,1,1)
#two_stage_power_OF_boundary_normal(c(0.2,0.3,0.4),0.2,1,1,0,1,c(1/3,2/3,1),0.025,777,1,2,1)
########################################################################

#################################
### Alpha -spending boundary#####
#################################

#' Phase 2/3 seamless combination Group sequential design: Normal distribution, Sample size calculation, Alpha-spending boundary
#'
#' @param mu_t Means of the test arm
#' @param mu_c Means of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger mu is better, 0 = smaller mu is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'test arm sample size', 'control arm sample size',
#' 'sample size for each dose in a two arm study with same critical boudary',
#' 'sample size for control arm in a two arm study with same critical boudary', 'type I error control', 'exit probabilities'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_alsp_boundary_normal(mu_t = 2.6, mu_c = 0.6, std_t = 9.43, std_c = 9.86,
#' margin = 0.1, direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1,
#' dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_alsp_boundary_normal<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given mu_t,mu_c,std_t,std_c,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #mu_t is the vector of means for the test arms
  #length(mu_t)>=2
  #mu_c is the mean of the control arm
  #std_t>0 is the common standard deviation of the test arms
  #std_c>0 is the standard deviation of the control arm
  #alpha_2 is the two-sided tye I error
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
  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################

  c_bry=as.numeric(unlist(bry_out[1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  ##########################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))

  #### saple size for specified distribution ####
  ratio=std_c^2+std_t^2/rand_ratio
  ########## two stage sample size###
  n_c=ceiling(total_time*ratio)
  n_t=n_c*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ###################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,info_time_00)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm in a two arm study with same critical boudary","type I error control","exit probabilities")
  return(out)
}

####################################################################
#two_stage_Sample_size_alsp_boundary_normal(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){

#two_stage_Sample_size_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#two_stage_Sample_size_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)

#two_stage_Sample_size_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,2,1)
#two_stage_Sample_size_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,2,2)
######################################################################
#' Phase 2/3 seamless combination Group sequential design: Normal distribution, Power calculation, Alpha-spending boundary
#'
#' @param mu_t Mean of the test arm
#' @param mu_c Mean of the control arm
#' @param std_t Standard deviation of the test arm
#' @param std_c Standard deviation of the control arm
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 0 = superiority test, 1 = non-inferiority test
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control 'sampsz_control' as used by this function; see Examples.
#' @param rand_ratio Randomization ratio, sample size of the test arm over the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'type I error control'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_alsp_boundary_normal(mu_t = c(0.2,0.3,0.4), mu_c = 0, std_t = 1, std_c = 1,
#' margin = 0.01, direction = 1, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023),
#' sampsz_control = 180, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_alsp_boundary_normal<-function(mu_t,mu_c,std_t,std_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for given mu_t,mu_c,std_t,std_c,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #rand_ratio is the randomization ratio of test arm:control arm
  #rand_ratio>0
  #alpha is one-sided type I error
  #0<alpha<1
  #typically, alpha=0.025
  #mu_t is the vector of means for the test arms
  #length(mu_t)>=2
  #mu_c is the mean of the control arm
  #std_t>0 is the common standard deviation of the test arms
  #std_c>0 is the standard deviation of the control arm
  #alpha_2 is the two-sided tye I error
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
  if (direction==1){
    mu=mu_t-mu_c+margin
  }else{
    mu=-(mu_t-mu_c-margin)
  }
  comp=length(mu)
  T=sampsz_control/(std_c^2+std_t^2/rand_ratio)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################################
  out=list(power_out,bry_out)
  names(out)=c("power","type I error control")
  return(out)
}

####################################################################
#two_stage_power_alsp_boundary_normal(mu_t,mu_c,std_t,std_c,margin,direcction,s,alpha,sampsz_control,rand_ratio,dist1,dist2)
#two_stage_power_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.01,1,c(1/3,2/3,1),c(0.001,0.001,0.023),180,1,1,1)
#two_stage_power_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.01,2,c(1/3,2/3,1),c(0.001,0.001,0.023),180,1,1,1)
#two_stage_power_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),110,1,1,1)
#two_stage_power_alsp_boundary_normal(c(0.2,0.3,0.4),0,1,1,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),110,1,1,1)

##############################################################


#############################################################################
#### Part III.2 power and Sample size clculation: binary distributions   ####
#############################################################################


#####################################################
####    supeiority trials and non-inferiority #######
#####################################################


################################
### O'Brien_Fleming boundary ###
################################

#############################################################
#' Phase 2/3 seamless combination Group sequential design: Binary distribution, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param p_t Event (response) proportion in the test arm
#' @param p_c Event (response) proportion in the control arm
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger p is better, 0 = smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'test arm sample size', 'control arm sample size',
#' 'sample size for each dose in a two arm study with same critical boudary',
#' 'sample size for control arm in a two arm study with same critical boudary', 'type I error control', 'exit probabilities'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_OF_boundary_binary_diff(p_t = 0.4, p_c = 0.3, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_OF_boundary_binary_diff<-function(p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  #for superiority and non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the vector of rates of the test arm
  #length(p_t)>=2
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
  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  #############################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####

  ########## two stage sample size###
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  n_c=ceiling(total_time*max(ratio))
  n_t=n_c*rand_ratio

  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ###################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,info_time_00)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm in a two arm study with same critical boudary","type I error control","exit probabilities")
  return(out)
}

####################################################################
#two_stage_Sample_size_OF_boundary_binary_diff(p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2)
#two_stage_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.2,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.15,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.1,0.15,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.1,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
######################################################################

#' Phase 2/3 seamless combination Group sequential design: Binary distribution, Power calculation, O'Brien-Fleming boundary
#'
#' @param p_t Event (response) proportion in the test arm
#' @param p_c Event (response) proportion in the control arm
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger p is better, 0 = smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'type I error control'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_OF_boundary_binary_diff(p_t = c(0.2,0.3,0.4), p_c = 0.5, margin = 0.15,
#' direction = 1, s = c(1/3,2/3,1), alpha = 0.025,
#' sampsz_control = 82, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_OF_boundary_binary_diff<-function(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  #for superiority and non-inferiroity trial
  # the critical boundary is O'Brien-Fleming boundary
  #p_t is the vector of rates of the test arm
  #length(p_t)>=2
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
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################
  out=list(power_out,bry_out)
  names(out)=c("power","type I error control")
  return(out)
}

####################################################################
#two_stage_power_OF_boundary_binary_diff(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2)
#two_stage_power_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.15,1,c(1/3,2/3,1),0.025,82,1,1,1)
#two_stage_power_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),0.025,82,1,1,1)
#two_stage_power_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.2,1,c(1/3,2/3,1),0.025,642,1,1,1)
#two_stage_power_OF_boundary_binary_diff(c(0.2,0.3,0.4),0.1,0.15,2,c(1/3,2/3,1),0.025,1654,1,1,1)

######################################################################

#################################
### Alpha -spending boundary#####
#################################

##########################################
#' Phase 2/3 seamless combination Group sequential design: Binary distribution, Sample size calculation, Alpha-spending boundary
#'
#' @param p_t Event (response) proportion in the test arm
#' @param p_c Event (response) proportion in the control arm
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger p is better, 0 = smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'test arm sample size', 'control arm sample size',
#' 'sample size for each dose in a two arm study with same critical boudary',
#' 'sample size for control arm in a two arm study with same critical boudary', 'type I error control', 'exit probabilities'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_alsp_boundary_binary_diff(p_t = 0.4, p_c = 0.3, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025,
#' beta = 0.1, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_alsp_boundary_binary_diff<-function(p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #p_t is the vector of rates of the test arm
  #length(p_t)>=2
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

  if (direction==1){
    mu=p_t-p_c+margin
  }else{
    mu=-(p_t-p_c-margin)
  }
  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  #############################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))
  #### saple size for specified distribution ####

  ########## two stage sample size###
  ratio=(1-p_t)*p_t/rand_ratio+(1-p_c)*p_c
  n_c=ceiling(total_time*max(ratio))
  n_t=n_c*rand_ratio

  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio
  ###################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,info_time_00)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm in a two arm study with same critical boudary","type I error control","exit probabilities")
  return(out)
}

####################################################################
#two_stage_Sample_size_alsp_boundary_binary_diff(p_t,p_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2)
  #two_stage_Sample_size_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#two_stage_Sample_size_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
######################################################################

#' Phase 2/3 seamless combination Group sequential design: Binary distribution, Power calculation, Alpha-spending boundary
#'
#' @param p_t Event (response) proportion in the test arm
#' @param p_c Event (response) proportion in the control arm
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger p is better, 0 = smaller p is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'type I error control'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_alsp_boundary_binary_diff(p_t = c(0.2,0.3,0.4), p_c = 0.5, margin = 0.15,
#' direction = 1, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), sampsz_control = 82,
#' rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_alsp_boundary_binary_diff<-function(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #p_t is the vector of rates of the test arm
  #length(p_t)>=2
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
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################
  out=list(power_out,bry_out)
  names(out)=c("power","type I error control")
  return(out)
}

####################################################################
#two_stage_power_alsp_boundary_binary_diff(p_t,p_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2)
#two_stage_power_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),82,1,1,1)
#two_stage_power_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),82,1,1,1)
#two_stage_power_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.15,1,c(1/3,2/3,1),c(0.001,0.001,0.023),2438,1,1,1)
#two_stage_power_alsp_boundary_binary_diff(c(0.2,0.3,0.4),0.5,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),44,1,1,1)

##########################################


#############################################################################
#### Part III.3 power and Sample size clculation: survival analysis      ####
#############################################################################


#############################################################
####    supeiority trials and non-inferiority trials  #######
#############################################################


################################
### O'Brien_Fleming boundary ###
################################

########################################################
#' Phase 2/3 seamless combination Group sequential design: Survival analysis, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param HR Hazard ratio (test vs control)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Randomization ratio, sample size of the test arm over the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'Combined number of events between a dose and control',
#' 'total number of events from all doses and control',
#' 'combined number of events in each dose and control in a two arm study with same critical boudary',
#' 'type I error control', 'exit probabilities', 'reference'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR = 0.7, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1,
#' rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_OF_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #HR is the vector of hazards ratio of test arms vs. control arm
  #length(HR)>=2
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
  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  ##################################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))

  #### sample size for specified distribution ####

  ########## two stage sample size###
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  evt_2=ceiling(total_time/(p_c*p_t))
  # n_tot=ceiling((n_2/2)*(length(mu)+1))
  evt_tot=ceiling(evt_2*((s[1]/2)*(comp+1)+s[length(s)]-s[1]))
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  evt_individual_2=ceiling(info_time_ind/(p_c*p_t))
  ##################################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(evt_2,evt_tot,evt_individual_2,bry_out,info_time_00,ref)
  names(out)=c("Combined number of events between a dose and control","total number of events from all doses and control","combined number of events in each dose and control in a two arm study with same critical boudary","type I error control","exit probabilities","reference")
  return(out)
}

####################################################################
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2)
###########################################################
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(1.1,1.2,1.25),margin=1,direction=1,s=c(0.3,0.7,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
##########################################################
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.2,0.3,0.4),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),0.025,0.1,1,2,1)
#xx=two_stage_Sample_size_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),0.025,0.1,1,2,1)
#print(xx)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(c(0.9,0.8),1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(c(0.8,0.8),1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(c(0.6,0.8),1,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
### Jounce###
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.55,0.6),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.55,0.55),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.55,0.65),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.55,0.7),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.55,0.75),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.55,0.8),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.6,0.6),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.6,0.65),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.6,0.7),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.6,0.75),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.6,0.8),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.65,0.65),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.65,0.7),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.65,0.75),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.65,0.8),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.7,0.7),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.7,0.75),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.7,0.8),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.75,0.75),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.75,0.8),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
#two_stage_Sample_size_OF_boundary_survival_Schoenfeld(HR=c(0.8,0.8),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,beta=0.1,rand_ratio=1,dist1=2,dist2=1)
######################################################################

#' Phase 2/3 seamless combination Group sequential design: Survival analysis, Power calculation, O'Brien-Fleming boundary
#'
#' @param HR Hazard ratio (test vs control)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param events_tot Total number of events
#' @param rand_ratio Randomization ratio, sample size of the test arm over the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'type I error control', 'reference'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_OF_boundary_survival_Schoenfeld(HR = c(0.55,0.55), margin = 1, direction = 2,
#' s = c(1/3,2/3,1), alpha = 0.025,
#' events_tot = 137, rand_ratio = 1, dist1 = 2, dist2 = 1)
#'
#' @export
two_stage_power_OF_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,events_tot,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for given HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #HR is the vector of hazards ratio of test arms vs. control arm
  #length(HR)>=2
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
  events_2=events_tot/((s[1]/2)*(comp+1)+(s[length(s)]-s[1]))
  T=events_2*(p_c*p_t)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))

  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################################

  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power_out,bry_out,ref)
  names(out)=c("power","type I error control","reference")
  return(out)
}

####################################################################
#two_stage_power_OF_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,events_tot,rand_ratio,dist1,dist2)
#xx=two_stage_power_OF_boundary_survival_Schoenfeld(HR=c(0.2,0.3,0.4),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,events_tot=32,rand_ratio=1,dist1=1,dist2=1)
#two_stage_power_OF_boundary_survival_Schoenfeld(HR=c(0.55,0.55),margin=1,direction=2,s=c(1/3,2/3,1),alpha=0.025,events_tot=137,rand_ratio=1,dist1=2,dist2=1)
#xx=two_stage_power_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),0.025,23,1,1,1)
#print(xx)
#two_stage_power_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1,2,c(1/3,2/3,1),0.025,24,1,2,1)
#two_stage_power_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,1,c(1/3,2/3,1),0.025,23,1,1,1)
#two_stage_power_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),0.025,24,1,2,1)
#two_stage_power_OF_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,1,c(1/3,2/3,1),0.025,24,1,2,1)
#two_stage_power_OF_boundary_survival_Schoenfeld(HR=c(1.1,1.2,1.25),margin=1,direction=1,s=c(0.3,0.7,1),alpha=0.025,events_tot=1511,rand_ratio=1,dist1=2,dist2=1)
##########################################

#################################
### Alpha -spending boundary#####
#################################

######################################################################
#' Phase 2/3 seamless combination Group sequential design: Survival analysis, Sample size calculation, Alpha-spending boundary
#'
#' @param HR Hazard ratio (test vs control)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Randomization ratio, sample size of the test arm over the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'Combined number of events between a dose and control',
#' 'total number of events from all doses and control',
#' 'combined number of events in each dose and control in a two arm study with same critical boudary',
#' 'type I error control', 'exit probabilities', 'reference'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_alsp_boundary_survival_Schoenfeld(HR = 0.7, margin = 0.1, direction = 0,
#' s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_alsp_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #HR is the vector of hazards ratio of test arms vs. control arm
  #length(HR)>=2
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
  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  ##################################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))

  #### saple size for specified distribution ####

  ########## two stage sample size###
  p_c=1/(1+rand_ratio)
  p_t=1-p_c
  evt_2=ceiling(total_time/(p_c*p_t))
  #   n_tot=ceiling((n_2/2)*(length(mu)+1))
  evt_tot=ceiling(evt_2*((s[1]/2)*(comp+1)+s[length(s)]-s[1]))
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  evt_individual_2=ceiling(info_time_ind/(p_c*p_t))
  ##################################################
  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(evt_2,evt_tot,evt_individual_2,bry_out,info_time_00,ref)
  names(out)=c("Combined number of events between a dose and control","total number of events from all doses and control","combined number of events in each dose and control in a two arm study with same critical boudary","type I error control","exit probabilities","reference")
  return(out)
}

####################################################################
#two_stage_Sample_size_alsp_boundary_survival_Schoenfeld(HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2)
#two_stage_Sample_size_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#two_stage_Sample_size_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#two_stage_Sample_size_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,2,1)
#two_stage_Sample_size_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,2,1)
######################################################################

#' Phase 2/3 seamless combination Group sequential design: Survival analysis, Power calculation, Alpha-spending boundary
#'
#' @param HR Hazard ratio (test vs control)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger HR is better, 0 = smaller HR is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param events_tot Total number of events
#' @param rand_ratio Randomization ratio, sample size of the test arm over the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'type I error control', 'reference'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_alsp_boundary_survival_Schoenfeld(HR = c(0.2,0.3,0.4), margin = 1.1,
#' direction = 1, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), events_tot = 20,
#' rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_alsp_boundary_survival_Schoenfeld<-function(HR,margin,direction,s,alpha,events_tot,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for given HR,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #HR is the vector of hazards ratio of test arms vs. control arm
  #length(HR)>=2
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
  events_2=events_tot/((s[1]/2)*(comp+1)+(s[length(s)]-s[1]))
  T=events_2*(p_c*p_t)

  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))

  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################################

  ref="Schoenfeld DA. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics. 39:499-503. "
  out=list(power_out,bry_out,ref)
  names(out)=c("power","type I error control","reference")
  return(out)
}

####################################################################
#two_stage_power_alsp_boundary_survival_Schoenfeld(HR,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
#two_stage_power_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),20,1,1,1)
#two_stage_power_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),20,1,1,1)
#two_stage_power_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),21,1,2,1)
#two_stage_power_alsp_boundary_survival_Schoenfeld(c(0.2,0.3,0.4),1.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),21,1,2,1)
##########################################
#############################################################################
#### Part III.4 power and Sample size clculation: Poisson distributions  ####
#############################################################################

#############################################################
####    supeiority trials and non-inferiority trials  #######
#############################################################


################################
### O'Brien_Fleming boundary ###
################################

#######################################################
#' Phase 2/3 seamless combination Group sequential design: Poisson distribution, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param lambda_t Event rate in the test arm (Poisson / negative binomial endpoint)
#' @param lambda_c Event rate in the control arm (Poisson / negative binomial endpoint)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'test arm sample size', 'control arm sample size',
#' 'sample size for each dose in a two arm study with same critical boudary',
#' 'sample size for control arm in a two arm study with same critical boudary', 'type I error control', 'exit probabilities'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_OF_boundary_poisson(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025,
#' beta = 0.1, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_OF_boundary_poisson<-function(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #lambda_t is the vector of rates of the test arm
  #length(lambda_t)>=2
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
  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  ######################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))

  #### saple size for specified distribution ####

  ########## two stage sample size###

  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(total_time*max(ratio))
  n_t=n_c*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio

  ###################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,info_time_00)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm in a two arm study with same critical boudary","type I error control","exit probabilities")
  return(out)
}

####################################################################
#two_stage_Sample_size_OF_boundary_poisson(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2)
#two_stage_Sample_size_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.15,1,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.15,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
#two_stage_Sample_size_OF_boundary_poisson(c(0.2,0.3,0.4),0,0.15,1,c(1/3,2/3,1),0.025,0.1,1,1,1)

#two_stage_Sample_size_OF_boundary_poisson(c(0.2,0.3,0.4),0,0.3,2,c(1/3,2/3,1),0.025,0.1,1,1,1)
######################################################################
####################################################
#' Phase 2/3 seamless combination Group sequential design: Poisson distribution, Power calculation, O'Brien-Fleming boundary
#'
#' @param lambda_t Event rate in the test arm (Poisson / negative binomial endpoint)
#' @param lambda_c Event rate in the control arm (Poisson / negative binomial endpoint)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'type I error control'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_OF_boundary_poisson(lambda_t = c(0.2,0.3,0.4), lambda_c = 0.5, margin = 0.12,
#' direction = 1, s = c(1/3,2/3,1), alpha = 0.025, sampsz_control = 150, rand_ratio = 1,
#' dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_OF_boundary_poisson<-function(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #lambda_t is the vector of rates of the test arm
  #length(lambda_t)>=2
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
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################################
  out=list(power_out,bry_out)
  names(out)=c("power","type I error control")
  return(out)
}

####################################################################
#two_stage_power_OF_boundary_poisson(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
#two_stage_power_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.12,1,c(1/3,2/3,1),0.025,150,1,1,1)
#two_stage_power_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.12,2,c(1/3,2/3,1),0.025,150,1,1,1)
#two_stage_power_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.15,1,c(1/3,2/3,1),0.025,4521,1,1,1)
#two_stage_power_OF_boundary_poisson(c(0.2,0.3,0.4),0.5,0.15,2,c(1/3,2/3,1),0.025,63,1,1,1)
##########################################

#################################
### Alpha -spending boundary#####
#################################

######################################################################
#' Phase 2/3 seamless combination Group sequential design: Poisson distribution, Sample size calculation, Alpha-spending boundary
#'
#' @param lambda_t Event rate in the test arm (Poisson / negative binomial endpoint)
#' @param lambda_c Event rate in the control arm (Poisson / negative binomial endpoint)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'test arm sample size', 'control arm sample size',
#' 'sample size for each dose in a two arm study with same critical boudary',
#' 'sample size for control arm in a two arm study with same critical boudary', 'type I error control', 'exit probabilities'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_alsp_boundary_poisson(lambda_t = 1.1, lambda_c = 1.0, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025,
#' beta = 0.1, rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_alsp_boundary_poisson<-function(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and saple size for lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #lambda_t is the vector of rates of the test arm
  #length(lambda_t)>=2
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
  if (max(mu)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))

  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  ######################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))

  #### saple size for specified distribution ####

  ########## two stage sample size###

  ratio=lambda_c+lambda_t/rand_ratio
  n_c=ceiling(total_time*max(ratio))
  n_t=n_c*rand_ratio
  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  n_c_individual=ceiling(info_time_ind*ratio)
  n_t_individual=n_c_individual*rand_ratio

  ###################################################
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,info_time_00)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm in a two arm study with same critical boudary","type I error control","exit probabilities")
  return(out)
}

####################################################################
#two_stage_Sample_size_alsp_boundary_poisson(lambda_t,lambda_c,margin,direction,s,alpha,beta,rand_ratio,dist1,dist2)
#two_stage_Sample_size_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.3,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
#two_stage_Sample_size_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.15,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1)
##########################################################################
#' Phase 2/3 seamless combination Group sequential design: Poisson distribution, Power calculation, Alpha-spending boundary
#'
#' @param lambda_t Event rate in the test arm (Poisson / negative binomial endpoint)
#' @param lambda_c Event rate in the control arm (Poisson / negative binomial endpoint)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger lambda is better, 0 = smaller lambda is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'type I error control'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_alsp_boundary_poisson(lambda_t = c(0.2,0.3,0.4), lambda_c = 0.5, margin = 0.12,
#' direction = 1, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), sampsz_control = 150,
#' rand_ratio = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_alsp_boundary_poisson<-function(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
  ### This function calculates the critical boundary and power for lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #lambda_t is the vector of rates of the test arm
  #length(lambda_t)>=2
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
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))
  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################################
  out=list(power_out,bry_out)
  names(out)=c("power","type I error control")
  return(out)
}

####################################################################
#two_stage_power_alsp_boundary_poisson(lambda_t,lambda_c,margin,direction,s,alpha,sampsz_control,rand_ratio,dist1,dist2){
#two_stage_power_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.12,1,c(1/3,2/3,1),c(0.001,0.001,0.023),150,1,1,1)
#two_stage_power_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.12,2,c(1/3,2/3,1),c(0.001,0.001,0.023),150,1,1,1)
#two_stage_power_alsp_boundary_poisson(c(0.2,0.3,0.4),0.5,0.3,1,c(1/3,2/3,1),c(0.001,0.001,0.023),336,1,1,1)

#######################################################################################
#### Part III.5 power and Sample size clculation: negative binomial distributions  ####
#######################################################################################

#############################################################
####    supeiority trials and non-inferiority trials  #######
#############################################################

################################
### O'Brien_Fleming boundary ###
################################

#########################################################
#' Phase 2/3 seamless combination Group sequential design: Negative binomial distribution, Sample size calculation, O'Brien-Fleming boundary
#'
#' @param r_t Rate parameter of the test arm (negative binomial endpoint)
#' @param r_c Rate parameter of the control arm (negative binomial endpoint)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger r is better, 0 = smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t Follow-up duration of the test arm (negative binomial endpoint)
#' @param kappa Dispersion (size) parameter of the negative binomial distribution
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'test arm sample size', 'control arm sample size',
#' 'sample size for each dose in a two arm study with same critical boudary',
#' 'sample size for control arm in a two arm study with same critical boudary', 'type I error control', 'exit probabilities', 'Reference'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' \donttest{
#' two_stage_Sample_size_OF_boundary_neg_binomial(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, beta = 0.1, rand_ratio = 1, nu_t = 6,
#' kappa = 1, dist1 = 1, dist2 = 1)
#' }
#'
#' @export
two_stage_Sample_size_OF_boundary_neg_binomial<-function(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2 for
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the vector of the rates of the test arm
  #length(r_t)>=2
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
  if (max(mu)<1e-4|| min(r_1)<1e-4){return(NA)}

  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))

  #c("Critical boundary","1-sided alpha spending at each interim analysis","1-sided accumulative alpha","2-sided alpha spending at each interim analysis","2-sided accumulative alpha","Critical boundary-phase 3 only","nominal p-value for phase 3")

  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  #############################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))

  #### saple size for specified distribution ####
  ########## two stage sample size###
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(max(V))
  n_c=ceiling(total_time*max(V))
  n_t=n_c*rand_ratio

  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  n_c_individual=ceiling(info_time_ind*max(V))
  n_t_individual=n_c_individual*rand_ratio
  ###################################################
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,info_time_00,ref)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm in a two arm study with same critical boudary","type I error control","exit probabilities","Reference")
  return(out)
}

####################################################################
#two_stage_Sample_size_OF_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2)
#two_stage_Sample_size_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,1,c(1/3,2/3,1),0.025,0.1,1,1,1,1,1)
#two_stage_Sample_size_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1,1,1)
#two_stage_Sample_size_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),0.025,0.1,1,1,1,2,1)

######################################################################

#' Phase 2/3 seamless combination Group sequential design: Negative binomial distribution, Power calculation, O'Brien-Fleming boundary
#'
#' @param r_t Rate parameter of the test arm (negative binomial endpoint)
#' @param r_c Rate parameter of the control arm (negative binomial endpoint)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger r is better, 0 = smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t Follow-up duration of the test arm (negative binomial endpoint)
#' @param kappa Dispersion (size) parameter of the negative binomial distribution
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'Critical boundaries', 'Critical boundaries for phase 3 only',
#' '1-sided alpha spending', '1-sided cummulative alpha', '2-sided alpha spending', '2-sided cummulative alpha',
#' 'selecting prob at end of phase 2', 'total selection prob at end of phase 2', 'Rejection rate for ech dose', 'total rejection rate',
#' 'conditional rejection rate for each selected dose', 'critical boundary for phase 3 stage', 'Nomimal p-value for phase 3 stage'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_OF_boundary_neg_binomial(r_t = c(1.1, 1.1, 1.1), r_c = 1.1, margin = 0.1,
#' direction = 0, s = c(1/3, 2/3, 1), alpha = 0.025, sampsz_control = 100, rand_ratio = 1,
#' nu_t = 6, kappa = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_OF_boundary_neg_binomial<-function(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2){
  ### This function calculates the critical boundary and power for given r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the vector of the rates of the test arm
  #length(r_t)>=2
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
  bry_out=two_stage_OF_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))

  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################################

  # out=list(power_tot,c_bry,c_bry_1,alpha_s_1_side,cumu_alpha_1_side,alpha_s_2_side,cumu_alpha_2_side,select_prob_0,tot_select_prob_0,exit_prob_vector,exit_prob_tot,cond_exit_prob_vector, c_bry_1, nom_p_3)
  #  names(out)=c("power","Critical boundaries","Critical boundaries for phase 3 only","1-sided alpha spending","1-sided cummulative alpha","2-sided alpha spending","2-sided cummulative alpha","selecting prob at end of phase 2","total selection prob at end of phase 2","Rejection rate for ech dose","total rejection rate","conditional rejection rate for each selected dose","critical boundary for phase 3 stage","Nomimal p-value for phase 3 stage")
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power_out,bry_out,ref)
  names(out)=c("power","type I error control","Reference")
  return(out)
}

####################################################################
#two_stage_power_OF_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2){
#two_stage_power_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,1,c(1/3,2/3,1),0.025,111,1,1,1,1,1)
#two_stage_power_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),0.025,111,1,1,1,1,1)
#two_stage_power_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),0.025,67,1,1,1,1,1)
#two_stage_power_OF_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),0.025,69,1,1,1,2,1)

######################################################################

#################################
### Alpha -spending boundary#####
#################################

###########################################################
#' Phase 2/3 seamless combination Group sequential design: Negative binomial distribution, Sample size calculation, Alpha-spending boundary
#'
#' @param r_t Rate parameter of the test arm (negative binomial endpoint)
#' @param r_c Rate parameter of the control arm (negative binomial endpoint)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger r is better, 0 = smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param beta Type II error rate, i.e. 1 - power
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t Follow-up duration of the test arm (negative binomial endpoint)
#' @param kappa Dispersion (size) parameter of the negative binomial distribution
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'test arm sample size', 'control arm sample size',
#' 'sample size for each dose in a two arm study with same critical boudary', 'Critical boundaries', 'Critical boundaries for phase 3 only',
#' '1-sided alpha spending', '1-sided cummulative alpha', '2-sided alpha spending', '2-sided cummulative alpha',
#' 'selecting prob at end of phase 2', 'total selection prob at end of phase 2', 'Rejection rate for ech dose',
#' 'total rejection rate', 'conditional rejection rate for each selected dose', 'critical boundary for phase 3 stage',
#' 'Nomimal p-value for phase 3 stage'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_Sample_size_alsp_boundary_neg_binomial(r_t = c(0.2,0.3,0.4), r_c = 0.6, margin = 0.1,
#' direction = 1, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.023), beta = 0.1, rand_ratio = 1,
#' nu_t = 1, kappa = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_Sample_size_alsp_boundary_neg_binomial<-function(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2){
  ### This function calculates the critical boundary and sample size for given r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2
  # the critical boundary is O'Brien-Fleming boundary
  #r_t is the vector of the rates of the test arm
  #length(r_t)>=2
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
  r=rand_ratio
  if (direction==1){
    r_1=r_t+margin
    mu=log((r_1)/r_c)
  }else{
    r_1=r_t-margin
    mu=-log((r_1)/r_c)
  }
  if (max(mu)<1e-4|| min(r_1)<1e-4){return(NA)}

  theta=rand_ratio
  comp=length(mu)
  ##### Step 1.1: boundary calcualtion -OF boundary###
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))

  #c("Critical boundary","1-sided alpha spending at each interim analysis","1-sided accumulative alpha","2-sided alpha spending at each interim analysis","2-sided accumulative alpha","Critical boundary-phase 3 only","nominal p-value for phase 3")

  #######Step 1.2: Common to all boundary types ######
  info_time_out=two_stage_info_time(mu,s,c_bry,beta,alpha,dist2)
  #############################################################
  info_time_00=info_time_out[-c(1:2)]
  total_time=as.numeric(unlist(info_time_out[1]))
  info_time_ind=as.numeric(unlist(info_time_out[2]))

  #### saple size for specified distribution ####
  ########## two stage sample size###
  V=(1/nu_t)*(1/r_0+1/(theta*r_1))+(1+theta)*kappa/theta
  sigma=sqrt(max(V))
  n_c=ceiling(total_time*max(V))
  n_t=n_c*rand_ratio

  ################################################################
  ### individual sape size per arm using same critical boundary ##
  ###          in a two arm sequential design                   ##
  ################################################################
  n_c_individual=ceiling(info_time_ind*max(V))
  n_t_individual=n_c_individual*rand_ratio

  #  out=list(n_t,n_c,n_t_individual,c_bry,c_bry_1,alpha_s_1_side,cumu_alpha_1_side,alpha_s_2_side,cumu_alpha_2_side,select_prob_0,tot_select_prob_0,exit_prob_vector,exit_prob_tot,cond_exit_prob_vector, c_bry_1, nom_p_3)
  #  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","Critical boundaries","Critical boundaries for phase 3 only","1-sided alpha spending","1-sided cummulative alpha","2-sided alpha spending","2-sided cummulative alpha","selecting prob at end of phase 2","total selection prob at end of phase 2","Rejection rate for ech dose","total rejection rate","conditional rejection rate for each selected dose","critical boundary for phase 3 stage","Nomimal p-value for phase 3 stage")
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(n_t,n_c,n_t_individual,n_c_individual,bry_out,info_time_00,ref)
  names(out)=c("test arm sample size","control arm sample size","sample size for each dose in a two arm study with same critical boudary","sample size for control arm in a two arm study with same critical boudary","type I error control","exit probabilities","Reference")
  return(out)
}

####################################################################
#two_stage_Sample_size_alsp_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,beta,rand_ratio,nu_t,kappa,dist1,dist2)
#two_stage_Sample_size_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1,1,1)
#two_stage_Sample_size_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),0.1,1,1,1,1,1)

#########################################################################

#' Phase 2/3 seamless combination Group sequential design: Negative binomial distribution, Power calculation, Alpha-spending boundary
#'
#' @param r_t Rate parameter of the test arm (negative binomial endpoint)
#' @param r_c Rate parameter of the control arm (negative binomial endpoint)
#' @param margin Non-inferiority / superiority margin
#' @param direction Test direction: 1 = larger r is better, 0 = smaller r is better
#' @param s Vector of information fractions of each analysis, e.g. c(1/3, 2/3, 1); the last element must be 1
#' @param alpha One-sided significance level (type I error rate)
#' @param sampsz_control Control arm sample size
#' @param rand_ratio Ratio of the number of subjects in the test arm vs. the number of subjects in the control arm
#' @param nu_t Follow-up duration of the test arm (negative binomial endpoint)
#' @param kappa Dispersion (size) parameter of the negative binomial distribution
#' @param dist1 Endpoint / distribution indicator for the first stage
#' @param dist2 Endpoint / distribution indicator for the second stage
#'
#' @return A named list containing the following elements: 'power', 'Critical boundaries', 'Critical boundaries for phase 3 only',
#' '1-sided alpha spending', '1-sided cummulative alpha', '2-sided alpha spending', '2-sided cummulative alpha',
#' 'selecting prob at end of phase 2', 'total selection prob at end of phase 2', 'Rejection rate for ech dose', 'total rejection rate',
#' 'conditional rejection rate for each selected dose', 'critical boundary for phase 3 stage', 'Nomimal p-value for phase 3 stage'.
#'
#' @references P. Gao, Y. Li (05 May 2024): Adaptive two-stage seamless sequential design for clinical trials, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2342518. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2342518}
#' @references Gao, P., Zhang, W. (2024). A systematic approach to adaptive sequential design for clinical trials: using simulations to select a design with desired operating characteristics. Journal of Biopharmaceutical Statistics,2024 Aug;34(5):737-752. doi:10.1080/10543406.2024.2358796. Epub 2024 May 30. To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2358796}
#' @references H.H. Müller, H. Sch?fer. Adaptive Group Sequential Designs for Clinical Trials: Combining the Advantages of Adaptive and of Classical Group Sequential Approaches. BIOMETRICs 57, 886-891. September 2001.
#' @references P. Gao, J. H. Ware, C. Mehta. Sample size re-estimation for adaptive sequential design in clinical trials. J Biopharm Stat. 2008;18(6):1184-96.
#'
#' @examples
#' two_stage_power_alsp_boundary_neg_binomial(r_t = c(0.2,0.3,0.4), r_c = 0.6, margin = 0.1,
#' direction = 1, s = c(1/3,2/3,1), alpha = c(0.001,0.001,0.048), sampsz_control = 111,
#' rand_ratio = 1, nu_t = 1, kappa = 1, dist1 = 1, dist2 = 1)
#'
#' @export
two_stage_power_alsp_boundary_neg_binomial<-function(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2){
  ### This function calculates the critical boundary and power for given r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2
  # the critical boundary is alpha-spending boundary
  #r_t is the vector of the rates of the test arm
  #length(r_t)>=2
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
  bry_out=two_stage_alsp_bundry(s,alpha,comp,dist1)
  ################################################
  c_bry=as.numeric(unlist(bry_out[1]))

  ###############Power calculation ##########
  power_out=two_stage_time_power(mu,s,c_bry,T,alpha,dist2)
  ##########################################################

  # out=list(power_tot,c_bry,c_bry_1,alpha_s_1_side,cumu_alpha_1_side,alpha_s_2_side,cumu_alpha_2_side,select_prob_0,tot_select_prob_0,exit_prob_vector,exit_prob_tot,cond_exit_prob_vector, c_bry_1, nom_p_3)
  #  names(out)=c("power","Critical boundaries","Critical boundaries for phase 3 only","1-sided alpha spending","1-sided cummulative alpha","2-sided alpha spending","2-sided cummulative alpha","selecting prob at end of phase 2","total selection prob at end of phase 2","Rejection rate for ech dose","total rejection rate","conditional rejection rate for each selected dose","critical boundary for phase 3 stage","Nomimal p-value for phase 3 stage")
  ref="H. Zhu,H. Lakkis.Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33(3):376-387  February 2014"
  out=list(power_out,bry_out,ref)
  names(out)=c("power","type I error control","Reference")
  return(out)
}

  ####################################################################
#two_stage_power_alsp_boundary_neg_binomial(r_t,r_c,margin,direction,s,alpha,sampsz_control,rand_ratio,nu_t,kappa,dist1,dist2){
#two_stage_power_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,1,c(1/3,2/3,1),c(0.001,0.001,0.048),111,1,1,1,1,1)
#two_stage_power_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.048),111,1,1,1,1,1)
#two_stage_power_alsp_boundary_neg_binomial(c(0.2,0.3,0.4),0.6,0.1,2,c(1/3,2/3,1),c(0.001,0.001,0.023),66,1,1,1,1,1)
