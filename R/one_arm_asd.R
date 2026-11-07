#rm(list = ls())

#library(clinfun)

###################################################################
####### part 1: Sample size for Simon's two-stage optimal design###
###################################################################

ph2simon <- function(pu, pa, ep1, ep2, nmax = 100) {
  if(nmax > 1000) stop("nmax cannot exceed 1000")
  nmax1 = nmax + 1
  m = (nmax * (nmax + 3))/2
  n = rep(1:nmax, (1:nmax) + 1)
  x = n * 0
  j = 1
  for(i in 1:nmax) {
    x[j + (0:i)] = 0:i
    j = j + i + 1
  }
  p0 = stats::dbinom(x, n, pu)
  p1 = stats::dbinom(x, n, pa)
  cp0 = 1 - stats::pbinom(x - 1, n, pu)
  cp1 = 1 - stats::pbinom(x - 1, n, pa)
  sym_f2bdry = getNativeSymbolInfo("f2bdry",PACKAGE="clinfun")
  zz = .Fortran(sym_f2bdry,
                 as.integer(m),
                 as.integer(nmax),
                 ep1,
                 ep2,
                 p0,
                 p1,
                 cp0,
                 cp1,
                 bdry = as.integer(rep(0, nmax*4)),
                 peten = rep(0, nmax*2),
                 as.integer(nmax1),
                 rep(0,nmax1),
                 rep(0,nmax1))
  ph2 = list()
  ph2out = cbind(matrix(zz$bdry,nmax,4),matrix(zz$peten,nmax,2))
  ph2out = ph2out[ph2out[,5]!=0,]
  if(nrow(ph2out)==0) {
    errmesg = paste("  No feasible solution found. \n\tIncrease maximum sample size.  Current nmax value = ",nmax,".",sep="")
    stop(message=errmesg)
  }
  dimnames(ph2out) = list(NULL,c("r1", "n1", "r", "n", "EN(p0)", "PET(p0)"))
  ph2$pu = pu
  ph2$pa = pa
  ph2$alpha = ep1
  ph2$beta = ep2
  ph2$out = ph2out
  ph2$nmax = nmax
  class(ph2) = "ph2simon"
  return(ph2)
}

#' Print method for a Simon's two-stage phase II design
#'
#' @param x An object of class `ph2simon`, as returned by the internal `ph2simon` helper
#' @param ... Further arguments passed to `print`
#'
#' @return `NULL`, invisibly, or the warning message, invisibly, if the optimal sample size
#'   is too close to `nmax`. Called for the side effect of printing the optimal and minimax
#'   designs.
#'
#' @export
print.ph2simon <- function(x, ...) {
  xout = x$out
  nmax = x$nmax
  n = nrow(xout)
  nopt = ((1:n)[xout[,5]==min(xout[,5])])[1]
  xopt = xout[c(nopt,1),]
  dimnames(xopt)[[1]] = c("Optimal","Minimax")
  cat("\n Simon 2-stage Phase II design \n\n")
  cat("Unacceptable response rate: ",x$pu,"\n")
  cat("Desirable response rate: ",x$pa,"\n")
  cat("Error rates: alpha = ",x$alpha,"; beta = ",x$beta,"\n\n")
  print(xopt, digits = 4, ...)
  cat("\n")
  if(xopt[1,4]>nmax-10) {
    warning(paste("  Optimal sample size too close to nmax. \n  Try increasing nmax (current value = ",nmax,")\n",sep=""))
  }
}

#' Plot method for a Simon's two-stage phase II design
#'
#' @param x An object of class `ph2simon`, as returned by the internal `ph2simon` helper
#' @param ... Further graphical arguments passed to `plot`
#'
#' @return `NULL`, invisibly. Called for the side effect of drawing the expected-sample-size curve.
#'
#' @export
plot.ph2simon <- function(x, ...) {
  xout = x$out
  n = nrow(xout)
  nopt = ((1:n)[xout[,5]==min(xout[,5])])[1]
  nopt1 = min(nopt+5,n)
  plot(xout[1:nopt1,4],xout[1:nopt1,5],type="l",xlab="Maximum number of patients",ylab="Expected trial size", ...)
  graphics::points(xout[1,4],xout[1,5],pch="M")
  graphics::points(xout[nopt,4],xout[nopt,5],pch="O")
}

oc.twostage.bdry <- function(pu, pa, r1, n1, r, n){
  pet = err0 = stats::pbinom(r1,n1,pu)
  ess = n1 + (n-n1)*(1-pet)
  err1 = stats::pbinom(r1,n1,pa)
  for(i in (r1+1):r) {
    err0 = err0 + stats::dbinom(i,n1,pu)*stats::pbinom(r-i,n-n1,pu)
    err1 = err1 + stats::dbinom(i,n1,pa)*stats::pbinom(r-i,n-n1,pa)
  }
  out = c(1-err0, 1-err1, pet, ess)
  names(out) = c("P(reject H0 | p0)","P(reject H0 | p1)","PET(p0)","EN(p0)")
  out
}

twostage.inference <- function(x, r1, n1, n, pu, alpha=0.05) {
  out = list()
  n2 = n - n1
  # UMVUE from Jung & Kim (2004) Stats in Medicine
  xlo = max(r1 + 1, x - n2)
  xhi = min(x, n1)
  x1 = xlo:xhi
  if (x <= r1) pumvue = x/n1
  else pumvue = sum(choose(n1-1, x1-1)*choose(n2, x-x1))/sum(choose(n1, x1)*choose(n2, x-x1))
  # p-value and CI from Koyama & Chen (2008) Stats in Medicine
  x1 = (r1+1):n1
  if (x <= r1) p.value = 1 - stats::pbinom(x-1, n1, pu)
  else p.value = sum(stats::dbinom(x1, n1, pu)*(1-stats::pbinom(x-x1-1,n2,pu)))
  # CI steps: first bracket the LCL and UCL
  pp = seq(0, 1, by=0.01)
  pval = rep(0, 101)
  if (x <= r1) {
    pval = 1 - stats::pbinom(x-1, n1, pp)
  } else {
    x2 = x - x1 - 1
    pval = sapply(pp, function(p, x1, x2, n1, n2) {
      sum(stats::dbinom(x1, n1, p)*(1-stats::pbinom(x2,n2,p)))
    }, x1, x2, n1, n2)
  }
  # LCL & UCL indices
  ii = which(pval > alpha)[1]-1
  jj = which(pval >= 1-alpha)[1]-1
  # LCL refinement
  pp0 = pp[ii] + pp/100
  pval = rep(0, 101)
  if (x <= r1) {
    pval = 1 - stats::pbinom(x-1, n1, pp0)
  } else {
    x2 = x - x1 - 1
    pval = sapply(pp0, function(p, x1, x2, n1, n2) {
      sum(stats::dbinom(x1, n1, p)*(1-stats::pbinom(x2,n2,p)))
    }, x1, x2, n1, n2)
  }
  LCL = pp0[which(pval > alpha)[1]]
  # UCL refinement
  pp0 = pp[jj] + pp/100
  pval = rep(0, 101)
  if (x <= r1) {
    pval = 1 - stats::pbinom(x-1, n1, pp0)
  } else {
    x2 = x - x1 - 1
    pval = sapply(pp0, function(p, x1, x2, n1, n2) {
      sum(stats::dbinom(x1, n1, p)*(1-stats::pbinom(x2,n2,p)))
    }, x1, x2, n1, n2)
  }
  UCL = pp0[which(pval >= 1-alpha)[1]-1]
  out = c(pumvue, p.value, LCL, UCL)
  names(out) = c("pumvue", "p.value", paste(100*alpha, "% LCL", sep=""), paste(100*(1-alpha), "% UCL", sep=""))
  out
}

# function to get all the admisssible 2 stage designs
twostage.admissible  <- function(x) {
  xout = x$out
  nmax = x$nmax
  nopt = which(xout[,5]==min(xout[,5]))[1]
  # admissible indicator for the rows
  adm = rep(0, nopt)
  adm[1] = 1 # minimax
  adm[nopt] = 1 # optimal
  mss = xout[1:nopt, 4] # max sample size
  ess = xout[1:nopt, 5] # expected sample size
  # search admissible points only if search space is not a singleton
  istart = 1
  while (nopt > istart) {
    # slope of each point from the minimax
    slopes = (ess[(istart+1):nopt] - ess[istart])/(mss[(istart+1):nopt]-mss[istart])
    imin = which.min(slopes)[1]
    adm[istart+imin] = 1
    istart = istart+imin
  }
  xout = xout[1:nopt,][adm==1,]
  # find weight probability at which each rule is admissible
  dmss = diff(xout[,4])
  dess = diff(xout[,5])
  qwt = round(dess/(dess - dmss),3)
  # admissible designs
  xout = cbind(xout, c(qwt, 0), c(1, qwt))
  colnames(xout)[7:8] = c("qLo", "qHi")
  rownames(xout) = c("Minimax", rep("Admissible", nrow(xout)-2), "Optimal")
  return(xout)
}


####################################
#ph2simon(pu, pa, ep1, ep2, nmax = 100)
#ph2simon(pu=0.2, pa=0.4, ep1=0.025, ep2=0.1)
#ph2simon(pu=0.2, pa=0.5, ep1=0.025, ep2=0.1)
#ph2simon(pu=0.2, pa=0.35, ep1=0.05, ep2=0.05, nmax=150)
#ph2simon(pu=0.2, pa=0.35, ep1=0.05, ep2=0.05, nmax=300)
#######################################################

##############################################################################
#### part 2: Sample size calculation for two stage expanded Simon's Design ###
#####                The hybrid design and midpoint design                 ###
##############################################################################
###############################################

#' Single arm phase 2 oncology trial: Two-stage design, Simon's design
#'
#' @param alpha Type I error for the one-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param p_0 Response rate indicating low activity with insufficient clinical benefit
#' @param p_1 Assumed response rate
#'
#' @return A matrix with one row for each of the minimax design (smallest n_2), the optimal
#'   design (smallest EN(p_0)) and the average design (an augmentation of the original
#'   Simon's design, Gao and Zhang 2024), and one column per design quantity:
#' \describe{
#'   \item{r1}{If no more than r1 responses are observed at the first look, then the trial
#'     would be stopped for futility}
#'   \item{n1}{Number of patients at first look}
#'   \item{r_2}{The null hypothesis will be rejected if at least r_2 responses are observed
#'     at the end of the trial}
#'   \item{n_2}{Number of patients at second look; n_2 is also the total sample size}
#'   \item{EN(p_0)}{Expected sample size under the null hypothesis}
#'   \item{PET(p_0)}{Probability of early termination under the null hypothesis}
#'   \item{t1_err}{Type I error, the probability of rejecting the null hypothesis under
#'     p <= p_0}
#'   \item{PET_p_1}{Probability of early termination under p = p_1}
#'   \item{power_p_1}{Power under p = p_1}
#' }
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials:
#' an expansion of Simon's design, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#'
#' @examples
#' simon_ph2_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.33)
#'
#' @export
simon_ph2_sz<-function(alpha,beta,p_0,p_1){
  ##############
  a=stats::qnorm(1-sum(alpha))
  b=stats::qnorm(1-beta)

  ##################################################################################
  ####### obtaning n_1 and r_1 by using the fucntion ph2simon (package clinfun) ####
  ##################################################################################
  p_ave=(p_0+p_1)/2
  n_approx=ceiling(p_ave*(1-p_ave)*(a+b)^2*(p_1-p_0)^(-2))-1
  n_max=floor(1.5*n_approx)

  ph2simon_out =ph2simon(pu=p_0, pa=p_1, ep1=sum(alpha), ep2=beta, nmax =n_max)
  ph2simon_out_mat=matrix(unlist(ph2simon_out[5]),ncol=6)
  colnames(ph2simon_out_mat)=c("r1", "n1","r_2","n_2", "EN(p_0)","PET(p_0)")

  ### Minimax ###
  min_max=ph2simon_out_mat[ph2simon_out_mat[,4]==min(ph2simon_out_mat[,4]),]
  if (length((as.numeric(min_max)))==6){
    min_max_opt=min_max
  }else{
    min_max_opt=min_max[min_max[,5]==min(min_max[,5]),]
  }

  ### Optimal ####
  optim=ph2simon_out_mat[ph2simon_out_mat[,5]==min(ph2simon_out_mat[,5]),]
  ##### average   #####################

  r_1_av=ceiling(mean(c(min_max_opt[1],optim[1])))
  n_1_av=ceiling(mean(c(min_max_opt[2],optim[2])))
  av_0=Stg2_sz(n_1=n_1_av,r_1=r_1_av,d_12_cut=0,alpha=alpha,beta=beta,p_0=p_0,p_1=p_1)
  av_n_tot=sum(c(av_0[1],av_0[4]))
  #Stg2_sz(n_1=30,r_1=7,d_12_cut=2,alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4)
  av_PET_p0=stats::pbinom(r_1_av,size=n_1_av, prob=p_0)
  av_EN_P0=n_1_av+(1-av_PET_p0)*av_0[4]
  av_1=c(av_0[2],av_0[1],av_0[5],av_n_tot,av_EN_P0,av_PET_p0)

  ###############################
  out_mat_f=rbind(min_max_opt,optim,av_1)
  rownames(out_mat_f)=c("minimax","optimal","average")
  r_1=as.numeric(out_mat_f[,1])
  n_1=as.numeric(out_mat_f[,2])
  r_2=as.numeric(out_mat_f[,3])
  d_12=as.numeric(out_mat_f[,4])-n_1
  ######
  PET_p0=stats::pbinom(r_1,size=n_1, prob=p_0)
  PET_p_1=stats::pbinom(r_1,size=n_1,prob=p_1)
  power_p_1=NULL
  t1_err=NULL
  for (i in 1:3){
    power_p_1[i]=rej_prob_1(p_00=p_1,r_1=r_1[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
    t1_err[i]=rej_prob_1(p_00=p_0,r_1=r_1[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
  }
  En_p0=n_1+(1-stats::pbinom(r_1,size=n_1, prob=p_0))*d_12
  out=cbind(out_mat_f,t1_err,PET_p_1,power_p_1)
  return(out)
}
##############################################
#simon_ph2_sz(alpha,beta,p_0,p_1)
#simon_ph2_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4)
#simon_ph2_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.33)


########################################################
## Function fro various steps for two stage  design  ###
########################################################

##############################################################


##################################################################################
####### rej_prob_1=1-f(r_1,n_1,min(r_1_e,n_1),r_2,n_2,p_00) see manuscript ####
####### rej_prob_1=1-g(r_1,n_1,min(r_1_e,n_1),r_2,n_2,p_00) see manuscript ####
##################################################################################
rej_prob_1<-function(p_00,r_1,r_1_e=2^8,r_2,n_1,d_12){
  a_0=NULL
  for (i in ((r_1+1):min(c(r_1_e,n_1,r_2)))){
    a_0[i]=stats::dbinom(i, size=n_1, prob=p_00)*stats::pbinom(r_2-i, size=d_12, prob=p_00)
  }
  a=a_0[(r_1+1):length(a_0)]
  b=stats::pbinom(r_1, size=n_1, prob=p_00)+sum(a)
  e=1-b
  return(e)
}
####################################################
#rej_prob_1(p_00=0.4,r_1=6,r_1_e=2^8,r_2=17,n_1=22,d_12=20)
# rej_prob_0(p_00=0.5,r_1=2,r_2=9,n_1=11,n_2=26)
##############################################################
## choose crticical boundary at n_1 when n_1> n_beta ####

r_1_e_fun<-function(alpha,r_1,n_1,p_00){
  r_1_e_0=NULL
  for (i in r_1:n_1){
    r_1_e_0[i]=stats::pbinom(i,size=n_1, prob=p_00)
  }
  r_1_e_1=cbind(r_1:n_1,r_1_e_0[-c(1:(r_1-1))])
  r_1_e_2=r_1_e_1[(r_1_e_1[,2]>1-alpha[1]),]
  r_1_e=r_1_e_2[1,1]
  return(r_1_e)
}
###############################
#r_1_e_fun(alpha=c(0.001,0.024),r_1=2,n_1=15,p_00=0.2)

##########################

#############################
##### seraching for n_2, r_2, after n_1 and r_1 have been previously selected

Stg2_sz<-function(n_1,r_1,d_12_cut,alpha,beta,p_0,p_1){
  a_1=stats::qnorm(1-alpha)
  a_tot=stats::qnorm(1-sum(alpha))
  b=stats::qnorm(1-beta)
  n_1_beta=(1-p_1)*p_1*(a_tot+b)^2*(p_1-p_0)^(-2)
  if (n_1>=n_1_beta){
    ## choose crticical boundary at n_1
    r_1_e=r_1_e_fun(alpha=alpha,r_1=r_1,n_1=n_1,p_00=p_0)

    ######## into loop 1 to identify n_2 and r_2#####
    keep_max=0
    d_12=max(d_12_cut-1,0)
    out_1=NULL
    while (keep_max==0){
      d_12=d_12+1
      n_2=n_1+d_12
      t1_err=NULL
      power=NULL
      keep=NULL
      out_0=NULL
      PET_p1=NULL
      for (i in (r_1+1):n_2){
        t1_err[i]=round(rej_prob_1(p_00=p_0,r_1=r_1,r_1_e=r_1_e,r_2=i,n_1=n_1,d_12=d_12),4)
        PET_p1[i]=round(stats::pbinom(r_1,size=n_1, prob=p_1),4)
        power[i]=round(rej_prob_1(p_00=p_1,r_1=r_1,r_1_e=r_1_e, r_2=i,n_1=n_1,d_12=d_12),4)

        if ( (t1_err[i]<=sum(alpha))&&(power[i]>1-beta)){
          keep[i]=1
        }else {keep[i]=0}

      }
      ######################################
      t1_err=t1_err[(r_1+1):n_2]
      power=power[(r_1+1):n_2]
      PET_p1=PET_p1[(r_1+1):n_2]
      keep=keep[(r_1+1):n_2]
      r_2=c((r_1+1):n_2)
      ##########################################
      out_0a=cbind(n_1,r_1, r_1_e,d_12,r_2,n_2,t1_err,power,PET_p1,keep)
      #     out_00=out_0a[out_0a[,10]==1,]

      out_0=out_0a[out_0a[,10]==1,]
      out_1=rbind(out_1,out_0)

      keep_max=max(keep)
    } #end of while loop
    cul_num=length(out_1[1,])
    out_01=out_1[out_1[,cul_num]==1,]
    if (length(as.numeric(out_01))>cul_num){
      out=out_01[1,]
    }else{
      out=out_01
    }
    #  names(out)=c("n_1","r_1","r_1_e","d_12","r_2","n_2","t1_err","power","PET_p1")
  } else {
    r_1_e=2^15

    ######## into loop 2 to identify n_2 and r_2#####
    keep_max=0
    d_12=max(d_12_cut-1,0)
    out_1=NULL
    while (keep_max==0){
      d_12=d_12+1
      t1_err=NULL
      power=NULL
      keep=NULL
      out_0=NULL
      PET_p1=NULL
      r_2=NULL
      n_2=n_1+d_12
      for (i in ((r_1+1):n_2)){
        r_1_e=2^15
        t1_err[i]=round(rej_prob_1(p_00=p_0,r_1=r_1,r_1_e=r_1_e,r_2=i,n_1=n_1,d_12=d_12),4)
        PET_p1[i]=round(stats::pbinom(r_1,size=n_1, prob=p_1),4)
        power[i]=round(rej_prob_1(p_00=p_1,r_1=r_1,r_1_e=r_1_e,r_2=i,n_1=n_1,d_12=d_12),4)
        #  En_p0=n_1+(1-stats::pbinom(r_1,size=n_1, prob=p_0))*n_2
        #  PET_p0=stats::pbinom(r_1,size=n_1, prob=p_0)

        if ( (t1_err[i]<=sum(alpha))&&(power[i]>1-beta)){
          keep[i]=1
        }else {keep[i]=0}

        r_1_e=NA
        r_2[i]=i
      }
      out_0=cbind(n_1,r_1,r_1_e,d_12,r_2,n_2,t1_err,power,PET_p1,keep)[(r_1+1):n_2,]
      out_1=rbind(out_1,out_0)
      keep_max=max(keep[(r_1+1):n_2])
    } #end of while loop

    colnames(out_1)=c("n_1","r_1","r_1_e","d_12","r_2","n_2","t1_err","power_p","PET_p","keep")
    # out_0=out_1[out_1[,10]==1,1:9]
    # out_00=out_1[!is.na(out_1[,10]),]
    cul_num=length(out_1[1,])
    out_01=out_1[out_1[,cul_num]==1,]
    if (length(as.numeric(out_01))>cul_num){
      out=out_01[1,]
    }else{
      out=out_01
    }
  }
  # out=as.numeric(out_0)
  out=round(out[-cul_num],4)

  # rownames(out)=""
  return(out)
}
#########################################################
#Stg2_sz(n_1,r_1,d_12_cut,alpha,beta,p_0,p_1)
#Stg2_sz(n_1=11,r_1=2,d_12_cut=3,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.6)
#Stg2_sz(n_1=11,r_1=2,d_12_cut=3,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.55)
#Stg2_sz(n_1=30,r_1=2,d_12_cut=3,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.5)
#Stg2_sz(n_1=30,r_1=2,d_12_cut=2,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.5)
#Stg2_sz(n_1=12,r_1=2,d_12_cut=2,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.4)
#Stg2_sz(n_1=60,r_1=2,d_12_cut=2,alpha=0.025,beta=0.1,p_0=0.2,p_1=0.45)
#Stg2_sz(n_1=60,r_1=2,d_12_cut=2,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.45)
#Stg2_sz(n_1=100,r_1=25,d_12_cut=2,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.4)
#Stg2_sz(n_1=30,r_1=7,d_12_cut=2,alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4)
##########################

stg2_mid_mat<-function(r_01,n_1,r_02,n_2){
  r_011=ceiling(r_01*1.2)
  r_021=ceiling(r_02*1.2)
  mat_1=NULL
  for (i in (1:r_011)){
    r_1=i
    r_2=c(r_1:r_021)
    mat_01=cbind(r_1,r_2)
    mat_1=rbind(mat_1,mat_01)
  }
  mat_2=cbind(mat_1,n_1,n_2)
  return(mat_2)
}
###########
#stg2_mid_mat(r_01,n_1,r_02,n_2)
#stg2_mid_mat(r_01=2,n_1=11,r_02=4,n_2=22)
###########################

###########################
Stg2_mid_sz<-function(n_1,r_01,d_12_cut,alpha,beta,p_0,p_low,p_1){
  PET_p_low_0=round(stats::pbinom(1,size=n_1, prob=p_low),4)
  if (PET_p_low_0>beta){
    out="n_1 is too small. Impossible for PET_p_low<beta. Increase n_1"
  }else{
    a_1=stats::qnorm(1-alpha)
    a_tot=stats::qnorm(1-sum(alpha))
    b=stats::qnorm(1-beta)
    n_1_beta=(1-p_1)*p_1*(a_tot+b)^2*(p_1-p_0)^(-2)
    if (n_1>=n_1_beta){
      ## choose crticical boundary at n_1
      r_1_e_0=NULL
      for (i in r_01:n_1){
        r_1_e_0[i]=stats::pbinom(i,size=n_1, prob=p_0)
      }
      r_1_e_1=cbind(r_01:n_1,r_1_e_0[-c(1:(r_01-1))])
      r_1_e_2=r_1_e_1[(r_1_e_1[,2]>1-alpha[1]),]
      # cul_num=length(r_1_e_1[1,])
      if (length(as.numeric(r_1_e_2))>2){
        r_1_e=r_1_e_2[1,1]
      }else{ r_1_e=r_1_e_2[1]}
    }else {
      r_1_e=2^15
    }

    ######## into loop 1 to identify n_2 and r_2#####
    keep_max=0
    d_12=max(d_12_cut-1,0)
    out_1=NULL
    while (keep_max==0){
      d_12=d_12+1
      n_2=n_1+d_12
      r_02=floor(n_2*p_1)
      stg2_mat=stg2_mid_mat(r_01=r_01,n_1=n_1,r_02=r_02,n_2=n_2)
      mat_row_count=length(stg2_mat[,1])
      t1_err=NULL
      PET_p1=NULL
      power_p=NULL
      keep=NULL
      keep_0=NULL
      PET_p0=NULL
      PET_p_low=NULL
      EN_p0=NULL
      for (i in (1:mat_row_count)){
        t1_err[i]=round(rej_prob_1(p_00=p_0,r_1=stg2_mat[i,1],r_1_e=r_1_e,r_2=stg2_mat[i,2],n_1=n_1,d_12=d_12),4)
        PET_p_low[i]=round(stats::pbinom(stg2_mat[i,1],size=n_1, prob=p_low),4)
        power_p[i]=round(rej_prob_1(p_00=p_1,r_1=stg2_mat[i,1],r_1_e=r_1_e,r_2=stg2_mat[i,2],n_1=n_1,d_12=d_12),4)
        PET_p0[i]=round(stats::pbinom(stg2_mat[i,1],size=n_1, prob=p_0),4)
        EN_p0[i]=n_1+(1-stats::pbinom(stg2_mat[i,1],size=n_1, prob=p_0))*(n_2-n_1)
        if ( (t1_err[i]<=sum(alpha))&&(power_p[i]>1-beta)&&(PET_p_low[i]<beta)){
          keep[i]=1
        }else {keep[i]=0}

      }
      ######################################
      #     cul_num_s3_mat=length(s3_mat[1,])
      #   if (r_1_e==2^15){r_1_e=NA}
      out_00=cbind(stg2_mat,r_1_e,t1_err,PET_p_low,power_p,PET_p0,EN_p0,keep)

      cul_num=length(out_00[1,])
      keep_all=out_00[,cul_num]
      keep_max=max(keep_all)
    } #end of while loop

    cul_num=length(as.numeric(out_00[1,]))
    out_01=out_00[keep==1,]
    # out_01[-c(1:100),]

    if (length(as.numeric(out_01))>cul_num){
      out_1=out_01[out_01[,(cul_num-1)]==min(out_01[,(cul_num-1)]),-cul_num]
    }else{
      out_1=out_01[-cul_num]
    }

    ##################

    if (length(as.numeric(out_1))>cul_num){
      out=out_1[1,]
    }else{
      out=round(out_1,4)
    }
    out=round(out,4)
    if (out[5]>out[3]){out[5]=NA}
  }
  return(out)
}
#########################################################
#Stg2_mid_sz(n_1,r_01,d_12_cut,alpha,beta,p_0,p_low,p_1)
#Stg2_mid_sz(n_1=11,r_01=2,d_12_cut=3,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.6)
#Stg2_mid_sz(n_1=11,r_01=2,d_12_cut=3,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.55)
#Stg2_mid_sz(n_1=30,r_01=2,d_12_cut=3,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.5)
#Stg2_mid_sz(n_1=30,r_01=2,d_12_cut=2,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.5)
#Stg2_mid_sz(n_1=24,r_01=6,d_12_cut=0,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.55)
#Stg2_mid_sz(n_1=15,r_01=4,d_12_cut=0,alpha=c(0.001,0.024),beta=beta,p_0=0.2,p_low=0.33,p_1=0.55)
#Stg2_mid_sz(n_1=16,r_01=4,d_12_cut=3,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.65)
#Stg2_mid_sz(n_1=100,r_01=25,d_12_cut=5,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.4)
#Stg2_mid_sz(n_1=75,r_01=18,d_12_cut=5,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.4)
#Stg2_mid_sz(n_1=100,r_01=25,d_12_cut=5,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.4)
##################################################

###############################################
simon_sz=function(alpha,beta,p_0,p_1,p_low,d_12_cut){
  ##############
  a=stats::qnorm(1-sum(alpha))
  b=stats::qnorm(1-beta)

  ##################################################################################
  ####### obtaning n_1 and r_1 by using the fucntion ph2simon (package clinfun) ####
  ##################################################################################
  p_ave=(p_0+p_1)/2
  n_approx=ceiling(p_ave*(1-p_ave)*(a+b)^2*(p_1-p_0)^(-2))-1
  n_max=floor(1.5*n_approx)

  ph2simon_out =ph2simon(pu=p_0, pa=p_1, ep1=sum(alpha), ep2=beta, nmax =n_max)
  ph2simon_out_mat=matrix(unlist(ph2simon_out[5]),ncol=6)
  colnames(ph2simon_out_mat)=c("r1", "n1","r_2","n_2", "EN(p0)","PET(p0)")

  ### Minimax ###
  min_max=ph2simon_out_mat[ph2simon_out_mat[,4]==min(ph2simon_out_mat[,4]),]
  if (length((as.numeric(min_max)))==6){
    min_max_opt=min_max
  }else{
    min_max_opt=min_max[min_max[,5]==min(min_max[,5]),]
  }

  ### Optimal ####
  optim=ph2simon_out_mat[ph2simon_out_mat[,5]==min(ph2simon_out_mat[,5]),]
  ##### average   #####################

  r_1_av=ceiling(mean(c(min_max_opt[1],optim[1])))
  n_1_av=ceiling(mean(c(min_max_opt[2],optim[2])))
  av_0=Stg2_sz(n_1=n_1_av,r_1=r_1_av,d_12_cut=d_12_cut,alpha=alpha,beta=beta,p_0=p_0,p_1=p_1)
  av_n_tot=sum(c(av_0[1],av_0[4]))
  #Stg2_sz(n_1=30,r_1=7,d_12_cut=2,alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4)
  av_PET_p0=stats::pbinom(r_1_av,size=n_1_av, prob=p_0)
  av_EN_P0=n_1_av+(1-av_PET_p0)*av_0[4]
  av_1=c(av_0[2],av_0[1],av_0[5],av_n_tot,av_EN_P0,av_PET_p0)

  ###############################
  out_mat_f=rbind(min_max_opt,optim,av_1)
  rownames(out_mat_f)=c("minimax","optimal","average")
  r_1=as.numeric(out_mat_f[,1])
  n_1=as.numeric(out_mat_f[,2])
  r_2=as.numeric(out_mat_f[,3])
  d_12=as.numeric(out_mat_f[,4])-n_1
  ######
  PET_p0=stats::pbinom(r_1,size=n_1, prob=p_0)
  PET_p_low=stats::pbinom(r_1,size=n_1, prob=p_low)
  PET_p=stats::pbinom(r_1,size=n_1,prob=p_1)
  power_p_low=NULL
  power_p=NULL
  t1_err=NULL
  for (i in 1:3){
    power_p_low[i]=rej_prob_1(p_00=p_low,r_1=r_1[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
    power_p[i]=rej_prob_1(p_00=p_1,r_1=r_1[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
    t1_err[i]=rej_prob_1(p_00=p_0,r_1=r_1[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
  }
  En_p0=n_1+(1-stats::pbinom(r_1,size=n_1, prob=p_0))*d_12
  out=cbind(out_mat_f,t1_err,PET_p,PET_p_low,power_p,power_p_low)
  return(out)
}
##############################################
#simon_sz(alpha,beta,p_0,p_1,p_low,d_12_cut)
#simon_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4,p_low=0.33,d_12_cut=0)
#simon_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4,p_low=0.33,d_12_cut=2)
#simon_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4,p_low=0.34,d_12_cut=2)

################################

###############################################
mod_simon_sz<-function(alpha,beta,p_0,p_1,p_low,d_12_cut=0){
  ##############
  a=stats::qnorm(1-sum(alpha))
  b=stats::qnorm(1-beta)

  ##################################################################################
  ####### obtaning n_1 and r_1 by using the fucntion ph2simon (package clinfun) ####
  ##################################################################################
  # p_ave=(p_0+p_1)/2
  #  n_approx=ceiling(p_ave*(1-p_ave)*(a+b)^2*(p_1-p_0)^(-2))-1
  p_ave_p_low=(p_0+p_low)/2
  n_approx_p_low=ceiling(p_ave_p_low*(1-p_ave_p_low)*(a+b)^2*(p_low-p_0)^(-2))-1
  n_max_p_low=floor(1.5*n_approx_p_low)

  ph2simon_out =ph2simon(pu=p_0, pa=p_1, ep1=sum(alpha), ep2=beta, nmax =n_max_p_low)
  ph2simon_out_mat=matrix(unlist(ph2simon_out[5]),ncol=6)
  #colnames(ph2simon_out_mat)=c("r1", "n1","r_2","n_2", "EN(p0)","PET(p0)")

  r_1=ph2simon_out_mat[,1]
  n_1=ph2simon_out_mat[,2]
  PET_plow=round(stats::pbinom(r_1,size=n_1, prob=p_low),4)
  keep=NULL

  for (i in 1:length(r_1)){
    # PET_plow[i]=round(stats::pbinom(r_1[i],size=n_1[i], prob=p_low),4)
    if(PET_plow[i]<beta ){keep[i]=1}else(keep[i]=0)
  }


  ph2simon_out_plow_mat_0=cbind(ph2simon_out_mat,PET_plow,keep)
  colnames(ph2simon_out_plow_mat_0)=c("r1", "n1","r_2","n_2", "EN(p0)","PET(p0)","PET_plow","keep")
  ph2simon_out_plow_mat=ph2simon_out_plow_mat_0[keep=1,-8]

  ### Minimax ###
  min_max=ph2simon_out_plow_mat[ph2simon_out_plow_mat[,4]==min(ph2simon_out_plow_mat[,4]),]
  if (length((as.numeric(min_max)))==7){
    min_max_opt=min_max
  }else{
    min_max_opt=min_max[min_max[,5]==min(min_max[,5]),]
  }

  ### Optimal ####
  optim=ph2simon_out_plow_mat[ph2simon_out_plow_mat[,5]==min(ph2simon_out_plow_mat[,5]),]
  ##### average   #####################

  r_1_av=ceiling(mean(c(min_max_opt[1],optim[1])))
  n_1_av=ceiling(mean(c(min_max_opt[2],optim[2])))
  av_0=Stg2_sz(n_1=n_1_av,r_1=r_1_av,d_12_cut=d_12_cut,alpha=alpha,beta=beta,p_0=p_0,p_1=p_1)
  #av_n_tot=sum(c(av_0[1],av_0[4]))
  #Stg2_sz(n_1=30,r_1=7,d_12_cut=2,alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4)
  av_PET_p0=stats::pbinom(r_1_av,size=n_1_av, prob=p_0)
  av_EN_P0=n_1_av+(1-av_PET_p0)*av_0[4]
  av_PET_p_low=stats::pbinom(r_1_av,size=n_1_av, prob=p_low)
  av_1=c(av_0[2],av_0[1],av_0[5],av_0[6],av_EN_P0,av_PET_p0,av_PET_p_low)

  ###############################
  out_mat_f=rbind(min_max_opt,optim,av_1)
  rownames(out_mat_f)=c("minimax","optimal","average")
  r_1=as.numeric(out_mat_f[,1])
  n_1=as.numeric(out_mat_f[,2])
  r_2=as.numeric(out_mat_f[,3])
  d_12=as.numeric(out_mat_f[,4])-n_1
  ######
  # PET_p0=stats::pbinom(r_1,size=n_1, prob=p_0)
  # PET_p_low=stats::pbinom(r_1,size=n_1, prob=p_low)
  PET_p_1=stats::pbinom(r_1,size=n_1,prob=p_1)
  power_p_low=NULL
  power_p_1=NULL
  t1_err=NULL
  for (i in 1:3){
    # power_p_low[i]=rej_prob_1(p_00=p_low,r_1=r_1[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
    power_p[i]=rej_prob_1(p_00=p_1,r_1=r_1[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
    t1_err[i]=rej_prob_1(p_00=p_0,r_1=r_1[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
  }
  # En_p0=n_1+(1-stats::pbinom(r_1,size=n_1, prob=p_0))*d_12
  # out=cbind(out_mat_f,t1_err,PET_p,PET_p_low,power_p,power_p_low)
  out=cbind(out_mat_f,t1_err,PET_p_1,power_p_1,power_p_low)
  return(out)
}
##############################################
#mod_simon_sz(alpha,beta,p_0,p_1,p_low,d_12_cut=0)
#mod_simon_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4,p_low=0.33,d_12_cut=0)
#mod_simon_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4,p_low=0.33,d_12_cut=2)
#mod_simon_sz(alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4,p_low=0.34,d_12_cut=2)

#############################################
################################################################################
#### part 3: Sample size calculation for three stage expanded Simon's Design ###
################################################################################

########################################################
## Function fro various steps for three stage  design  ###
########################################################

##############################################################
####### rej_prob_2=1-h(r_1,n_1,r_2,min(r_2_e,n_2),r_3,n_3,p_00) see manuscript ####
####### rej_prob_2=1-k(r_1,n_1,r_2,min(r_2_e,n_2),r_3,n_3,p_00) see manuscript ####
##########################################################################
###################
rej_prob_2<-function(p_00,r_1,n_1,r_2,r_2_e=2^15,d_12,r_3,d_23){
  a_0=NULL
  c_0=NULL
  n_2=n_1+d_12
  for (i in ((r_1+1):min(n_1,r_2))){
    a_0[i]=stats::dbinom(i, size=n_1, prob=p_00)*stats::pbinom(r_2-i, size=d_12, prob=p_00)
  }
  a=a_0[(r_1+1):length(a_0)]
  for (i in ((r_1+1):min(n_1,r_2))){
    b_0=NULL
    for (j in ((r_2+1-i):(min(c(n_2,r_2_e,r_3))-i))){
      b_0[j]=stats::dbinom(i, size=n_1, prob=p_00)*stats::dbinom(j, size=d_12, prob=p_00)*stats::pbinom(r_3-i-j, size=d_23, prob=p_00)
    }

    b=b_0[(r_2+1-i):length(b_0)]
    c_0[i]=sum(b)
  }
  c=c_0[(r_1+1):length(c_0)]

  d=stats::pbinom(r_1, size=n_1, prob=p_00)+sum(a)+sum(c)
  e=1-d
  return(e)
}
#############################################
#rej_prob_2(p_00=0.2,r_1=4,n_1=22,r_2=10,r_2_e=2^15,d_12=28,r_3=36,d_23=37)
#round(rej_prob_2(p_00=0.4,r_1=0,n_1=22,r_2=11,r_2_e=2^15,d_12=28,r_3=17,d_23=6),5)

## optimal design, p_0=0.2, p_1=0.4
#round(rej_prob_2(p_00=0.4,r_1=4,n_1=22,r_2=9,r_2_e=2^15,d_12=28,r_3=17,d_23=6),5)
#round(rej_prob_2(p_00=0.2,r_1=4,n_1=22,r_2=9,r_2_e=2^15,d_12=28,r_3=17,d_23=6),5)
#EN_3stg(r_1=4,n_1=22,r_2=9,n_2=50,n_3=56,p_00=0.2)
## minimax design, p_0=0.2, p_1=0.4
#round(rej_prob_2(p_00=0.4,r_1=10,n_1=37,r_2=25,r_2_e=2^15,d_12=63,r_3=29,d_23=4),5)
#round(rej_prob_2(p_00=0.2,r_1=10,n_1=37,r_2=25,r_2_e=2^15,d_12=63,r_3=29,d_23=4),5)
#EN_3stg(r_1=10,n_1=37,r_2=25,n_2=50,n_3=56,p_00=0.2)
###########################

##################### Calculating the expected sample size under p_00 in a three stage design #########
EN_3stg=function(r_1,n_1,r_2,r_2_e=2^15,n_2,n_3,p_00){
  d_12=n_2-n_1
  PET_n_1=stats::pbinom(r_1, size=n_1, prob=p_00)
  if (r_2_e>n_2){
    PET_n_2=1-rej_prob_1(p_00=p_00,r_1=r_1,r_2=r_2,n_1=n_1,d_12=d_12)
  }else{
    t_1_f=stats::pbinom(r_1, size=n_1, prob=p_00)
    t_1_e=1-stats::pbinom(r_2_e-1, size=n_1, prob=p_00)
    t_2_f_0=NULL
    t_2_e_0=NULL
    for (i in ((r_1+1):min(c(r_2,n_1)))){
      t_2_f_0[i]=stats::dbinom(i, size=n_1, prob=p_00)*(1-stats::pbinom(r_2-i, size=d_12, prob=p_00))
    }

    t_2_f_1=t_2_f_0[(r_1+1):length(t_2_f_0)]
    t_2_f=sum(t_2_f_1)
    for (i in ((r_1+1):min(c(r_2_e,n_1)))){
      t_2_e_0[i]=stats::dbinom(i, size=n_1, prob=p_00)*(1-stats::pbinom(r_2_e-i-1, size=d_12, prob=p_00))
    }

    t_2_e_1=t_2_e_0[(r_1+1):length(t_2_e_0)]
    t_2_e=sum(t_2_e_1)
    PET_n_2=t_1_f+t_2_f+t_1_e+t_2_e
  }
  EN=n_1+(1-PET_n_1)*d_12+(1-PET_n_2)*(n_3-n_2)
  EN=round(EN,2)
  return(EN)
}
#############################################
#EN_3stg(r_1=2,n_1=10,r_2=6,r_2_e=2^15,n_2=20,n_3=30,p_00=0.2)
#EN_3stg(r_1=5,n_1=22,r_2=11,r_2_e=19,n_2=50,n_3=61,p_00=0.2)
#####################################################

##########################
### creating matrix for searching r_1, r_2, r_3, n_3, after n_1 and n_2 have been previoly selected #######
stg3_mat=function(r_01,n_1,r_02,n_2,r_03,d_23){
  n_3=n_2+d_23
  r_001=floor(r_01*0.4)
  r_002=floor(r_02*0.4)
  r_003=floor(r_03*0.4)
  r_011=ceiling(r_01*1.1)
  r_021=ceiling(r_02*1.1)
  mat_1=NULL
  for (i in (r_001:r_011)){
    mat_2=NULL
    r_1=i
    for (j in (r_1:r_021)){
      r_2=j
      r_3=c(r_003:r_03)
      mat_02=cbind(r_2,r_3,n_1,n_2,n_3)
      mat_2=rbind(mat_2,mat_02)
    }

    mat_01=cbind(r_1,mat_2)
    mat_1=rbind(mat_1,mat_01)
  }
  return(mat_1)
}
###########
#stg3_mat(r_01,n_1,r_02,n_2,r_03,d_23)
#stg3_mat(r_01=5,n_1=22,r_02=11,n_2=40,r_03=29,d_23=3)[-c(1:100),]
#length(stg3_mat(r_01=5,n_1=22,r_02=11,n_2=40,r_03=29,d_23=3)[,1] )
##############################################

#########################################
##### search for r_2_e after r_1,n_1,r_2,n_2 have been previousle chosen #####
r_2_e_fun=function(alpha,r_1,n_1,r_2,n_2,p_00){
  ## choose crticical boundary at n_2
  d_12=n_2-n_1
  r_2_e_0=NULL
  for (i in r_2:n_2){
    r_2_e_0[i]=rej_prob_1(p_00=p_00,r_1=r_1,r_1_e=2^8,r_2=i,n_1=n_1,d_12=d_12)
  }
  r_2_e_1=cbind(r_2:n_2,r_2_e_0[-c(1:(r_2-1))])
  r_2_e_2=r_2_e_1[(r_2_e_1[,2]<alpha[1]),]
  r_2_e=r_2_e_2[1,1]
  return(r_2_e)
}
############################
#n_beta(alpha=0.025,beta=0.1,p_1=0.4,p_0=0.2)
#r_2_e_fun(alpha=c(0.001,0.024),r_1=3,n_1=22,r_2=11,n_2=100,p_00=0.4)
#########################################
#rm(list = ls())

###### searching for r_1,r_2,r_3, n_3 after n_1 and n_2 have been prevoulsy selected #########
Stg3_sz_no_r2e=function(r_01,n_1,r_02,n_2,d_23_cut,alpha,beta,p_0,p_1,p_low){
  ##### This function searchs for r_1,r_2,r_3,n_3, without r_2_e ##########
  r_2_e=NA
  a_1=stats::qnorm(1-alpha)
  a=stats::qnorm(1-sum(alpha))
  b=stats::qnorm(1-beta)
  d_12=n_2-n_1
  n_1_beta=(1-p_1)*p_1*(a+b)^2*(p_1-p_0)^(-2)
  n_max=ceiling(1.5*n_1_beta)
  ######## into loop to identify r_1, r_2, d_23, and r_3 #####
  keep_max=0
  n_3=n_2+max(d_23_cut-1,0)
  # out_0=NULL
  while (keep_max==0){
    n_3=n_3+1
    d_23=n_3-n_2
    r_03=floor(n_3*p_1)
    s3_mat=stg3_mat(r_01,n_1,r_02,n_2,r_03,d_23)
    mat_row_count=length(s3_mat[,1])
    t1_err=NULL
    PET_p1=NULL
    power=NULL
    keep=NULL
    EN_p0=NULL
    PET_plow=NULL
    for (i in (1:mat_row_count)){
      t1_err[i]=round(rej_prob_2(p_00=p_0,r_1=s3_mat[i,1],n_1=n_1,r_2=s3_mat[i,2],r_2_e=2^15,d_12=d_12,r_3=s3_mat[i,3],d_23=d_23),4)
      PET_p1[i]=1-round(rej_prob_1(p_00=p_1,r_1=s3_mat[i,1],r_2=s3_mat[i,2],n_1=n_1,d_12=d_12),4)
      PET_plow[i]=1-round(rej_prob_1(p_00=p_low,r_1=s3_mat[i,1],r_2=s3_mat[i,2],n_1=n_1,d_12=d_12),4)
      power[i]=round(rej_prob_2(p_00=p_1,r_1=s3_mat[i,1],n_1=n_1,r_2=s3_mat[i,2],r_2_e=2^15,d_12=d_12,r_3=s3_mat[i,3],d_23=d_23),4)
      EN_p0[i]=EN_3stg(r_1=s3_mat[i,1],n_1=n_1,r_2=s3_mat[i,2],n_2=n_2,n_3=s3_mat[i,6],p_00=p_0)
      if ((t1_err[i]<=sum(alpha))&&(power[i]>1-beta)&&(PET_plow[i]<beta)){
        keep[i]=1
      }else {keep[i]=0}
    }

    out_00=cbind(s3_mat,r_2_e,t1_err,PET_p1,PET_plow,power,EN_p0,keep)

    cul_num=length(out_00[1,])
    keep_all=out_00[,cul_num]
    keep_max=max(keep_all)

  } #end of while loop

  #### keep rows with power >=1-beta, and type I errr<alpha
  cul_num=length(as.numeric(out_00[1,]))
  out_01=out_00[out_00[,cul_num]==1,]

  if (length(as.numeric(out_01))>cul_num){
    out_1=out_01[out_01[,(cul_num-1)]==min(out_01[,(cul_num-1)]),]
    if (length(as.numeric( out_1))>cul_num){
      out=out_1[1,-cul_num]
    }else{
      out=out_1[-cul_num]
    }
  }else{
    out=out_01[-cul_num]
  }
  out=round(out,4)
  return(out)
}
#########################################################
#Stg3_sz_no_r2e(r_01,n_1,r_02,n_2,d_23_cut,alpha,beta,p_0,p_1,p_low)
#Stg3_sz_no_r2e(r_01=5,n_1=22,r_02=19,n_2=50,d_23_cut=2,alpha=0.025,beta=0.1,p_0=0.2,p_1=0.4,p_low=0.33)
#Stg3_sz_no_r2e(r_01=5,n_1=22,r_02=11,n_2=50,d_23_cut=2,alpha=0.025,beta=0.1,p_0=0.2,p_1=0.5,p_low=0.33)

#########################################################################################


###  Search for r_2_e, r_3,n_3 after r_1,n_1,r_2,n_2 have been previous chosen ####
Stg3_sz_r2e=function(r_1,n_1,r_2,n_2,r_03,r_2_e,PET_plow,d_23_cut,alpha,beta,p_0,p_1){
  a_1=stats::qnorm(1-alpha)
  a_tot=stats::qnorm(1-sum(alpha))
  b=stats::qnorm(1-beta)
  r_003=floor(r_03*0.6)
  d_12=n_2-n_1

  ######## into loop to identify r_1, r_2, d_23, and r_3 #####
  keep_max=0
  d_23=max(d_23_cut-1,0)
  # out_0=NULL
  while (keep_max==0){
    d_23=d_23+1
    n_3=n_2+d_23

    r_003u=ceiling(1.2*p_1*n_3)
    t1_err_0=NULL
    PET_p1_0=NULL
    power_0=NULL
    keep_0=NULL
    EN_p0_0=NULL
    for (i in (r_003:r_003u)){
      t1_err_0[i]=rej_prob_2(p_00=p_0,r_1=r_1,n_1=n_1,r_2=r_2,r_2_e=r_2_e,d_12=d_12,r_3=i,d_23=d_23)
      PET_p1_0[i]=1-rej_prob_1(p_00=p_1,r_1=r_1,r_1_e=2^15,r_2=r_2,n_1=n_1,d_12=d_12)
      power_0[i]=rej_prob_2(p_00=p_1,r_1=r_1,n_1=n_1,r_2=r_2,r_2_e=r_2_e,d_12=d_12,r_3=i,d_23=d_23)
      EN_p0_0[i]=EN_3stg(r_1=r_1,n_1=n_1,r_2=r_2,r_2_e=r_2_e,n_2=n_2,n_3=n_3,p_00=0.2)

      if ( (t1_err_0[i]<=sum(alpha))&&(power_0[i]>1-beta)){
        keep_0[i]=1
      }else {keep_0[i]=0}
    }

    t1_err=t1_err_0[r_003:r_003u]
    PET_p1=PET_p1_0[r_003:r_003u]
    power=power_0[r_003:r_003u]
    keep=keep_0[r_003:r_003u]
    EN_p0=EN_p0_0[r_003:r_003u]
    r_3=c(r_003:r_003u)
    out_00=cbind(r_1,r_2,r_3,n_1,n_2,n_3,r_2_e,t1_err,PET_p1,PET_plow,power,EN_p0,keep)

    # out_0=rbind(out_0,out_00)
    cul_num=length(out_00[1,])
    keep_all=out_00[,cul_num]
    keep_max=max(keep_all)
  } #end of while loop

  # out_000=out_0[out_0[,1]==4,]
  # out_001=out_000[(out_000[,2]==9),]
  # out_002=out_001[(out_001[,3]==17),]

  #### keep rows with power >=1-beta, and type I errr<alpha
  cul_num=length(as.numeric(out_00[1,]))
  out_01=out_00[out_00[,cul_num]==1,]

  if (length(as.numeric(out_01))>cul_num){
    out_1=out_01[out_01[,(cul_num-1)]==min(out_01[,(cul_num-1)]),]
    if (length(as.numeric( out_1))>cul_num){
      out=out_1[1,-cul_num]
    }else{
      out=out_1[-cul_num]
    }
  }else{
    out=out_01[-cul_num]
  }
  out=round(out,4)
  return(out)
}
#########################################################
#Stg3_sz_r2e=function(r_1,n_1,r_2,n_2,r_03,r_2_e,PET_plow,d_23_cut,alpha,beta,p_0,p_1)
#Stg3_sz_r2e(r_1=5,n_1=22,r_2=11,n_2=50,r_03=17,r_2_e=19,PET_plow=0.09,d_23_cut=2,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.45)
#Stg3_sz_r2e(r_1=5,n_1=22,r_2=11,n_2=50,r_03=17,r_2_e=19,PET_plow=0.09,d_23_cut=2,alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_1=0.48)
###############################################################

############################################################################
############################################################################
#### part 4: Design function for both two stage  and three stage designs ###
############################################################################
############################################################################

#############################

#' Single arm phase 2 oncology trial: Two-stage design, the hybrid design
#'
#' @param alpha Type I error for the one-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param p_0 Response rate indicating low activity with insufficient clinical benefit
#' @param p_low Minimally clinically beneficial response rate
#' @param p_1 Assumed response rate
#' @param d_12_cut Sets the minimal difference between n_2 and n_1 for the design; min(n_2 - n_1)
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Simon's design: p=p_1}{A matrix with one row per design (minimax, optimal and
#'     average) for Simon's design using p_1 as the target response rate, with columns r1,
#'     n1, r_2, n_2, EN(p0) (expected sample size under p_0), PET(p0) (probability of early
#'     termination under p_0), t1_err, PET_p (PET under p_1), PET_p_low (PET under p_low),
#'     power_p (power under p_1) and power_p_low (power under p_low)}
#'   \item{Simon's design: p=p_low}{The corresponding matrix for Simon's design using p_low
#'     as the target response rate}
#'   \item{Hybrid design}{A matrix with one row per design (minimax, optimal and average)
#'     and columns n_1, r_1, r_1_e (the trial is stopped for superiority if at least r_1_e
#'     responses are observed at the first look), r_2, n_2, t1_err, power (power under p_1),
#'     PET_p1, PET_p_low, PET_p0, En_p0 (expected sample size under p_0) and PW_p_low (power
#'     under p_low)}
#' }
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials:
#' an expansion of Simon's design, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#'
#' @examples
#' hybrid_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,d_12_cut=0)
#'
#' @export
hybrid_sz<-function(alpha,beta,p_0,p_low,p_1,d_12_cut){
  ##############
  alpha_tot=sum(alpha)
  a=stats::qnorm(1-sum(alpha))
  b=stats::qnorm(1-beta)

  ##################################################################################
  ####### obtaning n_1 and r_1 by using the fucntion ph2simon (package clinfun) ####
  ##################################################################################
  simon_p1=simon_sz(alpha=alpha_tot,beta=beta,p_0=p_0,p_1=p_1,p_low=p_low,d_12_cut=d_12_cut)
  simon_p_low=simon_sz(alpha=alpha_tot,beta=beta,p_0=p_0,p_1=p_low,p_low=p_low,d_12_cut)

  ######
  n_1_p1=as.numeric(simon_p1[,2])
  r_1_p1=as.numeric(simon_p1[,1])
  n_1_p_low=as.numeric(simon_p_low[,2])
  r_1_p_low=as.numeric(simon_p_low[,1])

  hybrid_design=NULL
  for (i in 1:3){
   hybrid_design=rbind(hybrid_design, Stg2_sz(n_1=n_1_p_low[i],r_1=r_1_p_low[i],d_12_cut=d_12_cut,alpha=alpha,beta=beta,p_0=p_0,p_1=p_1))
  }
  ##################################################################################

  n_1=hybrid_design[,1]
  r_1=hybrid_design[,2]
  r_1_e=hybrid_design[,3]
  d_12=hybrid_design[,4]
  r_2=hybrid_design[,5]
  PW_p_low=NULL
  PET_p_low=NULL
  PET_p0=NULL
  for (i in 1:3){
    if (is.na(r_1_e[i])){
      PW_p_low[i]=rej_prob_1(p_00=p_low,r_1=r_1[i],r_1_e=2^15,r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])

    }else{
      PW_p_low[i]=rej_prob_1(p_00=p_low,r_1=r_1[i],r_1_e=r_1_e[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
    }
    PET_p_low[i]=stats::pbinom(r_1[i],size=n_1[i], prob=p_low)
    PET_p0[i]=stats::pbinom(r_1[i],size=n_1[i], prob=p_0)
  }

  En_p0=n_1+(1-stats::pbinom(r_1,size=n_1, prob=p_0))*d_12
  hybrid_design=cbind(hybrid_design[,-4],PET_p_low,PET_p0,En_p0,PW_p_low)
  rownames(hybrid_design)=c("Minimax","Optimal","Average")

  out=list(simon_p1,simon_p_low,hybrid_design)
  names(out)=c("Simon's design: p=p_1","Simon's design: p=p_low","Hybrid design")
  return(out)
  #return(exp_simon_design)
}
#######################################
#hybrid_sz(alpha,beta,p_0,p_low,p_1,d_12_cut)
#hybrid_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,d_12_cut=0)
#hybrid_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,d_12_cut=0)
#hybrid_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,d_12_cut=2)
#hybrid_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.65,d_12_cut=3)
#hybrid_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.4,d_12_cut=5)
################################

#############################

#' Single arm phase 2 oncology trial: Two-stage design, the midpoint design
#'
#' @param alpha Type I error for the one-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param p_0 Response rate indicating low activity with insufficient clinical benefit
#' @param p_low Minimally clinically beneficial response rate
#' @param p_1 Assumed response rate
#' @param q Weight between 0 and 1 used to interpolate the first-look sample size of the
#'   mid-point design from the two Simon's designs: n_1 = q * n_1(p_1) + (1 - q) * n_1(p_low),
#'   where n_1(p_1) and n_1(p_low) are the first-look sample sizes of Simon's designs with
#'   target response rates p_1 and p_low; 0 < q < 1
#' @param d_12_cut Sets the minimal difference between n_2 and n_1 for the design; min(n_2 - n_1)
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Simon's design: p=p_1}{A matrix with one row per design (minimax, optimal and
#'     average) for Simon's design using p_1 as the target response rate, with columns r1,
#'     n1, r_2, n_2, EN(p0) (expected sample size under p_0), PET(p0) (probability of early
#'     termination under p_0), t1_err, PET_p (PET under p_1), PET_p_low (PET under p_low),
#'     power_p (power under p_1) and power_p_low (power under p_low)}
#'   \item{Simon's design: p=p_low}{The corresponding matrix for Simon's design using p_low
#'     as the target response rate}
#'   \item{The mid-point design: p=p_1}{A matrix with one row per design (minimax, optimal
#'     and average) and columns r_1, r_2, n_1, n_2, r_1_e (the trial is stopped for
#'     superiority if at least r_1_e responses are observed at the first look), t1_err,
#'     PET_p_low, power_p (power under p_1), PET_p0, EN_p0 (expected sample size under p_0),
#'     PW_p_low (power under p_low) and q. If q is too large for the probability of early
#'     termination under p_low to be below beta, a character string with a warning message is
#'     returned instead}
#' }
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials:
#' an expansion of Simon's design, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#'
#' @examples
#' midpnt_sz(alpha = c(0.001,0.024), beta = 0.1, p_0 = 0.2, p_low = 0.33, p_1 = 0.5, q = 0,
#' d_12_cut = 5)
#'
#' @export
midpnt_sz<-function(alpha,beta,p_0,p_low,p_1,q,d_12_cut){
  ##############
  alpha_tot=sum(alpha)
  a=stats::qnorm(1-sum(alpha))
  b=stats::qnorm(1-beta)

  ##################################################################################
  ####### obtaning n_1 and r_1 by using the fucntion ph2simon (package clinfun) ####
  ##################################################################################
  simon_p1=simon_sz(alpha=alpha_tot,beta=beta,p_0=p_0,p_1=p_1,p_low=p_low,d_12_cut=d_12_cut)
  simon_p_low=simon_sz(alpha=alpha_tot,beta=beta,p_0=p_0,p_1=p_low,p_low=p_low,d_12_cut=d_12_cut)
  ######
  n_1_p1=as.numeric(simon_p1[,2])
  r_1_p1=as.numeric(simon_p1[,1])
  n_1_p_low=as.numeric(simon_p_low[,2])
  r_1_p_low=as.numeric(simon_p_low[,1])

  n_1_mid=ceiling(q*n_1_p1+(1-q)*n_1_p_low)
  r_1_mid=ceiling(q*r_1_p1+(1-q)*r_1_p_low)
  PET_p_low_0=round(stats::pbinom(1,size=n_1_mid, prob=p_low),4)

  if (max(PET_p_low_0)>beta){
    return("q is too large and n_1 is too small. Impossible for PET_p_low<beta. Reduce q ")
  }

  #Stg2_mid_sz(n_1,r_01,d_12_cut,alpha,beta,p_0,p_1)
  midpnt_design=NULL
  for (i in 1:3){
    midpnt_design=rbind(midpnt_design,round(Stg2_mid_sz(n_1=n_1_mid[i],r_01=r_1_mid[i],d_12_cut=d_12_cut,alpha=alpha,beta=beta,p_0=p_0,p_low=p_low,p_1=p_1),4))
  }
  ##################################################################################

  n_1=midpnt_design[,3]
  r_1=midpnt_design[,1]
  r_1_e=midpnt_design[,5]
  r_2=midpnt_design[,2]
  n_2=midpnt_design[,4]
  d_12=n_2-n_1
  PW_p_low=NULL
  for (i in 1:3){
    if (is.na(r_1_e[i])){
      PW_p_low[i]=rej_prob_1(p_00=p_low,r_1=r_1[i],r_1_e=2^15,r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
    }else{
      PW_p_low[i]=rej_prob_1(p_00=p_low,r_1=r_1[i],r_1_e=r_1_e[i],r_2=r_2[i],n_1=n_1[i],d_12=d_12[i])
    }
  }

  PET_p_low=stats::pbinom(r_1,size=n_1,prob=p_low)
  En_p0=n_1+(1-stats::pbinom(r_1,size=n_1,prob=p_0))*d_12
  midpnt_design=cbind(midpnt_design,PW_p_low,q)
  rownames(midpnt_design)=c("Minimax","optimal","average")

  out=list(simon_p1,simon_p_low,midpnt_design)
  names(out)=c("Simon's design: p=p_1","Simon's design: p=p_low","The mid-point design: p=p_1")
  return(out)
  #return(exp_simon_design)
}
#######################################
#midpnt_sz(alpha,beta,p_0,p_low,p_1,q,d_12_cut)
#midpnt_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.4,q=0.5,d_12_cut=5)
#midpnt_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.4,q=0,d_12_cut=5)
#midpnt_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,q=0.4,d_12_cut=0)
#midpnt_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,q=0.6,d_12_cut=2)
#midpnt_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.65,q=0.9,d_12_cut=3)

################################


###############################

#' Single arm phase 2 oncology trial: Three-stage design
#'
#' @param alpha Type I error for the one-sided test
#' @param beta Type II error rate, i.e. 1 - power
#' @param p_0 Response rate indicating low activity with insufficient clinical benefit
#' @param p_low Minimally clinically beneficial response rate
#' @param p_1 Assumed response rate
#' @param d_23_cut Sets the minimal difference between n_3 and n_2 for the design; min(n_3 - n_2)
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{Simon's design: p=p_1}{A matrix with one row per design (minimax, optimal and
#'     average) for Simon's design using p_1 as the target response rate, with columns r1,
#'     n1, r_2, n_2, EN(p0) (expected sample size under p_0), PET(p0) (probability of early
#'     termination under p_0), t1_err, PET_p (PET under p_1), PET_p_low (PET under p_low),
#'     power_p (power under p_1) and power_p_low (power under p_low)}
#'   \item{Simon's design: p=p_low}{The corresponding matrix for Simon's design using p_low
#'     as the target response rate}
#'   \item{Expanded 3 stage design}{A matrix with one row per design (minimax, optimal and
#'     average) and columns r_1, r_2, r_3, n_1, n_2, n_3 (n_3 is also the total sample size),
#'     r_2_e (the trial is stopped for superiority if at least r_2_e responses are observed at
#'     the second look), t1_err, PET_p1, PET_plow, power (power under p_1), EN_p0 (expected
#'     sample size under p_0) and PW_p_low (power under p_low)}
#' }
#'
#' @references R. Simon. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 10:1-10 (1989).
#' @references P. Gao & W. Zhang (15 Apr 2024): Adaptive sequential design for phase II single-arm oncology trials:
#' an expansion of Simon's design, Journal of Biopharmaceutical Statistics, DOI: 10.1080/10543406.2024.2341673.
#' To link to this article: \url{https://www.tandfonline.com/doi/full/10.1080/10543406.2024.2341673}.
#'
#' @examples
#' EXP_3_stg_sz( alpha = c(0.001,0.024), beta = 0.1, p_0 = 0.2, p_low = 0.4, p_1 = 0.55,
#'   d_23_cut = 0)
#'
#' @export
EXP_3_stg_sz=function(alpha,beta,p_0,p_low,p_1,d_23_cut){
  ##############
  alpha_tot=sum(alpha)

  a_tot=stats::qnorm(1-alpha_tot)
  b=stats::qnorm(1-beta)
  n_1_beta=(1-p_1)*p_1*(a_tot+b)^2*(p_1-p_0)^(-2)
  ##################################################################################
  ####### obtaning simon_p1 and simon_p_low by using the fucntion ph2simon (package clinfun) ####
  ##################################################################################

  simon_p1=simon_sz(alpha=alpha_tot,beta=beta,p_0=p_0,p_1=p_1,p_low=p_low,d_12_cut=0)
  simon_p_low=simon_sz(alpha=alpha_tot,beta=beta,p_0=p_0,p_1=p_low,p_low=p_low,d_12_cut=0)

  ####### average of minimax and optimal ####
  n_1_p1=as.integer(simon_p1[,2])
  n_1_p_low=as.integer(simon_p_low[,2])
  r_1_p1=as.integer(simon_p1[,1])
  r_1_p_low=as.integer(simon_p_low[,1])

  #####
  n_1=n_1_p1
  n_2=n_1_p_low
  r_01=r_1_p1
  r_02=r_1_p_low
  ######

  Exp_3_stg_design_0=NULL
  for (i in 1:3){
    Exp3_0=Stg3_sz_no_r2e(r_01=r_01[i],n_1=n_1[i],r_02=r_02[i],n_2=n_2[i],d_23_cut=d_23_cut,alpha=alpha,beta=beta,p_0=p_0,p_1=p_1,p_low=p_low)
    if (n_2[i]>n_1_beta){
      r_2_e=r_2_e_fun(alpha=alpha,r_1=Exp3_0[1],n_1=Exp3_0[4],r_2=Exp3_0[2],n_2=Exp3_0[5],p_00=p_0)
      Exp3_1=Stg3_sz_r2e(r_1=Exp3_0[1],n_1=Exp3_0[4],r_2=Exp3_0[2],n_2=Exp3_0[5],r_03=Exp3_0[3],r_2_e=r_2_e,PET_plow=Exp3_0[10],d_23_cut=d_23_cut,alpha=alpha,beta=beta,p_0=p_0,p_1=p_1)
    }else{
      Exp3_1= Exp3_0
    }
    Exp_3_stg_design_0=rbind(Exp_3_stg_design_0,Exp3_1)
  }
  r_1=Exp_3_stg_design_0[,1]
  r_2=Exp_3_stg_design_0[,2]
  r_3=Exp_3_stg_design_0[,3]
  n_1=Exp_3_stg_design_0[,4]
  n_2=Exp_3_stg_design_0[,5]
  n_3=Exp_3_stg_design_0[,6]
  r_2_e=Exp_3_stg_design_0[,7]
  d_12=n_2-n_1
  d_23=n_3-n_2
  PW_p_low=NULL
  for (i in 1:3){
    if (is.na(r_2_e[i])){
      PW_p_low[i]=rej_prob_2(p_00=p_low,r_1=r_1[i],n_1=n_1[i],r_2=r_2[i],r_2_e=2^15,d_12=d_12[i],r_3=r_3[i],d_23=d_23[i])
    }else{
      PW_p_low[i]=rej_prob_2(p_00=p_low,r_1=r_1[i],n_1=n_1[i],r_2=r_2[i],r_2_e=r_2_e[i],d_12=d_12[i],r_3=r_3[i],d_23=d_23[i])
    }
  }

  Exp_3_stg_design=cbind(Exp_3_stg_design_0,PW_p_low)
  rownames(Exp_3_stg_design)=c("Minimax","Optimal","Average")
  #simon_sz(alpha,beta,p_0,p_1,p_low,d_12_cut)

  out=list(simon_p1,simon_p_low,Exp_3_stg_design)
  names(out)=c("Simon's design: p=p_1","Simon's design: p=p_low","Expanded 3 stage design")
  return(out)
  #return(exp_simon_design)
}
#######################################
#warnings()
#EXP_3_stg_sz(alpha,beta,p_0,p_low,p_1,d_23_cut)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,d_23_cut=0)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,d_23_cut=0)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.55,d_23_cut=2)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.65,d_23_cut=3)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.45,p_1=0.65,d_23_cut=3)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.32,p_1=0.4,d_23_cut=3)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.45,d_23_cut=3)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.34,p_1=0.4,d_23_cut=3)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.4,p_1=0.5,d_23_cut=3)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.45,p_1=0.55,d_23_cut=3)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.34,p_1=0.4,d_23_cut=3)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.4,d_23_cut=4)
#EXP_3_stg_sz(alpha=c(0.001,0.024),beta=0.1,p_0=0.2,p_low=0.33,p_1=0.4,d_23_cut=5)
#EXP_3_stg_sz(alpha=0.025,beta=0.1,p_1=0.5,p_0=0.2)


