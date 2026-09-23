#library(foreach)
#library(iterators)
#library(parallel)
#library(doParallel)

#' @import foreach
#' @import parallel
#' @import doParallel

#debug_simu = TRUE
debug_simu = FALSE

free_version = TRUE
# free_version = FALSE
 # use is rd cores 72
#default_use_cores = 14
#default_use_cores = 18
# 2024/3/30 use 22 cores
if (Sys.info()[4]!="simucalc") {
  if (debug_simu) {
    default_use_cores = 60
  } else {
    default_use_cores = 50
  }
} else {  # simucalc
# 2025/1/16 added simulation computing server
  if (debug_simu) {
    default_use_cores = 140
  } else {
    default_use_cores = 50
  }
}
#default_use_cores = 10



#paraller_num = 0
is_checking = !identical(Sys.getenv("NOT_CRAN"), "true")
if (is_checking) {
  paraller_num = 2
} else {
  # devtools::check(), 使用paraller_num = 2, 避免错误
  #paraller_num = 2
  # devtools::check(env_vars = c(NOT_CRAN = "false")), 这样检查
  paraller_num = detectCores() - 4
}
paraller_register = FALSE

setParallerNum <- function(core_num)
{
  cores = detectCores()
  if (core_num>=cores) {
    core_num = cores-1
  }
  if(core_num<1) {
    core_num = 1
  }
  paraller_num <<- core_num
}

getParallerNum <- function(sim_num)
{
  if (isWindow()) {
    cores = detectCores()
    if (sim_num<2) {
        cores = 2
    } else if (sim_num<cores) {
        cores = sim_num
    }
    return(cores)
  }
  return(paraller_num)
}

registerISParallel <- function()
{
  if (paraller_register) {
    return()
  }
  if (paraller_num == 0) {
    setParallerNum(detectCores() - 4)
  }
  registerDoParallel(paraller_num)
#  paraller_register <<- TRUE
}

getParLoopNum <- function(sim_num)
{
  out = NULL
  cores_dopar = paraller_num
  if (sim_num <= cores_dopar) {
    parloop = rep(1, sim_num)
    out = list(sim_num, parloop)
  } else {
    loopcount = floor(sim_num / cores_dopar)         # cores
    parloop = rep(loopcount, cores_dopar)
    r = sim_num %% cores_dopar
    j = 1
    while (j <= r) {
      parloop[j] = parloop[j] + 1
      j = j + 1
    }
    out = list(cores_dopar, parloop)
  }
  names(out) = c("cores", "loop")
  return(out)
}

getParLoop <- function(sim_num)
{
  registerISParallel()
  return(getParLoopNum(sim_num))
}

getParCount <- function(sim_num)
{
  registerISParallel()
  out = NULL
  cores_dopar = paraller_num
  if (sim_num <= cores_dopar) {
    cores_dopar = sim_num
  } else {
    cores_dopar = paraller_num
  }
  loopcount = floor(sim_num / cores_dopar)         # cores
  parloop = rep(loopcount, cores_dopar)
  r = sim_num %% cores_dopar
  j = 1
  while (j <= r) {
    parloop[j] = parloop[j] + 1
    j = j + 1
  }
  parCount = rep(0, cores_dopar)
  j = 2
  while (j <= cores_dopar) {
    parCount[j] = parCount[j-1] + parloop[j-1]
    j = j + 1
  }
  out = list(cores_dopar, parloop, parCount)
  names(out) = c("cores", "loop", "count")
  return(out)
}

simuloop <- function(sim_num, type, func)
{
  if (debug_simu) {
    simuloop_debug(sim_num, type, func)
  }
  return(simuloop_success(sim_num, type, func))
}

simuloop_nopar <- function(sim_num, type, func)
{
  parloop = getParLoop(sim_num)
  simr = foreach(q = 1:parloop$cores, .combine = "rbind") %do% {
    foreach(m = 1:parloop$loop[q], .combine = "rbind") %do% {
      out = func()
    }
  }#(end of q loop)
  return(simr)
}

clusterObj = NULL

parallelStart <- function(sim_num)
{
  parloop = getParLoopNum(sim_num)
  if (isWindow()) {
    clusterObj = makePSOCKcluster(parloop$cores)
  } else {
    clusterObj = makeForkCluster(nnodes = parloop$cores)
    #clusterObj = makeCluster(getOption("cl.cores", parloop$cores))
  }
  registerDoParallel(clusterObj)
}

parallelStop <- function()
{
  stopCluster(clusterObj)
}

simuloop_windows <- function(sim_num, type, func)
{
  cores = detectCores()
  if (sim_num>=cores) {
    cores = cores-1
  }
  if(cores<1) {
    cores = 1
  }
  cl = makeCluster(cores)
  if (type=="ASD") {
    clusterExport(cl, "coxph")
    clusterExport(cl, "last")
    clusterExport(cl, "pmvnorm")
    clusterExport(cl, "Surv")
    clusterExport(cl, "asd_c_power_theta")
    clusterExport(cl, "new_c_bry_com")
    clusterExport(cl, "new_info_time_com")
  } else if (type=="2ASD") {
    clusterExport(cl, "coxph")
    clusterExport(cl, "last")
    clusterExport(cl, "pmvnorm")
    clusterExport(cl, "Surv")
    clusterExport(cl, "asd_c_power_theta")
    clusterExport(cl, "new_c_bry_com")
    clusterExport(cl, "new_info_time_com")
  } else if (type=="MASD") {
    clusterExport(cl, "coxph")
    clusterExport(cl, "pmvnorm")
    clusterExport(cl, "Surv")
    clusterExport(cl, "masd_c_power")
    clusterExport(cl, "masd_c_power_new_theta")
    clusterExport(cl, "masd_c_tyIer")
    clusterExport(cl, "masd_new_info_time")
    clusterExport(cl, "asd_c_power_new")
    clusterExport(cl, "asd_new_info_time")
    clusterExport(cl, "new_OF_bundry")
    clusterExport(cl, "new_OF_bundry_s")
  } else if (type == "SINGLE") {
    #      clusterExport(cl, "pmvnorm")
  }
  registerDoParallel(cl)

  simr = foreach(q = 1:sim_num, .combine = "rbind") %dopar% {
    out = func()
  }#(end of q loop)

  stopCluster(cl)
  return(simr)
}

simuloop_linux_dopar <- function(sim_num, type, func)
{
  parloop = getParLoop(sim_num)
  simr = foreach(q = 1:parloop$cores, .combine = "rbind") %dopar% {
    foreach(m = 1:parloop$loop[q], .combine = "rbind") %do% {
      out = func()
    }
  }#(end of q loop)
  return(simr)
}

simuloop_success <- function(sim_num, type, func)
{
  if (isWindow()) {
    simr = simuloop_windows(sim_num, type, func)
  } else {
    simr = simuloop_linux_dopar(sim_num, type, func)
  }
  return(simr)
}

call2str_info <- function(call)
{
  str = deparse(expr=call, width.cutoff=500L)
  str = gsub(" ", "", str)
  str = paste(str, collapse="")
  return(str)
}

get_exec_info <- function(call, msg, info, row)
{
  k = 1
  while (k < row) {
    if (call==info[k,1] && msg==info[k,2]) {
      return(TRUE)
    }
    k = k+1
  }
  return(FALSE)
}

print_exec_info <- function(info, row, tips1, tips2)
{
  if (row>1) {
    print(paste0(tips1, row-1), quote=FALSE)
  }
  k = 1
  while (k < row) {
    print(paste0(tips2, k, " info:"), quote=FALSE)
    print(paste0("    Call: ", info[k,1]), quote=FALSE)
    print(paste0("    Message: ", info[k,2]), quote=FALSE)
    k = k+1
  }
}

simuloop_debug <- function(sim_num, type, func)
{
  warn_func <- function(w) {
    out = c(-99, "simuW",  w$call, w$message)
  }
  error_func <- function(e) {
    out = c(-199, "simuE", e$call, e$message)
  }

  parloop = getParLoop(sim_num)
  simr = foreach(q = 1:parloop$cores) %dopar% {
    foreach(m = 1:parloop$loop[q]) %do% {
      tryCatch({
        out = func()
      }, warning = function(w) {
        out = warn_func(w)
      }, error = function(e) {
        out = error_func(e)
      }, finally = {

      })
    }
  }#(end of q loop)

  mylist = unlist(simr)
  rowE_max = 20
  rowW_max = 20
  infoE = array(data="", dim=c(rowE_max, 2))
  infoW = array(data="", dim=c(rowW_max, 2))
  rowE = 1
  rowW = 1
  for (m in 1:length(mylist)) {
    if (is.numeric(mylist[[m]]) && ((mylist[[m]]==-99  && mylist[[m + 1]]=="simuW") ||
                                    (mylist[[m]]==-199 && mylist[[m + 1]]=="simuE"))) {
      call = call2str_info(mylist[[m + 2]])
      msg = as.character(mylist[[m + 3]])
      if (mylist[[m]] == -99 && rowW<=rowW_max && !get_exec_info(call, msg, infoW, rowW)) {
        infoW[rowW,1] = call
        infoW[rowW,2] = msg
        rowW = rowW+1
      } else if (mylist[[m]] == -199 && rowE<=rowE_max && !get_exec_info(call, msg, infoE, rowE)) {
        infoE[rowE,1] = call
        infoE[rowE,2] = msg
        rowE = rowE+1
      }
    }
  }

  print_exec_info(infoW, rowW, "Total number of warnings: ", "Warning")
  print_exec_info(infoE, rowE, "Total number of errors: ", "Error")

  if (rowE>1) {
    stop(paste0("\n\nDebug info start ==>>\nSimuError: ", "",
                "\nCall: ", infoE[1,1],
                "\nMessage: ", infoE[1,2],
                "\n<<== Debug info end\n"))
  }
}

call_uniroot <- function(func, interval, extendInt, maxiter=1000L, tol=1e-6) {
  PRINT(paste0("call_uniroot begin           time:", Sys.time()))
  priv_uniroot <- function(lower, upper) {
    #tryCatch({
    out = stats::uniroot(function(x) func(x),
                  interval = c(lower, upper),
                  extendInt = "no",
                  maxiter = maxiter,
                  tol = tol)[1]
    #print(paste0("out=",out))
    #myout[1]=out
    #stop("success out", domain =TRUE)
    #})
  }

  # left>1
  priv_low_high <- function(left, right, corecount) {
    num = left
    low = log2(left)
    high = log2(right)
    delt = (high-low)/corecount
    fact = 2^delt
    lower1 = NULL
    for(i in 1:corecount) {
      lower1 = rbind(lower1, num)
      num = num * fact
    }
    lower = lower1
    upper = lower1 * fact
    return(list(lower,upper))
  }

  # right<-1
  priv_low_high2 <- function(left, right, corecount) {
    num = left
    low = log2(-left)
    high = log2(-right)
    delt = (low-high)/corecount
    fact = 2^delt
    lower1 = NULL
    for(i in 1:corecount) {
      lower1 = rbind(lower1, num)
      num = num / fact
    }
    lower = lower1
    upper = lower1 / fact
    return(list(lower,upper))
  }

  # left<-1, right>1
  priv_low_high3 <- function(left, right, corecount) {
    low = log2(-left)
    high = log2(right)
    delt = (low+high+1)/corecount
    leftcore = trunc(low*delt)
    rightcore = corecount-leftcore-1

    ret=priv_low_high2(left, -1, leftcore)
    lower=rbind(ret[[1]], -1)
    upper=rbind(ret[[2]], 1)

    ret=priv_low_high(1, right, rightcore)
    lower=rbind(lower, ret[[1]])
    upper=rbind(upper, ret[[2]])
    return(list(lower,upper))
  }

  registerISParallel()
  corecount = getParallerNum()
  num=1
  if (interval[1] >= 1) {         # interval[1] >= 1
    ret=priv_low_high(interval[1], interval[2], corecount)
    lower=ret[[1]]
    upper=ret[[2]]
    upper[corecount] = interval[2]
  } else if (interval[1] >= 0) {  # 0 <= interval[1] < 1
    ret=priv_low_high(1, interval[2], corecount-1)
    lower=rbind(interval[1], ret[[1]])
    upper=rbind(1, ret[[2]])
    upper[corecount] = interval[2]
  } else if (interval[2] <= -1) { # interval[1] < 0 && interval[2] <= -1
    ret=priv_low_high2(interval[1], interval[2], corecount)
    lower=ret[[1]]
    upper=ret[[2]]
    upper[corecount] = interval[2]
  } else if (interval[2] <= 1) {  # interval[1] < 0 && -1 < interval[2] <= 1
    ret=priv_low_high2(interval[1], -1, corecount-1)
    lower=rbind(ret[[1]], -1)
    upper=rbind(ret[[2]], interval[2])
  } else if (interval[1] <= -1) {  # interval[1] < -1 && interval[2] > 1
    ret=priv_low_high3(interval[1], interval[2], corecount)
    lower=ret[[1]]
    upper=ret[[2]]
    upper[corecount] = interval[2]
  } else {                        # -1 < interval[1] < 0 || interval[2] > 1
    ret=priv_low_high(1, interval[2], corecount-1)
    lower=rbind(interval[1], ret[[1]])
    upper=rbind(1, ret[[2]])
    upper[corecount] = interval[2]
  }

  simr = foreach(q = 1:corecount, .combine = "rbind", .errorhandling = "remove") %dopar% {
    priv_uniroot(lower[q], upper[q])
  }#(end of q loop)

  #  print(geterrmessage())
  #  print(paste0("simr=",simr))
  time_tot=as.numeric(simr)

  PRINT(paste0("call_uniroot end             time:", Sys.time()))

  return(time_tot)
}

if (free_version) {  # nodename
#  && (grep("is-rd", Sys.info()[4]) == 0)
#  use_cores = default_use_cores
#  setParallerNum(use_cores)
}
