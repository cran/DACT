#library(jsonlite)
#library(dplyr)
#' @import jsonlite
#' @import dplyr
#' @import clinfun
#' @import mvtnorm
#' @import parallel
#' @import foreach
#' @import doParallel
#' @rawNamespace import(MASS, except = c(select))
#' @import survival

isWindow <- function()
{
  return(substr(utils::osVersion,1,7)=="Windows")
}

PRINT <- function(info) {
}

