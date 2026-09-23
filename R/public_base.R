#library(jsonlite)
#library(dplyr)
#' @import jsonlite
#' @import dplyr

#' @import clinfun
#' @import mvtnorm

isWindow <- function()
{
  return(substr(utils::osVersion,1,7)=="Windows")
}

PRINT <- function(info) {
}

