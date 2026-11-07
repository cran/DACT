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

LABEL <- function() {
  return(list(
    mu_t = "mu_t",
    mu_t1 = "mu_t1,...",
    mu_t_group = "mu_t group",
    mu_c = "mu_c",
    lambda_t = "lambda_t",
    lambda_t1 = "lambda_t1,...",
    lambda_t_group = "lambda_t group",
    lambda_c = "lambda_c",
    empty = ""
  ))
}

# LABEL <- function() {
#   return(list(
#     mu_t = "μ<sub>t</sub>",
#     mu_t1 = "μ<sub>t1</sub>,...",
#     mu_t_group = "μ<sub>t</sub> group",
#     mu_c = "μ<sub>c</sub>",
#     lambda_t = "λ<sub>t</sub>",
#     lambda_t1 = "λ<sub>t1</sub>,...",
#     lambda_t_group = "λ<sub>t</sub> group",
#     lambda_c = "λ<sub>c</sub>",
#     empty = ""
#   ))
# }


