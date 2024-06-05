#' @import odin
#' @importFrom utils modifyList
NULL

#' Create Malaria Model
#'
#' This function creates a malaria model using user-defined parameters and returns an object of class `odin_model`.
#'
#' @param params A list of parameters to pass to the model, including:
#' - S0, Ds0, As0, Ts0, DR0, AR0, TR0: Initial proportions of human population in different compartments
#' - Sv0, Ev_s0, Iv_s0, Ev_r0, Iv_r0: Initial proportions of mosquito population in different compartments
#' - m: Ratio of vectors to humans
#' - a: Biting rate (bites on humans per mosquito per day)
#' - b: Probability of transmission from vectors to humans
#' - Phi: Proportion of symptomatic infections
#' - fT: Proportion of symptomatic infections treated
#' - rD, rA, rTs, rTR_true: Recovery rates for different human compartments
#' - e: Mosquito emergence rate
#' - mu: Mosquito death rate
#' - n: Extrinsic incubation period
#' - c_A, c_D, c_T: Probabilities of transmission from different human compartments to vectors
#' - ton, toff: Time points for changing the recovery rate for treated resistant infections
#' - res_time: Time when resistance starts to appear
#' - res_start: Proportion of resistant infections at res_time
#' @return An object of class `odin_model`.
#' @export
#' @examples
#' params <- list(S0 = 0.9, Ds0 = 0.05, As0 = 0.03, Ts0 = 0.02, DR0 = 0, AR0 = 0,
#'                TR0 = 0, Sv0 = 0.9, Ev_s0 = 0.05, Iv_s0 = 0.05, Ev_r0 = 0, Iv_r0 = 0,
#'                m = 0.7, a = 0.3, b = 0.1, Phi = 0.7, fT = 0.2, rD = 0.05,
#'                rA = 0.03, rTs = 0.02, rTR_true = 0.01, e = 0.1, mu = 0.1,
#'                n = 10, c_A = 0.05, c_D = 0.06, c_T = 0.07, ton = 10000,
#'                toff = 20000, res_time = 50, res_start = 0.5)
#' model <- malaria_model(params)
#' results <- model$run(0:500)
#' print(results)
malaria_model <- function(params) {
  model_file <- system.file("odin", "model.R", package = "AMRSpreadModel")
  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
  }
  model <- odin::odin(model_file, verbose = FALSE)
  model_instance <- model$new(user = form_pars(params))
  return(model_instance)
}

#' Format parameters for the malaria model
#'
#' @param pars Parameter List
#'
#' @return Formatted parameter list
form_pars <- function(pars) {
  default_pars <- list(
    S0 = 0.9, Ds0 = 0.05, As0 = 0.03, Ts0 = 0.02, DR0 = 0, AR0 = 0, TR0 = 0,
    Sv0 = 0.9, Ev_s0 = 0.05, Iv_s0 = 0.05, Ev_r0 = 0, Iv_r0 = 0,
    m = 0.7, a = 0.3, b = 0.1, Phi = 0.7, fT = 0.2, rD = 0.05,
    rA = 0.03, rTs = 0.02, rTR_true = 0.01, e = 0.1, mu = 0.1,
    n = 10, c_A = 0.05, c_D = 0.06, c_T = 0.07, ton = 10000,
    toff = 20000, res_time = 50, res_start = 0.5
  )

  pars <- modifyList(default_pars, pars)

  total_D <- pars$Ds0 + pars$DR0
  total_A <- pars$As0 + pars$AR0
  total_T <- pars$Ts0 + pars$TR0
  total_Ev <- pars$Ev_s0 + pars$Ev_r0
  total_Iv <- pars$Iv_s0 + pars$Iv_r0

  pars$Ds0 <- total_D * (1 - pars$res_start)
  pars$As0 <- total_A * (1 - pars$res_start)
  pars$Ts0 <- total_T * (1 - pars$res_start)
  pars$Ev_s0 <- total_Ev * (1 - pars$res_start)
  pars$Iv_s0 <- total_Iv * (1 - pars$res_start)

  pars$DR0 <- total_D * pars$res_start
  pars$AR0 <- total_A * pars$res_start
  pars$TR0 <- total_T * pars$res_start
  pars$Ev_r0 <- total_Ev * pars$res_start
  pars$Iv_r0 <- total_Iv * pars$res_start

  return(pars)
}
