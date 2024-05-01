#' @import odin
NULL

#' Create Malaria Model
#'
#' This function creates a malaria model using user-defined parameters and returns an `odin` model object that can be used to run simulations.
#'
#' @param S0 Initial proportion of susceptible human population
#' @param Ds0 Initial proportion of diseased human population (sensitive strain)
#' @param As0 Initial proportion of asymptomatic human population (sensitive strain)
#' @param Ts0 Initial proportion of treated human population (sensitive strain)
#' @param DR0 Initial proportion of diseased human population (resistant strain)
#' @param AR0 Initial proportion of asymptomatic human population (resistant strain)
#' @param TR0 Initial proportion of treated human population (resistant strain)
#' @param Sv0 Initial proportion of susceptible mosquito population
#' @param Ev_s0 Initial proportion of exposed mosquito population (sensitive strain)
#' @param Iv_s0 Initial proportion of infectious mosquito population (sensitive strain)
#' @param Ev_r0 Initial proportion of exposed mosquito population (resistant strain)
#' @param Iv_r0 Initial proportion of infectious mosquito population (resistant strain)
#' @param mu Mosquito death rate
#' @param e Mosquito emergence rate
#' @param n Extrinsic incubation period
#' @param Lambda_s Force of infection (sensitive strain)
#' @param Lambda_R Force of infection (resistant strain)
#' @param Lambda_v_s Force of infection from humans to mosquitoes (sensitive strain)
#' @param Lambda_v_r Force of infection from humans to mosquitoes (resistant strain)
#' @param Phi Proportion of symptomatic infections
#' @param fT Proportion of symptomatic infections treated
#' @param rD Recovery rate for diseased individuals
#' @param rA Recovery rate for asymptomatic individuals
#' @param rTs Recovery rate for treated individuals (sensitive strain)
#' @param rTR Recovery rate for treated individuals (resistant strain)
#' @return An object of class `odin_model`.
#' @export
#' @examples
#' model <- create_malaria_model(S0 = 0.9, Ds0 = 0.05, As0 = 0.025, Ts0 = 0.0125,
#'                               DR0 = 0.005, AR0 = 0.0025, TR0 = 0.001,
#'                               Lambda_s = 0.3, Lambda_R = 0.1, Phi = 0.7,
#'                               fT = 0.2, rD = 0.05, rA = 0.03, rTs = 0.02,
#'                               rTR = 0.01, Sv0 = 0.9, Ev_s0 = 0, Iv_s0 = 0,
#'                               Ev_r0 = 0, Iv_r0 = 0,
#'                               e = 0.1, mu = 0.1, n = 10,
#'                               Lambda_v_s = 0.005, Lambda_v_r = 0.003)
#' results <- model$run(0:500)
#' print(results)
create_malaria_model <- function(S0, Ds0, As0, Ts0, DR0, AR0, TR0, Lambda_s, Lambda_R, Phi, fT, rD, rA, rTs, rTR,
                                 Sv0, Ev_s0, Iv_s0, Ev_r0, Iv_r0, e, mu, n, Lambda_v_s, Lambda_v_r) {
  model_code <- readLines(system.file("odin/model.R", package = "AMRSpreadModel"))
  model <- odin::odin(model_code)
  model$new(user = list(S0 = S0, Ds0 = Ds0, As0 = As0, Ts0 = Ts0, DR0 = DR0, AR0 = AR0, TR0 = TR0,
                        Lambda_s = Lambda_s, Lambda_R = Lambda_R, Phi = Phi, fT = fT,
                        rD = rD, rA = rA, rTs = rTs, rTR = rTR,
                        Sv0 = Sv0, Ev_s0 = Ev_s0, Iv_s0 = Iv_s0, Ev_r0 = Ev_r0, Iv_r0 = Iv_r0,
                        e = e, mu = mu, n = n, Lambda_v_s = Lambda_v_s, Lambda_v_r = Lambda_v_r))
}
