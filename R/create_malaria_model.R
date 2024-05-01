#' @import odin
NULL

#' Create Malaria Model
#'
#' This function creates a malaria model using user-defined parameters and returns an `odin` model object that can be used to run simulations.
#'
#' @param S0 Initial susceptible population
#' @param Ds0 Initial Ds population
#' @param As0 Initial As population
#' @param Ts0 Initial Ts population
#' @param DR0 Initial DR population
#' @param AR0 Initial AR population
#' @param TR0 Initial TR population
#' @param Lambda_s Force of infection from vectors to humans (general)
#' @param Lambda_R Force of infection from vectors to humans (Resistant)
#' @param Phi Probability of developing symptoms after infection
#' @param fT Frequency of treatment after symptom onset
#' @param rD Rate of recovery from symptomatic to asymptomatic disease
#' @param rA Rate of recovery from asymptomatic to susceptible status
#' @param rTs Rate of recovery from treatment to susceptible status
#' @param rTR Rate of recovery from treatment to resistant status
#' @param Sv0 Initial susceptible mosquito population
#' @param Ev_s0 Initial exposed mosquito population (susceptible strain)
#' @param Iv_s0 Initial infectious mosquito population (susceptible strain)
#' @param Ev_r0 Initial exposed mosquito population (resistant strain)
#' @param Iv_r0 Initial infectious mosquito population (resistant strain)
#' @param e Egg laying rate of mosquitoes
#' @param mu Mortality rate of mosquitoes
#' @param n Delay in days for mosquito development
#' @param Lambda_v_s Infection rate for susceptible strain
#' @param Lambda_v_r Infection rate for resistant strain
#' @return An object of class `odin_model`.
#' @export
#' @examples
#' model <- create_malaria_model(S0 = 1000, Ds0 = 100, As0 = 50, Ts0 = 25, DR0 = 10, AR0 = 5, TR0 = 2,
#'                               Lambda_s = 0.3, Lambda_R = 0.1, Phi = 0.7, fT = 0.2,
#'                               rD = 0.05, rA = 0.03, rTs = 0.02, rTR = 0.01,
#'                               Sv0 = 500, Ev_s0 = 0, Iv_s0 = 0, Ev_r0 = 0, Iv_r0 = 0,
#'                               e = 100, mu = 0.1, n = 10, Lambda_v_s = 0.005, Lambda_v_r = 0.003)
#' results <- model$run(0:50)
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
