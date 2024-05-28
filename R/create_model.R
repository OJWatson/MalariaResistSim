#' @import odin
NULL

#' Create Malaria Model
#'
#' This function creates a malaria model using user-defined parameters and returns a list containing the `odin` model object, FOI, EIR_s, EIR_R, Prevalence, and Prevalence_R.
#'
#' @param params A list of parameters for the model including:
#'   \itemize{
#'     \item S0: Initial proportion of susceptible human population
#'     \item Ds0: Initial proportion of diseased human population (sensitive strain)
#'     \item As0: Initial proportion of asymptomatic human population (sensitive strain)
#'     \item Ts0: Initial proportion of treated human population (sensitive strain)
#'     \item DR0: Initial proportion of diseased human population (resistant strain)
#'     \item AR0: Initial proportion of asymptomatic human population (resistant strain)
#'     \item TR0: Initial proportion of treated human population (resistant strain)
#'     \item Sv0: Initial proportion of susceptible mosquito population
#'     \item Ev_s0: Initial proportion of exposed mosquito population (sensitive strain)
#'     \item Iv_s0: Initial proportion of infectious mosquito population (sensitive strain)
#'     \item Ev_r0: Initial proportion of exposed mosquito population (resistant strain)
#'     \item Iv_r0: Initial proportion of infectious mosquito population (resistant strain)
#'     \item mu: Mosquito death rate
#'     \item e: Mosquito emergence rate
#'     \item n: Extrinsic incubation period
#'     \item a: Biting rate (bites on humans per mosquito per day) (optional)
#'     \item b_h: Probability of transmission from vectors to humans
#'     \item m: Ratio of vectors to humans (optional)
#'     \item Phi: Proportion of symptomatic infections
#'     \item fT: Proportion of symptomatic infections treated
#'     \item rD: Recovery rate for diseased individuals
#'     \item rA: Recovery rate for asymptomatic individuals
#'     \item rTs: Recovery rate for treated individuals (sensitive strain)
#'     \item rTR: Recovery rate for treated individuals (resistant strain)
#'     \item Lambda_v_s: Force of infection from humans to mosquitoes (sensitive strain)
#'     \item Lambda_v_r: Force of infection from humans to mosquitoes (resistant strain)
#'     \item EIR_s: Entomological Inoculation Rate for sensitive strain (optional)
#'     \item EIR_R: Entomological Inoculation Rate for resistant strain (optional)
#'   }
#' @return A list containing the `odin` model object, FOI (EIR_s and EIR_R), EIR_s, EIR_R, Prevalence, and Prevalence_R.
#' @export
#' @examples
#' params <- list(S0 = 0.9, Ds0 = 0.05, As0 = 0.025, Ts0 = 0.0125,
#'                DR0 = 0.005, AR0 = 0.0025, TR0 = 0.001,
#'                a = 0.3, b_h = 0.1, m = 0.7, Phi = 0.7,
#'                fT = 0.2, rD = 0.05, rA = 0.03, rTs = 0.02,
#'                rTR = 0.01, Sv0 = 0.9, Ev_s0 = 0, Iv_s0 = 0,
#'                Ev_r0 = 0, Iv_r0 = 0,
#'                e = 0.1, mu = 0.1, n = 10,
#'                Lambda_v_s = 0.005, Lambda_v_r = 0.003,
#'                EIR_s = NULL, EIR_R = NULL)
#' results <- create_model(params)
#' model_results <- results$model$run(0:500)
#' prevalence <- rowSums(model_results[, c("AR", "TR", "DR", "As", "Ds", "Ts")])
#' print(prevalence)
#' print(results$FOI)
#' print(results$EIR_s)
#' print(results$EIR_R)
create_model <- function(params) {
  with(as.list(params), {
    # If EIR_s is NULL, calculate EIR_s
    if (is.null(EIR_s)) {
      if (is.null(m) || is.null(a)) {
        stop("Either EIR_s or both m and a must be provided.")
      }
      EIR_s <- m * a * Iv_s0 * 365
    }

    # If EIR_R is NULL, calculate EIR_R
    if (is.null(EIR_R)) {
      if (is.null(m) || is.null(a)) {
        stop("Either EIR_R or both m and a must be provided.")
      }
      EIR_R <- m * a * Iv_r0 * 365
    }

    # Read the model code
    model_code <- '
    deriv(S) <- -S * EIR_s * b_h * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) -
                S * EIR_R * b_h * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) +
                Ts * rTs + As * rA + AR * rA + TR * rTR

    deriv(Ds) <- S * EIR_s * b_h * Phi * (1 - fT) +
                 EIR_s * b_h * As * Phi * (1 - fT) -
                 Ds * rD

    deriv(As) <- S * EIR_s * b_h * (1 - Phi) +
                 Ds * rD -
                 EIR_s * b_h * As * Phi * (1 - fT) -
                 EIR_s * b_h * As * Phi * fT -
                 As * rA

    deriv(Ts) <- S * EIR_s * b_h * Phi * fT +
                 EIR_s * b_h * As * Phi * fT -
                 Ts * rTs

    deriv(DR) <- S * EIR_R * b_h * Phi * (1 - fT) +
                 EIR_R * b_h * AR * Phi * (1 - fT) -
                 DR * rD

    deriv(AR) <- S * EIR_R * b_h * (1 - Phi) +
                 DR * rD -
                 EIR_R * b_h * AR * Phi * (1 - fT) -
                 EIR_R * b_h * AR * Phi * fT -
                 AR * rA

    deriv(TR) <- S * EIR_R * b_h * Phi * fT +
                 EIR_R * b_h * AR * Phi * fT -
                 TR * rTR

    ## Mosquito Equations
    deriv(Sv) <- e - (Lambda_v_s + Lambda_v_r) * Sv - mu * Sv

    delayed_Lambda_v_s_Sv_raw <- delay(Lambda_v_s * Sv, n)
    delayed_Lambda_v_s_Sv <- delayed_Lambda_v_s_Sv_raw * exp(-mu * n)

    deriv(Ev_s) <- Lambda_v_s * Sv - delayed_Lambda_v_s_Sv - mu * Ev_s
    deriv(Iv_s) <- delayed_Lambda_v_s_Sv - mu * Iv_s

    delayed_Lambda_v_r_Sv_raw <- delay(Lambda_v_r * Sv, n)
    delayed_Lambda_v_r_Sv <- delayed_Lambda_v_r_Sv_raw * exp(-mu * n)

    deriv(Ev_r) <- Lambda_v_r * Sv - delayed_Lambda_v_r_Sv - mu * Ev_r
    deriv(Iv_r) <- delayed_Lambda_v_r_Sv - mu * Iv_r

    ## Initial conditions
    initial(S) <- S0
    initial(Ds) <- Ds0
    initial(As) <- As0
    initial(Ts) <- Ts0
    initial(DR) <- DR0
    initial(AR) <- AR0
    initial(TR) <- TR0
    initial(Sv) <- Sv0
    initial(Ev_s) <- Ev_s0
    initial(Iv_s) <- Iv_s0
    initial(Ev_r) <- Ev_r0
    initial(Iv_r) <- Iv_r0

    ## User-defined parameters
    S0 <- user()
    Ds0 <- user()
    As0 <- user()
    Ts0 <- user()
    DR0 <- user()
    AR0 <- user()
    TR0 <- user()
    EIR_s <- user()
    EIR_R <- user()
    Phi <- user()
    fT <- user()
    rD <- user()
    rA <- user()
    rTs <- user()
    rTR <- user()
    Sv0 <- user()
    Ev_s0 <- user()
    Iv_s0 <- user()
    Ev_r0 <- user()
    Iv_r0 <- user()
    e <- user()
    mu <- user()
    n <- user()
    Lambda_v_s <- user()
    Lambda_v_r <- user()
    b_h <- user()
    '

    model <- odin::odin(model_code)

    # Initialize the model instance
    model_instance <- model$new(user = list(S0 = S0, Ds0 = Ds0, As0 = As0, Ts0 = Ts0, DR0 = DR0, AR0 = AR0, TR0 = TR0,
                                            EIR_s = EIR_s, EIR_R = EIR_R, Phi = Phi, fT = fT,
                                            rD = rD, rA = rA, rTs = rTs, rTR = rTR,
                                            Sv0 = Sv0, Ev_s0 = Ev_s0, Iv_s0 = Iv_s0, Ev_r0 = Ev_r0, Iv_r0 = Iv_r0,
                                            e = e, mu = mu, n = n, Lambda_v_s = Lambda_v_s, Lambda_v_r = Lambda_v_r, b_h = b_h))

    # Run the model and calculate Prevalence and Prevalence_R
    model_results <- model_instance$run(0:500)
    prevalence <- rowSums(model_results[, c("AR", "TR", "DR", "As", "Ds", "Ts")])
    final_prevalence <- prevalence[length(prevalence)]
    final_prevalence_R <- (model_results[length(prevalence), "AR"] + model_results[length(prevalence), "TR"] + model_results[length(prevalence), "DR"]) /
      sum(model_results[length(prevalence), c("AR", "TR", "DR", "As", "Ds", "Ts")])

    return(list(model = model_instance, FOI = list(EIR_s = EIR_s, EIR_R = EIR_R), EIR_s = EIR_s, EIR_R = EIR_R,
                Prevalence = final_prevalence, Prevalence_R = final_prevalence_R))
  })
}
