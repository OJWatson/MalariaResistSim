#' Run Malaria Model for Ranges of fT and EIR Values
#'
#' This function runs the malaria model for ranges of fT and EIR values and returns a data frame of results.
#'
#' @param fT_range A vector of fT values to run the model over
#' @param EIR_s_range A vector of EIR_s values to run the model over
#' @param EIR_R_range A vector of EIR_R values to run the model over
#' @param params A list of other parameters to pass to the `create_model` function
#' @return A data frame of results for each combination of fT and EIR values
#' @export
#' @examples
#' fT_range <- seq(0.1, 0.9, by = 0.1)
#' EIR_s_range <- seq(10, 200, by = 10)
#' EIR_R_range <- rep(0, length(EIR_s_range))
#' params <- list(S0 = 0.9, Ds0 = 0.05, As0 = 0.025, Ts0 = 0.0125,
#'                DR0 = 0, AR0 = 0, TR0 = 0,
#'                b_h = 0.1, Phi = 0.7,
#'                rD = 0.05, rA = 0.03, rTs = 0.02, rTR = 0.01,
#'                Sv0 = 0.9, Ev_s0 = 0, Iv_s0 = 0, Ev_r0 = 0, Iv_r0 = 0,
#'                e = 0.1, mu = 0.1, n = 10, Lambda_v_s = 0.005, Lambda_v_r = 0)
#' results <- run_model_for_ranges(fT_range, EIR_s_range, EIR_R_range, params)
#' print(results)
run_model_for_ranges <- function(fT_range, EIR_s_range, EIR_R_range, params) {
  results_list <- expand.grid(fT = fT_range, EIR_s = EIR_s_range, EIR_R = EIR_R_range)
  results_list <- apply(results_list, 1, function(p) {
    param_list <- params
    param_list$fT <- p["fT"]
    param_list$EIR_s <- p["EIR_s"]
    param_list$EIR_R <- p["EIR_R"]
    res <- create_model(param_list)
    data.frame(fT = p["fT"], EIR_s = p["EIR_s"], Prevalence = res$Prevalence)
  })
  results_df <- do.call(rbind, results_list)
  return(results_df)
}
