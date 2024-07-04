#' @import odin
#' @importFrom Rcpp sourceCpp
#' @useDynLib AMRSpreadModel, .registration = TRUE
NULL

#' Create Malaria Model
#'
#' @param params A list of parameters to pass to the model. If NULL, parameters will be generated based on EIR and ft.
#' @param EIR Entomological Inoculation Rate (used if params is NULL)
#' @param ft Treatment rate (used if params is NULL)
#' @param ton Time at which treatment is turned on
#' @param toff Time at which treatment is turned off
#' @param init_res Initial resistance level
#' @param res_time Time at which resistance is introduced
#' @param rTR_true True treatment rate for resistant parasites
#' @param verbose Logical. If TRUE, prints detailed logs. Default is FALSE.
#' @return An object of class `odin_model`.
#' @export
malaria_model <- function(params = NULL, EIR = NULL, ft = NULL,
                          ton = 5000, toff = 50000, init_res = 0.01, res_time = 3000, rTR_true = 0.1,
                          verbose = FALSE) {
  tryCatch({
    model_file <- system.file("odin", "model.R", package = "AMRSpreadModel")
    if (!file.exists(model_file)) {
      stop("Model file not found: ", model_file)
    }

    if (is.null(params) && (is.null(EIR) || is.null(ft))) {
      stop("Either params or both EIR and ft must be provided")
    }

    if (is.null(params)) {
      params <- phi_eir_rel(EIR, ft, ton, toff, init_res, res_time, rTR_true)
    } else {
      # If params are provided, update with user-specified values
      params$ton <- ton
      params$toff <- toff
      params$init_res <- init_res
      params$res_time <- res_time
      params$rTR_true <- rTR_true
    }

    # Ensure all parameters are scalar
    params <- lapply(params, function(x) if(length(x) > 1) x[1] else x)

    # Check if all required parameters are present and numeric
    required_params <- c("S0", "Ds0", "As0", "Ts0", "DR0", "AR0", "TR0",
                         "Sv0", "Ev_s0", "Iv_s0", "Ev_r0", "Iv_r0",
                         "m", "a", "b", "Phi", "fT", "rD", "rA", "rTs", "rTR_true",
                         "e", "mu", "n", "c_A", "c_D", "c_T",
                         "ton", "toff", "res_time", "init_res")
    missing_params <- setdiff(required_params, names(params))
    if (length(missing_params) > 0) {
      stop("Missing required parameters: ", paste(missing_params, collapse = ", "))
    }

    non_numeric_params <- names(params)[!sapply(params, is.numeric)]
    if (length(non_numeric_params) > 0) {
      stop("Non-numeric parameters: ", paste(non_numeric_params, collapse = ", "))
    }

    if (verbose) {
      cat("Creating odin model...\n")
    }
    model <- odin::odin(model_file, verbose = FALSE)
    if (verbose) {
      cat("odin model created successfully.\n")
    }

    if (verbose) {
      cat("Initializing model with parameters...\n")
    }
    model_instance <- model$new(user = params, unused_user_action = "ignore")
    if (verbose) {
      cat("Model initialized successfully.\n")
    }

    return(model_instance)
  }, error = function(e) {
    message("Error in malaria_model:")
    message(e$message)
    return(NULL)
  })
}

#' Format parameters for the malaria model
#'
#' @param pars A list of parameters to format
#' @return A formatted list of parameters
#' @keywords internal
form_pars <- function(pars) {
  default_pars <- list(
    S0 = 0.9, Ds0 = 0.05, As0 = 0.02, Ts0 = 0.03, DR0 = 0, AR0 = 0, TR0 = 0,
    Sv0 = 0.9, Ev_s0 = 0.08, Iv_s0 = 0.02, Ev_r0 = 0, Iv_r0 = 0,
    m = 10, a = 0.3, b = 0.5876259, Phi = 0.7, fT = 0.1, rD = 0.2,
    rA = 0.01, rTs = 0.2, rTR_true = 0.01, e = 0.132, mu = 0.132,
    n = 10, c_A = 0.05, c_D = 0.06, c_T = 0.02, ton = 8000,
    toff = 99999999, res_time = 3000, init_res = 0.2
  )
  pars <- modifyList(default_pars, pars)
  if (pars$res_time == 0) {
    initial_pars <- list(
      Ds0 = pars$Ds0 * (1 - pars$init_res),
      DR0 = pars$Ds0 * pars$init_res,
      As0 = pars$As0 * (1 - pars$init_res),
      AR0 = pars$As0 * pars$init_res,
      Ts0 = pars$Ts0 * (1 - pars$init_res),
      TR0 = pars$Ts0 * pars$init_res,
      Ev_s0 = pars$Ev_s0 * (1 - pars$init_res),
      Iv_s0 = pars$Iv_s0 * (1 - pars$init_res),
      Ev_r0 = pars$Ev_s0 * pars$init_res,
      Iv_r0 = pars$Iv_s0 * pars$init_res
    )
    pars <- modifyList(pars, initial_pars)
  }
  return(pars)
}
