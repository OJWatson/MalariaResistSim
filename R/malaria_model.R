#' @import odin
#' @importFrom Rcpp sourceCpp
#' @useDynLib AMRSpreadModel, .registration = TRUE
#' @importFrom AMRSpreadModel phi_eir_rel
NULL

#' Create Malaria Model
#'
#' @param params A list of parameters to pass to the model. If NULL, parameters will be generated based on EIR and ft.
#' @param EIR Entomological Inoculation Rate (used if params is NULL)
#' @param ft Treatment rate (used if params is NULL)
#' @param ton Time at which treatment is turned on
#' @param toff Time at which treatment is turned off
#' @param init_res Initial resistance level at res_time
#' @param day0_res Resistant at Day 0. Default = 0.01
#' @param res_time Time at which resistance is introduced
#' @param rTR_true True treatment rate for resistant parasites
#' @param verbose Logical. If TRUE, prints detailed logs. Default is FALSE.
#' @return An object of class `odin_model`.
#' @export
malaria_model <- function(params = NULL, EIR = NULL, ft = NULL,
                          ton = 5000, toff = 50000, day0_res = 0,
                          init_res = 0.01, res_time = 3000, rTR_true = 0.1,
                          verbose = FALSE) {
  tryCatch({
    model_file <- system.file("odin", "model.R", package = "AMRSpreadModel")
    if (!file.exists(model_file)) {
      stop("Model file not found: ", model_file)
    }

    if (is.null(params)) {
      if (is.null(EIR) || is.null(ft)) {
        stop("Either params or both EIR and ft must be provided")
      }
      params <- phi_eir_rel(EIR, ft, ton, toff, init_res, res_time, rTR_true)
    }

    # Update parameters
    params$ton <- ton
    params$toff <- toff
    params$init_res <- init_res
    params$res_time <- res_time
    params$rTR_true <- rTR_true

    # Generate initial parameters if not already present
    if (!"S0" %in% names(params)) {
      params$S0 <- params$S
      params$Ds0 <- params$D * (1 - day0_res)
      params$DR0 <- params$D * day0_res
      params$As0 <- params$A * (1 - day0_res)
      params$AR0 <- params$A * day0_res
      params$Ts0 <- params$T * (1 - day0_res)
      params$TR0 <- params$T * day0_res
      params$Sv0 <- params$Sv
      params$Ev_s0 <- params$Ev * (1 - day0_res)
      params$Iv_s0 <- params$Iv * (1 - day0_res)
      params$Ev_r0 <- params$Ev * day0_res
      params$Iv_r0 <- params$Iv * day0_res
    }

    # Rename and assign parameters
    params$Phi <- params$Phi %||% params$phi
    params$fT <- params$fT %||% params$ft
    params$rTs <- params$rTs %||% params$rT_S
    params$rTR <- params$rTR_true
    params$e <- params$e %||% params$mu
    params$c_A <- params$c_A %||% params$cA
    params$c_D <- params$c_D %||% params$cD
    params$c_T <- params$c_T %||% params$cT

    # Ensure all parameters are scalar
    params <- lapply(params, function(x) if(length(x) > 1) x[1] else x)

    # Check if all required parameters are present and numeric
    required_params <- c("S0", "Ds0", "As0", "Ts0", "DR0", "AR0", "TR0",
                         "Sv0", "Ev_s0", "Iv_s0", "Ev_r0", "Iv_r0",
                         "m", "a", "b", "Phi", "fT", "rD", "rA", "rTs", "rTR",
                         "e", "mu", "n", "c_A", "c_D", "c_T",
                         "ton", "toff", "res_time", "init_res")
    missing_params <- setdiff(required_params, names(params))
    if (length(missing_params) > 0) {
      stop("Missing required parameters: ", paste(missing_params, collapse = ", "))
    }

    # Calculate EIR if not provided
    if (!"EIR" %in% names(params)) {
      params$EIR <- params$m * params$a * params$Iv_s0 * 365  # Annual EIR
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

# Helper function to provide a default value
`%||%` <- function(x, y) if (is.null(x)) y else x
