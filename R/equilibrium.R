#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom tidyr expand_grid
NULL

#' @import ICDMM
NULL

#' Generate equilibrium parameters
#'
#' @param EIR Numeric. The Entomological Inoculation Rate.
#' @param ft Numeric. The treatment rate.
#' @param ton Time at which treatment is turned on
#' @param toff Time at which treatment is turned off
#' @param init_res Initial resistance level
#' @param res_time Time at which resistance is introduced
#' @param rTR_true True treatment rate for resistant parasites
#' @return A list of generated parameters.
#' @export
#' @importFrom ICDMM model_param_list_create equilibrium_init_create
#' @importFrom stats weighted.mean
phi_eir_rel <- function(EIR, ft, ton = 5000, toff = 50000, init_res = 0.01, res_time = 3000, rTR_true = 0.1) {
  tryCatch({
    mpl <- ICDMM::model_param_list_create(rho = 0, rA = 1/250, rU = Inf, rP = Inf, sigma2 = 0)

    age_vector <- c(0, 0.25, 0.5, 0.75, 1, 2, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

    eq <- ICDMM::equilibrium_init_create(
      age_vector = age_vector,
      EIR = EIR, ft = ft,
      model_param_list = mpl, het_brackets = 2,
      country = NULL,
      admin_unit = NULL
    )

    if (is.null(eq$phi0) || length(eq$phi0) == 0) {
      stop("eq$phi0 is NULL or empty")
    }
    Phi <- eq$phi0

    if (is.null(eq$cA) || length(eq$cA) == 0) {
      stop("eq$cA is NULL or empty")
    }
    c_A <- eq$cA

    if (is.null(eq$b0) || length(eq$b0) == 0) {
      stop("eq$b0 is NULL or empty")
    }
    b <- eq$b0

    S <- sum(eq$init_S) + sum(eq$init_P)
    D <- sum(eq$init_D)
    A <- sum(eq$init_A + eq$init_U)
    T <- sum(eq$init_T)

    lambda_v_scale <- ((eq$av0 * (c_A*A + eq$cD*D + eq$cT*T))/eq$FOIv_eq)

    params <- list(
      EIR = EIR, ft = ft,
      S0 = S, Ds0 = D, As0 = A, Ts0 = T,
      DR0 = 0, AR0 = 0, TR0 = 0,
      Phi = Phi, b = b,
      m = eq$mv0, Sv0 = eq$init_Sv, Ev_s0 = eq$init_Ev, Iv_s0 = eq$init_Iv,
      Ev_r0 = 0, Iv_r0 = 0,
      a = eq$av0,
      c_A = c_A, c_D = eq$cD, c_T = eq$cT,
      n = eq$delayMos,
      mu = eq$mu0,
      rD = eq$rD,
      rU = eq$rU,
      rA = eq$rA,
      rTs = eq$rT,
      rTR_true = rTR_true,
      lambda_v_scale = lambda_v_scale,
      e = eq$mu0,
      fT = ft,
      ton = ton,
      toff = toff,
      res_time = res_time,
      init_res = init_res
    )

    return(params)
  }, error = function(e) {
    message("Error in phi_eir_rel:")
    message(e$message)
    return(NULL)
  })
}

#' Generate starting parameters for a range of EIR and ft values
#'
#' @param EIRs Numeric vector. The range of Entomological Inoculation Rates.
#' @param fts Numeric vector. The range of treatment rates.
#' @return A data frame of starting parameters for each combination of EIR and ft.
#' @export
starting_params <- function(EIRs = c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100, 200),
                            fts = seq(0.1, 0.9, 0.1)) {
  pars <- tidyr::expand_grid("EIR" = EIRs, "ft" = fts)
  pars$n <- seq_along(pars$EIR)

  starting_params <- lapply(split(pars, pars$n),
                            function(x) {
                              phi_eir_rel(x$EIR, x$ft)
                            }) %>% dplyr::bind_rows()

  return(starting_params)
}
