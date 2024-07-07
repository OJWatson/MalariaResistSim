#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom tidyr expand_grid
#' @importFrom stats weighted.mean
NULL

#' @import ICDMM
NULL

#' Create equilibrium initial conditions
#'
#' @param par A list of parameters
#' @return A data frame of equilibrium conditions
#' @export
equilibrium_init_create <- function(par) {
  ## EIR
  EIRY_eq <- par$EIR  # initial annual EIR
  EIRd_eq <- EIR_eq <- EIRY_eq/365

  # FOI
  FOI_eq <- EIR_eq * par$b

  # FOI to T and D
  aT <- FOI_eq * par$phi * par$ft/par$rT_S
  aD <- FOI_eq * par$phi * (1 - par$ft)/par$rD

  Z_eq <- rep(NA, 3)
  Z_eq[1] <- 1/(1 + aT + aD)
  Z_eq[2] <- aT * Z_eq[1]
  Z_eq[3] <- aD * Z_eq[1]

  Y_eq <- Z_eq[1]
  T_eq <- Z_eq[2]
  D_eq <- Z_eq[3]

  betaS <- FOI_eq
  betaA <- FOI_eq * par$phi + par$rA

  A_eq <- (FOI_eq * (1 - par$phi) * Y_eq + par$rD * D_eq)/(betaA + FOI_eq * (1 - par$phi))
  S_eq <- Y_eq - A_eq

  FOIv_eq <- par$a * (par$cT*T_eq + par$cD*D_eq + par$cA*A_eq)

  # mosquito states
  Iv_eq <- FOIv_eq * exp(-par$mu * par$n)/(FOIv_eq + par$mu)
  Sv_eq <- par$mu * Iv_eq/(FOIv_eq * exp(-par$mu * par$n))
  Ev_eq <- 1 - Sv_eq - Iv_eq

  # mosquito density needed to give this EIR
  mv0 <- EIRd_eq/(Iv_eq * par$a)

  ## collate init
  list(
    EIR = par$EIR, ft = par$ft,
    S = S_eq, D = D_eq, A = A_eq, T = T_eq,
    phi = par$phi, b = par$b,
    m = mv0, Sv = Sv_eq, Ev = Ev_eq, Iv = Iv_eq, a = par$a,
    cA = par$cA, cD = par$cD, cT = par$cT,
    n = par$n,
    mu = par$mu,
    rD = par$rD,
    rA = 1/(1/par$rA +1/par$rU),
    rT_S = par$rT_S,
    rT_R = par$rT_R
  ) %>%
    as.data.frame()
}

#' Generate equilibrium parameters based on EIR and ft
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
phi_eir_rel <- function(EIR, ft, ton = 5000, toff = 50000, init_res = 0.01, res_time = 3000, rTR_true = 0.1) {
  mpl <- ICDMM::model_param_list_create(rho=0, rA = 1/(250), rU = Inf, rP = Inf, sigma2 = 0)
  eq <- ICDMM::equilibrium_init_create(
    age_vector=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
    EIR=EIR, ft=ft,
    model_param_list = mpl, het_brackets=1,
    country = NULL,
    admin_unit = NULL)

  # Safe function to calculate weighted mean
  safe_weighted_mean <- function(x, w) {
    if (is.null(x) || length(x) == 0) {
      return(NA)
    }
    if (is.null(dim(x))) {
      if (length(x) != length(w)) {
        return(mean(x, na.rm = TRUE))
      }
      return(weighted.mean(x, w, na.rm = TRUE))
    }
    if (nrow(x) != length(w)) {
      return(mean(x, na.rm = TRUE))
    }
    return(weighted.mean(rowMeans(x, na.rm = TRUE), w, na.rm = TRUE))
  }

  phi <- safe_weighted_mean(eq$phi_eq, eq$het_wt)
  if (is.na(phi)) phi <- mean(eq$phi0, na.rm = TRUE)

  c_A <- safe_weighted_mean(eq$cA_eq, eq$het_wt)
  if (is.na(c_A)) c_A <- mean(eq$cA, na.rm = TRUE)

  b <- safe_weighted_mean(eq$b0 * ((1 - eq$b1)/(1 + (eq$init_IB/eq$IB0)^eq$kB) + eq$b1), eq$den)
  if (is.na(b)) b <- mean(eq$b0, na.rm = TRUE)

  S <- sum(eq$init_S) + sum(eq$init_P)
  D <- sum(eq$init_D)
  A <- sum(eq$init_A + eq$init_U)
  T <- sum(eq$init_T)

  lambda_v_scale <- ((eq$av0 * (c_A*A + eq$cD*D + eq$cT*T))/eq$FOIv_eq)

  par <- list(
    EIR = EIR, ft = ft,
    S = S, D = D, A = A, T = T, phi = phi, b = b,
    m = eq$mv0, Sv = eq$init_Sv, Ev = eq$init_Ev, Iv = eq$init_Iv, a = eq$av0,
    cA = c_A, cD = mean(eq$cD, na.rm = TRUE), cT = mean(eq$cT, na.rm = TRUE),
    n = eq$delayMos,
    mu = eq$mu0,
    rD = eq$rD,
    rU = eq$rU,
    rA = eq$rA,
    rT_S = eq$rT,
    rT_R = rTR_true,
    lambda_v_scale = lambda_v_scale,
    ton = ton,
    toff = toff,
    res_time = res_time,
    init_res = init_res,
    rTR_true = rTR_true
  )

  equilibrium_init_create(par)
}

#' Generate starting parameters for a range of EIR and ft values
#'
#' @param EIRs Numeric vector. The range of Entomological Inoculation Rates.
#' @param fts Numeric vector. The range of treatment rates.
#' @return A data frame of starting parameters for each combination of EIR and ft.
#' @export
starting_params <- function(EIRs = c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100, 200),
                            fts = seq(0.1, 0.9, 0.1)) {
  pars <- expand_grid("EIR" = EIRs, "ft" = fts)
  pars$n <- seq_along(pars$EIR)

  starting_params <- lapply(split(pars, pars$n),
                            function(x) {
                              phi_eir_rel(x$EIR, x$ft)
                            }) %>% bind_rows()

  return(starting_params)
}
