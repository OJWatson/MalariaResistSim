#' Calculate EIR
#'
#' This function calculates the Entomological Inoculation Rate (EIR) based on the given parameters.
#'
#' @param m Ratio of vectors to humans
#' @param a Biting rate (bites on humans per mosquito per day)
#' @param Iv Proportion of infectious mosquito population
#' @return EIR The Entomological Inoculation Rate
#' @export
#' @examples
#' EIR <- calculate_eir(m = 0.7, a = 0.3, Iv = 0.05)
#' print(EIR)
calculate_eir <- function(m, a, Iv) {
  EIR <- m * a * Iv * 365
  return(EIR)
}
