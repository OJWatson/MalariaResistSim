% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_malaria_model.R
\name{create_malaria_model}
\alias{create_malaria_model}
\title{Create Malaria Model}
\usage{
create_malaria_model(
  S0,
  Ds0,
  As0,
  Ts0,
  DR0,
  AR0,
  TR0,
  Lambda_s,
  Lambda_R,
  Phi,
  fT,
  rD,
  rA,
  rTs,
  rTR,
  Sv0,
  Ev_s0,
  Iv_s0,
  Ev_r0,
  Iv_r0,
  e,
  mu,
  n,
  Lambda_v_s,
  Lambda_v_r
)
}
\arguments{
\item{S0}{Initial proportion of susceptible human population}

\item{Ds0}{Initial proportion of diseased human population (sensitive strain)}

\item{As0}{Initial proportion of asymptomatic human population (sensitive strain)}

\item{Ts0}{Initial proportion of treated human population (sensitive strain)}

\item{DR0}{Initial proportion of diseased human population (resistant strain)}

\item{AR0}{Initial proportion of asymptomatic human population (resistant strain)}

\item{TR0}{Initial proportion of treated human population (resistant strain)}

\item{Lambda_s}{Force of infection (sensitive strain)}

\item{Lambda_R}{Force of infection (resistant strain)}

\item{Phi}{Proportion of symptomatic infections}

\item{fT}{Proportion of symptomatic infections treated}

\item{rD}{Recovery rate for diseased individuals}

\item{rA}{Recovery rate for asymptomatic individuals}

\item{rTs}{Recovery rate for treated individuals (sensitive strain)}

\item{rTR}{Recovery rate for treated individuals (resistant strain)}

\item{Sv0}{Initial proportion of susceptible mosquito population}

\item{Ev_s0}{Initial proportion of exposed mosquito population (sensitive strain)}

\item{Iv_s0}{Initial proportion of infectious mosquito population (sensitive strain)}

\item{Ev_r0}{Initial proportion of exposed mosquito population (resistant strain)}

\item{Iv_r0}{Initial proportion of infectious mosquito population (resistant strain)}

\item{e}{Mosquito emergence rate}

\item{mu}{Mosquito death rate}

\item{n}{Extrinsic incubation period}

\item{Lambda_v_s}{Force of infection from humans to mosquitoes (sensitive strain)}

\item{Lambda_v_r}{Force of infection from humans to mosquitoes (resistant strain)}
}
\value{
An object of class \code{odin_model}.
}
\description{
This function creates a malaria model using user-defined parameters and returns an \code{odin} model object that can be used to run simulations.
}
\examples{
model <- create_malaria_model(S0 = 0.9, Ds0 = 0.05, As0 = 0.025, Ts0 = 0.0125,
                              DR0 = 0.005, AR0 = 0.0025, TR0 = 0.001,
                              Lambda_s = 0.3, Lambda_R = 0.1, Phi = 0.7,
                              fT = 0.2, rD = 0.05, rA = 0.03, rTs = 0.02,
                              rTR = 0.01, Sv0 = 0.9, Ev_s0 = 0, Iv_s0 = 0,
                              Ev_r0 = 0, Iv_r0 = 0,
                              e = 0.1, mu = 0.1, n = 10,
                              Lambda_v_s = 0.005, Lambda_v_r = 0.003)
results <- model$run(0:500)
print(results)
}
