## Equations
deriv(S) <- - S * Lambda_s * (Phi * fT + Phi * (1-fT) + (1-Phi)) -
  S * Lambda_R * (Phi * fT + Phi * (1-fT) + (1-Phi)) +
  Ts * rTs + As * rA + AR * rA + TR * rTR

deriv(Ds) <- S * Lambda_s * Phi * (1-fT) +
  Lambda_s * As * Phi * (1-fT) -
  Ds * rD

deriv(As) <- S * Lambda_s * (1-Phi) +
  Ds * rD -
  Lambda_s * As * Phi * (1-fT) -
  Lambda_s * As * Phi * fT -
  As * rA

deriv(Ts) <- S * Lambda_s * Phi * fT +
  Lambda_s * As * Phi * fT -
  Ts * rTs

deriv(DR) <- S * Lambda_R * Phi * (1-fT) +
  Lambda_R * AR * Phi * (1-fT) -
  DR * rD

deriv(AR) <- S * Lambda_R * (1-Phi) +
  DR * rD -
  Lambda_R * AR * Phi * (1-fT) -
  Lambda_R * AR * Phi * fT -
  AR * rA

deriv(TR) <- S * Lambda_R * Phi * fT +
  Lambda_R * AR * Phi * fT -
  TR * rTR

## Initial conditions
initial(S) <- S0
initial(Ds) <- Ds0
initial(As) <- As0
initial(Ts) <- Ts0
initial(DR) <- DR0
initial(AR) <- AR0
initial(TR) <- TR0

## User-defined parameters
S0 <- user()
Ds0 <- user()
As0 <- user()
Ts0 <- user()
DR0 <- user()
AR0 <- user()
TR0 <- user()
Lambda_s <- user()
Lambda_R <- user()
Phi <- user()
fT <- user()
rD <- user()
rA <- user()
rTs <- user()
rTR <- user()
