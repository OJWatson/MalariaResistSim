## Human Equations
deriv(S) <- -S * Lambda_s * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) -
  if (t >= res_time) (S * Lambda_R * (Phi * fT + Phi * (1 - fT) + (1 - Phi))) else 0 +
  Ts * rTs + As * rA + if (t >= res_time) (AR * rA + TR * rTR) else 0

deriv(Ds) <- S * Lambda_s * Phi * (1 - fT) +
  Lambda_s * As * Phi * (1 - fT) -
  Ds * rD

deriv(As) <- S * Lambda_s * (1 - Phi) +
  Ds * rD -
  Lambda_s * As * Phi * (1 - fT) -
  Lambda_s * As * Phi * fT -
  As * rA

deriv(Ts) <- S * Lambda_s * Phi * fT +
  Lambda_s * As * Phi * fT -
  Ts * rTs

deriv(DR) <- if (t >= res_time) (S * Lambda_R * Phi * (1 - fT) +
                                   Lambda_R * AR * Phi * (1 - fT) -
                                   DR * rD) else - DR * rD

deriv(AR) <- if (t >= res_time) (S * Lambda_R * (1 - Phi) +
                                   DR * rD -
                                   Lambda_R * AR * Phi * (1 - fT) -
                                   Lambda_R * AR * Phi * fT -
                                   AR * rA) else - AR * rA

deriv(TR) <- if (t >= res_time) (S * Lambda_R * Phi * fT +
                                   Lambda_R * AR * Phi * fT -
                                   TR * rTR) else - TR * rTR

## Mosquito Equations
deriv(Sv) <- e - (Lambda_v_s + if (t >= res_time) Lambda_v_r else 0) * Sv - mu * Sv

delayed_Lambda_v_s_Sv <- delay(Lambda_v_s * Sv * exp(-mu * n), n)

deriv(Ev_s) <- Lambda_v_s * Sv - delayed_Lambda_v_s_Sv - mu * Ev_s
deriv(Iv_s) <- delayed_Lambda_v_s_Sv - mu * Iv_s

Lambda_v_r_delayed <- if (t >= res_time) Lambda_v_r else 0
delayed_Lambda_v_r_Sv <- delay(Lambda_v_r_delayed * Sv * exp(-mu * n), n)

deriv(Ev_r) <- if (t >= res_time) (Lambda_v_r * Sv - delayed_Lambda_v_r_Sv - mu * Ev_r) else - mu * Ev_r
deriv(Iv_r) <- if (t >= res_time) (delayed_Lambda_v_r_Sv - mu * Iv_r) else - mu * Iv_r

output(prevalence) <- S + As + Ds + Ts + if (t >= res_time) (AR + DR + TR) else 0
output(prevalence_res) <- if (t >= res_time && (AR + DR + TR + As + Ds + Ts) > 0) ((AR + DR + TR) / (As + Ds + Ts + AR + DR + TR)) else 0

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
m <- user()
a <- user()
b <- user()
Lambda_s <- m * a * b * Iv_s
Lambda_R <- if (t >= res_time) (m * a * b * Iv_r) else 0
Phi <- user()
fT <- user()
rD <- user()
rA <- user()
rTs <- user()
rTR_true <- user()
rTR <- if (t > ton && t < toff) rTR_true else rTs
Sv0 <- user()
Ev_s0 <- user()
Iv_s0 <- user()
Ev_r0 <- user()
Iv_r0 <- user()
e <- user()
mu <- user()
n <- user()
c_A <- user()
c_D <- user()
c_T <- user()
Lambda_v_s <- a * (c_A * As + c_D * Ds + c_T * Ts)
Lambda_v_r <- if (t >= res_time) (a * (c_A * AR + c_D * DR + c_T * TR)) else 0
ton <- user()
toff <- user()
res_time <- user()
res_start <- user()
