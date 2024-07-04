# AMRSpreadModel

`AMRSpreadModel` is an R package designed to simulate the spread of malaria with a focus on antimalarial resistance. This package allows users to model the impact of different treatment rates, infection rates, and other parameters on malaria prevalence and resistance.

## Installation

You can install the development version of `AMRSpreadModel` from GitHub:

### Install the devtools package if you haven't already

```{r}
install.packages("devtools")
```

### Install AMRSpreadModel from GitHub

```{r}
devtools::install_github("ChaokunHong/MalariaResistSim")
```

## Usage

### Creating a Malaria Model

To create a malaria model, use the create_model function. You need to specify various initial conditions and parameters for the model.

```{r}
library(AMRSpreadModel)

params <- list(
  S0 = 0.9, Ds0 = 0.05, As0 = 0.025, Ts0 = 0.0125,
  DR0 = 0.005, AR0 = 0.0025, TR0 = 0.001,
  a = 0.3, b_h = 0.1, m = 0.7, Phi = 0.7,
  fT = 0.2, rD = 0.05, rA = 0.03, rTs = 0.02,
  rTR = 0.01, Sv0 = 0.9, Ev_s0 = 0, Iv_s0 = 0,
  Ev_r0 = 0, Iv_r0 = 0,
  e = 0.1, mu = 0.1, n = 10,
  Lambda_v_s = 0.005, Lambda_v_r = 0.003,
  EIR_s = NULL, EIR_R = NULL
)

results <- create_model(params)

# Run the model for 500 time steps
model_output <- results$model$run(0:500)
print(model_output)
print(results$FOI)
print(results$EIR_s)
print(results$EIR_R)
print(results$Prevalence)
print(results$Prevalence_R)
```

### Running the Model for Ranges of Parameters

You can also run the model over a range of treatment rates (fT) and entomological inoculation rates (EIR_s and EIR_R).

```{r}
fT_range <- seq(0.1, 0.9, by = 0.1)
EIR_s_range <- seq(10, 200, by = 10)
EIR_R_range <- rep(0, length(EIR_s_range))

params <- list(
  S0 = 0.9, Ds0 = 0.05, As0 = 0.025, Ts0 = 0.0125,
  DR0 = 0.005, AR0 = 0.0025, TR0 = 0.001,
  b_h = 0.1, Phi = 0.7,
  rD = 0.05, rA = 0.03, rTs = 0.02, rTR = 0.01,
  Sv0 = 0.9, Ev_s0 = 0, Iv_s0 = 0, Ev_r0 = 0, Iv_r0 = 0,
  e = 0.1, mu = 0.1, n = 10, Lambda_v_s = 0.005, Lambda_v_r = 0
)

results <- run_model_for_ranges(fT_range, EIR_s_range, EIR_R_range, params)
print(results)
```

## New Features

### Adding Antimalarial Resistance

To simulate the appearance of antimalarial resistance at a specific time, you can specify the res_time parameter along with initial conditions for resistant compartments (DR0, AR0, TR0). The following example demonstrates how to set these parameters:

```{r}
params <- list(
  S0 = 0.9, Ds0 = 0.05, As0 = 0.02, Ts0 = 0.03,
  DR0 = 0, AR0 = 0, TR0 = 0,
  Sv0 = 0.9, Ev_s0 = 0.08, Iv_s0 = 0.02, Ev_r0 = 0, Iv_r0 = 0,
  m = 1.271483, a = 0.3, b = 0.5876259, Phi = 0.7, fT = 0.1,
  rD = 0.2, rA = 0.005, rTs = 0.2, rTR_true = 0.01, e = 0.132,
  mu = 0.132, n = 10, c_A = 0.05, c_D = 0.06, c_T = 0.02,
  ton = 10000, toff = 20000, res_time = 100, res_start = 0.5
)

results <- create_model(params)
```
In this example, resistance compartments (DR, AR, TR) will remain zero until res_time is reached.

### Ensuring Total Population Normalization

The model automatically normalizes the total population to ensure it remains constant throughout the simulation. This is achieved by calculating a correction factor applied to each compartment at every time step.


### Run the model for 500 time steps
model_output <- results$model$run(0:500)
print(model_output)


## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## License

This package is licensed under the MIT License.
