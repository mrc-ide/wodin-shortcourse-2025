## Initial conditions

N_init                  <- user(1e6)  # Initial total population size
I_init                  <- user(1)     # Initial number infected
theta                   <- user(0.4)   # Proportion at-risk (in S class)

initial(S_f)                     <- (N_init - I_init) * theta / 2
initial(S_m)                     <- (N_init - I_init) * theta / 2
initial(I_f)                     <- I_init / 2
initial(I_m)                     <- I_init / 2
initial(P_f)                     <- (1 - theta) * (N_init - I_init) / 2
initial(P_m)                     <- (1 - theta) * (N_init - I_init) / 2
initial(cumulative_infections)   <- 0
initial(cumulative_hiv_deaths)   <- 0


## Parameter values

c_m          <- user(11)       # Partner change rate, males
#c_f          <- user(11)       # Partner change rate, females
p_mtof       <- user(0.1)      # Per partner HIV transmission probability, male to female
p_ftom       <- user(0.1)      # Per partner HIV transmission probability, female to male
sigma        <- user(1/12)   # Mortality rate per person per year due to HIV/AIDS (1 / mean duration in years)
mu           <- user(0.008)    # Crude mortality rate due to causes other than AIDS, scaled to rate per person
alpha        <- user(0.0332)   # Birth rate scaled to per person

## Other equations

N_f      <- S_f + I_f + P_f
N_m      <- S_m + I_m + P_m
N        <- N_f + N_m

c_f <- c_m * (S_m + I_m) / (S_f + I_f)

beta_mtof <- p_mtof * c_f
beta_ftom <- p_ftom * c_m

foi_mtof <- beta_mtof * I_m / N_m  # Force of infection: male → female
foi_ftom <- beta_ftom * I_f / N_f  # Force of infection: female → male

B <- alpha * N             # Entry rate, exponentially growing population
# B <- mu * N + sigma * (I_f + I_m)  # Entry rate, constant population size - useful for model checking

## Differential equations

deriv(S_f) <- theta * B / 2 - foi_mtof * S_f - mu * S_f
deriv(S_m) <- theta * B / 2 - foi_ftom * S_m - mu * S_m

deriv(I_f) <- foi_mtof * S_f - (mu + sigma) * I_f
deriv(I_m) <- foi_ftom * S_m - (mu + sigma) * I_m

deriv(P_f) <- (1 - theta) * B / 2 - mu * P_f
deriv(P_m) <- (1 - theta) * B / 2 - mu * P_m

deriv(cumulative_infections)     <- foi_mtof * S_f + foi_ftom * S_m
deriv(cumulative_hiv_deaths)     <- sigma * (I_f + I_m)

## Additional outputs

output(prevalence_f)       <- I_f / N_f
output(prevalence_m)       <- I_m / N_m
output(prevalence)     <- (I_f + I_m) / N

output(partnerships_f)     <- c_f * (S_f + I_f)
output(partnerships_m)     <- c_m * (S_m + I_m)

output(N_m) <- TRUE
output(N_f) <- TRUE
output(N) <- TRUE
output(S) <- S_m + S_f
output(I) <- I_m + I_f
