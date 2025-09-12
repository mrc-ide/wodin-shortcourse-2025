## Initial conditions

initial(S) <- N_init - I_init
initial(I) <- I_init
initial(cumulative_infections) <- 0
initial(cumulative_hiv_deaths) <- 0


## Parameter values
N_init <- user(1e6)    # Initial population size
I_init <- user(1)       # Initial number infected
c      <- user(2.4)     # Partner change rate
p      <- user(0.1)     # Per partner HIV transmission probability
sigma  <- user(1/12)  # Mortality rate per person per year due to HIV/AIDS (1 / mean duration in years)
mu     <- user(0.008)   # Crude mortality rate due to causes other than AIDS, scaled to rate per person
alpha  <- user(0.0332)  # Birth rate scaled to per person


## Differential equations

deriv(S) <- B - foi * S - mu * S
deriv(I) <- foi * S - (mu + sigma) * I
deriv(cumulative_infections) <- foi * S
deriv(cumulative_hiv_deaths) <- sigma * I


## Other equations

N <- S + I
B <- alpha * N                # Entry rate, exponentially growing population
# B <- mu * N + sigma * I     # Entry rate, constant population size - useful for model checking

beta <- p * c                 # Population size
foi <- beta * I / N           # Force of infection per capita

R0 <- beta / (mu + sigma)     # Basic reproduction number


## Additional outputs
output(R0)                     <- TRUE
output(prevalence)             <- I / N
output(N)                      <- TRUE
