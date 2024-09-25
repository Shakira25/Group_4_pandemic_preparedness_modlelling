
#title: "Group 4 SIR Model"
#author: "Group 4 members"
#date: "2024-09-25"

cases<-read.csv("C:/Users/bbrsh/Desktop/Shakira/Shakira/Infectious diseases modelling opportunity Kwame Nkuruma/Hands On/Group_4/MPox_docs/Mpox_data/mpox_2024.csv")
observed_data<-cases$new_cases #Extract the weekly case counts
# Load the necessary libraries
library(deSolve)
library(ggplot2)
library(reshape2)
library(MASS)


# Define base parameters
base_params <- c(
  beta = 0.05,    # Transmission rate
  gamma_I = 0.7,    # Recovery rate for general infectious population
  sigma = 0.05,       # Detection rate
  gamma_Id = 0.14   # Recovery rate for detected population
)

# Define initial state
initial_state <- c(
  S = 105e6,    # Susceptible
  I = 20,     # Infectious
  Id = 2,     # Detected population
  R = 0       # Recovered
)

# Define the time sequence (57 weeks)
times <- seq(0, 57, by = 1)  # Replace with appropriate time frame

# Define the ODE system
mpox_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # ODEs based on the flowchart parameters
    dS <- -beta * S * I
    dI <- beta * S * I - gamma_I * I - sigma * I
    dId <- sigma * I - gamma_Id * Id
    dR <- gamma_I * I + gamma_Id * Id
    
    return(list(c(dS, dI, dId, dR)))
  })
}

# Run the model
output <- ode(y = initial_state, times = times, func = mpox_model, parms = base_params)

# Rename columns for clarity
colnames(output) <- c("time", "Susceptible", "Infectious", "Detected", "Recovered")

# Check the output to ensure it's being generated (Optional)
print(head(output))

# Plot the results
ggplot(data = output, aes(x = time)) +
  geom_line(aes(y = Susceptible, color = "Susceptible"), size = 1) +
  geom_line(aes(y = Infectious, color = "Infectious"), size = 1) +
  geom_line(aes(y = Detected, color = "Detected"), size = 1) +
  geom_line(aes(y = Recovered, color = "Recovered"), size = 1) +
  labs(title = "SIR Model Simulation", x = "Days", y = "Population") +
  scale_color_manual("", 
                     breaks = c("Susceptible", "Infectious", "Detected", "Recovered"),
                     values = c("blue", "red", "brown", "green")) +
  theme_minimal() +
  theme(legend.position = "bottom")


# Define cost function for parameter estimation
cost_function <- function(pars, observed_data) {
  # Create a copy of base_params and update sigma
  params <- base_params
  params["sigma"] <- pars[1]
  
  # Simulate the model
  simulation <- ode(y = initial_state, times = times, func = mpox_model, parms = params)
  
  # Extract the number of detected individuals (Id)
  simulated_Id <- simulation[, "Id"]
  
  # Calculate the sum of squared errors
  return(sum((observed_data - simulated_Id)^2))
}
  
# Define lower and upper bounds for the parameters
lower_bounds <- c(sigma = 0.01)
upper_bounds <- c(sigma = 1)

# Initial guesses for the parameters to be estimated
initial_guesses <- lower_bounds

# Perform parameter estimation using optim with bounds
fit <- optim(par = initial_guesses, fn = cost_function, observed_data = observed_data, 
             method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, hessian = TRUE)

# Extract the estimated parameters
estimated <- fit$par
##names(estimated) <- c("sigma")
print(estimated)

# Compute confidence intervals
hessian_matrix <- fit$hessian
cov_matrix <- solve(hessian_matrix)  # Invert the Hessian matrix to get the covariance matrix
std_errors <- sqrt(diag(cov_matrix))  # Standard errors of the parameters
conf_intervals <- cbind(estimated - 1.96 * std_errors, estimated + 1.96 * std_errors)
print(conf_intervals)

params["sigma"] <- estimated["sigma"]