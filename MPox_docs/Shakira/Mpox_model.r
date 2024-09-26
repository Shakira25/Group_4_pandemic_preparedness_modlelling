
#title: "Group 4 SIR Model"
#author: "Group 4 members"
#date: "2024-09-25"

cases<-read.csv("mpox_2024.csv")
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
times <- seq(0, 56, by = 1)  # Replace with appropriate time frame

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

# Simulate again with new parameters
simulation <- ode(y = initial_state, times = times, func = mpox_model, parms = params)
simulation_df <- as.data.frame(simulation)

# Define start date for the timeline
start_date <- as.Date("2023-08-20")  # Example start date

# Create a date vector for the time points (based on weeks, as the time vector is weekly)
date_vector <- seq.Date(from = start_date, by = "week", length.out = length(times))
project_time <- seq(max(times),max(times)+52*1) # 1 year more (52 Weeks)
projection <- ode(y = initial_state, times = project_time, func = mpox_model, parms = params)
projection_df <- as.data.frame(projection)

# Create a date vector for the projection period
projection_date_vector <- seq.Date(from = max(date_vector), by = "week", length.out = length(project_time))

# Plot the observed and simulated results with dates on the x-axis
print(
  plot(date_vector, simulation_df$Ihd, type = "l", col = "red", ylab = "Detected Infected Humans (Id)", xlab = "Date", xaxt = "n")
)
# Overlay observed data (make sure 'observed_data' has the correct length)
points(date_vector, observed_data, col = "blue", pch = 16)

# Format the x-axis to show dates
axis.Date(1, at = date_vector, format = "%b %Y")  # Custom format: "Month Year"
legend("topleft", legend = c("Simulated Id", "Observed Data"), col = c("red", "blue"), lty = c(1, NA), pch=c(NA,16))

# Plot the projection with dates on the x-axis
print(
  plot(projection_date_vector, projection_df$Id, type = "l", col = "red", ylab = "Projected detected Infected Humans (Id)", xlab = "Date", xaxt = "n")
)
axis.Date(1, at = projection_date_vector, format = "%b %Y")


#################### Projection ( 2023 to 2025)
#################### 

project_time <- seq(0,max(times)+52*1) # 1 year more (52 Weeks)
projection <- ode(y = initial_state, times = project_time, func = mpox_model, parms = params)
projection_df <- as.data.frame(projection)

# Create a date vector for the projection period
projection_date_vector <- seq.Date(from = start_date, by = "week", length.out = length(project_time))


#Plot the projection with dates on the x-axis
print(
  plot(projection_date_vector, projection_df$Id, type = "l", col = "red", ylab = "Projected Detected Infected Humans (Id)", xlab = "Date", xaxt = "n")
)
# Overlay observed data on the projection graph (observed data only goes to the current time)
points(date_vector, observed_data, col = "blue", pch = 16)  # Observed data only goes up to current time
axis.Date(1, at = projection_date_vector, format = "%b %Y")
legend("topleft", legend = c("Projected MPOX Id", "Observed Data"), col = c("red", "blue"), lty = c(1, NA), pch=c(NA,16))