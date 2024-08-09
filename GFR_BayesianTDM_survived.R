require(R2jags)
require(matrixStats) ## For rowVars function
devtools::install_github("rpruim/CalvinBayes", force = TRUE) #to get posterior samples
require(CalvinBayes)

# Load the CKD data
load("ckd.rdata")


# Filter data for patients without graft failure (failure == 0)
success_data <- subset(ckd, failure == 0)

# Extract variables
id <- success_data$id
years <- success_data$years
proteinuria <- success_data$proteinuria
gfr <- success_data$gfr
hematocrit <- success_data$hematocrit
weight <- success_data$weight
age <- success_data$age

# B-spline basis function for time-dynamic model
tpower <- function(x, t, p) {
  return((x - t)^p * (x > t))
}

bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 2) {
  dx <- (xr - xl) / nseg
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  return(B)
}

nknots <- 10
thi <- max(years)
tlo <- min(years)
nseg <- nknots - 2
deg <- 2

bases <- bbase(years, tlo, thi, nseg, deg)
Penalty <- crossprod(diff(diag(ncol(bases)), diff = 2)) + 1e-06 * diag(ncol(bases))

# Data preparation for JAGS
X <- cbind(rep(1, length(gfr)), proteinuria, hematocrit, weight, age)
nX <- dim(X)[2]
N <- length(gfr)

jagsdata <- list(y = gfr, X = X, nX = nX, N = N, nknots = nknots, bases = bases, Penalty_bases = Penalty)

# JAGS model
pspline_jags <- function() {
  for (k in 1:nX) {
    beta[1:N, k] <- bases[, 1:nknots] %*% alpha[k, 1:nknots]
  }
  
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta[i, 1:nX] %*% X[i, ]
  }
  
  for (i in 1:nX) {
    alpha[i, 1:nknots] ~ dmnorm(rep(0, nknots), tau.alpha[i] * Penalty_bases[, ])
    tau.alpha[i] ~ dgamma(1, .005)
  }
  
  tau ~ dgamma(1, .005)
  sigma <- 1 / tau
}

# JAGS fit
params <- c("alpha", "sigma")

fit_pspline <- jags(data = jagsdata, parameters.to.save = params, model.file = pspline_jags,
                    n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 5, DIC = F)

traceplot(fit_pspline)
print(fit_pspline)

# Accessing the simulated values correctly
alpha_samples <- fit_pspline$BUGSoutput$sims.list$alpha

# Plotting the results
gridPoints <- seq(min(years), max(years), len = 20)
bases_grid <- bbase(gridPoints, tlo, thi, nseg, deg)

alpha_0 <- colMeans(alpha_samples[, 1, 1:nknots])
alpha_1 <- colMeans(alpha_samples[, 2, 1:nknots])
alpha_2 <- colMeans(alpha_samples[, 3, 1:nknots])
alpha_3 <- colMeans(alpha_samples[, 4, 1:nknots])
alpha_4 <- colMeans(alpha_samples[, 5, 1:nknots])

beta_0 <- alpha_0 %*% t(bases_grid)
beta_1 <- alpha_1 %*% t(bases_grid)
beta_2 <- alpha_2 %*% t(bases_grid)
beta_3 <- alpha_3 %*% t(bases_grid)
beta_4 <- alpha_4 %*% t(bases_grid)

par(mfrow = c(2, 3))
plot(gridPoints, beta_0, type = "l", main = "beta_0 (Intercept)")


plot(gridPoints, beta_1, type = "l", main = "beta_1 (Proteinuria)")
plot(gridPoints, beta_2, type = "l", main = "beta_2 (Hematocrit)")
plot(gridPoints, beta_3, type = "l", main = "beta_3 (Weight)")
plot(gridPoints, beta_4, type = "l", main = "beta_4 (Age)")


#---------------------------PLOT--------------------
# Function to calculate credible intervals
credible_interval <- function(samples, prob = 0.95) {
  alpha <- (1 - prob) / 2
  lower <- apply(samples, 2, quantile, probs = alpha)
  upper <- apply(samples, 2, quantile, probs = 1 - alpha)
  list(lower = lower, upper = upper)
}

alpha_0 <- colMeans(alpha_samples[, 1, 1:nknots])
alpha_1 <- colMeans(alpha_samples[, 2, 1:nknots])
alpha_2 <- colMeans(alpha_samples[, 3, 1:nknots])
alpha_3 <- colMeans(alpha_samples[, 4, 1:nknots])
alpha_4 <- colMeans(alpha_samples[, 5, 1:nknots])

ci_alpha_0 <- credible_interval(alpha_samples[, 1, 1:nknots])
ci_alpha_1 <- credible_interval(alpha_samples[, 2, 1:nknots])
ci_alpha_2 <- credible_interval(alpha_samples[, 3, 1:nknots])
ci_alpha_3 <- credible_interval(alpha_samples[, 4, 1:nknots])
ci_alpha_4 <- credible_interval(alpha_samples[, 5, 1:nknots])

beta_0 <- alpha_0 %*% t(bases_grid)
beta_1 <- alpha_1 %*% t(bases_grid)
beta_2 <- alpha_2 %*% t(bases_grid)
beta_3 <- alpha_3 %*% t(bases_grid)
beta_4 <- alpha_4 %*% t(bases_grid)

ci_beta_0_lower <- ci_alpha_0$lower %*% t(bases_grid)
ci_beta_0_upper <- ci_alpha_0$upper %*% t(bases_grid)

ci_beta_1_lower <- ci_alpha_1$lower %*% t(bases_grid)
ci_beta_1_upper <- ci_alpha_1$upper %*% t(bases_grid)

ci_beta_2_lower <- ci_alpha_2$lower %*% t(bases_grid)
ci_beta_2_upper <- ci_alpha_2$upper %*% t(bases_grid)

ci_beta_3_lower <- ci_alpha_3$lower %*% t(bases_grid)
ci_beta_3_upper <- ci_alpha_3$upper %*% t(bases_grid)

ci_beta_4_lower <- ci_alpha_4$lower %*% t(bases_grid)
ci_beta_4_upper <- ci_alpha_4$upper %*% t(bases_grid)

# Plotting
par(mfrow = c(2, 3))

plot(gridPoints, beta_0, type = "l", main = "beta_0 (Intercept)", ylim = range(ci_beta_0_lower, ci_beta_0_upper))
lines(gridPoints, ci_beta_0_lower, lty = 2)
lines(gridPoints, ci_beta_0_upper, lty = 2)

plot(gridPoints, beta_1, type = "l", main = "beta_1 (Proteinuria)", ylim = range(ci_beta_1_lower, ci_beta_1_upper))
lines(gridPoints, ci_beta_1_lower, lty = 2)
lines(gridPoints, ci_beta_1_upper, lty = 2)

plot(gridPoints, beta_2, type = "l", main = "beta_2 (Hematocrit)", ylim = range(ci_beta_2_lower, ci_beta_2_upper))
lines(gridPoints, ci_beta_2_lower, lty = 2)
lines(gridPoints, ci_beta_2_upper, lty = 2)

plot(gridPoints, beta_3, type = "l", main = "beta_3 (Weight)", ylim = range(ci_beta_3_lower, ci_beta_3_upper))
lines(gridPoints, ci_beta_3_lower, lty = 2)
lines(gridPoints, ci_beta_3_upper, lty = 2)

plot(gridPoints, beta_4, type = "l", main = "beta_4 (Age)", ylim = range(ci_beta_4_lower, ci_beta_4_upper))
lines(gridPoints, ci_beta_4_lower, lty = 2)
lines(gridPoints, ci_beta_4_upper, lty = 2)

# ___________ PLOT 2 _______________________

# Define the grid points for plotting
gridPoints <- seq(min(years), max(years), len = 20)
bases <- bbase(gridPoints, tlo, thi, nseg, deg) ## basis function. 

# Extracting the alpha values for each predictor
alpha_0 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 1, 1:nknots])
alpha_1 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 2, 1:nknots])
alpha_2 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 3, 1:nknots])
alpha_3 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 4, 1:nknots])
alpha_4 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 5, 1:nknots])

# Calculating beta values
beta_0 <- alpha_0 %*% t(bases)
beta_1 <- alpha_1 %*% t(bases)
beta_2 <- alpha_2 %*% t(bases)
beta_3 <- alpha_3 %*% t(bases)
beta_4 <- alpha_4 %*% t(bases)

# Standard deviations for credible intervals
std_beta_0 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 1, 1:nknots]) %*% t(bases))
std_beta_1 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 2, 1:nknots]) %*% t(bases))
std_beta_2 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 3, 1:nknots]) %*% t(bases))
std_beta_3 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 4, 1:nknots]) %*% t(bases))
std_beta_4 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 5, 1:nknots]) %*% t(bases))

# Plotting the results
par(mfrow = c(2, 3))
plot_names <- c("beta_0 (Intercept)", "beta_1 (Proteinuria)", "beta_2 (Hematocrit)", "beta_3 (Weight)", "beta_4 (Age)")

beta_list <- list(beta_0, beta_1, beta_2, beta_3, beta_4)
std_beta_list <- list(std_beta_0, std_beta_1, std_beta_2, std_beta_3, std_beta_4)

for (i in 1:5) {
  plot(gridPoints, beta_list[[i]], type = "l", lwd = 4, ylab = plot_names[i], xlab = "Time (years)", xaxs = "i", yaxs = "i", ylim = range(c(beta_list[[i]] - 1.96 * std_beta_list[[i]], beta_list[[i]] + 1.96 * std_beta_list[[i]])))
  polygon(c(rev(gridPoints), gridPoints), c(rev(beta_list[[i]] - 1.96 * std_beta_list[[i]]), beta_list[[i]] + 1.96 * std_beta_list[[i]]), col = "cornflowerblue", border = NA)
  lines(gridPoints, beta_list[[i]], lty = 1, col = "black", lwd = 4)
  abline(h = 0, lwd = 2)
}

# Remove the empty plot in the 6th subplot
plot.new()

# ------------- PLOT 3 --------------------
# Define the grid points for plotting
gridPoints <- seq(min(years), max(years), len = 20)
bases <- bbase(gridPoints, tlo, thi, nseg, deg) ## basis function.

# Extracting the alpha values for each predictor
alpha_0 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 1, 1:nknots])
alpha_1 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 2, 1:nknots])
alpha_2 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 3, 1:nknots])
alpha_3 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 4, 1:nknots])
alpha_4 <- colMeans(fit_pspline$BUGSoutput$sims.list$alpha[, 5, 1:nknots])

# Calculating beta values
beta_0 <- alpha_0 %*% t(bases)
beta_1 <- alpha_1 %*% t(bases)
beta_2 <- alpha_2 %*% t(bases)
beta_3 <- alpha_3 %*% t(bases)
beta_4 <- alpha_4 %*% t(bases)

# Standard deviations for credible intervals
std_beta_0 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 1, 1:nknots]) %*% t(bases))
std_beta_1 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 2, 1:nknots]) %*% t(bases))
std_beta_2 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 3, 1:nknots]) %*% t(bases))
std_beta_3 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 4, 1:nknots]) %*% t(bases))
std_beta_4 <- sqrt(colVars(fit_pspline$BUGSoutput$sims.list$alpha[, 5, 1:nknots]) %*% t(bases))

# Define y-axis limits for each beta
ylims <- list(
  c(-52, 10), # ylim for beta_0 (Intercept)
  c(-14, 7.5), # ylim for beta_1 (Proteinuria)
  c(-0.2, 1.7), # ylim for beta_2 (Hematocrit)
  c(-0.2, 1), # ylim for beta_3 (Weight)
  c(-0.7, 0.3) # ylim for beta_4 (Age)
)

# Plotting the results
par(mfrow = c(2, 3))
covariate_names <- c("Intercept", "Proteinuria", "Hematocrit", "Weight", "Age")

beta_list <- list(beta_0, beta_1, beta_2, beta_3, beta_4)
std_beta_list <- list(std_beta_0, std_beta_1, std_beta_2, std_beta_3, std_beta_4)

for (i in 1:5) {
  plot(gridPoints, beta_list[[i]], type = "l", lwd = 2, 
       ylab = expression(hat(beta)(t)), xlab = "Years since the transplant", 
       xaxs = "i", yaxs = "i", ylim = ylims[[i]], main = covariate_names[i])
  polygon(c(rev(gridPoints), gridPoints), 
          c(rev(beta_list[[i]] - 1.96 * std_beta_list[[i]]), beta_list[[i]] + 1.96 * std_beta_list[[i]]), 
          col = "cornflowerblue", border = NA)
  lines(gridPoints, beta_list[[i]], lty = 1, col = "black", lwd = 2)
  abline(h = 0, lwd = 1)
}

# Remove the empty plot in the 6th subplot
plot.new()

