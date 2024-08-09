# Load necessary libraries
require(R2jags)
require(matrixStats) ## For rowVars function
devtools::install_github("rpruim/CalvinBayes", force = TRUE) # to get posterior samples
require(CalvinBayes)

# Load the CKD data
load("ckd.rdata")

# Filter data for patients with graft failure (failure == 1)
failed_data <- subset(ckd, failure == 1)

# Extract variables
id <- failed_data$id
years <- failed_data$years
proteinuria <- failed_data$proteinuria
gfr <- failed_data$gfr
hematocrit <- failed_data$hematocrit
weight <- failed_data$weight
age <- failed_data$age
gender <- failed_data$gender

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
X <- cbind(rep(1, length(gfr)), proteinuria, hematocrit, weight, age, gender)
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

credible_interval <- function(samples, prob = 0.95) {
  alpha <- (1 - prob) / 2
  lower <- apply(samples, 2, quantile, probs = alpha)
  upper <- apply(samples, 2, quantile, probs = 1 - alpha)
  list(lower = lower, upper = upper)
}

ci_alpha_0 <- credible_interval(alpha_samples[, 1, 1:nknots])
ci_alpha_1 <- credible_interval(alpha_samples[, 2, 1:nknots])
ci_alpha_2 <- credible_interval(alpha_samples[, 3, 1:nknots])
ci_alpha_3 <- credible_interval(alpha_samples[, 4, 1:nknots])
ci_alpha_4 <- credible_interval(alpha_samples[, 5, 1:nknots])
ci_alpha_5 <- credible_interval(alpha_samples[, 6, 1:nknots])

beta_0 <- colMeans(alpha_samples[, 1, 1:nknots]) %*% t(bases_grid)
beta_1 <- colMeans(alpha_samples[, 2, 1:nknots]) %*% t(bases_grid)
beta_2 <- colMeans(alpha_samples[, 3, 1:nknots]) %*% t(bases_grid)
beta_3 <- colMeans(alpha_samples[, 4, 1:nknots]) %*% t(bases_grid)
beta_4 <- colMeans(alpha_samples[, 5, 1:nknots]) %*% t(bases_grid)
beta_5 <- colMeans(alpha_samples[, 6, 1:nknots]) %*% t(bases_grid)

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
ci_beta_5_lower <- ci_alpha_5$lower %*% t(bases_grid)
ci_beta_5_upper <- ci_alpha_5$upper %*% t(bases_grid)

par(mfrow = c(2, 3))

plot(gridPoints, beta_0, type = "l", lwd = 4, ylab = expression(hat(beta)[0](t)), xlab = "Years since transplant", xaxs = "i", yaxs = "i", ylim = range(ci_beta_0_lower, ci_beta_0_upper))
polygon(c(rev(gridPoints), gridPoints), c(rev(ci_beta_0_lower), ci_beta_0_upper), col = "coral", border = NA)
lines(gridPoints, beta_0, lty = 1, col = "black", lwd = 4)
abline(h = 0)
title("Intercept")

plot(gridPoints, beta_1, type = "l", lwd = 4, ylab = expression(hat(beta)[1](t)), xlab = "Years since transplant", xaxs = "i", yaxs = "i", ylim = range(ci_beta_1_lower, ci_beta_1_upper))
polygon(c(rev(gridPoints), gridPoints), c(rev(ci_beta_1_lower), ci_beta_1_upper), col = "coral", border = NA)
lines(gridPoints, beta_1, lty = 1, col = "black", lwd = 4)
abline(h = 0)
title("Proteinuria")

plot(gridPoints, beta_2, type = "l", lwd = 4, ylab = expression(hat(beta)[2](t)), xlab = "Years since transplant", xaxs = "i", yaxs = "i", ylim = range(ci_beta_2_lower, ci_beta_2_upper))
polygon(c(rev(gridPoints), gridPoints), c(rev(ci_beta_2_lower), ci_beta_2_upper), col = "coral", border = NA)
lines(gridPoints, beta_2, lty = 1, col = "black", lwd = 4)
abline(h = 0)
title("Hematocrit")

plot(gridPoints, beta_3, type = "l", lwd = 4, ylab = expression(hat(beta)[3](t)), xlab = "Years since transplant", xaxs = "i", yaxs = "i", ylim = range(ci_beta_3_lower, ci_beta_3_upper))
polygon(c(rev(gridPoints), gridPoints), c(rev(ci_beta_3_lower), ci_beta_3_upper), col = "coral", border = NA)
lines(gridPoints, beta_3, lty = 1, col = "black", lwd = 4)
abline(h = 0)
title("Weight")

plot(gridPoints, beta_4, type = "l", lwd = 4, ylab = expression(hat(beta)[4](t)), xlab = "Years since transplant", xaxs = "i", yaxs = "i", ylim = range(ci_beta_4_lower, ci_beta_4_upper))
polygon(c(rev(gridPoints), gridPoints), c(rev(ci_beta_4_lower), ci_beta_4_upper), col = "coral", border = NA)
lines(gridPoints, beta_4, lty = 1, col = "black", lwd = 4)
abline(h = 0)
title("Age")

plot(gridPoints, beta_5, type = "l", lwd = 4, ylab = expression(hat(beta)[5](t)), xlab = "Years since transplant", xaxs = "i", yaxs = "i", ylim = range(ci_beta_5_lower, ci_beta_5_upper))
polygon(c(rev(gridPoints), gridPoints), c(rev(ci_beta_5_lower), ci_beta_5_upper), col = "coral", border = NA)
lines(gridPoints, beta_5, lty = 1, col = "black", lwd = 4)
abline(h = 0)
title("Gender")
