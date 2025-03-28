# calculate nagelkerke R2
nagelkerke_r2 <- function(fit) {
  tidy_fit <- broom::glance(fit)
  n <- tidy_fit$nobs
  dev <- tidy_fit$deviance
  nul <- tidy_fit$null.deviance
  r2 <- (1 - exp((dev - nul) / n)) / (1 - exp(-nul / n))
  return(r2)
}

# pool vector of R2 values
pool_r2 <- function(r2_values, n_obs) {
  r <- sqrt(r2_values)
  z <- 0.5 * log((r + 1) / (1 - r))
  df <- 1 / (n_obs - 3)
  pooling <- mice::pool.scalar(z, df)
  out <- array(((exp(2 * pooling$qbar) - 1) / (1 + exp(2 * pooling$qbar)))^2,
               dim = c(1, 4)
  )
  dimnames(out) <- list("R^2", c("est", "lo 95", "hi 95", "fmi"))
  out[, 2] <- ((exp(2 * (pooling$qbar - 1.96 * sqrt(pooling$t))) - 1) / (1 + exp(2 * (pooling$qbar - 1.96 * sqrt(pooling$t)))))^2
  out[, 3] <- ((exp(2 * (pooling$qbar + 1.96 * sqrt(pooling$t))) - 1) / (1 + exp(2 * (pooling$qbar + 1.96 * sqrt(pooling$t)))))^2
  out[, 4] <- pooling$fmi
  return(out)
}

# pool nagelkerke R2 values
pool_nagelkerke_r2 <- function(mira) {
  r2_values <- sapply(mira$analyses, nagelkerke_r2)
  n <- broom::glance(mira$analyses[[1]])$nobs
  pooled <- pool_r2(r2_values, n_obs = n)
  return(pooled)
}

