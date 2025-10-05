library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Set seed
set.seed(1)

# Generate data
sim_data <- function(n = 1000, a = 2, b = 2, sd = 1) {
  X <- rbeta(n, a, b)
  m <- function(x) sin((x/3 + 0.1)^(-1))
  Y <- m(X) + rnorm(n, 0, sd)
  data.frame(X = X, Y = Y)
}

# Create blocks with equal number of points
make_blocks <- function(d, N) {
  d <- d[order(d$X), ]
  q <- quantile(d$X, probs = seq(0, 1, length.out = N + 1))
  d$block <- cut(d$X, breaks = q, include.lowest = TRUE, labels = FALSE)
  return(d)
}

# Fit polynomial models in blocks and compute estimates
get_est <- function(d, N) {
  n <- nrow(d)
  d <- make_blocks(d, N)
  
  preds <- numeric(n)
  d2 <- numeric(n)
  
  for (j in 1:N) {
    bd <- d[d$block == j, ]
    
    bd$X2 <- bd$X^2
    bd$X3 <- bd$X^3
    bd$X4 <- bd$X^4
    
    mod <- lm(Y ~ X + X2 + X3 + X4, data = bd)
    
    idx <- which(d$block == j)
    pdata <- data.frame(
      X = d$X[idx],
      X2 = d$X[idx]^2,
      X3 = d$X[idx]^3,
      X4 = d$X[idx]^4
    )
    
    preds[idx] <- predict(mod, newdata = pdata)
    
    cfs <- coef(mod)
    b2 <- cfs[3]; b3 <- cfs[4]; b4 <- cfs[5]
    d2[idx] <- 2*b2 + 6*b3*d$X[idx] + 12*b4*d$X[idx]^2
  }
  
  t22 <- mean(d2^2, na.rm = TRUE)
  s2 <- sum((d$Y - preds)^2, na.rm = TRUE) / (n - 5*N)
  
  list(t22 = t22, s2 = s2)
}

# Compute optimal h
h_opt <- function(n, s2, t22, L = 1) {
  if (t22 <= 0 || is.na(t22)) return(NA)
  n^(-1/5) * (35 * s2 * L / t22)^(1/5)
}

# Compute RSS for a given block size N
rss <- function(d, N) {
  d <- make_blocks(d, N)
  r <- 0
  
  for (j in 1:N) {
    bd <- d[d$block == j, ]
    bd$X2 <- bd$X^2
    bd$X3 <- bd$X^3
    bd$X4 <- bd$X^4
    
    mod <- lm(Y ~ X + X2 + X3 + X4, data = bd)
    
    idx <- which(d$block == j)
    pdata <- data.frame(
      X = d$X[idx],
      X2 = d$X[idx]^2,
      X3 = d$X[idx]^3,
      X4 = d$X[idx]^4
    )
    
    bp <- predict(mod, newdata = pdata)
    r <- r + sum((d$Y[idx] - bp)^2, na.rm = TRUE)
  }
  return(r)
}

# Compute Mallow's Cp
cp <- function(d, N, N_max) {
  n <- nrow(d)
  r_N <- rss(d, N)
  r_max <- rss(d, N_max)
  cp_val <- (r_N / (r_max / (n - 5 * N_max))) - (n - 10 * N)
  return(cp_val)
}

# Minimizes Mallow's Cp
opt_blocks <- function(d) {
  n <- nrow(d)
  N_max <- max(min(floor(n / 20), 5), 1)
  N_vals <- 1:N_max
  cp_vals <- map_dbl(N_vals, ~ cp(d, ., N_max))
  N_opt <- N_vals[which.min(cp_vals)]
  return(N_opt)
}

# Beta parameters
beta_pars <- list(
  c(0.5, 0.5),  
  c(1, 1),      
  c(2, 2),      
  c(5, 1),    
  c(1, 5)     
)

# 1. Effect of sample size n
n_vals <- c(200, 500, 1000, 2000, 5000,15000)
res_n <- data.frame()

for (n in n_vals) {
  for (par in beta_pars) {
    a <- par[1]
    b <- par[2]
    
    d <- sim_data(n = n, a = a, b = b)
    
    N_opt <- opt_blocks(d)
    est <- get_est(d, N_opt)
    h <- h_opt(n, est$s2, est$t22)
    
    res_n <- rbind(res_n, data.frame(
      n = n,
      a = a,
      b = b,
      dist = paste0("B(", a, ",", b, ")"),
      N_opt = N_opt,
      h_opt = h,
      t22 = est$t22,
      s2 = est$s2
    ))
  }
}

# Plot 1: Effect of sample size
p1 <- ggplot(res_n, aes(x = n, y = h_opt, color = dist, group = dist)) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_log10() +
  labs(title = "Sample Size vs Optimal Bandwidth",
       x = "n",
       y = "h_AMISE",
       color = "Distribution") +
  theme_minimal() +
  theme(legend.position = "bottom")

p1
ggsave("plots1/n_vs_h.png", p1, width = 9, height = 6, dpi = 300)

# 2. Effect of block size N
N_vals <- 1:10
res_N <- data.frame()

for (par in beta_pars) {
  a <- par[1]
  b <- par[2]
  
  d <- sim_data(n = 1000, a = a, b = b)
  
  N_opt <- opt_blocks(d)
  
  for (N in N_vals) {
    est <- get_est(d, N)
    h <- h_opt(1000, est$s2, est$t22)
    
    res_N <- rbind(res_N, data.frame(
      a = a,
      b = b,
      dist = paste0("B(", a, ",", b, ")"),
      N = N,
      h_opt = h,
      t22 = est$t22,
      s2 = est$s2,
      opt = (N == N_opt)
    ))
  }
}

# Plot 2: Effect of block size
p2 <- ggplot(res_N, aes(x = N, y = h_opt, color = dist, group = dist)) +
  geom_point(aes(shape = opt), size = 2) +
  geom_line(alpha = 0.7) +
  scale_shape_manual(values = c(16, 17), guide = "none") +
  labs(title = "Block Size vs Optimal Bandwidth",
       x = "N",
       y = "h_AMISE",
       color = "Distribution") +
  theme_minimal() +
  theme(legend.position = "bottom")

p2
ggsave("plots1/N_vs_h.png", p2, width = 9, height = 6, dpi = 300)

# 3. Effect of Beta distribution param
res_dist <- res_n %>%
  filter(n == 1000) %>%
  group_by(dist) %>%
  slice(1)

p3 <- ggplot(res_dist, aes(x = h_opt, y = reorder(dist, h_opt))) +
  geom_point(aes(color = dist), size = 4) +
  geom_segment(aes(xend = 0, yend = dist, color = dist), linewidth = 1) +
  geom_text(aes(label = paste("N =", N_opt)), 
            hjust = -0.2, size = 3, color = "black") +
  labs(title = "Optimal Bandwidth by Distribution",
       x = "h_AMISE",
       y = "Beta Distribution",
       subtitle = "Ordered by bandwidth value\nN shows optimal block count") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank()) +
  expand_limits(x = max(res_dist$h_opt) * 1.1)

p3
ggsave("plots1/distr_vs_h.png", p3, width = 8, height = 6, dpi = 300)




