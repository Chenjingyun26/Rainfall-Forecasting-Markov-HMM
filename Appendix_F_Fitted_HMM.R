# Appendix F: Fitted two-state Bernoulli HMM for Wilmington
# Base R only. Training: 2016-2023. Held-out evaluation: 2024-2025.
path <- file.choose()
d <- read.csv(path, stringsAsFactors = FALSE)
stopifnot(all(c('DATE', 'PRCP', 'STATION') %in% names(d)))
stopifnot(all(d$STATION == 'USW00013748'))
d$DATE <- as.Date(d$DATE)
d <- d[order(d$DATE), ]
stopifnot(identical(as.character(d$DATE), as.character(seq(as.Date('2016-01-01'), as.Date('2025-12-31'), by='day'))))
stopifnot(all(is.finite(d$PRCP)), all(d$PRCP >= 0))
y <- as.integer(d$PRCP > 0)
n <- sum(d$DATE < as.Date('2024-01-01'))
y_train <- y[seq_len(n)]
# theta = (P(1->2), P(2->1), P(rain|1), P(rain|2)).
# Use the stationary initial distribution rather than estimating a fifth parameter.
forward <- function(theta, obs, keep = FALSE) {
  a <- theta[1]; b <- theta[2]; e1 <- theta[3]; e2 <- theta[4]
  p <- a/(a+b); ll <- 0
  if (keep) { forecast <- filtered <- numeric(length(obs)) }
  for (i in seq_along(obs)) {
    if (i > 1L) p <- (1-p)*a + p*(1-b)
    if (keep) forecast[i] <- (1-p)*e1 + p*e2
    if (obs[i] == 1L) { w1 <- (1-p)*e1; w2 <- p*e2 }
    else { w1 <- (1-p)*(1-e1); w2 <- p*(1-e2) }
    s <- w1+w2
    ll <- ll + log(s); p <- w2/s
    if (keep) filtered[i] <- p
  }
  if (keep) return(list(nll=-ll, forecast=forecast, filtered=filtered))
  -ll
}
starts <- rbind(c(.15,.30,.10,.80), c(.05,.10,.15,.60),
                c(.30,.30,.05,.90), c(.40,.15,.10,.50),
                c(.10,.40,.20,.80), c(.50,.50,.20,.70),
                c(.20,.20,.25,.45), c(.70,.30,.05,.60))
eps <- 1e-6
fits <- lapply(seq_len(nrow(starts)), function(i) {
  optim(starts[i,], forward, obs=y_train, method='L-BFGS-B',
        lower=rep(eps,4), upper=rep(1-eps,4),
        control=list(maxit=2000, factr=1e5, pgtol=1e-7, ndeps=rep(1e-6,4)))
})
fit_checks <- data.frame(start=seq_along(fits),
  logLik=sapply(fits,function(f)-f$value),
  convergence=sapply(fits,function(f)f$convergence))
ok <- which(fit_checks$convergence==0 & is.finite(fit_checks$logLik))
if (!length(ok)) stop('No successful optimization; inspect fit diagnostics.')
best <- fits[[ok[which.max(fit_checks$logLik[ok])]]]
theta <- best$par
# Relabel so state 1 is the lower-rain-probability state.
if(theta[3] > theta[4]) theta <- theta[c(2,1,4,3)]
a <- theta[1]; b <- theta[2]
A <- matrix(c(1-a,a,b,1-b),2,byrow=TRUE,
            dimnames=list(c('State_1','State_2'),c('State_1','State_2')))
params <- data.frame(state=c('State_1','State_2'),
  initial_probability=c(b/(a+b),a/(a+b)),
  to_state_1=A[,1],to_state_2=A[,2],rain_emission=theta[3:4])
# Forecast on target day uses observations through previous day only.
f <- forward(theta,y,keep=TRUE)
idx <- seq.int(n+1L,length(y))
from <- y_train[-n]; to <- y_train[-1]
q <- sapply(0:1,function(k) mean(to[from==k]))
clim <- mean(y_train)
pred <- data.frame(DATE=d$DATE[idx], Actual=y[idx],
  HMM=f$forecast[idx], Binary_Markov=q[y[idx-1L]+1L],
  Historical_Frequency=rep(clim,length(idx)))
loss <- (as.matrix(pred[,3:5])-pred$Actual)^2
score <- colMeans(loss)
# Paired moving-block bootstrap for fixed fitted models, 7 consecutive days.
set.seed(2026); B <- 2000L; L <- 7L; N <- nrow(loss)
boot <- replicate(B, {
  st <- sample.int(N-L+1L,ceiling(N/L),replace=TRUE)
  ii <- as.vector(outer(0:(L-1L),st,'+'))[seq_len(N)]
  colMeans(loss[ii,,drop=FALSE])
})
ci <- t(apply(boot,1,quantile,probs=c(.025,.975)))
scores <- data.frame(model=names(score),N=N,Brier=as.numeric(score),
  lower95=ci[,1],upper95=ci[,2],row.names=NULL)
delta <- rbind(boot[1,]-boot[2,],boot[1,]-boot[3,])
dci <- t(apply(delta,1,quantile,probs=c(.025,.975)))
comparison <- data.frame(comparison=c('HMM - Binary_Markov','HMM - Historical_Frequency'),
  difference=c(score[1]-score[2],score[1]-score[3]),lower95=dci[,1],upper95=dci[,2],row.names=NULL)
out <- file.path(dirname(path),'hmm_fitted_results'); dir.create(out,showWarnings=FALSE)
write.csv(params,file.path(out,'hmm_parameters.csv'),row.names=FALSE)
write.csv(fit_checks,file.path(out,'hmm_fit_checks.csv'),row.names=FALSE)
write.csv(pred,file.path(out,'hmm_test_predictions.csv'),row.names=FALSE)
write.csv(scores,file.path(out,'hmm_brier_scores.csv'),row.names=FALSE)
write.csv(comparison,file.path(out,'hmm_brier_comparison.csv'),row.names=FALSE)
write.csv(data.frame(training_days=n,testing_days=N,training_rain_frequency=clim,
  markov_rain_given_no_rain=q[1],markov_rain_given_rain=q[2],training_logLik=-best$value),
  file.path(out,'hmm_data_summary.csv'),row.names=FALSE)
print(params); print(fit_checks); print(scores); print(comparison)
if(any(theta <= eps*1.01 | theta >= 1-eps*1.01))
  message('At least one estimate is near the numerical boundary. Avoid ordinary Wald intervals and strong latent-regime interpretations.')
message('Results saved to: ',out)

# Parameter bootstrap confidence intervals
set.seed(2026)
B_param <- 500L

simulate_hmm <- function(theta, n) {
  a <- theta[1]
  b <- theta[2]
  e <- theta[3:4]
  z <- integer(n)
  z[1] <- sample(1:2, 1, prob = c(b, a)/(a+b))
  
  for (t in 2:n) {
    prob <- if (z[t-1] == 1) a else 1-b
    z[t] <- 1L + rbinom(1, 1, prob)
  }
  
  rbinom(n, 1, e[z])
}

parameter_vector <- function(th) {
  if (th[3] > th[4]) th <- th[c(2,1,4,3)]
  
  a <- th[1]
  b <- th[2]
  
  c(Initial_State_1 = b/(a+b),
    Initial_State_2 = a/(a+b),
    State_1_to_1 = 1-a,
    State_1_to_2 = a,
    State_2_to_1 = b,
    State_2_to_2 = 1-b,
    Rain_State_1 = th[3],
    Rain_State_2 = th[4])
}

point_params <- parameter_vector(theta)

boot_params <- matrix(
  NA_real_, nrow = B_param, ncol = 8
)

boundary_fit <- rep(NA, B_param)

for (i in seq_len(B_param)) {
  ys <- simulate_hmm(theta, n)
  
  initial_values <- rbind(
    theta, c(0.15, 0.30, 0.10, 0.80)
  )
  
  bf <- lapply(1:2, function(j) {
    tryCatch(
      optim(initial_values[j,], forward,
            obs = ys,
            method = "L-BFGS-B",
            lower = rep(eps, 4),
            upper = rep(1-eps, 4),
            control = list(maxit = 2000,
                           factr = 1e5,
                           ndeps = rep(1e-6, 4))),
      error = function(e) NULL
    )
  })
  
  good <- Filter(function(f) {
    !is.null(f) && f$convergence == 0 &&
      is.finite(f$value)
  }, bf)
  
  if (length(good) > 0) {
    best_boot <- good[[which.min(
      sapply(good, function(f) f$value)
    )]]
    
    th <- best_boot$par
    boot_params[i,] <- parameter_vector(th)
    
    boundary_fit[i] <- any(
      th <= eps*1.01 |
        th >= 1-eps*1.01
    )
  }
  
  if (i %% 50 == 0) {
    cat("Bootstrap:", i, "/", B_param, "\n")
  }
}

valid_boot <- complete.cases(boot_params)
n_success <- sum(valid_boot)

if (n_success < 0.9*B_param) {
  stop("Too many unsuccessful bootstrap fits.")
}

param_ci <- t(apply(
  boot_params[valid_boot,,drop=FALSE],
  2, quantile, probs=c(0.025,0.975)
))

parameter_intervals <- data.frame(
  Parameter = names(point_params),
  Estimate = as.numeric(point_params),
  Lower95 = param_ci[,1],
  Upper95 = param_ci[,2],
  Successful_Bootstrap = n_success
)

write.csv(
  parameter_intervals,
  file.path(out, "hmm_parameter_intervals.csv"),
  row.names = FALSE
)

bootstrap_diagnostics <- data.frame(
  Requested = B_param,
  Successful = n_success,
  Boundary_Fits = sum(boundary_fit, na.rm=TRUE)
)

write.csv(
  bootstrap_diagnostics,
  file.path(out, "hmm_bootstrap_diagnostics.csv"),
  row.names = FALSE
)

print(parameter_intervals)
print(bootstrap_diagnostics)

