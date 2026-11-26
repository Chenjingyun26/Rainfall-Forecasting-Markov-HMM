# Appendix I: empirical two-day-ahead rainfall forecasts, Wilmington
# Uses Appendix F fitted parameters; no test-period refitting.
# Base R only. First choose Wilmington CSV, then hmm_parameters.csv.
path <- file.choose()
d <- read.csv(path, stringsAsFactors=FALSE)
stopifnot(all(c('DATE','PRCP','STATION') %in% names(d)),
          all(d$STATION == 'USW00013748'))
d$DATE <- as.Date(d$DATE)
d <- d[order(d$DATE), ]
stopifnot(identical(as.character(d$DATE), as.character(seq(as.Date('2016-01-01'),as.Date('2025-12-31'),by='day'))),
          all(is.finite(d$PRCP)), all(d$PRCP>=0))
p <- read.csv(file.choose(), stringsAsFactors=FALSE)
stopifnot(nrow(p)==2L, setequal(p$state,c('State_1','State_2')))
p <- p[match(c('State_1','State_2'),p$state), ]
A <- as.matrix(p[,c('to_state_1','to_state_2')])
e <- p$rain_emission
pi0 <- p$initial_probability
stopifnot(all(is.finite(c(A,e,pi0))), all(c(A,e,pi0)>=0), all(c(A,e,pi0)<=1),
          max(abs(rowSums(A)-1))<1e-8, abs(sum(pi0)-1)<1e-8)
y <- as.integer(d$PRCP>0)
n <- sum(d$DATE<as.Date('2024-01-01'))
y_train <- y[seq_len(n)]
# Binary Markov transition probabilities estimated only from training pairs.
q <- sapply(0:1,function(k) mean(y_train[-1][y_train[-n]==k]))
M <- cbind(1-q,q)
clim <- mean(y_train)
# Filter sequentially: alpha[t,] depends only on observations through t.
alpha <- matrix(NA_real_,length(y),2)
for(t in seq_along(y)) {
  prior <- if(t==1L) pi0 else as.numeric(alpha[t-1L,] %*% A)
  w <- prior * if(y[t]==1L) e else (1-e)
  alpha[t,] <- w/sum(w)
}
# Common 731 target dates at both horizons, including Jan 1, 2024.
idx <- seq.int(n+1L,length(y))
out <- file.path(dirname(path),'hmm_two_day_results')
dir.create(out,showWarnings=FALSE)
set.seed(2026); B <- 2000L; L <- 7L; N <- length(idx)
# Paired moving blocks; same resamples for all models and both horizons.
boot_idx <- replicate(B, {
  st <- sample.int(N-L+1L,ceiling(N/L),replace=TRUE)
  as.vector(outer(0:(L-1L),st,'+'))[seq_len(N)]
})
all_pred <- all_scores <- all_diff <- list()
for(h in 1:2) {
  origin <- idx-h
  Ah <- if(h==1L) A else A %*% A
  Mh <- if(h==1L) M else M %*% M
  pred <- data.frame(Horizon=h, Origin_Date=d$DATE[origin], Target_Date=d$DATE[idx],
                     Actual=y[idx], HMM=as.numeric(alpha[origin,,drop=FALSE] %*% Ah %*% e),
                     Binary_Markov=Mh[y[origin]+1L,2], Historical_Frequency=rep(clim,N))
  loss <- (as.matrix(pred[,c('HMM','Binary_Markov','Historical_Frequency')])-pred$Actual)^2
  score <- colMeans(loss)
  boot <- apply(boot_idx,2,function(ii) colMeans(loss[ii,,drop=FALSE]))
  ci <- t(apply(boot,1,quantile,probs=c(.025,.975)))
  scores <- data.frame(Horizon=h,Model=names(score),N=N,Brier=as.numeric(score),Lower95=ci[,1],Upper95=ci[,2],row.names=NULL)
  delta <- rbind(boot[1,]-boot[2,],boot[1,]-boot[3,])
  dci <- t(apply(delta,1,quantile,probs=c(.025,.975)))
  diff <- data.frame(Horizon=h,Comparison=c('HMM - Binary_Markov','HMM - Historical_Frequency'),Difference=as.numeric(score[1]-score[2:3]),Lower95=dci[,1],Upper95=dci[,2],row.names=NULL)
  all_pred[[h]] <- pred; all_scores[[h]] <- scores; all_diff[[h]] <- diff
}
# Intervals quantify test performance conditional on fixed fitted models.
write.csv(do.call(rbind,all_pred),file.path(out,'hmm_horizon_predictions.csv'),row.names=FALSE)
write.csv(do.call(rbind,all_scores),file.path(out,'hmm_horizon_scores.csv'),row.names=FALSE)
write.csv(do.call(rbind,all_diff),file.path(out,'hmm_horizon_comparison.csv'),row.names=FALSE)
print(do.call(rbind,all_scores),row.names=FALSE)
print(do.call(rbind,all_diff),row.names=FALSE)
message('Results saved to: ',out)
