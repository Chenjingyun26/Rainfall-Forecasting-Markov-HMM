# Training-period uncertainty for thresholds and four-state transitions.
# Run the entire script with Source; select any of the five original CSVs.
folder <- dirname(file.choose())
out <- file.path(folder, 'transition_interval_results')
dir.create(out, showWarnings = FALSE)
cities <- c('Wilmington','Denver','Los Angeles','Seattle','Miami')
files <- c('wilmington','denver','los_angeles','seattle','miami')
states <- c('Dry','Rain','Heavy Rain','Thunder-Reported Day')
set.seed(2026)
B <- 2000L; L <- 7L
thresholds <- list(); transitions <- list()
classify <- function(r, th, q) {
  s <- rep(NA_integer_, length(r))
  s[th] <- 4L
  s[!th & !is.na(r) & r == 0] <- 1L
  s[!th & !is.na(r) & r > 0 & r < q] <- 2L
  s[!th & !is.na(r) & r >= q] <- 3L
  s
}
get_q <- function(r, p) as.numeric(quantile(r[!is.na(r) & r > 0], p, type=7))
get_counts <- function(s, idx, consecutive) {
  idx <- idx[idx < length(s)]
  idx <- idx[consecutive[idx] & !is.na(s[idx]) & !is.na(s[idx+1L])]
  matrix(tabulate(s[idx] + 4L*(s[idx+1L]-1L), nbins=16L), 4, 4)
}
get_P <- function(counts) {
  den <- rowSums(counts)
  P <- counts/den
  P[den == 0,] <- NA_real_
  P
}
interval <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < .95*B) return(c(NA_real_,NA_real_,n))
  ci <- as.numeric(quantile(x,c(.025,.975)))
  if (ci[1] == ci[2]) ci[] <- NA_real_
  c(ci,n)
}
for (i in seq_along(cities)) {
  path <- file.path(folder,paste0(files[i],'.csv'))
  if (!file.exists(path)) path <- file.path(folder,paste0(files[i],'(1).csv'))
  d <- read.csv(path,stringsAsFactors=FALSE)
  d$DATE <- as.Date(d$DATE)
  stopifnot(!anyNA(d$DATE), !anyDuplicated(d$DATE))
  d <- d[d$DATE >= as.Date('2016-01-01') & d$DATE < as.Date('2024-01-01'),]
  d <- d[order(d$DATE),]
  n <- nrow(d); r <- d$PRCP
  stopifnot(n == 2922L, all(diff(d$DATE)==1))
  th <- !is.na(d$WT03) & d$WT03==1
  consecutive <- as.integer(diff(d$DATE))==1L
  # Resample calendar days and their original outgoing pairs together.
  # No artificial transition across block boundaries or the training endpoint.
  idx_boot <- replicate(B, {
    starts <- sample.int(n, ceiling(n/L), replace=TRUE)
    idx <- as.vector(outer(0:(L-1L), starts, '+'))
    ((idx-1L) %% n + 1L)[seq_len(n)]
  })
  for (p in c(.80,.85)) {
    q <- get_q(r,p)
    counts <- get_counts(classify(r,th,q),seq_len(n),consecutive)
    P <- get_P(counts)
    boot <- vapply(seq_len(B), function(j) {
      idx <- idx_boot[,j]
      qb <- get_q(r[idx],p)
      Pb <- get_P(get_counts(classify(r,th,qb),idx,consecutive))
      c(qb,as.vector(Pb))
    }, numeric(17))
    ci <- t(apply(boot,1,interval))
    tag <- paste0(files[i],'_',round(100*p))
    thresholds[[tag]] <- data.frame(City=cities[i],Percentile=p,Estimate=q,Lower95=ci[1,1],Upper95=ci[1,2],Valid_Bootstrap=ci[1,3])
    transitions[[tag]] <- data.frame(City=cities[i],Percentile=p,
      Current=rep(states,4),Next_State=rep(states,each=4),
      Count=as.vector(counts),Origin_Count=rep(rowSums(counts),4),
      Estimate=as.vector(P),Lower95=ci[-1,1],Upper95=ci[-1,2],Valid_Bootstrap=ci[-1,3])
    cat(cities[i],p,': threshold =',q,'; valid transitions =',sum(counts),'\n')
  }
}
threshold_results <- do.call(rbind,thresholds)
transition_results <- do.call(rbind,transitions)
write.csv(threshold_results,file.path(out,'threshold_intervals.csv'),row.names=FALSE)
write.csv(transition_results,file.path(out,'transition_intervals.csv'),row.names=FALSE)
print(threshold_results,row.names=FALSE)
cat('Transition records:',nrow(transition_results),'\n')
cat('Intervals unavailable:',sum(is.na(transition_results$Lower95)), '\n')
message('Results saved to: ',out)
