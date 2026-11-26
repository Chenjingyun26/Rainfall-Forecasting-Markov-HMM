# Appendix G: Fitted HMM regime diagnostics for Wilmington
# Select the original Wilmington CSV, then the Appendix F parameter CSV.
path <- file.choose()
d <- read.csv(path, stringsAsFactors=FALSE)
stopifnot(all(c('DATE','STATION','PRCP','TMAX','TMIN','AWND') %in% names(d)), all(d$STATION=='USW00013748'))
d$DATE <- as.Date(d$DATE)
d <- d[order(d$DATE),]
stopifnot(identical(as.character(d$DATE),as.character(seq(as.Date('2016-01-01'),as.Date('2025-12-31'),by='day'))), all(is.finite(d$PRCP)), all(d$PRCP>=0))
par_path <- file.choose()
p <- read.csv(par_path, stringsAsFactors=FALSE)
p <- p[match(c('State_1','State_2'),p$state),]
A <- as.matrix(p[,c('to_state_1','to_state_2')]); e <- p$rain_emission
pi <- p$initial_probability
stopifnot(all(is.finite(c(A,e,pi))),all(abs(rowSums(A)-1)<1e-8),abs(sum(pi)-1)<1e-8,all(c(A,e,pi)>=0),all(c(A,e,pi)<=1),e[1]<=e[2])
y <- as.integer(d$PRCP>0)
filt <- forecast <- numeric(nrow(d)); q <- pi
for(i in seq_len(nrow(d))) {
  if(i>1L) q <- as.numeric(q %*% A)
  forecast[i] <- sum(q*e)
  w <- q * if(y[i]==1L) e else (1-e)
  q <- w/sum(w)
  filt[i] <- q[2]
}
daily <- data.frame(DATE=d$DATE,Period=ifelse(d$DATE<as.Date('2024-01-01'),'Training','Testing'), State=ifelse(filt>0.5,'State_2','State_1'), Filtered_State_2=filt, Rain_Forecast=forecast,Rain=y,Temperature=(d$TMAX+d$TMIN)/2,Wind=d$AWND,Precipitation=d$PRCP)
vars <- c('Temperature','Wind','Precipitation','Rain'); states <- c('State_1','State_2')
# Conditional on the fitted HMM, resample seven-day blocks within each period.
# Resample joint daily records, not individual state subsets.
stat <- function(z) {
  ans <- sapply(vars,function(v) sapply(states,function(s) {
    x <- z[[v]][z$State==s]; x <- x[is.finite(x)]
    if(length(x)) mean(x) else NA_real_
  }))
  c(as.vector(ans),ans[2,]-ans[1,])
}
set.seed(2026); B <- 2000L; L <- 7L
summaries <- intervals <- matches <- list()
for(period in c('Training','Testing')) {
  z <- daily[daily$Period==period,]; N <- nrow(z)
  for(s in states) for(v in vars) {
    x <- z[[v]][z$State==s]; x <- x[is.finite(x)]
    summaries[[length(summaries)+1L]] <- data.frame(Period=period,State=s,Variable=v,N=length(x),Mean=mean(x),SD=sd(x),Q1=unname(quantile(x,.25)),Median=median(x),Q3=unname(quantile(x,.75)))
  }
  point <- stat(z)
  boot <- replicate(B,{
    starts <- sample.int(N-L+1L,ceiling(N/L),replace=TRUE)
    idx <- as.vector(outer(0:(L-1L),starts,'+'))[seq_len(N)]
    stat(z[idx,])
  })
  for(j in seq_along(point)) {
    x <- boot[j,is.finite(boot[j,])]; ci <- if(length(x)>=.9*B) unname(quantile(x,c(.025,.975))) else c(NA_real_,NA_real_)
    status <- if(anyNA(ci)) 'Insufficient valid replicates' else if(abs(diff(ci))<1e-12) 'Degenerate: not informative' else 'Approximate conditional interval'
    if(status=='Degenerate: not informative') ci <- c(NA_real_,NA_real_)
    intervals[[length(intervals)+1L]] <- data.frame(Period=period,Group=c(rep(states,4),rep('State_2_minus_State_1',4))[j],Variable=c(rep(vars,each=2),vars)[j],Estimate=point[j],Lower95=ci[1],Upper95=ci[2],Valid_Bootstrap=length(x),Status=status)
  }
  tab <- table(State=factor(z$State,levels=states),Observed_Rain=factor(z$Rain,levels=0:1))
  matches[[period]] <- data.frame(Period=period,as.data.frame(tab))
}
out <- file.path(dirname(path),'hmm_regime_results'); dir.create(out,showWarnings=FALSE)
write.csv(daily,file.path(out,'regime_daily.csv'),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(out,'regime_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,intervals),file.path(out,'regime_intervals.csv'),row.names=FALSE)
write.csv(do.call(rbind,matches),file.path(out,'regime_rain_agreement.csv'),row.names=FALSE)
for(period in c('Training','Testing')) {
  z <- daily[daily$Period==period,]
  pdf(file.path(out,paste0('regime_boxplots_',tolower(period),'.pdf')),width=10,height=4)
  par(mfrow=c(1,3),mar=c(4,4,3,1))
  boxplot(Temperature~State,data=z,main='Daily mean temperature',ylab='Temperature (F)',xlab='Filtered state')
  boxplot(Wind~State,data=z,main='Daily mean wind speed',ylab='Wind speed (mph)',xlab='Filtered state')
  boxplot(Precipitation~State,data=z,main='Daily precipitation',ylab='Precipitation (inches)',xlab='Filtered state')
  dev.off()
}
print(do.call(rbind,summaries),row.names=FALSE)
print(do.call(rbind,matches),row.names=FALSE)
print(do.call(rbind,intervals),row.names=FALSE)
message('Results saved to: ',out)
