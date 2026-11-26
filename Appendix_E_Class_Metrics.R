# Run in RStudio with the five CSV files in one folder.
folder <- dirname(file.choose())
cities <- c('Wilmington','Denver','Los Angeles','Seattle','Miami')
files <- c('wilmington','denver','los_angeles','seattle','miami')
states <- c('Dry','Rain','Heavy Rain','Thunder-Reported Day')
output <- file.path(folder,'class_metrics_results')
dir.create(output,showWarnings=FALSE)
metrics <- list(); summaries <- list()

# Conditional test-period uncertainty: circular blocks of 7 calendar days.
# Missing transitions remain missing, rather than joining separated dates.
set.seed(2026)
B <- 2000L
L <- 7L
intervals <- list()
metric_vector <- function(x, y) {
  cm <- table(factor(x, levels=states), factor(y, levels=states))
  tp <- diag(cm); n <- rowSums(cm); k <- colSums(cm)
  c(sum(tp)/sum(cm),
    ifelse(k>0,tp/k,NA_real_),
    ifelse(n>0,tp/n,NA_real_),
    ifelse(n+k>0,2*tp/(n+k),NA_real_))
}

for(i in seq_along(cities)) {
  path <- file.path(folder,paste0(files[i],'.csv'))
  if(!file.exists(path)) path <- file.path(folder,paste0(files[i],'(1).csv'))
  d <- read.csv(path,stringsAsFactors=FALSE)
  d$DATE <- as.Date(d$DATE); d <- d[order(d$DATE),]
  train <- d$DATE >= as.Date('2016-01-01') & d$DATE < as.Date('2024-01-01')
  test <- d$DATE >= as.Date('2024-01-01') & d$DATE <= as.Date('2025-12-31')
  for(p in c(.80,.85)) {
    q <- as.numeric(quantile(d$PRCP[train & !is.na(d$PRCP) & d$PRCP>0],p,type=7))
    s <- rep(NA_character_,nrow(d)); th <- !is.na(d$WT03) & d$WT03==1
    s[th] <- states[4]
    s[!th & !is.na(d$PRCP) & d$PRCP==0] <- states[1]
    s[!th & !is.na(d$PRCP) & d$PRCP>0 & d$PRCP<q] <- states[2]
    s[!th & !is.na(d$PRCP) & d$PRCP>=q] <- states[3]
    k <- seq_len(nrow(d)-1L)
    valid <- !is.na(s[k]) & !is.na(s[k+1L]) & as.integer(d$DATE[k+1L]-d$DATE[k])==1L
    a <- k[valid & train[k] & train[k+1L]]
    b <- k[valid & test[k] & test[k+1L]]
    counts <- table(factor(s[a],levels=states),factor(s[a+1L],levels=states))
    P <- counts/rowSums(counts)
    b <- b[!is.na(P[s[b],1])]
    pred <- vapply(s[b],function(x) states[which.max(P[x,])],character(1))
    actual <- s[b+1L]
    cm <- table(Actual=factor(actual,levels=states),Predicted=factor(pred,levels=states))
    tp <- diag(cm); support <- rowSums(cm); predicted <- colSums(cm)
    precision <- ifelse(predicted>0,tp/predicted,NA_real_)
    recall <- ifelse(support>0,tp/support,NA_real_)
    f1 <- ifelse(support+predicted>0,2*tp/(support+predicted),NA_real_)
    tag <- paste0(files[i],'_',round(p*100))
    metrics[[tag]] <- data.frame(City=cities[i],Percentile=p,State=states,Support=support,Predicted=predicted,Precision=precision,Recall=recall,F1=f1)
    summaries[[tag]] <- data.frame(City=cities[i],Percentile=p,Threshold=q,N=sum(cm),Accuracy=sum(tp)/sum(cm))
    write.csv(cm,file.path(output,paste0('confusion_',tag,'.csv')))

    # Keep the complete test-origin calendar, including invalid transitions.
    origins <- k[test[k] & test[k+1L]]
    x <- s[origins+1L]
    y <- rep(NA_character_,length(origins))
    y[match(b,origins)] <- pred
    x[is.na(y)] <- NA_character_
    n_orig <- length(origins)
    point <- metric_vector(x,y)
    boot <- replicate(B, {
      starts <- sample.int(n_orig,ceiling(n_orig/L),replace=TRUE)
      idx <- as.vector(outer(0:(L-1L),starts,'+'))
      idx <- ((idx-1L) %% n_orig)+1L
      idx <- idx[seq_len(n_orig)]
      metric_vector(x[idx],y[idx])
    })
    lower <- upper <- rep(NA_real_,length(point))
    status <- rep("OK",length(point))
    for(j in seq_along(point)) {
      v <- boot[j,is.finite(boot[j,])]
      if(!is.finite(point[j])) {
        status[j] <- "Undefined: zero denominator"
      } else if(length(v)<0.95*B) {
        status[j] <- "Too few defined bootstrap replicates"
      } else if(length(unique(v))<2L) {
        status[j] <- "Degenerate bootstrap: interval not informative"
      } else {
        ci <- quantile(v,c(.025,.975),names=FALSE)
        lower[j] <- ci[1]; upper[j] <- ci[2]
      }
    }
    intervals[[tag]] <- data.frame(City=cities[i],Percentile=p,
      State=c("Overall",rep(states,3)),
      Metric=c("Accuracy",rep(c("Precision","Recall","F1"),each=4)),
      Estimate=point,Lower95=lower,Upper95=upper,
      ValidReplicates=rowSums(is.finite(boot)),Status=status)

    print(cm)
  }
}
write.csv(do.call(rbind,metrics),file.path(output,'class_metrics.csv'),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(output,'accuracy_check.csv'),row.names=FALSE)
print(do.call(rbind,summaries),row.names=FALSE)
print(do.call(rbind,metrics),row.names=FALSE)

write.csv(do.call(rbind,intervals),file.path(output,"confidence_intervals.csv"),row.names=FALSE)
print(do.call(rbind,intervals),row.names=FALSE)
