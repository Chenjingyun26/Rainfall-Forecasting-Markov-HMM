# Supply five NOAA CSV paths as command-line arguments, ordered:
# Wilmington, Denver, Los Angeles, Seattle, Miami.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L) stop('Provide five station CSV paths in city order.')
cities <- c('Wilmington','Denver','Los Angeles','Seattle','Miami')
states <- c('Dry','Rain','Heavy Rain','Thunder-Reported Day')
classify <- function(prcp,thunder,cutoff) {
  out <- rep(NA_character_,length(prcp))
  th <- !is.na(thunder) & thunder==1
  out[th] <- states[4]
  out[!th & !is.na(prcp) & prcp==0] <- states[1]
  out[!th & !is.na(prcp) & prcp>0 & prcp<cutoff] <- states[2]
  out[!th & !is.na(prcp) & prcp>=cutoff] <- states[3]
  factor(out,levels=states)
}
fit_city <- function(file,city,p) {
 d <- read.csv(file,stringsAsFactors=FALSE)
 stopifnot(all(c('DATE','PRCP','WT03') %in% names(d)))
 d$DATE <- as.Date(d$DATE)
 d <- d[order(d$DATE),]
 train <- d$DATE < as.Date('2024-01-01') & d$DATE >= as.Date('2016-01-01')
 test <- d$DATE >= as.Date('2024-01-01') & d$DATE <= as.Date('2025-12-31')
 positive <- d$PRCP[train & !is.na(d$PRCP) & d$PRCP>0]
 threshold <- as.numeric(quantile(positive,p,na.rm=TRUE,type=7))
 state <- classify(d$PRCP,d$WT03,threshold)
 train_idx <- which(train)
 valid_train <- train_idx[-length(train_idx)]
 valid_train <- valid_train[valid_train+1L %in% train_idx &
    as.integer(d$DATE[valid_train+1L]-d$DATE[valid_train])==1L &
    !is.na(state[valid_train]) & !is.na(state[valid_train+1L])]
 counts <- table(factor(state[valid_train],levels=states),factor(state[valid_train+1L],levels=states))
 probs <- counts / rowSums(counts)
 probs[rowSums(counts)==0,] <- NA_real_
 test_idx <- which(test)
 origins <- test_idx[-length(test_idx)]
 origins <- origins[origins+1L %in% test_idx &
   as.integer(d$DATE[origins+1L]-d$DATE[origins])==1L &
   !is.na(state[origins]) & !is.na(state[origins+1L])]
 origins <- origins[!is.na(probs[as.character(state[origins]),1])]
 pred <- vapply(as.character(state[origins]),function(x) states[which.max(probs[x,])],character(1))
 actual <- as.character(state[origins+1L]); current <- as.character(state[origins])
 common <- names(which.max(table(state[train])))
 recall <- vapply(states,function(x) if(sum(actual==x)>0) mean(pred[actual==x]==x) else NA_real_,numeric(1))
 summary <- data.frame(city=city,percentile=p,threshold_inches=threshold,
   n_test_transitions=length(origins),accuracy=mean(pred==actual),
   persistence=mean(current==actual),most_common=mean(common==actual),
   heavy_rain_recall=recall[3],heavy_rain_test_count=sum(actual==states[3]))
 list(summary=summary,transitions=probs,train_states=table(state[train]),test_states=table(state[test]))
}
results <- list(); summaries <- list()
for(i in seq_along(cities)) for(p in c(.80,.85)) {
  z <- fit_city(args[i],cities[i],p)
  tag <- paste0(gsub(' ','_',cities[i]),'_',round(p*100))
  results[[tag]] <- z
  summaries[[tag]] <- z$summary
  write.csv(z$transitions,paste0('transition_',tag,'.csv'))
}
write.csv(do.call(rbind,summaries),'sensitivity_summary.csv',row.names=FALSE)
print(do.call(rbind,summaries),row.names=FALSE)
