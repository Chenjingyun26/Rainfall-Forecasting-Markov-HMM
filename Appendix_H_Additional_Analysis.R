# Additional analysis: run the entire script using Source in RStudio.
folder <- dirname(file.choose())
out <- file.path(folder, 'additional_analysis_results')
dir.create(out, showWarnings=FALSE)
cities <- c('Wilmington','Denver','Los Angeles','Seattle','Miami')
files <- c('wilmington','denver','los_angeles','seattle','miami')
states <- c('Dry','Rain','Heavy Rain','Thunder-Reported Day')
short <- c('Dry','Rain','Heavy Rain','Thunder')
models <- c('Markov','Persistence','Most Common')
set.seed(2026)
B <- 2000L; L <- 7L
scores <- list(); matrices <- list(); audit <- list(); excluded <- list(); missing <- list()
for (i in seq_along(cities)) {
  path <- file.path(folder,paste0(files[i],'.csv'))
  if (!file.exists(path)) path <- file.path(folder,paste0(files[i],'(1).csv'))
  d <- read.csv(path,stringsAsFactors=FALSE)
  d$DATE <- as.Date(d$DATE); d <- d[order(d$DATE),]
  stopifnot(!anyNA(d$DATE), !anyDuplicated(d$DATE))
  train <- d$DATE >= as.Date('2016-01-01') & d$DATE < as.Date('2024-01-01')
  test <- d$DATE >= as.Date('2024-01-01') & d$DATE <= as.Date('2025-12-31')
  k <- seq_len(nrow(d)-1L)
  consecutive <- as.integer(d$DATE[k+1L]-d$DATE[k])==1L
  for (p in c(.80,.85)) {
    q <- as.numeric(quantile(d$PRCP[train & !is.na(d$PRCP) & d$PRCP>0],p,type=7))
    s <- rep(NA_character_,nrow(d)); th <- !is.na(d$WT03) & d$WT03==1
    s[th] <- states[4]
    s[!th & !is.na(d$PRCP) & d$PRCP==0] <- states[1]
    s[!th & !is.na(d$PRCP) & d$PRCP>0 & d$PRCP<q] <- states[2]
    s[!th & !is.na(d$PRCP) & d$PRCP>=q] <- states[3]
    valid <- consecutive & !is.na(s[k]) & !is.na(s[k+1L])
    tr <- k[train[k] & train[k+1L]]
    te <- k[test[k] & test[k+1L]]
    a <- tr[valid[tr]]; b <- te[valid[te]]
    counts <- table(factor(s[a],levels=states),factor(s[a+1L],levels=states))
    P <- counts/rowSums(counts)
    b <- b[!is.na(P[s[b],1])]
    actual <- s[b+1L]
    pred <- vapply(s[b],function(x) states[which.max(P[x,])],character(1))
    common <- states[which.max(table(factor(s[train],levels=states)))]
    correct <- cbind(Markov=pred==actual, Persistence=s[b]==actual, Most_Common=common==actual)
    z <- matrix(NA_real_,length(te),3); z[match(b,te),] <- correct
    point <- colMeans(z,na.rm=TRUE)
    boot <- replicate(B, {
      starts <- sample.int(nrow(z),ceiling(nrow(z)/L),replace=TRUE)
      idx <- as.vector(outer(0:(L-1L),starts,'+'))
      idx <- ((idx-1L) %% nrow(z))+1L
      colMeans(z[idx[seq_len(nrow(z))],,drop=FALSE],na.rm=TRUE)
    })
    ci <- t(apply(boot,1,quantile,c(.025,.975),na.rm=TRUE))
    tag <- paste0(files[i],'_',round(p*100))
    scores[[tag]] <- data.frame(City=cities[i],Percentile=p,Model=models,N=length(b),Accuracy=as.numeric(point),Lower95=ci[,1],Upper95=ci[,2])
    matrices[[tag]] <- P
    write.csv(P,file.path(out,paste0('transition_',tag,'.csv')))
    if (files[i]=='seattle') {
      for (period in c('Training','Testing')) {
        mask <- if(period=='Training') train else test
        candidate <- if(period=='Training') tr else te
        kept <- if(period=='Training') a else b
        key <- paste(tag,period)
        audit[[key]] <- data.frame(Period=period,Percentile=p,Days=sum(mask),Missing_PRCP=sum(mask & is.na(d$PRCP)),Missing_State=sum(mask & is.na(s)),Candidate_Transitions=length(candidate),Valid_Transitions=length(kept),Excluded_Transitions=length(setdiff(candidate,kept)),Heavy_Rain_Days=sum(s[mask]==states[3],na.rm=TRUE),Heavy_Rain_Targets=sum(s[kept+1L]==states[3],na.rm=TRUE))
        bad <- setdiff(candidate,kept)
        excluded[[key]] <- data.frame(Period=period,Percentile=p,Origin_Date=d$DATE[bad],Target_Date=d$DATE[bad+1L],Origin_PRCP=d$PRCP[bad],Target_PRCP=d$PRCP[bad+1L],Origin_State=s[bad],Target_State=s[bad+1L])
      }
      not_scored <- setdiff(which(test & !is.na(s) & s==states[3]),b+1L)
      missing[[tag]] <- data.frame(Percentile=rep(p,length(not_scored)),Target_Date=d$DATE[not_scored],PRCP=d$PRCP[not_scored],Previous_Date=d$DATE[not_scored-1L],Previous_PRCP=d$PRCP[not_scored-1L],Previous_State=s[not_scored-1L])
    }
  }
}
results <- do.call(rbind,scores)
write.csv(results,file.path(out,'accuracy_comparison.csv'),row.names=FALSE)
write.csv(do.call(rbind,audit),file.path(out,'seattle_transition_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,excluded),file.path(out,'seattle_excluded_transitions.csv'),row.names=FALSE)
write.csv(do.call(rbind,missing),file.path(out,'seattle_heavy_rain_excluded.csv'),row.names=FALSE)
# Heatmaps: one page per percentile, common probability scale.
pal <- colorRampPalette(c('white','#2166AC'))(101)
for (p in c(80,85)) {
  pdf(file.path(out,paste0('transition_heatmaps_',p,'.pdf')),width=12,height=8)
  par(mfrow=c(2,3),mar=c(5,5,3,1))
  for (i in seq_along(cities)) {
    P <- matrices[[paste0(files[i],'_',p)]]
    plot(NA,xlim=c(.5,4.5),ylim=c(4.5,.5),xlab='',ylab='',axes=FALSE,main=cities[i],asp=1)
    for (r in 1:4) for (j in 1:4) {
      v <- P[r,j]
      rect(j-.5,r-.5,j+.5,r+.5,col=pal[1+round(100*v)],border='white')
      text(j,r,sprintf('%.3f',v),col=if(v>.55) 'white' else 'black',cex=.9)
    }
    axis(1,1:4,short,las=2,cex.axis=.75); axis(2,1:4,short,las=2,cex.axis=.75)
    mtext('Next state',1,4,cex=.8); mtext('Current state',2,4,cex=.8)
  }
  plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
  for(j in 1:100) rect((j-1)/100,.45,j/100,.58,col=pal[j],border=NA)
  text(c(0,.5,1),.38,c('0','0.5','1')); text(.5,.7,'Transition probability'); text(.5,.2,paste0(p,'th percentile'))
  dev.off()
}
# Primary accuracy figure with conditional block-bootstrap intervals.
for (p in c(.80,.85)) {
  pdf(file.path(out,paste0('accuracy_comparison_',round(p*100),'.pdf')),width=10,height=5.5)
  par(mar=c(5,5,4,1))
  r <- results[results$Percentile==p,]
  mat <- matrix(r$Accuracy,nrow=3); lo <- matrix(r$Lower95,nrow=3); hi <- matrix(r$Upper95,nrow=3)
  pos <- barplot(mat,beside=TRUE,names.arg=cities,ylim=c(0,1.13),col=c('#2166AC','#999999','#D1DCE5'),ylab='Test accuracy',main=paste0('Model comparison (',round(p*100),'th percentile)'))
  arrows(pos,lo,pos,hi,angle=90,code=3,length=.035)
  legend('top',models,fill=c('#2166AC','#999999','#D1DCE5'),horiz=TRUE,bty='n',cex=.85)
  dev.off()
}
print(results,row.names=FALSE)
print(do.call(rbind,audit),row.names=FALSE)
print(do.call(rbind,missing),row.names=FALSE)
message('Results saved to: ',out)
