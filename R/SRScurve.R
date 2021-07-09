SRScurve <- function(data, metric = "richness", step = 50, sample = 0, max.sample.size = 0,
                     rarefy.comparison = FALSE, rarefy.repeats = 10, 
                     rarefy.comparison.legend = FALSE, xlab = "sample size", 
                     ylab = "richness", label = FALSE, col, lty, ...) {
  #function that calculates diversity indices in a series of sequencing depth/sample size
  #currently the diversity indices are limited by the available options in vegan::diversity
  #data: feature-table
  #Cminseq: vector with sequencing depths
  #metric: richness/shannon/simpson/invsimpson
  SRSdiversity<-function(data,Cminseq,metric){
    nspec=as.data.frame(matrix(data=NA,nrow = ncol(data),ncol=length(Cminseq))) #create output dataframe
    colnames(nspec) <- paste("N", Cminseq, sep="")
    rownames(nspec) <- colnames(data)
    if (metric=="richness"){  #if the chosen metric was richness, we simply count the number of observed species
      for (i in seq(1,ncol(data),1)){
        j=1
        for (cutoff in Cminseq){
          nspec[i,j]=vegan::specnumber(t(SRS(data = as.data.frame(data[,i]), Cmin = cutoff)))
          j=j+1
        }
      }
    }
    else{ #if the chosen metric was any other, we use diversity function from vegan package
      for (i in seq(1,ncol(data),1)){
        j=1
        for (cutoff in Cminseq){
          nspec[i,j]=vegan::diversity(SRS(data = as.data.frame(data[,i]), Cmin = cutoff),index=metric)
          j=j+1
        }
      }
    }
    attr(nspec, "Subsample") <- Cminseq
    return(nspec)
  }
  
  #the function RSdiversity does the same than SRS diversity but after random subsampling instead of SRS
  RSdiversity<-function(data,Cminseq,metric){
    nspec=as.data.frame(matrix(data=NA,nrow = ncol(data),ncol=length(Cminseq))) #create output dataframe
    colnames(nspec) <- paste("N", Cminseq, sep="")
    rownames(nspec) <- colnames(data)
    if (metric=="richness"){  #if the chosen metric was richness, we simply count the number of observed species
      for (i in seq(1,ncol(data),1)){
        j=1
        for (cutoff in Cminseq){
          nspec[i,j]=median(replicate(vegan::specnumber(rrarefy(t(data[,i]),sample=cutoff)),n=rarefy.repeats))
          j=j+1
        }
      }
    }
    else{ #if the chosen metric was any other, we use diversity function from vegan package
      for (i in seq(1,ncol(data),1)){
        j=1
        for (cutoff in Cminseq){
          nspec[i,j]=median(replicate(vegan::diversity(rrarefy(t(data[,i]),sample=cutoff),index=metric),n=rarefy.repeats))#repeat the estimate 100 times
          j=j+1
        }
      }
    }
    attr(nspec, "Subsample") <- Cminseq
    return(nspec)
  }
  
  ## matrix is faster than data.frame
  x <- as.matrix(t(data))
  ## check input data: must be counts
  if (!identical(all.equal(x, round(x)), TRUE))
    stop("function accepts only integers (counts)")
  ## sort out col and lty
  if (missing(col)){
    col <- par("col")
    ncolors <- 0}
  else{
    ncolors <- length(col)
  }
  if (missing(lty)){
    lty <- par("lty")
    nltypes <- 0}
  else{ 
    nltypes <- length(lty)
    }
  tot <- rowSums(x)
  S <- specnumber(x)
  ## remove empty rows or we fail
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0,, drop =FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  ncolors=length(col)
  ## rep col and lty to appropriate length
  if(rarefy.comparison==T & length(col)<2*nr){
    col <- rep(col, length.out = 2*nr)}
  if(rarefy.comparison==T & length(lty)<2*nr){
    lty <- rep(lty, length.out = 2*nr)}
  if(rarefy.comparison==F){
    col <- rep(col, length.out = nr)
    lty <- rep(lty, length.out = nr)}
  ## Normalize by SRS
  outsrs <- lapply(seq_len(nr), function(i) {
    if (max.sample.size!=0 & tot[i] > max.sample.size){
      tot[i] <- max.sample.size 
    }
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i])
      n <- c(n, tot[i])
    drop(SRSdiversity(as.data.frame(x[i,]),Cminseq = n, metric = metric)) #using the function declared above
  })
  if(rarefy.comparison==T){ #Normalize by RS
    outrs <- lapply(seq_len(nr), function(i) {
      if (max.sample.size!=0 & tot[i] > max.sample.size){
        tot[i] <- max.sample.size 
      }
      n <- seq(1, tot[i], by = step)
      if (n[length(n)] != tot[i])
        n <- c(n, tot[i])
      drop(RSdiversity(as.data.frame(x[i,]),Cminseq = n, metric = metric)) #using the other function declared above
    })}
  Nmax <- sapply(outsrs, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(outsrs, max)
  ## set up plot
  plot(c(1, max(Nmax)), c(0, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
  ## rarefied richnesses for given 'sample'
  if (sample>0) {
    col1<-col
    lty1<-lty
    abline(v = sample)
    for (ln in seq_along(outsrs)) {
      rare <- (approx(x = attr(outsrs[[ln]], "Subsample"), y = outsrs[[ln]],xout = sample, rule = 1)$y)
      abline(h = rare, lwd=0.5,col = col1[ln], lty = lty1[ln])
      if (rarefy.comparison==T){ #if the user chose to compare with random subsampling
        rare <- (approx(x = attr(outrs[[ln]], "Subsample"), y = outrs[[ln]],xout = sample, rule = 1)$y)
        abline(h = rare, lwd=0.5,col = col1[ln+1], lty = lty1[ln+1])
        col1<-col1[-1]
        lty1<-lty1[-1]}
    }
  }
  ## rarefaction curves
  for (ln in seq_along(outsrs)) {
    N <- attr(outsrs[[ln]], "Subsample")
    outln<-as.numeric(outsrs[[ln]][1,])
    lines(N,outln,col = col[ln], lty = lty[ln], ...)
    if (rarefy.comparison==T){ #if the user chose to compare with random subsampling
      outln<-as.numeric(outrs[[ln]][1,])
      lines(N,outln,col = col[ln+1], lty = lty[ln+1], ...)
      col<-col[-1]
      lty<-lty[-1]}
    if (label){
      ordilabel(cbind(max(N), tail(outln,n=1)), labels=rownames(x)[ln], ...)}
  }
  if (rarefy.comparison.legend==T){
    if (!missing(col)&ncolors==2&(nltypes!=2|missing(lty))){
      if(ncol(data)%%2!=0){
      legend("bottomright", legend=c("SRS", "rarefy"),
             col=c(col[2], col[1]),lty="solid", cex=0.8)}
      else{
        legend("bottomright", legend=c("SRS", "rarefy"),
               col=c(col[1], col[2]),lty="solid", cex=0.8)
      }
    }
    if (!missing(lty)&nltypes==2&(ncolors!=2|missing(col))){
      if(ncol(data)%%2!=0){
      legend("bottomright", legend=c("SRS", "rarefy"),
             col="black", lty=c(lty[2], lty[1]), cex=0.8)}
      else{
      legend("bottomright", legend=c("SRS", "rarefy"),
               col="black", lty=c(lty[1], lty[2]), cex=0.8)  
      }
    }
    if (!missing(lty)&!missing(col)&nltypes==2&ncolors==2){
      if(ncol(data)%%2!=0){
      legend("bottomright", legend=c("SRS", "rarefy"),
             col=c(col[2], col[1]), lty=c(lty[2], lty[1]), cex=0.8) }
      else{
        legend("bottomright", legend=c("SRS", "rarefy"),
               col=c(col[1], col[2]), lty=c(lty[1], lty[2]), cex=0.8) 
      }
  }
  if((nltypes!=2|missing(lty))&(ncolors!=2|missing(col))){
  warning("rarefy.comparison.legend only works when exactly 2 colors ('col') and/or 2 line types ('lty') were chosen")
    }
  }
  plot.base <- recordPlot()
  return (plot.base)
}