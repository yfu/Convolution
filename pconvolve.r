
# setwd("~/Challenge/results/2013-02-11/")
p.convolve <- function(fn, tau.amp=0.5, tau.del=-0.5, p.value.cutoff=0.05, fdr.cutoff=0.25, plot.pdf=TRUE){
# p.value.cutoff = 0.05
# fdr.value.cutoff = 0.25

# rat represents the data frame that contains log2 ratios
# file = "genome_data.tsv"
# fn = "chrom7"
original.tab <- read.table(file=fn, header=TRUE)
rat.table = original.tab
rat.table[is.na(rat.table)] <- 0 # Handle the NAs
head(rat.table)
sample.num<- length(rat.table)
# sample.num = sample.num - 2
# rat.part <- rat.table[, 3:(sample.num + 2)]
rat.part <- rat.table
head(rat.part, n=10)


thres.amp = tau.amp

thres.del = tau.del


# TODO: Maybe we need to normalize the log2ratios across samples?

# TODO: add a option to let users choose whether do 1[cond] or not
# See it here:
# https://projects.zoho.com/portal/montilab#wiki/465555000000032025/Scoring-SCNAs.html
rat.amp <- rat.part
rat.amp[rat.amp < thres.amp] <- 0
rat.del <- rat.part
rat.del[rat.del >thres.del] <- 0
rats = list(rat.amp, rat.del)
# As convolve() expects the probability density distribution as the
# input, I need to generate this here. It is something like this:
# convolve(c(0.1, 0.2, 0.3, 0.5, 0.3, 0.2, 0.1), rev(c(0.1, 0.2, 0.3, 0.5,
# 0.3, 0.2, 0.1)), type="open"). The values in the vector represent the pdf
# at each corresponding x

cur.status= 1 # Flag whether amp (1) or del(2) is being processes
final.results = list()
for (rat in rats) {
  # When cur.status == 1, the current loop is dealing with amp only
  # When cur.status == 2, the current loop is dealing with del only
  if (cur.status==1) {status = "amp"} else {status = "del"}
  max.vector <- sapply(rat, max)
  rat.max <- max(max.vector)
  min.vector <- sapply(rat, min)
  rat.min <- min(min.vector)

  ## step = 0.001
  step = 0.01 # The bigger the value is, the more accurate pdf (by convolution)
            # is. However, execution time will increase dramatically.
# Given the accuracy of the step, I need to round rat.max and rat.min
# to th nearest value so that I can put them in bins.
  rat.max <- ceiling(rat.max * 1 / step) * step
  rat.min <- floor(rat.min * 1 / step) * step
  breaks <- seq(rat.min, rat.max, step)

# Prepare for the next loop
  h <- hist(rat[ , 1], breaks=breaks, plot=FALSE) # Store the results of convolution
  conv.result = h$counts / sum(h$counts)

  for(i in 2:(sample.num)) {
  # num.zero.pad <- length(conv.result) - length(rat[,1])
  # One weird thing: You cannot use plot=FALSE and freq=FALSE at
  # the same time
    h <- hist(rat[,i], breaks=breaks, plot=FALSE)
    sample.pdf <- h$counts / sum(h$counts)
  # sample.pdf <- c(rep(0, times=num.zero.pad/2), sample.pdf, rep(0,times=num.zero.pad/2)) # Add zeroes surrounding
    conv.result <- convolve(conv.result, rev(sample.pdf), type="o")
}
  breaks.conv <- seq(rat.min * sample.num, rat.max * sample.num - step * sample.num, step)
# Round them as there are some nasty accuracy things.
  breaks.conv <- round(breaks.conv / step) * step
# I substract step * sample.num so that the lengths of conv.results and breaks.conv are equal.
# plot(breaks.conv, conv.result)
# A closer look at the distribution
  if (plot.pdf == TRUE) {
  # If true, plot it to file, else print it to screen (can be controlled
  # from outside)
    pdf(paste(fn, "_cdf.pdf", sep=""))
    plot(breaks.conv, conv.result, xlim=c(-5,5))
    dev.off()} else {
      plot(breaks.conv, conv.result, xlim=c(-5,5))
    }
    
# Get the sum for each locus
  sum.cols <- apply(rat, 1, sum)

  sum.pdf.df <- data.frame(x=breaks.conv, dens=conv.result)
  cdf.sum <- cumsum(conv.result)
  sum.cdf.df <- data.frame(x=breaks.conv, cumulative=cdf.sum)

# plot(sum.pdf.df, main="pdf of sum of ratios for loci")
# plot for cdf
# plot(sum.cdf.df, main="cdf of sum of ratios for loci")

# Get significant sites
  median.sum <- with(sum.cdf.df, min(x[cumulative>0.5]))
  sig.vector = c()
  flag.vector = c()

  sum.cols.round <- (round(sum.cols * 1 / step)) * step


  for ( i in 1:(length(sum.cols)) ) {
  # Just a p-value for one locus
    significance = NA
    significance <- with(sum.cdf.df, cumulative[x==sum.cols.round[i]])

    if (status == "amp") {
      # When a correspoding cumulative probability cannot be found
      # Give it a p-value of 1. This is because the whole row only
      # contains 0s
      if (length(significance) == 0) {sig.vector[i] = 1}
      else { sig.vector[i] <- 1 - significance }
    }
    else {
      if (length(significance) == 0) {sig.vector[i] = 1}
      else{ sig.vector[i] <- significance }
    }
  }

  results.df <- data.frame(sum=sum.cols, pval=sig.vector)

  results.df <- cbind(original.tab, results.df)
  ## head(original.tab)
  fdr <- p.adjust(results.df$pval, method="fdr")

  results.df <- cbind(results.df, fdr)
  final.results[[cur.status]] = results.df

  new.file <- paste(fn, status, "results.tsv", sep="_")
  write.table(results.df, file=new.file, quote=F, row.names=TRUE, sep="\t")
  
  sig.results.df <- with(results.df, results.df[pval < p.value.cutoff & fdr < fdr.cutoff, ])
  new.file <- paste(fn, status, "significant_results.tsv", sep="_")
  write.table(sig.results.df, file=new.file, quote=F, row.names=TRUE, sep="\t")
  cur.status <- cur.status + 1
}
# TODO: Change the following line
return(final.results)
}
