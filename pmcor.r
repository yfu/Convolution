## Generate the SCNA significance plots of two types:
## whole-genome based: if a SCNA region spans exome boundaries, interpolate
## between the two exomes (by showing CN level in the intron segment at the
## same level as in the two flanking exons).
## "coverage-based": i.e., plotting only the region covered by the sequencing
## reads (in our case, the MCOR regions). This will be achieved by
## constructing the x-axis as the joining of all the regions.
## for both, we discussed having different options on how to join
## non-adjacent regions ('zero', 'interpolation', 'min').

## Function prototype:
## pmcor(file, base.on="coverage", join.method)
## file is the location of the file containing the mcor output +
## their genomic coordiantes. base.on should be "coverage" or "genome"
## joni.method should be "zero", "interpolate" or "min"


pmcor <- function(file, base.on="coverage", joining.method="interp"){
  library(ggplot2)
  # setwd("~/Challenge/results/2013-03-18/data/")
  mcor.with.coor.fn = file
  mcor <- read.table(mcor.with.coor.fn, header = FALSE)
  loc.pq <- mcor[,c(1, 2, 3, 5, 6)] # This is what I used to plot
  loc.pq <- mcor[,c(15, 16, 17, 13, 14)] # This is what I used to plot

  names(loc.pq) = c("chr", "start", "end", "pval", "fdr")
  # Fuzzy matching
  base.on = pmatch(base.on, c("coverage", "genome"))

  fdr.cutoff = c(-log(0.05), -log(0.25))
  fdr.color = c("blue", "red")
  
  if (base.on == 1){
  # Concatenate all the mcor coordiantes
    cur.pos = 1                             # 1-based coordinates
    loc.pq.new = loc.pq
    for ( i in seq(1, length(loc.pq[,1]))) {
      offset <- loc.pq.new[i, 2] - cur.pos
      loc.pq.new[i, 2] = loc.pq.new[i,2] - offset
      loc.pq.new[i, 3] = loc.pq.new[i,3] - offset
      cur.pos <- loc.pq.new[i,3] + 1
    }
    names(loc.pq.new) = c("chr", "start", "end", "pval", "fdr")
    temp <- rbind(data.frame(chr=loc.pq.new$chr, coor=loc.pq.new$start, fdr=loc.pq.new$fdr), data.frame(chr=loc.pq.new$chr, coor=loc.pq.new$start, fdr=loc.pq.new$fdr))

    # Due to the finite accuracy of the convolution method:
    temp$fdr[temp$fdr<0] <- min(abs(temp$fdr))

    # Generate a dataframe to plot the two lines for cutoffs
    df.cutoff=data.frame(x=c(-Inf, Inf, -Inf, Inf), y=rep(fdr.cutoff, each=2), cutoff=c("0.05", "0.05", "0.25", "0.25"))
      
    pdf(paste(file, "_coverage_based.pdf", sep=""))
    print(ggplot(temp, aes(coor, -log(fdr))) + geom_line() + scale_x_continuous("MCOR coordinate") + scale_y_continuous("-log(fdr)") 
          ## + ggtitle("Coverage-based Plot")
          + geom_line(data=df.cutoff, aes(x=x, y=y, color=cutoff)))
    dev.off()
  }
    else if (base.on == 2){
      joining.method = pmatch(joining.method, c("interp", "zero", "min", "max"))
      if (joining.method == 1){
        temp2 <- rbind(data.frame(chr=loc.pq$chr, coor=loc.pq$start, fdr=loc.pq$fdr), data.frame(chr=loc.pq$chr, coor=loc.pq$start, fdr=loc.pq$fdr))
        temp2$fdr[temp2$fdr < 0] = min(abs(temp2$fdr))
        pdf(paste(file, "_genome_based_interp.pdf", sep=""))
        print(ggplot(temp2, aes(coor,-log(fdr))) + geom_line() + ggtitle("Genome-baesd Plot (Interpolation)") + geom_hline(y=fdr.cutoff))
        dev.off()
      }
      else if (joining.method == 2) {
        # Assign 0 to those regions with no values
        loc.pq.2 <- data.frame()
        previous.pos = 0 # For the simplicity of the loop
        for ( i in seq(1, length(loc.pq[,1]))) {
          cur.chr <- loc.pq$chr[i]
          # print(cur.chr)
          cur.start <- loc.pq$start[i]
          cur.end <- loc.pq$end[i]
          if (cur.start - previous.pos != 1) {
            # This means that the two segments are not adjacent
            loc.pq.2 = rbind(loc.pq.2, list(cur.chr, previous.pos + 1, cur.start - 1, 0.99, 0.99))
            # Inefficient, but I cannot figure out a way to initialize
            # an empty dataframe with colnames
            names(loc.pq.2) <- c("chr", "start", "end", "pval", "fdr")
          }
          # print(loc.pq[i,])
          # print(loc.pq.2)
          loc.pq.2 = rbind(loc.pq.2, loc.pq[i,])
          previous.pos = cur.end
        }
        temp.2 <- rbind(data.frame(chr=loc.pq.2$chr, coor=loc.pq.2$start, fdr=loc.pq.2$fdr), data.frame(chr=loc.pq.2$chr, coor=loc.pq.2$end, fdr=loc.pq.2$fdr))
        temp.2$fdr[temp.2$fdr < 0] = min(abs(temp.2$fdr))
        pdf(paste(file, "_genome_based_zero.pdf", sep=""))
        print(ggplot(temp.2, aes(coor, -log(fdr))) + geom_line() + ggtitle("Genome-based Plot (Zero)") + geom_hline(y=fdr.cutoff))
        dev.off()
      }
      else if (joining.method == 3) {
        # Assign the smaller value of two adjacent peaks to the regions with
        # no values
        loc.pq.3 <- data.frame()
        previous.pos <- 0 # For the simplicity of the loop
        previous.pval <- 0.99 # Same as the above
        previous.fdr <- 0.99 # Same as the above
        
        for ( i in seq(1, length(loc.pq[,1]))) {
          cur.chr <- loc.pq$chr[i]
          # print(cur.chr)
          cur.start <- loc.pq$start[i]
          cur.end <- loc.pq$end[i]
          cur.pval <- loc.pq$pval[i]
          cur.fdr <- loc.pq$fdr[i]
          if (cur.start - previous.pos != 1) {
            # This means that the two segments are not adjacent
            loc.pq.3 = rbind(loc.pq.3, list(cur.chr, previous.pos + 1, cur.start - 1, max(previous.pval, cur.pval), max(previous.fdr, cur.fdr)))
            # Inefficient, but I cannot figure out a way to initialize
            # an empty dataframe with colnames
            names(loc.pq.3) <- c("chr", "start", "end", "pval", "fdr")
          }
          # print(loc.pq[i,])
          # print(loc.pq.2)
          loc.pq.3 = rbind(loc.pq.3, loc.pq[i,])
          previous.pos = cur.end
          previous.pval <- cur.pval
          previous.fdr <- cur.fdr
        }
        temp.3 <- rbind(data.frame(chr=loc.pq.3$chr, coor=loc.pq.3$start, fdr=loc.pq.3$fdr), data.frame(chr=loc.pq.3$chr, coor=loc.pq.3$end, fdr=loc.pq.3$fdr))
        temp.3$fdr[temp.3$fdr < 0] = min(abs(temp.3$fdr))
        pdf(paste(file, "_genome_based_min.pdf", sep=""))
        print(ggplot(temp.3, aes(coor, -log(fdr))) + geom_line() + ggtitle("Genome-based Plot (Minimum)") + geom_hline(y=fdr.cutoff))
        dev.off()
      }
      else if (joining.method == 4) {
              # Assign the smaller value of two adjacent peaks to the regions with
        # no values
        loc.pq.4 <- data.frame()
        previous.pos <- 0 # For the simplicity of the loop
        previous.pval <- 0.99 # Same as the above
        previous.fdr <- 0.99 # Same as the above
      
        for ( i in seq(1, length(loc.pq[,1]))) {
          cur.chr <- loc.pq$chr[i]
                                        # print(cur.chr)
          cur.start <- loc.pq$start[i]
          cur.end <- loc.pq$end[i]
          cur.pval <- loc.pq$pval[i]
          cur.fdr <- loc.pq$fdr[i]
          if (cur.start - previous.pos != 1) {
          # This means that the two segments are not adjacent
            loc.pq.4 = rbind(loc.pq.4, list(cur.chr, previous.pos + 1, cur.start - 1, min(previous.pval, cur.pval), min(previous.fdr, cur.fdr)))
          # Inefficient, but I cannot figure out a way to initialize
          # an empty dataframe with colnames
            names(loc.pq.4) <- c("chr", "start", "end", "pval", "fdr")
          }
        # print(loc.pq[i,])
        # print(loc.pq.2)
          loc.pq.4 = rbind(loc.pq.4, loc.pq[i,])
          previous.pos = cur.end
          previous.pval <- cur.pval
          previous.fdr <- cur.fdr
        }
        temp.4 <- rbind(data.frame(chr=loc.pq.4$chr, coor=loc.pq.4$start, fdr=loc.pq.4$fdr), data.frame(chr=loc.pq.4$chr, coor=loc.pq.4$end, fdr=loc.pq.4$fdr))
        temp.4$fdr[temp.4$fdr < 0] = min(abs(temp.4$fdr))
        pdf(paste(file, "_genome_based_max.pdf", sep=""))
        print(ggplot(temp.4, aes(coor, -log(fdr))) + geom_line() + ggtitle("Genome-based Plot (Maximum)") + geom_hline(y=fdr.cutoff))
        dev.off()
      }
      else {warning("Wong parameter!")}
    }
}
