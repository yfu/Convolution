\name{pconvolve}
\alias{pconvolve}
\title{Calculate p-values for Copynumber Variations across Samples}
\description{
  The amplification and deletion scores for a given loci i are computed
  as equations which can be found in the 1/30/2013 meeting
  minutes. Conceptually, the computation of the p-value is based on
  permutation, where by each "column" of the matrix is shuffled multiple
  times, abd a distribution of permuted score is computed for each locus
  and compared to the observed score for that locus. See the equation in
  the meeting minutes. In practice, the p-value can be computed in a
  closed form using the concept of convolution. In R, this function uses
  \code{convolve} to calculate the p-value analytically.
}
\usage{
pconvolve(fn, tau.amp=0.5, tau.del=-0.5, p.value.cutoff=0.05,
fdr.cutoff=0.25, plot.pdf=TRUE)
}
\arguments{
  \item{fn}{The file containing the log2 ratios for all segments across
    samples. Each row represents a segment and each column is from one
    sample.}
  \item{tau.amp}{log2 ratios above this value are kept.}

  \item{tau.del}{log2 ratios below this value are kept}
  \item{p.value.cutoff}{Cutoff for p-values.}
  \item{fdr.cutoff}{Cutoff for FDRs}
  \item{plot.pdf}{If this is set to \code{TRUE}, then this function will
  skip plotting on the default device and print it to a PDF file.}
}
\seealso{
  \code{\link{convolve}}, \code{\link{p.adjust}}
}
\examples{
  ## Calculate the p-value and FDR using the default parameters
  pconvolve("matrix10_10.T.txt")

  ## Calculate the p-value and FDR using user defined parameters
  pconvolve("matrix", tau.amp=1, tau.del=-1, p.value.cutoff=0.02, fdr.c
  utoff=0.05, plot.pdf=TRUE)
}
\keyword{file}
\references{
  \url{https://projects.zoho.com/portal/montilab#wiki/465555000000032025/Scoring-SCNAs.html}
  
  \url{http://sas.uwaterloo.ca/~dlmcleis/s901/chapt5.pdf}
}
\author{Yu Fu \email{f@yfu.me}}