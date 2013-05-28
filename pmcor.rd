\name{pmcor}
\alias{pmcor}
\title{Generate the SCNA Significance Plots of Two Types}
\description{
  Generate the SCNA significance plots of two types:

  whole-genome based:
  
  if a SCNA region spans exome boundaries, interpolate between the two
  exomes (by showing CN level in the intron segment at the same level as
  in the two flanking exons).

  coverage-based: i.e., plotting only the region covered by the
  sequencing reads (in our case, the MCOR regions). This will be
  achieved by constructing the x-axis as the joining of all the regions.

  For both, we discussed having different options on how to join
  non-adjacent regions ('zero', 'interpolation', 'min').
}
\usage{
pmcor(file, base.on="coverage", join.method)
}
\arguments{
  \item{file}{file is the location of the file containing the mcor
    output their genomic coordiantes.}
  \item{base.on}{base.on should be "coverage" or "genome"}
  \item{join.method}{join.method should be "zero", "interpolate" or "min"}
}
\seealso{
  \code{\link{pmatch}}, \code{\link{ggplot}}
}
\examples{
## Using the default parameters
pmcor(file, base.on="coverage", join.method)

## Using other "bases" and joining method
pmcor(file, base.on="coverage", joining.method="min")
}
\keyword{file}
\references{
  \url{https://projects.zoho.com/portal/montilab#wiki/465555000000032025/Scoring-SCNAs.html}
}
\author{Yu Fu \email{f@yfu.me}}