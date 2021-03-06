% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ORF_features.R
\name{isOverlapping}
\alias{isOverlapping}
\title{find frame for each orf relative to cds}
\usage{
isOverlapping(dists)
}
\arguments{
\item{dists}{a vector of distances between ORF and cds}
}
\value{
a logical vector
}
\description{
Input of this function, is the output of the function ORFDistToCds
See article:  10.1074/jbc.R116.733899
}
\examples{
grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1,10), "+"))
fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1,20), "+"))
dist <- ORFDistToCds(grl, fiveUTRs, extension = 0)
isOverlapping <- isOverlapping(dist)
}
\seealso{
Other features: \code{\link{ORFDistToCds}},
  \code{\link{ORFScores}}, \code{\link{OrfRankOrder}},
  \code{\link{RibosomeReleaseScore}},
  \code{\link{RibosomeStallingScore}},
  \code{\link{computeFeatures}},
  \code{\link{disengagementScore}}, \code{\link{entropy}},
  \code{\link{floss}}, \code{\link{fpkm_calc}},
  \code{\link{fpkm}}, \code{\link{fractionLength}},
  \code{\link{insideOutsideORF}}, \code{\link{isInFrame}},
  \code{\link{kozakSequenceScore}},
  \code{\link{subset_coverage}}, \code{\link{te}}
}
