% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ORF_features.R
\name{OrfRankOrder}
\alias{OrfRankOrder}
\title{Get the orf rank in transcripts}
\usage{
OrfRankOrder(grl)
}
\arguments{
\item{grl}{a \code{\link[GenomicRanges]{GRangesList}} object with ORFs}
}
\value{
a numeric vector of integers
}
\description{
ig. second orf _2 -> 2
}
\details{
See article:  10.1074/jbc.R116.733899
}
\seealso{
Other features: \code{\link{ORFDistToCds}},
  \code{\link{ORFScores}},
  \code{\link{RibosomeReleaseScore}},
  \code{\link{RibosomeStallingScore}},
  \code{\link{computeFeatures}},
  \code{\link{disengagementScore}}, \code{\link{entropy}},
  \code{\link{floss}}, \code{\link{fpkm_calc}},
  \code{\link{fpkm}}, \code{\link{fractionLength}},
  \code{\link{insideOutsideORF}}, \code{\link{isInFrame}},
  \code{\link{isOverlapping}},
  \code{\link{kozakSequenceScore}},
  \code{\link{subset_coverage}}, \code{\link{te}}
}
