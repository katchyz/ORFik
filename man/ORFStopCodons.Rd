% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ORFs_helpers.R
\name{ORFStopCodons}
\alias{ORFStopCodons}
\title{Get the Stop codons(3 bases) from a GRangesList of orfs grouped by orfs}
\usage{
ORFStopCodons(grl, is.sorted = FALSE)
}
\arguments{
\item{grl}{a \code{\link[GenomicRanges]{GRangesList}} object}

\item{is.sorted}{a boolean, a speedup if you know the ranges are sorted}
}
\value{
a GRangesList of stop codons, since they might be split on exons
}
\description{
In ATGTTTTGC, get the positions TGC.
It takes care of exons boundaries, with exons < 3 length.
}
