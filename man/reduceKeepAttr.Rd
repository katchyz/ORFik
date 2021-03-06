% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRanges_helpers.R
\name{reduceKeepAttr}
\alias{reduceKeepAttr}
\title{Reduce GRanges / GRangesList}
\usage{
reduceKeepAttr(grl, keep.names = FALSE, drop.empty.ranges = FALSE,
  min.gapwidth = 1L, with.revmap = FALSE, with.inframe.attrib = FALSE,
  ignore.strand = FALSE)
}
\arguments{
\item{grl}{a \code{\link[GenomicRanges]{GRangesList}} or GRanges object}

\item{keep.names}{(FALSE) keep the names and meta columns of the GRangesList}

\item{drop.empty.ranges}{(FALSE) if a group is empty (width 0), delete it.}

\item{min.gapwidth}{(1L) how long gap can it be to say they belong together}

\item{with.revmap}{(FALSE) return info on which mapped to which}

\item{with.inframe.attrib}{(FALSE) For internal use.}

\item{ignore.strand}{(FALSE), can different strands be reduced together.}
}
\value{
A reduced GRangesList
}
\description{
Extends function \code{\link[GenomicRanges]{reduce}}
by trying to keep names and meta columns if it is a
GRangesList. It also does not loose sorting for GRangesList,
since original reduce sorts all by ascending.
If keep.names == F, it's just the normal GenomicRanges::reduce
with sorting negative strands descending for GRangesList
}
\details{
Only tested for orfik, might not work for other naming conventions.
}
\examples{
ORF <- GRanges(seqnames = "1",
ranges = IRanges(start = c(1, 2, 3),
                end = c(1, 2, 3)),
strand = "+")
##For GRanges
reduceKeepAttr(ORF, keep.names = TRUE)
## For GRangesList
grl <- GRangesList(tx1_1 = ORF)
reduceKeepAttr(grl, keep.names = TRUE)
##Only 1 GRanges object in GRangesList returned
}
