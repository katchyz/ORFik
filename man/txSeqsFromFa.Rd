% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRanges_helpers.R
\name{txSeqsFromFa}
\alias{txSeqsFromFa}
\title{Get transcript sequence from a GrangesList and a faFile}
\usage{
txSeqsFromFa(grl, faFile, is.sorted = F)
}
\arguments{
\item{grl}{a GRangesList object}

\item{faFile}{FaFile used to find the transcripts,}

\item{is.sorted}{a speedup, if you know the ranges are sorted}
}
\value{
a DNAStringSet of the transcript sequences
}
\description{
Get transcript sequence from a GrangesList and a faFile
}
